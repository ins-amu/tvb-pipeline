#!/usr/bin/env python


import datetime
import json
import logging
import os
import re
import sys
import zipfile

import nibabel as nib
import numpy as np
import pandas as pd

from .elecs import Contacts
from . import nifti
from . import read_eeg

DATADIR = os.path.join(os.path.dirname(__file__), "data")

def get_sec(time):
    if pd.isna(time):
        return None
    elif type(time) == float:
        # Already in seconds
        return time
    elif type(time) == datetime.time:
        return datetime.timedelta(hours=time.hour,
                                 minutes=time.minute,
                                 seconds=time.second,
                                 microseconds=time.microsecond).total_seconds()
    elif type(time) == str:
        h, m, s = time.split(':')
        return int(h)*3600 + int(m)*60 + float(s)
    else:
        raise ValueError("Unexpected time type: %s" % type(time))



def add_same_occurence_index(df, column):
    df['_%s_repeated' % column] = False
    df['_%s_index' % column] = 1

    for key in pd.unique(df[column]):
        if pd.isna(key):
            continue

        subdf = df[df[column] == key]
        if len(subdf) > 1:
            for i, (index, row) in enumerate(subdf.iterrows()):
                df.loc[index, '_%s_repeated' % column] = True
                df.loc[index, '_%s_index' % column] = i + 1


def expand_channels(ch_list):
    ch_list = [a.replace("’", "'") for a in ch_list]

    new_list = []
    for string in ch_list:
        if not string.strip():
            continue

        # A'1
        match = re.match("^([A-Za-z]+[']*)([0-9]+)$", string)
        if match:
            new_list.append(string)
            continue

        # A'1-10
        match = re.match("^([A-Za-z]+[']*)([0-9]+)-([0-9]+)$", string)
        if match:
            name, fst_idx, last_idx = match.groups()
            new_list.extend([name + str(i) for i in range(int(fst_idx), int(last_idx) + 1)])
            continue

        # A'1-A10
        match = re.match("^([A-Za-z]+[']*)([0-9]+)-([A-Za-z]+[']*)([0-9]+)$", string)
        if match:
            name1, fst_idx, name2, last_idx = match.groups()
            if name1 == name2:
                new_list.extend([name1 + str(i) for i in range(int(fst_idx), int(last_idx) + 1)])
                continue

        print("expand_channels: Cannot parse this: %s" % string)

    return new_list


def get_bad_channels(cell_value):
    if pd.isna(cell_value) or cell_value == 0:
        return []
    else:
        return expand_channels([a.strip() for a in re.split("[,;.]", cell_value)])


def get_converted_filename(filenames):
    TARGET_FMT = '.raw.fif'

    if type(filenames) == str:
        filenames = [filenames]

    roots = [os.path.splitext(filename.strip())[0] for filename in filenames]
    return "_".join(roots) + TARGET_FMT


def get_sidecar_name(filename, is_repeated, file_index):
    known_extensions = ['.eeg', '.raw.fif']

    basename = os.path.splitext(filename)[0]
    for ext in known_extensions:
        if filename[-len(ext):] == ext:
            basename = filename[:-len(ext)]

    if not is_repeated:
        return basename + ".json"
    else:
        return basename + "_" + str(file_index) + ".json"


def convert_recordings(xlsx_file, seeg_rec_dir, contacts_file, output_direc):
    df = pd.read_excel(xlsx_file, sheet_name="Recordings")
    add_same_occurence_index(df, 'File')

    contact_names = np.genfromtxt(contacts_file, usecols=(0,), dtype=str)

    iterrows = df.iterrows()
    for index, row in iterrows:
        if pd.notna(row['File']):
            if row['Termination'] == '>':
                # Merge two files
                index2, row2 = next(iterrows)
                assert row2['Onset'] == '<'

                onset, termination = get_sec(row['Onset']), get_sec(row2['Termination'])
                bad_channels1 = get_bad_channels(row['Bad channels'])
                bad_channels2 = get_bad_channels(row2['Bad channels'])
                orig_filename1 = row['File']
                orig_filename2 = row2['File']
                conv_filename = get_converted_filename([orig_filename1, orig_filename2])
                jsonname = get_sidecar_name(conv_filename, False, None)

                eeg = read_eeg.EEG(os.path.join(seeg_rec_dir, orig_filename1)).to_fif()
                eeg2 = read_eeg.EEG(os.path.join(seeg_rec_dir, orig_filename2)).to_fif()

                assert len(eeg.ch_names) == len(eeg2.ch_names)
                assert all([ch1 == ch2 for ch1, ch2 in zip(eeg.ch_names, eeg2.ch_names)])
                assert eeg.info['sfreq'] == eeg2.info['sfreq']
                assert row['Seizure type'] == row2['Seizure type']

                bad_channels = sorted(list(set(bad_channels1 + bad_channels2)))
                termination += (eeg.n_times - 1) * (1./eeg.info['sfreq'])

                eeg.append(eeg2)


            else:
                onset, termination = get_sec(row['Onset']), get_sec(row['Termination'])
                bad_channels = get_bad_channels(row['Bad channels'])
                orig_filename = row['File']
                conv_filename = get_converted_filename(row['File'])
                jsonname = get_sidecar_name(row['File'], row['_File_repeated'], row['_File_index'])
                eeg = read_eeg.EEG(os.path.join(seeg_rec_dir, orig_filename)).to_fif()

            eeg.save(os.path.join(output_direc, conv_filename), overwrite=True)

            data = {
                'filename': conv_filename,
                'onset': onset,
                'termination': termination,
                'bad_channels': bad_channels,
                'non_seeg_channels': sorted(list(set(eeg.ch_names) - set(contact_names))),
                'type': row['Seizure type']
            }

            with open(os.path.join(output_direc, jsonname), 'w') as outfile:
                json.dump(data, outfile, indent=4)


def get_ez_from_regions(xlsx_file, region_names):
    """Return list of indices of EZ regions given in the patient spreadsheet"""

    LH_NAMES_IND = 9
    LH_EZ_IND = 10
    RH_NAMES_IND = 12
    RH_EZ_IND = 13

    df = pd.read_excel(xlsx_file, sheet_name="EZ hypothesis and EI", header=1)

    ez_names = []
    for names_ind, ez_ind in [(LH_NAMES_IND, LH_EZ_IND), (RH_NAMES_IND, RH_EZ_IND)]:
        names_col = df.iloc[:, names_ind]
        mask = names_col.notnull()
        names = names_col[mask]
        ez_mask = df.iloc[:, ez_ind][mask].astype(str) == 'YES'
        ez_names.extend(names[ez_mask])

    return [region_names.index(name) for name in ez_names]


def get_ez_from_contacts(xlsx_file, contacts_file, label_volume_file):
    """Return list of indices of EZ regions given by the EZ contacts in the patient spreadsheet"""

    CONTACTS_IND = 6
    EZ_IND = 7

    df = pd.read_excel(xlsx_file, sheet_name="EZ hypothesis and EI", header=1)

    ez_contacts = []
    contacts_col = df.iloc[:, CONTACTS_IND]
    mask = contacts_col.notnull()
    contacts_names = contacts_col[mask]
    ez_mask = df.iloc[:, EZ_IND][mask] == 'YES'
    ez_contacts.extend(contacts_names[ez_mask])

    contacts = Contacts(contacts_file)
    label_vol = nib.load(label_volume_file)

    ez_inds = []
    for contact in ez_contacts:
        coords = contacts.get_coords(contact)
        region_ind = nifti.point_to_brain_region(coords, label_vol, tol=3.0) - 1   # Minus one to account for the shift
        if region_ind != -1:
            ez_inds.append(region_ind)

    return ez_inds


def get_nregions(tvb_zipfile):
    with zipfile.ZipFile(tvb_zipfile) as zf:
        with zf.open("centres.txt") as fl:
            region_names = list(np.genfromtxt(fl, usecols=(0,), dtype=str))
    return len(region_names)


def save_ez_hypothesis(xlsx_file, tvb_zipfile, contacts_file, label_volume_file_dk, output_file,
                       label_volume_file_trg=None):
    """
    Extract the EZ hypothesis from the xlsx file and save it to plain text file.

    Args:
    xlsx_file (str):              Path to the patient excel file.
    tvb_zipfile (str):            Path to the TVB zipfile for target parcellation.
    contacts_file (str):          Path to the text file with contact coordinates.
    label_volume_file_dk (str):   Path to the nifti label volume file for Desikan-Killiany parcellation.
                                  DK parcellation is used for the EZ specification and the parcellation is thus needed
                                  even if the target parcellation is different.
    output_file (str):            Path to the generated text file with EZ hypothesis.
    label_volume_file_trg (str):  (Optional) Path to the nifti label volume file for the target parcellation.
                                  If absent, Desikan-Killiany is thought to be the desired target parcellation.
    """
    region_names_dk = list(np.genfromtxt(os.path.join(DATADIR, "region_names.dk.txt"), usecols=(0,), dtype=str))
    nreg_dk = len(region_names_dk)

    nreg_trg = get_nregions(tvb_zipfile)

    # Epileptogenic regions in DK parcellation
    ez_inds_dk_from_regions = get_ez_from_regions(xlsx_file, region_names_dk)

    ez_hyp_dk = np.zeros(nreg_dk, dtype=int)
    ez_hyp_dk[ez_inds_dk_from_regions] = 1

    # Translate to the target parcellation if needed
    if label_volume_file_trg is not None:
        ez_hyp_trg = nifti.translate_ez_hypothesis(label_volume_file_dk, label_volume_file_trg, ez_hyp_dk, nreg_trg)
    else:
        ez_hyp_trg = ez_hyp_dk
        label_volume_file_trg = label_volume_file_dk

    # Epileptogenic regions from contact specification (parcellation independent)
    ez_inds_trg_from_contacts = get_ez_from_contacts(xlsx_file, contacts_file, label_volume_file_trg)
    ez_hyp_trg[ez_inds_trg_from_contacts] = 1

    np.savetxt(output_file, ez_hyp_trg, fmt='%i')



if __name__ == '__main__':
    loglevel = logging.INFO
    if os.environ.get('VERBOSE', False):
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel)
    cmd = sys.argv[1]
    eval(cmd)(*sys.argv[2:])
