#!/usr/bin/env python

import json
import os
import re
import sys

import pandas as pd


def get_sec(time_str):
    if type(time_str) == float:
        # Already in seconds
        return time_str

    h, m, s = time_str.split(':')
    return int(h)*3600 + int(m)*60 + float(s)


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

        print(string)

        # A'1
        match = re.match("^([A-Za-z]+[']*)([0-9]+)$", string)
        if match:
            new_list.append(string)
            continue

        # A'1-10
        match = re.match("^([A-Za-z]+[']*)([0-9]+)-([0-9]+)$", string)
        print(match)
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


def gen_sidecar_files(xlsx_file, output_direc, convert_format=None):
    df = pd.read_excel(xlsx_file, sheet_name="Recordings")
    add_same_occurence_index(df, 'File')

    # TODO: check that everything is where expected

    for index, row in df.iterrows():
        if pd.notna(row['File']):
            onset = get_sec(row['Onset']) if pd.notna(row['Onset']) else None
            termination = get_sec(row['Termination']) if pd.notna(row['Termination']) else None

            bad_channels = row['Bad channels']
            if pd.isna(bad_channels) or bad_channels == 0:
                bad_channels = []
            else:
                bad_channels = expand_channels([a.strip() for a in bad_channels.split(",")])

            filename = row['File']
            if convert_format is not None:
                root, ext = os.path.splitext(filename)
                if ext in convert_format:
                    filename = root + convert_format[ext]

            data = {
                'filename': filename,
                'onset': onset,
                'termination': termination,
                'bad_channels': bad_channels,
                'type': row['Seizure type']
            }

            basename = os.path.splitext(row['File'])[0]

            if not row['_File_repeated']:
                jsonfile = basename + ".json"
            else:
                jsonfile = basename + "_" + str(row['_File_index']) + ".json"

            sidecar_name = os.path.join(output_direc, jsonfile)
            with open(sidecar_name, 'w') as outfile:
                json.dump(data, outfile, indent=4)


if __name__ == '__main__':
    gen_sidecar_files(sys.argv[1], sys.argv[2], {'.EEG': '.raw.fif', '.eeg': '.raw.fif'})
