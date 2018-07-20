"""
Prepares dataset from reconstruction steps (FreeSurfer, MRtrix, TVB-ification)
for statistical models with Stan.

"""

import os
import glob
import json
import logging
import pickle
import zipfile

import numpy as np
import scipy.signal
import mne

try:
    import matplotlib as mpl
    mpl.use('agg')
except:
    pass
import pycmdstan as pcs


def compute_raw_slp(raw, cfg):
    nperseg: int = cfg.get('nperseg', 4 * int(raw.info['sfreq']))
    hpf: float = cfg.get('hpf', 10.0)
    lpf: float = cfg.get('lpf', 100.0)
    Cs = []
    for y in raw._data:
        F, T, C = scipy.signal.spectrogram(
            y, raw.info['sfreq'], nperseg=nperseg)
        fmask = np.ones(F.shape, 'bool')
        if hpf:
            fmask *= F > hpf
        if lpf:
            fmask *= F < lpf
        Cs.append(np.log(C[fmask].sum(axis=0)))
    Cs = np.array(Cs)
    return Cs


def read_gain(subj_proc_dir):
    np_fname = os.path.join(subj_proc_dir,
                            'elec/gain_inv-square.destrieux.txt')
    return np.loadtxt(np_fname)


def _process_one_fif(js, cfg):
    exclude = js['bad_channels'] + js['non_seeg_channels']
    fif_fname = os.path.join(os.path.dirname(js['_source']), js['filename'])
    raw = mne.io.Raw(fif_fname, verbose='WARNING')
    picks = set(raw.ch_names) - set(exclude)
    assert js['onset'] is not None and js['termination'] is not None
    raw.crop(tmin=js['onset'], tmax=js['termination'])
    raw.load_data()
    raw.pick_channels(picks)
    slp = compute_raw_slp(raw, cfg)
    return picks, slp, raw.ch_names


def _process_one_bids_vhdr(vhdrname, cfg):
    raw = mne.io.read_raw_brainvision(vhdrname, preload=True)
    raw = raw.pick_types(meg=False, eeg=True)
    slp = compute_raw_slp(raw, cfg)
    return set(raw.ch_names), slp, raw.ch_names


def _load_js(js_fname: os.PathLike) -> dict:
    with open(js_fname, 'r') as fd:
        js = json.load(fd)
    js['_source'] = js_fname
    return js


def _is_seizure(js) -> bool:
    return js['type'] == 'Spontaneous seizure'


def _read_all_jsons(subj_proc_dir):
    pattern: str = os.path.join(subj_proc_dir, 'seeg/fif/*.json')
    matches: [os.PathLike] = glob.glob(pattern)
    for match in matches:
        yield _load_js(match)


def _many_picks_intersection(many_picks: [set]):
    first, *rest = many_picks
    for next in rest:
        first = first.intersection(next)
    return first


def _is_bids(subj_proc_dir):
    # return True if no fifs found (not great, but..)
    return len(list(_read_all_jsons(subj_proc_dir))) == 0


def _find_vhdrs(subj_proc_dir):
    subj_id = os.path.basename(subj_proc_dir)
    raw_path = os.path.join(subj_proc_dir, '..', '..', '0-Raw', subj_id)
    vhdr_pattern = os.path.join(raw_path, '*/ieeg/*run-01*.vhdr')
    return glob.glob(vhdr_pattern)


def read_all_seeg_data(subj_proc_dir, cfg: dict):
    # read all datasets
    data = []
    if _is_bids(subj_proc_dir):
        for vhdr in _find_vhdrs(subj_proc_dir):
            data.append(_process_one_bids_vhdr(vhdr, cfg))
    else:
        for js in _read_all_jsons(subj_proc_dir):
            if _is_seizure(js):
                data.append(_process_one_fif(js, cfg))
    picks, slps, chs = zip(*data)
    # find intersection of channels across datasets
    picks = _many_picks_intersection(picks)
    # remap data consistently
    slps_ = []
    is_first_ = []  # 1 if first samp of seizure
    for slp, ch in zip(slps, chs):
        im = np.array([i for i, c in enumerate(ch) if c in picks])
        is_first = np.zeros(slp.shape[1], 'i')
        is_first[0] = 1
        slps_.append(slp[im])
        is_first_.append(is_first)
    slp = np.concatenate(slps_, axis=1)
    is_first = np.concatenate(is_first_)
    chs = [_ for _ in chs[0] if _ in picks]
    return picks, slp, chs, is_first


def read_seeg_xyz(subj_proc_dir):
    lines = []
    fname = os.path.join(subj_proc_dir, 'elec/seeg.xyz')
    with open(fname, 'r') as fd:
        for line in fd.readlines():
            name, *sxyz = line.strip().split()
            xyz = [float(_) for _ in sxyz]
            lines.append((name, xyz))
    return lines


def read_weights(subj_proc_dir):
    roi_names = []
    fname = os.path.join(subj_proc_dir, 'tvb/connectivity.destrieux.zip')
    with zipfile.ZipFile(fname) as zf:
        with zf.open('weights.txt') as fd:
            weights = np.loadtxt(fd)
        with zf.open('centres.txt', 'r') as fd:
            for line in fd.readlines():
                roi_name, *_ = line.decode('ascii').strip().split(' ')
                roi_names.append(roi_name)
    weights_triu = weights[np.triu_indices(weights.shape[0], 1)]
    return weights_triu, roi_names


def build_data(subj_proc_dir, cfg=None):
    cfg = cfg or {}
    counts_triu, roi_names = read_weights(subj_proc_dir)
    gain = read_gain(subj_proc_dir)
    seeg_xyz = read_seeg_xyz(subj_proc_dir)
    picks, slp, ch_names, is_first = read_all_fif(subj_proc_dir, cfg)
    gain_pick = np.array(
        [i for i, (label, _) in enumerate(seeg_xyz) if label in picks])
    gain = gain[gain_pick]
    data = dict(
        nn=gain.shape[1],
        ns=gain.shape[0],
        nt=slp.shape[1],
        gain=gain,
        seeg=slp,
        is_first=is_first,
        counts_triu=counts_triu)
    return data


def retro_proc_dir(id):
    return os.path.join('/home/vep/RetrospectivePatients/1-Processed', id)


def retro_ids():
    return [
        os.path.basename(_)
        for _ in glob.glob(retro_proc_dir('id*'))
    ]


def ensure_vep_topic_dir(subj_proc_dir):
    path = os.path.join(subj_proc_dir, 'vep')
    if not os.path.exists(path):
        os.mkdir(path)


def build_and_save_one(subj_proc_dir, cfg=None):
    ensure_vep_topic_dir(subj_proc_dir)
    cfg_ = ''.join([f'-{k}_{v}' for k, v in cfg.items()]) if cfg else ''
    Rfname = os.path.join(subj_proc_dir, 'vep', f'data{cfg_}.R')
    if os.path.exists(Rfname):
        print(f'skipping existing {Rfname}')
        return
    data = build_data(subj_proc_dir, cfg)
    plot_dataset(subj_proc_dir, data)
    pcs.rdump(Rfname, data)


def build_and_save_all_retro():
    for id in retro_ids():
        print(id)
        build_and_save_one(retro_proc_dir(id))


def plot_dataset(subj_proc_dir, data):
    import matplotlib as mpl
    mpl.use('agg')
    import pylab as pl
    # TODO more informative ticks/labels
    pl.figure(figsize=(10, 6))
    pl.subplot(121)
    pl.imshow(data['seeg'], vmin=8, vmax=14, aspect='auto', cmap='binary')
    for ft in np.argwhere(data['is_first'])[:, 0]:
        pl.axvline(ft, color='r')
    pl.title('Concatenated Seizures')
    pl.xlabel('Time Window (2 s)')
    pl.ylabel('sEEG Monopolar Contact')
    pl.grid(1)
    pl.subplot(222)
    pl.imshow(data['gain'], aspect='auto', cmap='binary')
    pl.title('Gain Matrix (Destrieux -> sEEG)')
    pl.xlabel('Destrieux ROI')
    pl.ylabel('sEEG Monopolar Contact')
    pl.grid(1)
    pl.subplot(224)
    pl.hist(np.log(data['counts_triu'] + 1), 30, color='k')
    pl.title('log Track Count Histogram')
    pl.xlabel('log Track Count')
    pl.ylabel('# log Track Count')
    pl.grid(1)
    pl.tight_layout()
    pl.savefig(os.path.join(
        subj_proc_dir, 'vep', 'data.png'))

# testing
# spd = retro_proc_dir('id023_br')
# data = build_data(spd, {'hpf': 20.0, 'lpf': None})
# plot_dataset(spd, data)
# print('done')

if __name__ == '__main__':
    subj_proc_dir, = sys.argv[1:]
    build_and_save_one(subj_proc_dir)
    
