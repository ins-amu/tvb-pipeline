#!/usr/bin/env python3


from seegrecording import SeegRecording

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from math import floor, log10

import sys
import os
import mne
import glob
import json


def pow10floor(x):
    return 10**floor(log10(x))

def plot_traces(seeg, filename, representation='monopolar', interval=(None, None)):
    matplotlib.rcParams.update({'font.size': 24, 'font.family': 'monospace'})

    plt.figure(figsize=(48, 48))

    if representation in ['monopolar', 'avgref']:
        data = seeg.get_data()
        names = seeg.get_channel_names()
    elif representation == 'bipolar':
        data = seeg.get_data_bipolar()
        names = seeg.get_channel_names_bipolar()
    else:
        raise ValueError("Unexpected representation choice: %s" % representation)

    tid1, tid2 = seeg.interval_to_index(interval)
    data = data[:, tid1:tid2]
    time = seeg.t[tid1:tid2]

    maxrange = 0
    for i in range(data.shape[0]):
        data[i, :] -= np.mean(data[i, :])
        contact_range = np.max(data[i, :]) - np.min(data[i, :])
        maxrange = max(maxrange, contact_range)
    data /= maxrange

    if representation == 'avgref':
        data -= np.mean(data, axis=0)
        representation = 'monopolar'

    nchannels = data.shape[0]
    for i in range(nchannels):
        plt.plot(time, 3*data[nchannels - i - 1, :] + i, 'b', lw=0.4)
    plt.yticks(np.r_[:nchannels], reversed(["%-9s" % name for name in names]))

    plt.gca().xaxis.set_tick_params(labeltop='on')

    plt.xlabel("t [s]")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_spectrogram(seeg, filename, representation='monopolar', interval=(None, None)):
    matplotlib.rcParams.update({'font.size': 14, 'font.family': 'monospace'})

    if representation in ['monopolar', 'avgref']:
        data = seeg.get_data()
        electrodes = seeg.electrodes
        ch_names = seeg.get_channel_names()
    elif representation == 'bipolar':
        data = seeg.get_data_bipolar()
        electrodes = seeg.electrodes_bipolar
        ch_names = seeg.get_channel_names_bipolar()

    if representation == 'avgref':
        for i in range(data.shape[0]):
            data[i, :] -= np.mean(data[i, :])
        data -= np.mean(data, axis=0)
        representation = 'monopolar'


    tid1, tid2 = seeg.interval_to_index(interval)
    data = data[:, tid1:tid2]

    nelec = len(electrodes.keys())
    max_nchannels = max([len(chs) for chs in electrodes.values()])
    plt.figure(figsize=(5*max_nchannels, 2*nelec))

    caxs = []
    maxval = 0
    freqs = np.linspace(1., 60., 30)
    for i, (elec, inds) in enumerate(electrodes.items()):
        tfr = mne.time_frequency.tfr_array_multitaper(np.expand_dims(data[inds, :], 0),
                                                      seeg.sampling_rate,
                                                      freqs,
                                                      zero_mean=False,
                                                      output='avg_power')

        maxval = max(maxval, np.max(tfr[:, :, :]))

        for j, ind in enumerate(inds):
            ax = plt.subplot(nelec, max_nchannels, i*max_nchannels + j + 1)
            cax = ax.imshow(tfr[j, :, :], aspect='auto', origin='lower',
                            extent=[seeg.t[tid1], seeg.t[tid2], freqs[0], freqs[-1]],
                            norm=matplotlib.colors.LogNorm(1, 1e4),
                            cmap='jet'
            )
            plt.gcf().colorbar(cax)

            caxs.append(cax)
            ax.set_xlabel("t [s]")
            ax.set_ylabel("f [Hz]")
            ax.set_title(ch_names[ind])

    for cax in caxs:
        vmax = pow10floor(maxval)
        vmin = vmax*1e-4
        cax.set_clim([vmin, vmax])

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def process_seizure(name, seeg, interval, out_dir):
    if len(interval) == 0:
        interval = (None, None)

    plot_traces(seeg, os.path.join(out_dir, "traces_bipolar_%s.png" % name), representation='bipolar', interval=interval)
    plot_traces(seeg, os.path.join(out_dir, "traces_avgref_%s.png" % name), representation='avgref', interval=interval)





def plot_seeg_recording(jsonname, representation, imgname):
    TEMPORAL_MARGIN = 20

    with open(jsonname, 'r') as fl:
        metadata = json.load(fl)

    filename = os.path.join(os.path.dirname(jsonname), metadata['filename'])

    ext = os.path.splitext(filename)[1].lower()
    if ext == ".fif":
        onset, termination = metadata['onset'], metadata['termination']
        if onset is not None:
            onset = onset - TEMPORAL_MARGIN
        if termination is not None:
            termination = termination + TEMPORAL_MARGIN

        seeg = SeegRecording.from_fif(filename, drop_channels=metadata['bad_channels'])
        seeg.trim(interval=(onset, termination))
    else:
        raise ValueError("Unexpected file format %s" % ext)

    if representation in ['avgref', 'bipolar']:
        plot_traces(seeg, imgname, representation=representation)
    elif representation == 'spectrogram':
        plot_spectrogram(seeg, imgname, 'avgref')
    else:
        raise ValueError("Unexpected representation '%s'" % representation)

if __name__ == '__main__':
    mne.set_log_level('WARNING')
    plot_seeg_recording(sys.argv[1], sys.argv[2], sys.argv[3])
