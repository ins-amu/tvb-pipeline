
import os
import re
import mne

import numpy as np


class SeegRecording():

    def __init__(self, contacts, data, sampling_rate):
        # Sort the contacts first
        contacts, indices = zip(*sorted((c, i) for i, c in enumerate(contacts)))

        self.contacts = contacts
        self.ncontacts = len(self.contacts)
        self.contact_names = [str(a) + str(b) for a, b in self.contacts]

        assert data.ndim == 2
        assert data.shape[0] == self.ncontacts
        self.data = data[indices, :]

        self.sampling_rate = sampling_rate
        nsamples = self.data.shape[1]
        self.t = np.linspace(0, (nsamples - 1)*(1./self.sampling_rate), nsamples)

        self.electrodes = {}
        for i, (name, number) in enumerate(self.contacts):
            if name not in self.electrodes:
                self.electrodes[name] = []
            self.electrodes[name].append(i)

        self.set_bipolar()


    @classmethod
    def from_ades(cls, filename):
        data_file = os.path.splitext(filename)[0] + ".dat"

        contact_names = []
        contacts = []
        sampling_rate = None
        nsamples = None
        seeg_idxs = []


        bad_channels = []
        bad_file = filename + ".bad"
        if os.path.isfile(bad_file):
            with open(bad_file, 'r') as fd:
                bad_channels = [ch.strip() for ch in fd.readlines() if ch.strip() != ""]

        with open(filename, 'r') as fd:
            fd.readline() # ADES header file

            kw, sampling_rate = [s.strip() for s in fd.readline().strip().split('=')]
            assert kw == 'samplingRate'
            sampling_rate = float(sampling_rate)

            kw, nsamples = [s.strip() for s in fd.readline().strip().split('=')]
            assert kw == 'numberOfSamples'
            nsamples = int(nsamples)

            channel_idx = 0
            for line in fd.readlines():
                if not line.strip():
                    continue

                parts = [p.strip() for p in line.strip().split('=')]
                if len(parts) > 1 and parts[1] == 'SEEG':
                    name, idx = re.match("([A-Za-z]+[']*)([0-9]+)", parts[0]).groups()
                    idx = int(idx)
                    if parts[0] not in bad_channels:
                        contacts.append((name, idx))
                        seeg_idxs.append(channel_idx)


                channel_idx += 1

        data = np.fromfile(data_file, dtype='f4')
        ncontacts = data.size // nsamples

        if data.size != nsamples*ncontacts:
            print("!! data.size != nsamples*ncontacts")
            print("!! %d != %d %d" % (data.size, nsamples, ncontacts))
            print("!! Ignoring nsamples")
            nsamples = int(data.size/ncontacts)

        data = data.reshape((nsamples, ncontacts)).T
        data = data[seeg_idxs, :]

        return cls(contacts, data, sampling_rate)


    @classmethod
    def from_fif(cls, filename, drop_channels=None, rename_channels=None):
        raw = mne.io.read_raw_fif(filename)

        if rename_channels is not None:
            raw.rename_channels(rename_channels)

        return cls._from_mne_raw(raw, drop_channels)

    @classmethod
    def from_edf(cls, filename, drop_channels=None, rename_channels=None):
        raw = mne.io.read_raw_edf(filename, preload=True)

        if rename_channels is not None:
            raw.rename_channels(rename_channels)

        return cls._from_mne_raw(raw, drop_channels)

    @classmethod
    def _from_mne_raw(cls, raw, drop_channels=None):
        contacts = []
        seeg_idxs = []

        if drop_channels is None:
            drop_channels = []

        for i, ch_name in enumerate(raw.ch_names):
            if ch_name in raw.info['bads'] or ch_name in drop_channels:
                continue

            match = re.match("^([A-Za-z]+[']*)([0-9]+)$", ch_name)
            if match is not None:
                name, idx = match.groups()
                contacts.append((name, int(idx)))
                seeg_idxs.append(i)

        return cls(contacts, raw.get_data()[seeg_idxs, :], raw.info['sfreq'])


    def set_bipolar(self):
        self.bipolar = []
        self.electrodes_bipolar = {}

        nchannels = 0
        for i1, (name, contact_ind) in enumerate(self.contacts):
            try:
                i2 = self.contacts.index((name, contact_ind + 1))
                bip_name = name + str(contact_ind + 1) + "-" + str(contact_ind)
                self.bipolar.append((bip_name, i1, i2))

                if name not in self.electrodes_bipolar:
                    self.electrodes_bipolar[name] = []
                self.electrodes_bipolar[name].append(nchannels)
                nchannels += 1

            except ValueError:
                pass

    def interval_to_index(self, interval):
        tid1, tid2 = 0, -1
        if interval[0] is not None:
            if interval[0] < self.t[0]:
                tid1 = 0
            else:
                tid1 = np.argmax(self.t >= interval[0])
        if interval[1] is not None:
            if interval[1] > self.t[-1]:
                tid2 = -1
            else:
                tid2 = np.argmax(self.t >= interval[1])

        return (tid1, tid2)

    def time(self, interval=(None, None)):
        tid1, tid2 = self.interval_to_index(interval)
        return self.t[tid1:tid2]

    def ts(self, name, interval=(None, None)):
        idx = self.contact_names.index(name)
        tid1, tid2 = self.interval_to_index(interval)
        return self.data[idx, tid1:tid2]

    def get_data(self):
        return np.copy(self.data)

    def get_data_bipolar(self):
        data_bipolar = np.zeros((len(self.bipolar), len(self.t)))
        for i, (_, i1, i2) in enumerate(self.bipolar):
            data_bipolar[i, :] = self.data[i1, :] - self.data[i2, :]
        return data_bipolar

    def get_channel_names(self):
        return self.contact_names

    def get_channel_names_bipolar(self):
        return [name for name, i1, i2 in self.bipolar]

    def remove_avg(self):
        self.data = self.data - np.mean(self.data, axis=0)

    def trim(self, interval=(None, None)):
        tid1, tid2 = self.interval_to_index(interval)

        self.data = self.data[:, tid1:tid2]
        self.t = self.t[tid1:tid2]
