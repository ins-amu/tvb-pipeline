#!/usr/bin/env python

import os
import logging
import matplotlib.pyplot as plt
import numpy as np
import re
import zipfile

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors


class NamedPoints():
    def __init__(self, fl):
        data = np.genfromtxt(fl, dtype=None)
        self.xyz = np.array([[l[1], l[2], l[3]] for l in data])
        self.names = [l[0].decode('ascii') for l in data]


class Contacts(NamedPoints):
    def __init__(self, filename):
        super().__init__(filename)
        self.electrodes = {}
        for i, name in enumerate(self.names):
            elec_name, _ = re.match("([A-Za-z]+[']*)([0-9]+)", name).groups()
            if elec_name in self.electrodes:
                self.electrodes[elec_name].append(i)
            else:
                self.electrodes[elec_name] = [i]


def seeg_elecs(tvbzip_file, seegxyz, out_fig):
    with zipfile.ZipFile(tvbzip_file) as tvbzip:
        with tvbzip.open("centres.txt") as centres_file:
            regions = NamedPoints(centres_file)
        with tvbzip.open("cortical.txt") as cortical_file:
            cortical = np.genfromtxt(cortical_file, dtype=int)
    contacts = Contacts(seegxyz)

    fig = plt.figure(figsize=(10, 10))
    labels = ['L --- R', 'P --- A', 'I --- S']

    reg_color = ['royalblue' if c == 1 else 'darkblue' for c in cortical]

    for pos, id1, id2 in [(221, 0, 1), (224, 1, 2), (223, 0, 2)]:
        ax = fig.add_subplot(pos)
        ax.scatter(regions.xyz[:, id1], regions.xyz[:, id2], color=reg_color, s=40)
        for name, idxs in contacts.electrodes.items():
            ax.scatter(contacts.xyz[idxs, id1], contacts.xyz[idxs, id2], s=10, label=name)
        ax.set_xlabel(labels[id1])
        ax.set_ylabel(labels[id2])
        ax.set_aspect('equal')

    ax = fig.add_subplot(222, projection='3d')
    ax.scatter(regions.xyz[:, 0], regions.xyz[:, 1], regions.xyz[:, 2], color=reg_color, s=40)
    for name, idxs in contacts.electrodes.items():
        ax.scatter(contacts.xyz[idxs, 0], contacts.xyz[idxs, 1], contacts.xyz[idxs, 2], s=10, label=name)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])
    ax.set_aspect('equal')

    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(out_fig)





if __name__ == '__main__':
    import sys
    loglevel = logging.INFO
    if os.environ.get('VERBOSE', False):
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel)
    cmd = sys.argv[1]
    eval(cmd)(*sys.argv[2:])
