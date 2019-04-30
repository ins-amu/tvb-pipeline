#!/usr/bin/env python3

import sys


import numpy as np
import pandas as pd

from convert_to_vep_parc import load_rules, expand_wildcards_hemisphere


def create_vep_fs_lut(fs_lut, vep_rules, vep_fs_lut):
    fsregs = list(np.genfromtxt(fs_lut, usecols=(1,), dtype=str))
    rules = load_rules("VepAtlasRules.txt")

    newregs = []
    for rule in rules:
        if rule[0] in ["merge", "rename"]:
            newregs.append(rule[2])
        elif rule[0] in ["split", "splitnl"]:
            newregs.extend(rule[2].split(","))

    # Filter temp regions
    newregs = [reg for reg in newregs if reg not in ["%%%d" % i for i in range(10)]]

    assert all(["%H" in reg for reg in newregs])

    np.random.seed(42)
    colors = np.random.choice(256, (1000, 3))

    # Filter existing regions
    newregs = [reg for reg in newregs if reg not in fsregs]


    with open(fs_lut) as fl:
        lines = fl.readlines()
    with open(vep_fs_lut, "w") as fl:
        fl.writelines(lines)
        fl.write("\n\n#\n# Labels for the VEP parcellation\n#\n\n")

        for hemi, hnum in [("Left", 70000), ("Right", 71000)]:
            for i, reg in enumerate(newregs):
                fl.write("%5d  %-60s %3d %3d %3d  0\n" % (hnum + i + 1, reg.replace("%H", hemi),
                                                          colors[i, 0], colors[i, 1], colors[i, 2]))


def create_vep_mrtrix_lut(vep_fs_lut, vep_regions_file, vep_mrtrix_lut):
    names = list(np.genfromtxt(vep_fs_lut, usecols=(1,), dtype=str))
    colors = np.genfromtxt(vep_fs_lut, usecols=(2,3,4,5), dtype=int)
    regions = np.genfromtxt(vep_regions_file, usecols=(0,), dtype=str)

    with open(vep_mrtrix_lut, 'w') as fl:
        fl.write("   0   %-60s  0   0   0   0\n" % ("Unknown") )
        i = 1
        for hemi in ["Left", "Right"]:
            for reg in regions:
                regname = hemi + "-" + reg
                ind = names.index(regname)
                fl.write("%4d   %-60s  %4d %4d %4d %4d\n" % (i, regname, *colors[ind]))
                i += 1




if __name__ == "__main__":
    create_vep_fs_lut("FreeSurferColorLUT.txt", "VepAtlasRules.txt", "VepFreeSurferColorLut.txt")
    create_vep_mrtrix_lut("VepFreeSurferColorLut.txt", "VepRegions.txt", "lut.vep.txt")
