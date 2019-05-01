#!/usr/bin/env python3

import sys


import numpy as np
import pandas as pd

from convert_to_vep_parc import load_rules, expand_wildcards_hemisphere

SHIFT_LH = 70000
SHIFT_RH = 71000

def create_vep_fs_lut(fs_lut, vep_rules, vep_fs_lut):
    fsregs = list(np.genfromtxt(fs_lut, usecols=(1,), dtype=str))
    rules = load_rules(vep_rules)

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

    regsl = [reg.replace("%H", "Left") for reg in newregs]
    regsr = [reg.replace("%H", "Right") for reg in newregs]

    # Filter existing regions
    regsl = [reg for reg in regsl if reg not in fsregs]
    regsr = [reg for reg in regsr if reg not in fsregs]

    with open(fs_lut) as fl:
        lines = fl.readlines()
    with open(vep_fs_lut, "w") as fl:
        fl.writelines(lines)
        fl.write("\n\n#\n# Labels for the VEP parcellation\n#\n\n")

        for regs, hnum in [(regsl, SHIFT_LH), (regsr, SHIFT_RH)]:
            for i, reg in enumerate(regs):
                fl.write("%5d  %-60s %3d %3d %3d  0\n" % (hnum + i + 1, reg,
                                                          colors[i, 0], colors[i, 1], colors[i, 2]))


def create_vep_mrtrix_lut(vep_fs_lut, vep_regions_file, vep_mrtrix_lut):
    names = list(np.genfromtxt(vep_fs_lut, usecols=(1,), dtype=str))
    colors = np.genfromtxt(vep_fs_lut, usecols=(2,3,4,5), dtype=int)
    regions = np.genfromtxt(vep_regions_file, usecols=(0,), dtype=str)

    with open(vep_mrtrix_lut, 'w') as fl:
        fl.write("   0   %-60s  0   0   0   0\n" % ("Unknown") )
        i = 1
        for reg in regions:
            ind = names.index(reg)
            fl.write("%4d   %-60s  %4d %4d %4d %4d\n" % (i, reg, *colors[ind]))
            i += 1

def create_subcort_list(vep_fs_lut, vep_regions_file, vep_subcort_list):
    fs_names = list(np.genfromtxt(vep_fs_lut, usecols=(1,), dtype=str))
    fs_inds = np.genfromtxt(vep_fs_lut, usecols=(0,), dtype=int)

    vep_names = list(np.genfromtxt(vep_regions_file, usecols=(0,), dtype=str))
    vep_iscort = np.genfromtxt(vep_regions_file, usecols=(1,), dtype=int).astype(bool)

    with open(vep_subcort_list, 'w') as fl:
        for name, iscort in zip(vep_names, vep_iscort):
            if not iscort:
                fl.write("%d\n" % fs_inds[fs_names.index(name)])



def create_parc_lut(vep_fs_lut, vep_regions_file, vep_parc_lut):
    fs_names = list(np.genfromtxt(vep_fs_lut, usecols=(1,), dtype=str))
    fs_inds = np.genfromtxt(vep_fs_lut, usecols=(0,), dtype=int)

    vep_names = np.genfromtxt(vep_regions_file, usecols=(0,), dtype=str)
    vep_iscort = np.genfromtxt(vep_regions_file, usecols=(1,), dtype=int).astype(bool)

    # Without hemisphere
    names = [reg[5:] for reg in vep_names[vep_iscort] if reg[:5] == "Left-"]

    # Just check
    # assert all([fs_inds[fs_names.index("Left-%s" % names[i])] == SHIFT_LH + i + 1 for i in range(len(names))])
    # for i, name in enumerate(names):
    #    if fs_inds[fs_names.index("Left-%s" % name)] != SHIFT_LH + i + 1:
    #        print(name, i)





if __name__ == "__main__":
    # create_vep_fs_lut("data/FreeSurferColorLUT.txt", "data/VepAtlasRules.txt", "data/VepFreeSurferColorLut.txt")
    # create_vep_mrtrix_lut("data/VepFreeSurferColorLut.txt", "data/VepRegions.txt", "data/lut.vep.txt")
    # create_subcort_list("data/VepFreeSurferColorLut.txt", "data/VepRegions.txt", "data/subcort.vep.txt")

    create_parc_lut("data/VepFreeSurferColorLut.txt", "data/VepRegions.txt", "data/VepAparcColorLut.txt")
