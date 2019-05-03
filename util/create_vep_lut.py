#!/usr/bin/env python3

import sys


import numpy as np
import pandas as pd

from convert_to_vep_parc import load_rules, expand_wildcards_hemisphere
from convert_to_vep_parc import SHIFT_LH, SHIFT_RH


def create_luts(fs_lut_file, vep_rules_file, vep_regions_file,
                vep_fs_lut_file, vep_mrtrix_lut_file, vep_subcort_file, vep_aparc_lut_file):

    # Load original table
    fs_regs = list(np.genfromtxt(fs_lut_file, usecols=(1,), dtype=str))

    # Load rules
    rules = load_rules(vep_rules_file)
    newregs = []
    for rule in rules:
        if rule[0] in ["merge", "rename"]:
            newregs.append(rule[2])
        elif rule[0] in ["split", "splitnl"]:
            newregs.extend(rule[2].split(","))

    # Filter temp regions
    newregs = [reg for reg in newregs if reg not in ["%%%d" % i for i in range(10)]]
    assert all(["%H" in reg for reg in newregs])

    vep_regs = np.genfromtxt(vep_regions_file, usecols=(1,), dtype=str)
    vep_iscort = np.genfromtxt(vep_regions_file, usecols=(0,), dtype=int).astype(bool)

    # Make sure that every cortical region is in rules
    for reg in vep_regs[vep_iscort]:
        if ("%%H-%s" % reg) not in newregs:
            raise Exception("Rule for region '%s' is missing" % reg)

    # Make sure that every subcortical region is either in rules or in Freesurfer table
    for reg in vep_regs[~vep_iscort]:
        assert ("%H-"+reg in newregs) or (("Left-"+reg in fs_regs) and ("Right-"+reg in fs_regs))


    # Make sure all subcortical regions are at the end
    assert np.all(vep_iscort[:-1] >= vep_iscort[1:])

    create_vep_fs_lut(vep_regs, fs_regs, fs_lut_file, vep_fs_lut_file)
    create_vep_mrtrix_lut(vep_regs, vep_fs_lut_file, vep_mrtrix_lut_file)
    create_subcort_list(vep_regs[~vep_iscort], vep_fs_lut_file, vep_subcort_file)
    create_parc_lut(vep_regs[vep_iscort], vep_fs_lut_file, vep_aparc_lut_file)


def create_parc_lut(vep_regs, vep_fs_lut_file, vep_aparc_lut_file):
    names = list(np.genfromtxt(vep_fs_lut_file, usecols=(1,), dtype=str))
    colors = np.genfromtxt(vep_fs_lut_file, usecols=(2,3,4,5), dtype=int)

    with open(vep_aparc_lut_file, 'w') as fl:
        fl.write("  0 %-60s   0   0   0   0\n" % "Unknown")
        for i, reg in enumerate(vep_regs):
            ind = names.index("Left-" + reg)
            fl.write("%3d %-60s %3d %3d %3d %3d\n" % (i+1, reg, *colors[ind]))


def create_vep_fs_lut(vep_regs, fs_regs, fs_lut_file, vep_fs_lut_file):
    regsl = ["Left-"+reg  for reg in vep_regs]
    regsr = ["Right-"+reg for reg in vep_regs]

    # Filter existing regions
    regsl = [reg for reg in regsl if reg not in fs_regs]
    regsr = [reg for reg in regsr if reg not in fs_regs]

    # Get random, distinct colors
    np.random.seed(42)
    colors = [(a % 2**8, (a % 2**16) // 2**8, (a % 2**24) // 2**16)
              for a in np.random.choice(2**24, len(vep_regs), replace=False)]


    with open(fs_lut_file) as fl:
        lines = fl.readlines()
    with open(vep_fs_lut_file, "w") as fl:
        fl.writelines(lines)
        fl.write("\n\n#\n# Labels for the VEP parcellation\n#\n\n")

        for regs, hnum in [(regsl, SHIFT_LH), (regsr, SHIFT_RH)]:
            for i, reg in enumerate(regs):
                fl.write("%5d  %-60s %3d %3d %3d  0\n" % (hnum + i + 1, reg, *colors[i]))


def create_vep_mrtrix_lut(vep_regs, vep_fs_lut_file, vep_mrtrix_lut_file):
    names = list(np.genfromtxt(vep_fs_lut_file, usecols=(1,), dtype=str))
    colors = np.genfromtxt(vep_fs_lut_file, usecols=(2,3,4,5), dtype=int)

    with open(vep_mrtrix_lut_file, 'w') as fl:
        fl.write("   0   %-60s  0   0   0   0\n" % ("Unknown") )
        i = 1
        for hemi in ['Left', 'Right']:
            for reg in vep_regs:
                regname = hemi + "-" + reg
                ind = names.index(regname)
                fl.write("%4d   %-60s  %4d %4d %4d %4d\n" % (i, regname, *colors[ind]))
                i += 1


def create_subcort_list(vep_subcort_regions, vep_fs_lut, vep_subcort_list):
    fs_names = list(np.genfromtxt(vep_fs_lut, usecols=(1,), dtype=str))
    fs_inds = np.genfromtxt(vep_fs_lut, usecols=(0,), dtype=int)

    with open(vep_subcort_list, 'w') as fl:
        for hemi in ['Left', 'Right']:
            for reg in vep_subcort_regions:
                name = hemi + "-" + reg
                fl.write("%d\n" % fs_inds[fs_names.index(name)])




if __name__ == "__main__":
    create_luts("data/FreeSurferColorLUT.txt", "data/VepAtlasRules.txt", "data/VepRegions.txt",
                "data/VepFreeSurferColorLut.txt", "data/VepMrtrixLut.txt", "data/subcort.vep.txt",
                "data/VepAparcColorLut.txt")
