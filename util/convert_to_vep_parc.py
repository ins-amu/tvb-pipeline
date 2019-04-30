


import sys

import numpy as np
import nibabel as nib
from sklearn.manifold import Isomap


class ColorLut():
    """Color look-up table"""
    def __init__(self, filename):
        self.inds  = np.genfromtxt(filename, usecols=(0,), dtype=int)
        self.names = np.genfromtxt(filename, usecols=(1,), dtype=str)

        # Temporary regions
        self.inds = np.insert(self.inds, 0, [-10-i for i in range(10)])
        self.names = np.insert(self.names, 0, ['__temporary_region_%d' % i for i in range(10)])

        self.name2ind = {name: ind for name, ind in zip(self.names, self.inds)}
        self.ind2name = {ind: name for name, ind in zip(self.names, self.inds)}



def load_rules(filename):
    with open(filename, 'r') as fl:
        lines = [line.strip() for line in fl.readlines()]

    rules = [line.split() for line in lines if len(line) > 0 and line[0] != "#"]
    return rules


def expand_wildcards_hemisphere(rules):
    rules_l = []
    rules_r = []
    rules_n = []

    for rule in rules:
        if any([("%h" in elem) or ("%H" in elem) for elem in rule]):
            rules_l.append([elem.replace("%h", "lh").replace("%H", "Left")  for elem in rule])
            rules_r.append([elem.replace("%h", "rh").replace("%H", "Right") for elem in rule])
        else:
            rules_n.append(rule)
    return rules_n + rules_l + rules_r


def expand_wildcards_temp(rules):
    new_rules = []
    for rule in rules:
        for i in range(10):
            rule = [elem.replace("%%%d" % i, "__temporary_region_%d" % i) for elem in rule]
        new_rules.append(rule)
    return new_rules


def parc_merge(parc, lut, regions_in, region_out):
    """
    In-place merging of several regions.

    Formerly `update_parc_regroup`
    """
    ind_out = lut.name2ind[region_out]
    for reg in regions_in:
        ind_in = lut.name2ind[reg]
        mask = (parc == ind_in)
        parc[mask] = ind_out


def parc_rename(parc, lut, region_in, region_out):
    """
    In-place rename of a single region.

    Formerly `update_parc_rename`
    """
    parc_merge(parc, lut, [region_in], region_out)


def parc_split(parc, affine, lut, region_in, regions_out, factors=None):
    """
    In-place split of a single region along a anterior-posterior axis.

    Regions in `regions_out` should be ordered in the anterior-posterior direction.
    If `factors` are missing, equal length split is performed.

    Formerly `update_divid_multipul`
    """

    if factors is None:
        factors = np.ones(len(regions_out), dtype=int)

    assert len(regions_out) == len(factors)

    inds = np.argwhere(parc == lut.name2ind[region_in])
    indsl = np.nonzero(parc == lut.name2ind[region_in])  # Just a different format

    coords = (affine.dot(np.c_[inds, np.ones(inds.shape[0])].T).T)[:, :3]
    xcoords = project_on_principal_axis(coords)

    limits = np.hstack([0, np.cumsum(factors)])
    limits = limits/limits[-1]
    # Better be careful about the floating point comparison
    limits[0]  = -1
    limits[-1] =  2


    for region_out, xfr, xto in zip(reversed(regions_out), limits[:-1], limits[1:]):
        ind_out = lut.name2ind[region_out]
        imask = (xcoords >= xfr) * (xcoords < xto)
        mask = [idxs[imask] for idxs in indsl]
        parc[mask] = ind_out


def parc_splitnl(parc, affine, lut, region_in, regions_out, factors=None):
    """
    In-place split split of a single region using a nonlinear embedding.

    Regions in `regions_out` should be ordered in the anterior-posterior direction.
    If `factors` are missing, equal length split (on the manifold) is performed.
    """

    if factors is None:
        factors = np.ones(len(regions_out), dtype=int)

    assert len(regions_out) == len(factors)

    inds = np.argwhere(parc == lut.name2ind[region_in])
    indsl = np.nonzero(parc == lut.name2ind[region_in])  # Just a different format

    coords = (affine.dot(np.c_[inds, np.ones(inds.shape[0])].T).T)[:, :3]

    isomap = Isomap(n_components=1, n_neighbors=10)
    xcoords = isomap.fit_transform(coords)[:, 0]

    # Normalize
    xcoords -= np.min(xcoords)
    xcoords /= np.max(xcoords)

    # Orient them from posterior to anterior
    if np.mean(coords[:, 1][xcoords < 0.1]) > np.mean(coords[:, 1][xcoords > 0.9]):
        xcoords = 1 - xcoords


    limits = np.hstack([0, np.cumsum(factors)])
    limits = limits/limits[-1]
    # Better be careful about the floating point comparison
    limits[0]  = -1
    limits[-1] =  2

    for region_out, xfr, xto in zip(reversed(regions_out), limits[:-1], limits[1:]):
        ind_out = lut.name2ind[region_out]
        imask = (xcoords >= xfr) * (xcoords < xto)
        mask = [idxs[imask] for idxs in indsl]
        parc[mask] = ind_out



def project_on_principal_axis(points):
    """
    Project all `points` (Nx3 array) on its principal axis and normalize to [0; 1].
    The axis is set to be oriented to positive in its second component (anterior direction in RAS coordinates).
    """

    center = np.mean(points, axis=0)
    centered_points = points - center

    # find principal eigendirection
    m_cov = np.dot(centered_points.T, centered_points)
    w, vr = np.linalg.eig(m_cov)
    eigendir = vr[:, np.argmax(w)]

    # orient them from posterior to anterior
    if eigendir[1] < 0:
        eigendir *= -1

    proj = np.dot(points, eigendir)
    proj -= proj.min()
    proj /= proj.max()

    return proj



def convert_to_vep_parc(destrieux_file, lut_file, rules_file, vep_file):
    mgz_destrieux = nib.load(destrieux_file)
    labelvol = mgz_destrieux.get_data().copy()
    affine = mgz_destrieux.affine
    colorlut = ColorLut(lut_file)

    rules = load_rules(rules_file)
    rules = expand_wildcards_hemisphere(rules)
    rules = expand_wildcards_temp(rules)

    for rule in rules:
        if rule[0] == "merge":
            parc_merge(labelvol, colorlut, rule[1].split(","), rule[2])
        elif rule[0] == "rename":
            parc_rename(labelvol, colorlut, rule[1], rule[2])
        elif rule[0] == "split":
            factors = [int(a) for a in rule[3].split(",")] if (len(rule) == 4) else None
            parc_split(labelvol, affine, colorlut, rule[1], rule[2].split(","), factors=factors)
        elif rule[0] == "splitnl":
            factors = [int(a) for a in rule[3].split(",")] if (len(rule) == 4) else None
            parc_splitnl(labelvol, affine, colorlut, rule[1], rule[2].split(","), factors=factors)
        else:
            raise ValueError("Unknown rule %s" % rule[0])

    assert np.sum(labelvol < 0) == 0


    mgz_vep = nib.freesurfer.mghformat.MGHImage(labelvol, affine, mgz_destrieux.header)
    nib.save(mgz_vep, vep_file)



if __name__ == "__main__":
    convert_to_vep_parc(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
