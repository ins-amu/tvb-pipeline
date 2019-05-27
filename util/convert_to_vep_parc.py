


import sys
import warnings

import numpy as np
import nibabel as nib
from sklearn.manifold import Isomap


class ColorLut():
    """Color look-up table"""
    def __init__(self, filename):
        self.inds  = np.genfromtxt(filename, usecols=(0,), dtype=int)
        self.names = np.genfromtxt(filename, usecols=(1,), dtype=str)
        self.colors = np.genfromtxt(filename, usecols=(2,3,4,5), dtype=int)

        self.name2ind = {name: ind for name, ind in zip(self.names, self.inds)}
        self.ind2name = {ind: name for name, ind in zip(self.names, self.inds)}


def load_rules(filename, section=None):
    with open(filename, 'r') as fl:
        lines = [line.strip() for line in fl.readlines()]

    # Remove comments
    rules = [line.split() for line in lines if len(line) > 0 and line[0] != "#"]

    # Get section index
    index = [('__begin__', 0)]
    rules_ = []
    for rule in rules:
        if rule[0] == 'Section':
            if len(rule) != 2:
                raise ValueError(f"Unexpected Section: '{rule}'")
            index.append((rule[1], len(rules_)))
        else:
            rules_.append(rule)
    index.append(('__end__', len(rules_)))
    rules = rules_

    if section is None:
        return rules
    else:
        ind = [a[0] for a in index].index(section)
        return rules[index[ind][1]:index[ind+1][1]]


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


def project_on_principal_axis(points):
    """
    Project all `points` (Nx3 array) on its principal axis.
    The axis is set to be oriented to positive in its second component (anterior direction in RAS coordinates).
    """

    # find principal eigendirection
    m_cov = np.dot(points.T, points)
    w, vr = np.linalg.eig(m_cov)
    eigendir = vr[:, np.argmax(w)]

    # orient them from posterior to anterior
    if eigendir[1] < 0:
        eigendir *= -1

    proj = np.dot(points, eigendir)
    return proj, eigendir


def find_interface_voxels(parc, regs):
    """Return indices of all voxels on the interface of multiple regions

    Straightforward and slow version.
    """

    # Add boundary layers
    assert -1 not in regs
    bparc = -1 * np.ones((parc.shape[0] + 2, parc.shape[1] + 2, parc.shape[2] + 2))
    bparc[1:-1, 1:-1, 1:-1] = parc

    # Possible voxels
    mask = np.zeros_like(parc, dtype=bool)
    for reg in regs:
        mask[parc == reg] = True

    # Get those on the interface
    interface = []
    for i, j, k in np.argwhere(mask):
        neigh_regs = set([bparc[i+1,j+1,k+1],
                          bparc[i  ,j+1,k+1],
                          bparc[i+2,j+1,k+1],
                          bparc[i+1,j  ,k+1],
                          bparc[i+1,j+2,k+1],
                          bparc[i+1,j+1,k  ],
                          bparc[i+1,j+1,k+2]])
        if all([reg in neigh_regs for reg in regs]):
            interface.append((i, j, k))

    return np.array(interface)


def find_interface_verts(triangs, labels, reg1, reg2, reg3):
    """
    Not the most efficient version, but it is enough.
    """

    interface = []
    for v1, v2, v3 in triangs:
        neighs = set([labels[v1], labels[v2], labels[v3]])
        if all([reg in neighs for reg in [reg1, reg2, reg3]]):
            interface.extend([v1, v2, v3])

    return np.array(interface)


def op_merge(labels, labels_in, label_out):
    """In-place region merge"""
    for lab in labels_in:
        labels[labels == lab] = label_out


def op_rename(labels, label_in, label_out):
    """In-place region rename"""
    op_merge(labels, [label_in], label_out)


def op_split(labels, mode, geom, label_in, labels_out, method, factors=None):
    """
    In-place split of a single region along a anterior-posterior axis.

    Method can be either:
      'pca'     for split after a linear projection
      'isomap'  for split after a nonlinear projection

    Regions in `labels_out` should be ordered in the anterior-posterior direction.
    If `factors` are missing, equal length split is performed.
    """

    if factors is None:
        factors = np.ones(len(labels_out), dtype=int)

    assert len(labels_out) == len(factors)

    inds = np.argwhere(labels == label_in)
    indsl = np.nonzero(labels == label_in)  # Just a different format

    if mode == 'voxel':
        # Affine transformation
        coords = (geom.dot(np.c_[inds, np.ones(inds.shape[0])].T).T)[:, :3]
    elif mode == 'triang':
        # Just coordinates
        coords = geom[inds[:, 0]]

    if method == 'pca':
        center = np.mean(coords, axis=0)
        xcoords, _ = project_on_principal_axis(coords - center)
    elif method == 'isomap':
        isomap = Isomap(n_components=1, n_neighbors=20)
        xcoords = isomap.fit_transform(coords)[:, 0]

        isomap_ori = 1
        if np.mean(coords[:, 1][xcoords < np.median(xcoords)]) > np.mean(coords[:, 1][xcoords > np.median(xcoords)]):
            isomap_ori = -1
        xcoords *= isomap_ori
    else:
        raise ValueError("Unknown method %s." % method)

    # Normalize
    xcoords -= np.min(xcoords)
    xcoords /= np.max(xcoords)

    limits = np.hstack([-np.inf, np.cumsum(factors)])
    limits = limits/limits[-1]
    limits[-1] = np.inf

    for label_out, xfr, xto in zip(reversed(labels_out), limits[:-1], limits[1:]):
        imask = (xcoords >= xfr) * (xcoords < xto)
        mask = [idxs[imask] for idxs in indsl]
        labels[mask] = label_out


def op_splitto(labels, mode, geom, label_in, labels_out, method):
    """
    In-place split of a single region between multiple other regions.
    Method can be either:
      'pca'     for split after a linear projection
      'isomap'  for split after a nonlinear projection

    `labels_out` should be ordered in the anterior-posterior direction.
    """
    # Reorder posterior-anterior
    labels_out = list(reversed(labels_out))

    inds = np.argwhere(labels == label_in)
    indsl = np.nonzero(labels == label_in)  # Just a different format

    if mode == 'voxel':
        # Affine transformation
        coords = (geom.dot(np.c_[inds, np.ones(inds.shape[0])].T).T)[:, :3]
    elif mode == 'triang':
        # Just coordinates
        verts, triangs = geom
        coords = verts[inds[:, 0]]

    if method == 'pca':
        center = np.mean(coords, axis=0)
        xcoords, direc = project_on_principal_axis(coords - center)
    elif method == 'isomap':
        isomap = Isomap(n_components=1, n_neighbors=20)
        xcoords = isomap.fit_transform(coords)[:, 0]

        isomap_ori = 1
        if np.mean(coords[:, 1][xcoords < np.median(xcoords)]) > np.mean(coords[:, 1][xcoords > np.median(xcoords)]):
            isomap_ori = -1
        xcoords *= isomap_ori
    else:
        raise ValueError("Unknown method %s." % method)

    # For each pair of to-regions
    limits = [-np.inf]
    for lab1, lab2 in zip(labels_out[:-1], labels_out[1:]):
        if mode == 'voxel':
            inds = find_interface_voxels(labels, (label_in, lab1, lab2))
        elif mode == 'triang':
            inds = find_interface_verts(triangs, labels, label_in, lab1, lab2)

        if len(inds) > 0:
            if mode == 'voxel':
                posras = np.mean((geom.dot(np.c_[inds, np.ones(inds.shape[0])].T).T)[:, :3], axis=0)
            elif mode == 'triang':
                posras = np.mean(verts[inds], axis=0)

            if method == 'pca':
                limits.append(np.dot(posras - center, direc))
            elif method == 'isomap':
                limits.append(isomap_ori * isomap.transform(posras.reshape(1, -1))[0, 0])

        else:
            warnings.warn("No interface found between %s,%s,%s." % (label_in, lab1, lab2))
            limits.append(limits[-1])
    limits.append(np.inf)

    for lab, x_fr, x_to in zip(labels_out, limits[:-1], limits[1:]):
        imask = (xcoords >= x_fr) * (xcoords < x_to)
        mask = [idxs[imask] for idxs in indsl]
        labels[mask] = lab


def op_splitmes(labels, hemi, verts, triangs, label_in, labels_out):
    """
    In-place split of a regions to mesial and 'lateral' parts.

    Criterion for this split is a orientation of a normal on an inflated surface.
    The operation is available only for triangulated surface.
    """
    NORMAL_THRESHOLD = -0.5

    inds = np.where(labels == label_in)
    triang_coords = verts[triangs]

    triang_normals = np.cross(triang_coords[:, 1] - triang_coords[:, 0],
                              triang_coords[:, 2] - triang_coords[:, 0])
    vert_normals = np.zeros_like(verts)
    vert_normals[triangs[:, 0]] += triang_normals
    vert_normals[triangs[:, 1]] += triang_normals
    vert_normals[triangs[:, 2]] += triang_normals
    vert_normals /= np.linalg.norm(vert_normals, axis=1)[:, None]

    mask = (labels == label_in)
    direc = 1 if hemi in ['rh', 'Right'] else -1
    labels[mask] = labels_out[0]
    labels[mask * (direc * vert_normals[:, 0] > NORMAL_THRESHOLD)] = labels_out[1]


def dehemize_name(name):
    if name[:7] == "ctx_%h_":
        return name[7:]
    elif name[:3] == "%H-":
        return name[3:]
    else:
        return name


def convert_parc(destrieux_annot_file, pial_file, inflated_file, hemisphere, parc_lut_file, rules_file, vep_annot_file):
    labels, _, names = nib.freesurfer.io.read_annot(destrieux_annot_file)
    names = [n.decode('ascii').replace("&", "_and_") for n in names]
    verts_pial, triangs_pial = nib.freesurfer.io.read_geometry(pial_file)
    verts_infl, triangs_infl = nib.freesurfer.io.read_geometry(inflated_file)

    rules = load_rules(rules_file, section='Cortex')
    colorlut = ColorLut(parc_lut_file)

    # Create master region list with old regions, new regions, and temporary regions
    names.extend(colorlut.names)                       # Names from VEP parcellation
    names.extend(["%%%d" % i for i in range(10)])      # Temporary regions

    # Shorthands
    def n2i(n):    return names.index(dehemize_name(n))
    def ns2i(ns):  return [names.index(dehemize_name(n)) for n in ns.split(",")]

    # Perform the in-place operations
    for rule in rules:
        if rule[0] == "merge":
            op_merge(labels, ns2i(rule[1]), n2i(rule[2]))
        elif rule[0] == "rename":
            op_rename(labels, n2i(rule[1]), n2i(rule[2]))
        elif rule[0] == "split":
            factors = [int(a) for a in rule[3].split(",")] if (len(rule) == 4) else None
            op_split(labels, 'triang', verts_pial, n2i(rule[1]), ns2i(rule[2]), method='pca', factors=factors)
        elif rule[0] == "splitnl":
            factors = [int(a) for a in rule[3].split(",")] if (len(rule) == 4) else None
            op_split(labels, 'triang', verts_pial, n2i(rule[1]), ns2i(rule[2]), method='isomap', factors=factors)
        elif rule[0] == "splitto":
            op_splitto(labels, 'triang', (verts_pial, triangs_pial), n2i(rule[1]), ns2i(rule[2]), method='pca')
        elif rule[0] == "splittonl":
            op_splitto(labels, 'triang', (verts_pial, triangs_pial), n2i(rule[1]), ns2i(rule[2]), method='isomap')
        elif rule[0] == "splitmes":
            op_splitmes(labels, hemisphere, verts_infl, triangs_infl, n2i(rule[1]), ns2i(rule[2]))
        else:
            raise ValueError("Unknown rule %s" % rule[0])

    # Keep only those in VEP parcellation color table
    newlabels = -1 * np.ones_like(labels)
    for i, name in zip(colorlut.inds, colorlut.names):
        newlabels[labels == names.index(name)] = i

    # Save
    nib.freesurfer.io.write_annot(vep_annot_file, newlabels, colorlut.colors, colorlut.names, fill_ctab=True)


def convert_seg(orig_label_file, lut_file, rules_file, vep_label_file):
    mgz_orig = nib.load(orig_label_file)
    labelvol = mgz_orig.get_data().copy()
    affine = mgz_orig.affine
    colorlut = ColorLut(lut_file)

    rules = load_rules(rules_file, section='Subcortical')
    rules = expand_wildcards_hemisphere(rules)

    n2i_ = {name: ind for name, ind in zip(colorlut.names, colorlut.inds)}
    n2i_.update({("%%%d" % i): -10-i for i in range(10)})  # Temporary regions

    # Shorthands
    def n2i(n):   return n2i_[n]
    def ns2i(ns): return [n2i_[n] for n in ns.split(",")]

    for rule in rules:
        if rule[0] == "merge":
            op_merge(labelvol, ns2i(rule[1]), n2i(rule[2]))
        elif rule[0] == "rename":
            op_rename(labelvol, n2i(rule[1]), n2i(rule[2]))
        elif rule[0] == "split":
            factors = [int(a) for a in rule[3].split(",")] if (len(rule) == 4) else None
            op_split(labelvol, 'voxel', affine, n2i(rule[1]), ns2i(rule[2]), method='pca', factors=factors)
        elif rule[0] == "splitnl":
            factors = [int(a) for a in rule[3].split(",")] if (len(rule) == 4) else None
            op_split(labelvol, 'voxel', affine, n2i(rule[1]), ns2i(rule[2]), method='isomap', factors=factors)
        elif rule[0] == "splitto":
            op_splitto(labelvol, 'voxel', affine, n2i(rule[1]), ns2i(rule[2]), method='pca')
        elif rule[0] == "splittonl":
            op_splitto(labelvol, 'voxel', affine, n2i(rule[1]), ns2i(rule[2]), method='isomap')
        elif rule[0] == "splitmes":
            raise ValueError("Rule 'splitmes' is available only in surface mode.")
        else:
            raise ValueError("Unknown rule %s" % rule[0])

    assert np.sum(labelvol < 0) == 0

    mgz_vep = nib.freesurfer.mghformat.MGHImage(labelvol, affine, mgz_orig.header)
    nib.save(mgz_vep, vep_label_file)


if __name__ == '__main__':
    cmd = sys.argv[1]
    eval(cmd)(*sys.argv[2:])
