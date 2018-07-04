#!/usr/bin/env python

import os
import logging
import numpy as np
import re
import textwrap
import zipfile


import nibabel as nib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from .nifti import gen_volume_regions
from .elecs import NamedPoints, Contacts


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

ras_slice_indices = {
    'transversal': ([1, 0], 2),
    'coronal':     ([2, 0], 1),
    'sagittal':    ([2, 1], 0)
}



def get_slice(img, slice_type, ras_coord):
    eps = 1e-6

    affine = img.affine

    # Check the the transformation
    assert np.all(np.sum(np.abs(affine[:3, :3]) > eps, axis=1) == 1)

    slice_dims, slice_param = ras_slice_indices[slice_type]

    ras_coords = [ras_coord for _ in range(3)] + [1]
    ras_coords[slice_dims[0]] = ras_coords[slice_dims[1]] = 0
    vox_coords = [int(x) for x in np.linalg.solve(affine, ras_coords)[0:3]]

    ax0 = np.where(np.abs(affine[slice_dims[0], 0:3]) > eps)[0][0]
    ax1 = np.where(np.abs(affine[slice_dims[1], 0:3]) > eps)[0][0]

    slicer = vox_coords
    slicer[ax0] = slice(None, None, 1 if affine[slice_dims[0], ax0] > 0 else -1)
    slicer[ax1] = slice(None, None, 1 if affine[slice_dims[1], ax1] > 0 else -1)
    slicer = tuple(slicer)

    data = img.get_data()[slicer]

    # Extent
    corner_0 = np.dot(affine, [0, 0, 0, 1])[0:3]
    corner_1 = np.dot(affine, np.append(img.shape, 1))[0:3]
    extent_ax0 = sorted([corner_0[slice_dims[0]], corner_1[slice_dims[0]]])
    extent_ax1 = sorted([corner_0[slice_dims[1]], corner_1[slice_dims[1]]])

    # Transpose if necessary
    if ax0 > ax1:
        data = data.T

    extent = extent_ax1 + extent_ax0

    return data, extent



def plot_t1_elec(t1_img_file, elec_img_file, elec_pos_file, contact_name, filename):
    t1_img = nib.load(t1_img_file)
    elec_img = nib.load(elec_img_file)
    contacts = Contacts(elec_pos_file)

    labels = ["L --- R", "P --- A", "I --- S"]

    t1_lims = np.min(t1_img.get_data()), np.max(t1_img.get_data())
    elec_lims = np.min(elec_img.get_data()), np.max(elec_img.get_data())

    contact_ind = contacts.names.index(contact_name)
    elec_inds = contacts.electrodes[contacts.get_elec(contact_name)]

    fig = plt.figure(figsize=(18, 18))
    for i, slice_type in enumerate(['coronal', 'sagittal', 'transversal']):
        plt.subplot(2, 2, i + 1)

        slice_dims, slice_param = ras_slice_indices[slice_type]

        contact_xy = contacts.xyz[contact_ind][slice_dims]
        contact_z = contacts.xyz[contact_ind][slice_param]

        # Alpha colormap
        cmap = plt.cm.hot
        alpha_cmap = cmap(np.arange(cmap.N))
        alpha_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        alpha_cmap = matplotlib.colors.ListedColormap(alpha_cmap)

        elec_img_data, elec_img_extent = get_slice(elec_img, slice_type, contact_z)
        plt.imshow(elec_img_data, extent=elec_img_extent,
                   vmin=elec_lims[0], vmax=elec_lims[1],
                   cmap='hot',
                   origin='lower', zorder=0)

        t1_img_data, t1_img_extent = get_slice(t1_img, slice_type, contact_z)
        plt.imshow(t1_img_data, extent=t1_img_extent,
                   vmin=t1_lims[0], vmax=t1_lims[1],
                   cmap='Greys_r', origin='lower', zorder=1, alpha=0.5)


        plt.scatter([contact_xy[1]], [contact_xy[0]], facecolors='none', edgecolors='g',
                    s=80, zorder=2)

        for ind in elec_inds:
            if ind == contact_ind:
                continue
            plt.scatter([contacts.xyz[ind, slice_dims[1]]], [contacts.xyz[ind, slice_dims[0]]],
                        c='g', s=10, marker='x', zorder=2)

        plt.xlabel(labels[slice_dims[1]], fontsize=16)
        plt.ylabel(labels[slice_dims[0]], fontsize=16)

    # Text description
    ax = plt.subplot(2, 2, 4)
    ax.axis('off')
    ax.text(0.1, 0.9, contact_name, size=16, ha="left")

    desc = """
    In black-white: T1 scan
    In black-red-yellow: Scan with implanted electrodes
    In green: electrode position
    """
    ax.text(0.1, 0.7, textwrap.dedent(desc), size=12, ha="left")

    plt.tight_layout(pad=1.4)
    plt.savefig(filename)
    plt.close()


def plot_t1_plus_elecs(t1_img_file, elec_img_file, elec_pos_file, out_direc):
    if not os.path.isdir(out_direc):
        os.makedirs(out_direc)

    contacts = Contacts(elec_pos_file)
    for elec, elec_inds in contacts.electrodes.items():
        contact_name = contacts.names[elec_inds[0]]
        contact_name_safe = contact_name.replace("'", "p")
        filename = os.path.join(out_direc, "t1_elec-img_elec-pos_%s.png" % contact_name_safe)
        plot_t1_elec(t1_img_file, elec_img_file, elec_pos_file, contact_name, filename)


def plot_connectivity(tvbzip_file, out_img_file):
    from matplotlib.colors import LogNorm

    with zipfile.ZipFile(tvbzip_file) as tvbzip:
        with tvbzip.open("weights.txt") as weights_file:
            weights = np.genfromtxt(weights_file)
        with tvbzip.open("tract_lengths.txt") as lengths_file:
            lengths = np.genfromtxt(lengths_file)
        with tvbzip.open("centres.txt") as centres_file:
            names = np.genfromtxt(centres_file, usecols=(0,), dtype=str)

    nreg = len(names)

    plt.figure(figsize=(24, 12))

    plt.subplot(1, 2, 1)
    plt.title("Weights")
    plt.imshow(weights, cmap=plt.cm.gray_r, norm=LogNorm(1e-5*np.max(weights), np.max(weights)))
    plt.xticks(np.r_[:nreg], names, rotation='vertical', fontsize=8)
    plt.yticks(np.r_[:nreg], names, fontsize=8)
    plt.colorbar()

    plt.subplot(1, 2, 2)
    plt.title("Lengths")
    plt.imshow(lengths)
    plt.xticks(np.r_[:nreg], names, rotation='vertical', fontsize=8)
    plt.yticks(np.r_[:nreg], names, fontsize=8)
    plt.colorbar()

    plt.tight_layout()
    plt.savefig(out_img_file)
    plt.close()


def plot_ez_hypothesis(ez_hyp_file, tvb_zipfile, label_volume_file, t1_volume_file, out_img_file):
    label_vol = nib.load(label_volume_file)
    t1_vol = nib.load(t1_volume_file)
    ez_hyp_vector = np.genfromtxt(ez_hyp_file, dtype=int)
    ez_hyp_vol = nib.Nifti1Image(gen_volume_regions(ez_hyp_vector, label_vol),
                                 label_vol.affine)

    with zipfile.ZipFile(tvb_zipfile) as zf:
        with zf.open("centres.txt") as fl:
            region_names = np.genfromtxt(fl, usecols=(0,), dtype=str)
    ez_region_names = region_names[ez_hyp_vector == 1]


    lims = (np.dot(t1_vol.affine, [0, 0, 0, 1])[2],
            np.dot(t1_vol.affine, np.append(t1_vol.shape, 1))[2])

    t1_lims = np.min(t1_vol.get_data()), np.max(t1_vol.get_data())

    nrows = 5
    ncols = 6

    # Restrict the range a little bit
    # TODO: more rigorously? Or not at all?
    lims = (lims[0] + 0.2 * (lims[1] - lims[0]), lims[1] - 0.2 * (lims[1] - lims[0]))
    delta = (lims[1] - lims[0])/(nrows*ncols)

    plt.figure(figsize=(20, 20))

    for i in range(nrows * ncols):
        ax = plt.subplot(nrows + 1, ncols, ncols + i + 1)
        infsup_coord = lims[1] - i*delta

        img_data, img_extent = get_slice(t1_vol, 'transversal', infsup_coord)
        plt.imshow(img_data, extent=img_extent,
                   vmin=t1_lims[0], vmax=t1_lims[1], cmap='Greys_r', origin='lower', zorder=1)

        try:
            img_data, img_extent = get_slice(ez_hyp_vol, 'transversal', infsup_coord)
            plt.imshow(img_data, extent=img_extent,
                       vmin=0, vmax=1, cmap='bwr', origin='lower', zorder=2, alpha=0.5)
        except IndexError:
            pass

        plt.title("%.1f" % infsup_coord)
        plt.xlabel("L --- R")
        plt.ylabel("P --- A")

    plt.suptitle("Epileptogenic zone hypothesis", fontsize=20)
    plt.figtext(0.5, .92, "EZ: %s" % ", ".join(ez_region_names), fontsize=14, ha='center', wrap=True)
    plt.tight_layout()
    plt.savefig(out_img_file)
    plt.close()



if __name__ == '__main__':
    import sys
    loglevel = logging.INFO
    if os.environ.get('VERBOSE', False):
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel)
    cmd = sys.argv[1]
    eval(cmd)(*sys.argv[2:])
