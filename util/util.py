#!/usr/bin/env python

import os
import logging
import numpy as np
import itertools
import csv

import nibabel as nb

def cores():
    import multiprocessing
    print(multiprocessing.cpu_count()//2)


def vol_val_xyz(vol, aff, val):
    vox_idx = np.argwhere(vol == val)
    xyz = aff.dot(np.c_[vox_idx, np.ones(vox_idx.shape[0])].T)[:3].T
    return xyz


def compute_label_volume_centers(label_volume, affine=None):
    try:
        vol = label_volume.get_data()
        aff = label_volume.affine
    except:
        vol = label_volume
        aff = affine
    for val in np.unique(vol):
        xyz = vol_val_xyz(vol, aff, val)
        x, y, z = xyz.mean(axis=0)
        yield val, (x, y, z)


def build_fs_label_name_map(lut_path):
    lut = {}
    with open(lut_path, 'r') as fd:
        for line in fd.readlines():
            if not line[0] == '#' and line.strip():
                val, name, _, _, _, _ = line.strip().split()
                lut[int(val)] = name
    return lut


def label_volume_centers(label_volume, output_tsv):
    log = logging.getLogger('label_volume_centers')

    log.info('reading %r', label_volume)
    vol = nb.load(label_volume)

    log.info('computing centers')
    centers = list(compute_label_volume_centers(vol))

    log.info('loading FS LUT')
    lut_path = os.path.join(os.environ['FREESURFER_HOME'], 'FreeSurferColorLUT.txt')
    lut_map = build_fs_label_name_map(lut_path)

    log.info('writing results to %r', output_tsv)
    with open(output_tsv, 'w') as fd:
        for val, (x, y, z) in centers:
            val_ = lut_map[val] if lut_map else val
            fd.write('%f\t%f\t%f\t%s\n' % (x, y, z, val_))


def label_with_dilation(masked_CT_fname, dilated_CT_fname, label_CT_fname):
    import scipy.ndimage
    log = logging.getLogger('label_with_dilation')
    log.info('reading mask %r', masked_CT_fname)
    mask = nb.load(masked_CT_fname)
    mask_data = mask.get_data()
    log.info('reading dilated mask %r', dilated_CT_fname)
    dil_mask = nb.load(dilated_CT_fname)
    dil_mask_data = dil_mask.get_data()
    log.info('labeling dilated mask..')
    lab, n = scipy.ndimage.label(dil_mask_data)
    log.info('found %d objects', n)
    lab_xyz = list(compute_label_volume_centers(lab, dil_mask.affine))
    lab_sort = np.r_[:n+1]
    # sort labels along AP axis
    for i, (val, _) in enumerate(sorted(lab_xyz, key=lambda t: t[1][1])):
        lab_sort[val] = i
    lab = lab_sort[lab]
    mask_data *= lab
    log.info('saving labeled mask to %r', label_CT_fname)
    label_CT = nb.nifti1.Nifti1Image(mask_data, mask.affine)
    nb.save(label_CT, label_CT_fname)

def periodic_xyz_for_object(lab, val, aff, bw=0.1, doplot=False):
    "Find blob centers for object in lab volume having value val."
    # TODO handle oblique with multiple spacing
    # TODO we know its 3.5 mm...
    # TODO constraints such as enough voxels, linear enough, etc
    # vox coords onto first mode
    log = logging.getLogger('periodic_xyz_for_object')
    vox_idx = np.argwhere(lab == val)
    xyz = aff.dot(np.c_[vox_idx, np.ones(vox_idx.shape[0])].T)[:3].T
    xyz_mean = xyz.mean(axis=0)
    xyz -= xyz_mean
    u, s, vt = np.linalg.svd(xyz, 0)
    xi = u[:, 0] * s[0]
    # histogram and ft to find spacing and offset
    bn, bxi_ = np.histogram(
        xi, np.r_[min(xi) - 0.5: max(xi) + 0.5: bw])
    bxi = bxi_[:-1] + bw / 2.0
    w = np.r_[3.0: 4.5: 1000j]
    f = (1.0 / w)[:, None]
    Bf = (np.exp(-2 * np.pi * 1j * bxi * f) * bn * bw).sum(axis=-1)
    i_peak = np.argmax(np.abs(Bf))
    theta = np.angle(Bf[i_peak])
    log.info('val=%d, dist=%f, theta=%f', val, 1 / f[i_peak][0], theta)
    xi_o = -theta / (2 * np.pi * f[i_peak])
    xi_pos = np.r_[xi_o: xi.max(): w[i_peak]]
    xi_neg = np.r_[-xi_o: -xi.min(): w[i_peak]]
    xi_pos = np.sort(np.r_[-xi_neg, xi_pos[1:]])
    xyz_pos = np.c_[xi_pos, np.zeros(
        (len(xi_pos), 2))].dot(vt) + xyz_mean
    return xyz_pos


def gen_seeg_xyz(labeled_CT_fname, schema_fname, seeg_xyz_fname):
    label_to_name = {}
    with open(schema_fname, "r") as fd:
        for line in fd.readlines():
            label_num, label_name = line.strip().split(' ')
            label_num = int(label_num)
            assert label_num not in label_to_name
            label_to_name[int(label_num)] = label_name
    nii = nb.load(labeled_CT_fname)
    lab_bin = nii.get_data()
    aff = nii.affine
    ulab = np.unique(lab_bin)
    ulab = ulab[ulab > 0]
    # TODO find closest voxel in parcellation, provide name
    fmt = '%s%d\t%f\t%f\t%f\n'
    with open(seeg_xyz_fname, "w") as fd:
        for ul in ulab[ulab > 0]:
            if ul not in label_to_name:
                print("skipping object label %d" % (ul, ))
                continue
            xyz_pos = periodic_xyz_for_object(lab_bin, ul, aff)
            xyz_pos = xyz_pos[np.argsort(np.abs(xyz_pos[:, 0]))]
            name = label_to_name[ul]
            for i, (x, y, z) in enumerate(xyz_pos):
                fd.write(fmt % (name, i, x, y, z))


def gen_contacts_on_electrode(name, target, entry, ncontacts, spacing_pattern):
    orientation = entry - target
    orientation /= np.linalg.norm(orientation)

    dists = itertools.chain([0.0], itertools.accumulate(itertools.cycle(spacing_pattern)))

    contacts = []
    for i in range(ncontacts):
        contact_name = name + str(i+1)
        position = target + next(dists)*orientation
        contacts.append((contact_name, position))

    return contacts

def transform(coords, src_img, dest_img, transform_mat):
    import subprocess
    coords_str = " ".join([str(x) for x in coords])

    cp = subprocess.run("echo %s | img2imgcoord -mm -src %s -dest %s -xfm %s" \
                            % (coords_str, src_img, dest_img, transform_mat),
                        shell=True, stdout=subprocess.PIPE)
    transformed_coords_str = cp.stdout.decode('ascii').strip().split('\n')[-1]
    return np.array([float(x) for x in transformed_coords_str.split(" ") if x])


def gen_seeg_xyz_from_endpoints(scheme_fname, out_fname, transform_mat=None,
                                src_img=None, dest_img=None):
    """
    Read the file with electrode endpoints (`scheme_fname`) and write the file with position of the electrode contacts
    (`out_fname`), possibly including linear transformation given either by the transformation matrix (`transform_mat`)
    or by source and destination image (`src_img` and `dest_img`).

    Each electrode in the schema file should be described by one line containing the following fields:

    Name Target_x Target_y Target_z Entry_x Entry_y Entry_z Num_contacts [Spacing_pattern]

    Spacing pattern should be a double quoted string with distances between the neighbouring contacts. If there are more
    contacts than elements in the spacing pattern, the pattern is repeated. If absent, default spacing "3.5" is used.
    All distances should be in mm.
    """

    DEFAULT_SPACING_PATTERN= "3.5"

    infile = open(scheme_fname, "r")
    outfile = open(out_fname, "w")

    for line in infile:
        line = line.strip()
        if not line or line[0] == '#':
            continue

        # Using csv.reader to allow for quoted strings
        items = next(csv.reader([line], delimiter=' ', quotechar='"'))
        # Skip empty fields created by multiple delimiters
        items = [item for item in items if item != ""]

        if len(items) == 8:
            # Using default spacing pattern of 3.5 mm
            name, tgx, tgy, tgz, enx, eny, enz, ncontacts = items
            spacing_pattern_str = DEFAULT_SPACING_PATTERN
        elif len(items) == 9:
            name, tgx, tgy, tgz, enx, eny, enz, ncontacts, spacing_pattern_str = items
        else:
            raise ValueError("Unexpected number of items:\n%s" % line)

        target = np.array([float(x) for x in [tgx, tgy, tgz]])
        entry  = np.array([float(x) for x in [enx, eny, enz]])
        ncontacts = int(ncontacts)
        spacing_pattern = [float(x) for x in spacing_pattern_str.split()]

        if transform_mat is not None:
            assert src_img is not None and dest_img is not None
            target = transform(target, src_img, dest_img, transform_mat)
            entry = transform(entry, src_img, dest_img, transform_mat)

        contacts = gen_contacts_on_electrode(name, target, entry, ncontacts, spacing_pattern)

        for contact_name, pos in contacts:
            outfile.write("%-6s %7.2f %7.2f %7.2f\n" % (contact_name, pos[0], pos[1], pos[2]))

    infile.close()
    outfile.close()


def transform_gardel_coords_to_tvb(gardel_file, src_img_file, target_img_file, transform_mat, tvb_coords_file):
    orig_coords = np.genfromtxt(gardel_file, names=True,
                                dtype=None,
                                usecols=(0, 1, 2, 3, 4))

    header = orig_coords.dtype.names
    if 'name' in header:
        names = orig_coords['name'].astype(str)
    elif 'Electrode' in header and 'Contact' in header:
        names = [name + str(ind) for name, ind in zip(orig_coords['Electrode'].astype(str), orig_coords['Contact'])]
    else:
        raise ValueError("Unexpected header in gardel file '%s'" % gardel_file)

    nsensors = len(names)

    src_img = nb.load(src_img_file)

    target_coords = np.zeros((nsensors, 3))
    for i in range(nsensors):
        vox_src = [orig_coords[i]['x'], orig_coords[i]['y'], orig_coords[i]['z']]
        ras_src = np.dot(src_img.affine, vox_src + [1])[0:3]
        target_coords[i, :] = transform(ras_src, src_img_file, target_img_file, transform_mat)

    with open(tvb_coords_file, 'w') as fl:
        for i in range(nsensors):
            fl.write("%-6s %7.2f %7.2f %7.2f\n" % (names[i],
                                                   target_coords[i, 0], target_coords[i, 1], target_coords[i, 2]))



def seeg_gain(seeg_xyz_fname, aa_xyz_fname, gain_mat_fname):
    # NB this is only a pseudo gain matrix
    seeg_xyz = np.loadtxt(seeg_xyz_fname, usecols=(1,2,3))  # (n, 3)
    aa_xyz = np.loadtxt(aa_xyz_fname, usecols=(0,1,2))  # (m, 3)
    dx = seeg_xyz[:, None] - aa_xyz  # (n, m, 3)
    r = np.sqrt((dx**2).sum(axis=-1))  # (n, m)
    gain = 1.0 / r
    np.savetxt(gain_mat_fname, gain)


def _label_objects_one(args):
    uv, labels, affine = args
    if uv == 0:
        return None
    xyz = vol_val_xyz(labels, affine, uv)
    if len(xyz) < 3:
        return None
    u, s, vt = np.linalg.svd(xyz, full_matrices=False)
    p1 = (u[:, 0] * s[0]).ptp()
    return s[0], uv, s, p1


def label_objects(masked_fname, labeled_fname, *args):
    import scipy.ndimage
    import multiprocessing
    nii = nb.load(masked_fname)
    labels, _ = scipy.ndimage.label(nii.get_data())
    log = logging.getLogger('label_objects')
    label_pcs = []
    uvs = np.unique(labels.reshape((-1, )))
    args = [(uv, labels, nii.affine) for uv in uvs]
    nthread = 8
    pool = multiprocessing.Pool(nthread)
    label_pcs = [r for r in pool.map(_label_objects_one, args) if r]
    label_map = np.zeros(5000)
    for i, (s0, uv, s, p1) in enumerate(sorted(label_pcs, key=lambda t: t[3])):
        label_map[uv] = i
        print(i, uv, s0, p1)
    labels = label_map[labels]
    nii_out = nb.nifti1.Nifti1Image(labels, nii.affine)
    nb.save(nii_out, labeled_fname)


def postprocess_connectome(in_mat_fname, out_mat_fname, *args):
    from numpy import loadtxt, fill_diagonal, savetxt
    mat = loadtxt(in_mat_fname)
    mat += mat.T
    fill_diagonal(mat, 0.0)
    savetxt(out_mat_fname, mat)


if __name__ == '__main__':
    import sys
    loglevel = logging.INFO
    if os.environ.get('VERBOSE', False):
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel)
    cmd = sys.argv[1]
    eval(cmd)(*sys.argv[2:])
