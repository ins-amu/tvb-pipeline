

import itertools

import nibabel as nib
import numpy as np


def add_min_max(volume):
    """Ugly hack to deal with the MRtrix bug (?) that causes MRview to crop min/max values"""

    volume[0, 0, 0] = np.nanmin(volume) - 1
    volume[0, 0, 1] = np.nanmax(volume) + 1


def gen_volume_regions(values, label_volume):
    new_volume = np.zeros(label_volume.shape)
    new_volume[:] = np.nan

    for i, value in enumerate(values):
        region = i + 1
        mask = label_volume.get_data() == region

        new_volume[mask] = value

    return new_volume


def gen_volume_points(values, positions, ref_volume, ref_aff, dist=0):
    new_volume = np.zeros(ref_volume.shape)
    new_volume[:, :] = np.nan

    kx, ky, kz = np.mgrid[-dist:dist+1, -dist:dist+1, -dist:dist+1]

    for val, pos in zip(values, positions):
        ix, iy, iz = np.linalg.solve(ref_aff, np.append(pos, 1.0))[0:3].astype(int)
        #new_volume[inds[0], inds[1], inds[2]] = val
        new_volume[ix + kx, iy + ky, iz + kz] = val

    return new_volume


def save_nifti_regions(values, label_volume_file, out_file):
    label_nii = nib.load(label_volume_file)
    label_volume = label_nii.get_data()
    new_volume = gen_volume_regions(values, label_volume)
    add_min_max(new_volume)
    new_nii = nib.Nifti1Image(new_volume, label_nii.affine)
    nib.save(new_nii, out_file)


def save_nifti_points(values, names, position_file, ref_volume_file, out_file, skip_missing=False, dist=0):
    ref_nii = nib.load(ref_volume_file)

    contact_names = list(np.genfromtxt(position_file, dtype=str, usecols=(0,)))
    contact_positions = np.genfromtxt(position_file, dtype=float, usecols=(1, 2, 3))

    positions = np.zeros((len(values), 3))
    for i, name in enumerate(names):
        try:
            contact_ind = contact_names.index(name)
            pos = contact_positions[contact_ind, :]
        except ValueError:
            pos = np.array([np.nan, np.nan, np.nan])

        positions[i, :] = pos

    missing_mask = np.isnan(positions[:, 0])
    if skip_missing:
        values = values[~missing_mask]
        positions = positions[~missing_mask, :]
    else:
        if np.any(missing_mask):
            raise ValueError("Missing contact position(s) for: %s." % ", ".join([
                name for name, missing in zip(names, missing_mask) if missing]))


        #np.array(names)[missing_mask]))

    new_volume = gen_volume_points(values, positions, ref_nii.get_data(), ref_nii.affine, dist=dist)
    add_min_max(new_volume)

    new_nii = nib.Nifti1Image(new_volume, ref_nii.affine)
    nib.save(new_nii, out_file)


# ----------- Brain region <-> point mapping ------------------- #

coord_sequences = []

def gen_coord_sequence(affine_transform):
    ncells = 20

    indx = np.arange(-ncells, ncells + 1, 1)
    indy = np.arange(-ncells, ncells + 1, 1)
    indz = np.arange(-ncells, ncells + 1, 1)

    coords = list(itertools.product(indx, indy, indz))
    dists = np.array([np.linalg.norm(affine_transform[0:3, 0:3] @ np.array(c))
                      for c in coords])

    sdists, scoords = zip(*sorted(zip(dists, coords)))

    return {'coords': np.array(scoords), 'dists': np.array(sdists)}


def get_coord_sequence(affine_transform):
    global coord_sequences

    for mat, sequence in coord_sequences:
        if np.all(mat == affine_transform):
            return sequence

    # Add a new one
    sequence = gen_coord_sequence(affine_transform)
    coord_sequences.append((affine_transform, sequence))
    return sequence


def point_to_brain_region(point, label_volume, outside_index=0, tol=0.0):
    """
    Return the label_value at a given point. If the label_value is equal to outside_index,
    return the closest label_value (if it is below tol).
    """

    seq = get_coord_sequence(label_volume.affine)

    coords0 = np.linalg.solve(label_volume.affine, np.append(point, 1.0))[0:3].astype(int)
    data = label_volume.get_data()

    for i, (coords_delta, dist) in enumerate(zip(seq['coords'], seq['dists'])):
        if dist > tol:
            return outside_index

        coords1 = coords0 + coords_delta
        if (coords1[0] < 0 or coords1[0] >= data.shape[0] or
            coords1[1] < 0 or coords1[1] >= data.shape[1] or
            coords1[2] < 0 or coords1[2] >= data.shape[2]):
            continue

        if data[tuple(coords1)] != outside_index:
            return data[tuple(coords1)]

    raise ValueError("Run out of coord sequence!")

# -------------------------------------------------------------------- #
