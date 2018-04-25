

import enum
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

def voxel_neighbours(voxel, shape, r=1):
    """
    >>> voxel_neighbours((0, 0, 0), (3, 3, 1))
    [(0, 1, 0), (1, 0, 0), (1, 1, 0)]
    """
    i, j, k = voxel
    neighs = [(ii, jj, kk) for ii in np.arange(i - r, i + r + 1)
                           for jj in np.arange(j - r, j + r + 1)
                           for kk in np.arange(k - r, k + r + 1)
              if ((ii != i) or (jj != j) or (kk != k))
                  and (ii >= 0) and (jj >= 0) and (kk >=0)
                  and (ii < shape[0]) and (jj < shape[1]) and (kk < shape[2])]

    return neighs


def separate_components(volume, r=1):
    """
    Separate a boolean volume into components

    >>> nc, comps = separate_components(np.array([[[False, True,  False], \
                                                   [False, False, False], \
                                                   [True,  True,  True]]]))
    >>> nc
    2
    >>> comps.tolist()
    [[[-1, 0, -1], [-1, -1, -1], [1, 1, 1]]]
    """

    components = -1 * np.ones_like(volume, dtype=int)
    current_comp = -1

    for init_voxel in zip(*np.where(volume)):
        if components[init_voxel] == -1:
            current_comp += 1
            lifo = [init_voxel]
            while len(lifo) > 0:
                voxel = lifo.pop()
                if components[voxel] == -1 and volume[voxel] == True:
                    components[voxel] = current_comp
                    lifo.extend(voxel_neighbours(voxel, volume.shape, r))


    return current_comp + 1, components


def translate_ez_hypothesis(label_file_src, label_file_trg, ez_hyp_src, nreg_trg):
    """
    Translate ez_hyp_src ({0, 1} vector of length nregions_src) to ez_hyp_trg.
    First refers to the regions in label_file_src, second to label_file_trg.
    """

    TRG_EZ_THRESHOLD = 0.5

    if label_file_src == label_file_trg:
        return ez_hyp_src

    # -1 to account for the shift
    label_src = nib.load(label_file_src).get_data().astype(int) - 1
    label_trg = nib.load(label_file_trg).get_data().astype(int) - 1
    assert label_src.shape == label_trg.shape

    # Translate hypothesis in source segmentation to a voxel-based hypothesis
    ez_volume =  np.zeros_like(label_src, dtype=int)
    nreg_src = len(ez_hyp_src)
    for i in range(nreg_src):
        ez_volume[label_src == i] = ez_hyp_src[i]

    # Generate target hypothesis:
    # First, set as epileptogenic the regions with more than TRG_EZ_THRESHOLD of epileptogenic voxels
    ez_hyp_trg = np.zeros(nreg_trg, dtype=int)
    for i in range(nreg_trg):
        region_mask = label_trg == i
        if np.mean(ez_volume[region_mask]) > TRG_EZ_THRESHOLD:
            ez_hyp_trg[i] = 1

    # Second, make sure that for each region in source ez_hyp there is at least one in target ez_hyp
    for i in range(nreg_src):
        if ez_hyp_src[i] == 1:
            reg_mask = label_src == i
            max_overlapping_region = np.argmax([np.sum(reg_mask * (label_trg == j)) for j in range(nreg_trg)])
            ez_hyp_trg[max_overlapping_region] = 1

    return ez_hyp_trg


if __name__ == "__main__":
    import sys
    translate_ez_hypothesis(*sys.argv[1:])
