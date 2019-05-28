#!/usr/bin/env python3

import enum
import os
import logging
import re
import subprocess
import sys
import tempfile
import time
import glob
from typing import List, Optional
import unittest
from zipfile import ZipFile

import numpy as np

import nibabel as nib
import nibabel.freesurfer.io


""" For each parcellation, following tuple is supplied:
(
  Left hemisphere shift in LUT,
  Right hemishere shift in LUT
)
"""
PARC_SHIFTS = {
    'dk':        ( 1000,  2000),
    'destrieux': (11100, 12100),
    'vep':       (71000, 72000)
}



class StructuralDataset:
    def __init__(self,
                 orientations: np.ndarray,
                 areas: np.ndarray,
                 centers: np.ndarray,
                 cortical: np.ndarray,
                 weights_ut: np.ndarray,
                 tract_lengths_ut: np.ndarray,
                 names: List[str],
                 volumes: np.ndarray=None):

        nregions = len(names)

        assert orientations.shape == (nregions, 3)
        assert areas.shape == (nregions,)
        assert centers.shape == (nregions, 3)
        assert cortical.shape == (nregions,)
        assert weights_ut.shape == (nregions, nregions)
        assert tract_lengths_ut.shape == (nregions, nregions)

        # Upper triangular -> symmetric matrices
        assert np.sum(np.tril(weights_ut, -1)) == 0
        assert np.sum(np.tril(tract_lengths_ut, -1)) == 0
        self.weights = weights_ut + weights_ut.transpose() - np.diag(np.diag(weights_ut))
        self.tract_lengths = tract_lengths_ut + tract_lengths_ut.transpose() - np.diag(np.diag(tract_lengths_ut))

        self.orientations = orientations
        self.areas = areas
        self.centers = centers
        self.cortical = cortical
        self.names = names
        self.volumes = volumes

    def save_to_txt_zip(self, filename: os.PathLike):

        tmpdir = tempfile.TemporaryDirectory()

        file_areas = os.path.join(tmpdir.name, 'areas.txt')
        file_orientations = os.path.join(tmpdir.name, 'average_orientations.txt')
        file_centres = os.path.join(tmpdir.name, 'centres.txt')
        file_cortical = os.path.join(tmpdir.name, 'cortical.txt')
        file_weights = os.path.join(tmpdir.name, 'weights.txt')
        file_tract_lengths = os.path.join(tmpdir.name, 'tract_lengths.txt')
        file_volumes = os.path.join(tmpdir.name, 'volumes.txt')

        np.savetxt(file_areas, self.areas, fmt='%.2f')
        np.savetxt(file_orientations, self.orientations, fmt='%.2f %.2f %.2f')
        np.savetxt(file_cortical, self.cortical, fmt='%d')
        np.savetxt(file_weights, self.weights, fmt='%d')
        np.savetxt(file_tract_lengths, self.tract_lengths, fmt='%.3f')
        if self.volumes is not None:
            np.savetxt(file_volumes, self.volumes, fmt='%.2f')

        with open(file_centres, 'w') as f:
            for i, name in enumerate(self.names):
                f.write('%s %.4f %.4f %.4f\n' % (name, self.centers[i, 0], self.centers[i, 1], self.centers[i, 2]))

        with ZipFile(filename, 'w') as zip_file:
            zip_file.write(file_areas, os.path.basename(file_areas))
            zip_file.write(file_orientations, os.path.basename(file_orientations))
            zip_file.write(file_centres, os.path.basename(file_centres))
            zip_file.write(file_cortical, os.path.basename(file_cortical))
            zip_file.write(file_weights, os.path.basename(file_weights))
            zip_file.write(file_tract_lengths, os.path.basename(file_tract_lengths))
            if self.volumes is not None:
                zip_file.write(file_volumes, os.path.basename(file_volumes))


class Hemisphere(enum.Enum):
    rh = 'rh'
    lh = 'lh'


class Surface:
    def __init__(self, vertices: np.array, triangles: np.array, region_mapping: np.array):
        assert vertices.ndim == 2
        assert triangles.ndim == 2
        assert region_mapping.ndim == 1

        assert vertices.shape[1] == 3
        assert triangles.shape[1] == 3
        assert region_mapping.shape[0] == vertices.shape[0]

        self.vertices = vertices
        self.triangles = triangles
        self.region_mapping = region_mapping

        self.nverts = self.vertices.shape[0]
        self.ntriangs = self.triangles.shape[0]

        self.vertex_triangles = compute_vertex_triangles(self.nverts, self.ntriangs, self.triangles)
        self.triangle_normals = compute_triangle_normals(self.triangles, self.vertices)
        self.triangle_angles = compute_triangle_angles(self.vertices, self.ntriangs, self.triangles)
        self.vertex_normals = compute_vertex_normals(self.nverts, self.vertex_triangles, self.triangles,
                                                     self.triangle_angles, self.triangle_normals, self.vertices)
        self.triangle_areas = compute_triangle_areas(self.vertices, self.triangles)

    def save_surf_zip(self, filename):
        tmpdir = tempfile.TemporaryDirectory()

        file_vertices = os.path.join(tmpdir.name, 'vertices.txt')
        file_triangles = os.path.join(tmpdir.name, 'triangles.txt')
        file_normals = os.path.join(tmpdir.name, 'normals.txt')

        np.savetxt(file_vertices, self.vertices, fmt='%.6f %.6f %.6f')
        np.savetxt(file_triangles, self.triangles, fmt='%d %d %d')
        np.savetxt(file_normals, self.vertex_normals, fmt='%.6f %.6f %.6f')

        with ZipFile(filename, 'w') as zip_file:
            zip_file.write(file_vertices, os.path.basename(file_vertices))
            zip_file.write(file_triangles, os.path.basename(file_triangles))
            zip_file.write(file_normals, os.path.basename(file_normals))

    def remap(self, remap_dict):
        for old_ind, new_ind in remap_dict.items():
            self.region_mapping[self.region_mapping == old_ind] = new_ind

    def save_region_mapping_txt(self, filename):
        np.savetxt(filename, self.region_mapping, fmt="%d")


class ColorLut:
    def __init__(self, filename: os.PathLike):
        table = np.genfromtxt(os.fspath(filename), dtype=None)

        if len(table.dtype) == 6:
            # id name R G B A
            self.inds = table[table.dtype.names[0]]
            self.names = table[table.dtype.names[1]].astype('U')
            self.r = table[table.dtype.names[2]]
            self.g = table[table.dtype.names[3]]
            self.b = table[table.dtype.names[4]]
            self.a = table[table.dtype.names[5]]
            self.shortnames = np.zeros(len(self.names), dtype='U')

        elif len(table.dtype) == 7:
            # id shortname name R G B A
            self.inds = table[table.dtype.names[0]]
            self.shortnames = table[table.dtype.names[1]].astype('U')
            self.names = table[table.dtype.names[2]].astype('U')
            self.r = table[table.dtype.names[3]]
            self.g = table[table.dtype.names[4]]
            self.b = table[table.dtype.names[5]]
            self.a = table[table.dtype.names[6]]


class RegionIndexMapping:

    def __init__(self, color_lut_src_file: os.PathLike, color_lut_trg_file: os.PathLike):
        self.src_table = ColorLut(color_lut_src_file)
        self.trg_table = ColorLut(color_lut_trg_file)

        names_to_trg = dict(zip(self.trg_table.names, self.trg_table.inds))

        self.src_to_trg = dict()
        for src_ind, src_name in zip(self.src_table.inds, self.src_table.names):
            trg_ind = names_to_trg.get(src_name, None)
            if trg_ind is not None:
                self.src_to_trg[src_ind] = trg_ind

        self.unknown_ind = names_to_trg.get('Unknown', 0)   # zero as the default unknown area

    def source_to_target(self, index):
        return self.src_to_trg.get(index, self.unknown_ind)


def merge_surfaces(surfaces: List[Surface]) -> Surface:
    offsets = np.cumsum([0] + [vs.shape[0] for vs in [surf.vertices for surf in surfaces]][:-1])
    vertices = np.vstack([surf.vertices for surf in surfaces])
    triangles = np.vstack([ts + offset for ts, offset in zip([surf.triangles for surf in surfaces], offsets)])
    region_mappings = np.hstack([surf.region_mapping for surf in surfaces])
    return Surface(vertices, triangles, region_mappings)


def compute_vertex_triangles(number_of_vertices, number_of_triangles, triangles):
    vertex_triangles = [[] for _ in range(number_of_vertices)]
    for k in range(number_of_triangles):
        vertex_triangles[triangles[k, 0]].append(k)
        vertex_triangles[triangles[k, 1]].append(k)
        vertex_triangles[triangles[k, 2]].append(k)
    return vertex_triangles


def compute_triangle_normals(triangles, vertices):
    """Calculates triangle normals."""
    tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
    tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
    tri_norm = np.cross(tri_u, tri_v)

    try:
        triangle_normals = tri_norm / np.sqrt(np.sum(tri_norm ** 2, axis=1))[:, np.newaxis]
    except FloatingPointError:
        # TODO: NaN generation would stop execution, however for normals this case could maybe be
        #  handled in a better way.
        triangle_normals = tri_norm
    return triangle_normals


def compute_triangle_angles(vertices, number_of_triangles, triangles):
    """
    Calculates the inner angles of all the triangles which make up a surface
    """
    verts = vertices
    # TODO: Should be possible with arrays, ie not nested loops...
    # A short profile indicates this function takes 95% of the time to compute normals
    # (this was a direct translation of some old matlab code)
    angles = np.zeros((number_of_triangles, 3))
    for tt in range(number_of_triangles):
        triangle = triangles[tt, :]
        for ta in range(3):
            ang = np.roll(triangle, -ta)
            angles[tt, ta] = np.arccos(np.dot(
                (verts[ang[1], :] - verts[ang[0], :]) /
                np.sqrt(np.sum((verts[ang[1], :] - verts[ang[0], :]) ** 2, axis=0)),
                (verts[ang[2], :] - verts[ang[0], :]) /
                np.sqrt(np.sum((verts[ang[2], :] - verts[ang[0], :]) ** 2, axis=0))))

    return angles


def compute_vertex_normals(number_of_vertices, vertex_triangles, triangles,
                           triangle_angles, triangle_normals, vertices):
    """
    Estimates vertex normals, based on triangle normals weighted by the
    angle they subtend at each vertex...
    """
    vert_norms = np.zeros((number_of_vertices, 3))
    bad_normal_count = 0
    for k in range(number_of_vertices):
        try:
            tri_list = list(vertex_triangles[k])
            angle_mask = triangles[tri_list, :] == k
            angles = triangle_angles[tri_list, :]
            angles = angles[angle_mask][:, np.newaxis]
            angle_scaling = angles / np.sum(angles, axis=0)
            vert_norms[k, :] = np.mean(angle_scaling * triangle_normals[tri_list, :], axis=0)
            # Scale by angle subtended.
            vert_norms[k, :] = vert_norms[k, :] / np.sqrt(np.sum(vert_norms[k, :] ** 2, axis=0))
            # Normalise to unit vectors.
        except (ValueError, FloatingPointError):
            # If normals are bad, default to position vector
            # A nicer solution would be to detect degenerate triangles and ignore their
            # contribution to the vertex normal
            vert_norms[k, :] = vertices[k] / np.sqrt(vertices[k].dot(vertices[k]))
            bad_normal_count += 1
    if bad_normal_count:
        print(" %d vertices have bad normals" % bad_normal_count)
    return vert_norms


def compute_triangle_areas(vertices, triangles):
    """Calculates the area of triangles making up a surface."""
    tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
    tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
    tri_norm = np.cross(tri_u, tri_v)
    triangle_areas = np.sqrt(np.sum(tri_norm ** 2, axis=1)) / 2.0
    triangle_areas = triangle_areas[:, np.newaxis]
    return triangle_areas


def compute_region_orientations(regions, vertex_normals, region_mapping):
    """Compute the orientation of given regions from vertex_normals and region mapping"""

    average_orientation = np.zeros((regions.size, 3))
    # Average orientation of the region
    for i, k in enumerate(regions):
        orient = vertex_normals[region_mapping == k, :]
        if orient.shape[0] > 0:
            avg_orient = np.mean(orient, axis=0)
            average_orientation[i, :] = avg_orient / np.sqrt(np.sum(avg_orient ** 2))

    return average_orientation


def compute_region_areas(regions, triangle_areas, vertex_triangles, region_mapping):
    """Compute the areas of given regions"""

    region_surface_area = np.zeros(regions.size)
    avt = np.array(vertex_triangles)
    # NOTE: Slightly overestimates as it counts overlapping border triangles,
    #       but, not really a problem provided triangle-size << region-size.
    for i, k in enumerate(regions):
        regs = map(set, avt[region_mapping == k])
        region_triangles = set.union(*regs)
        if region_triangles:
            region_surface_area[i] = triangle_areas[list(region_triangles)].sum()

    return region_surface_area


def compute_region_centers(regions, vertices, region_mapping):

    region_centers = np.zeros((regions.size, 3))
    for i, k in enumerate(regions):
        vert = vertices[region_mapping == k, :]
        if vert.shape[0] > 0:
            region_centers[i, :] = np.mean(vert, axis=0)

    return region_centers


def compute_region_params(surface: Surface, subcortical: bool=False)\
        -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
    verts, triangs, region_mapping = surface.vertices, surface.triangles, surface.region_mapping

    regions = np.unique(region_mapping)
    areas = compute_region_areas(regions, surface.triangle_areas, surface.vertex_triangles, region_mapping)
    orientations = compute_region_orientations(regions, surface.vertex_normals, region_mapping)
    centers = compute_region_centers(regions, verts, region_mapping)

    return regions, areas, orientations, centers


def get_region_volumes(labelvol_file):
    labelvol = nib.load(labelvol_file)
    labels = labelvol.get_data()
    nreg = np.max(labels)

    dx, dy, dz = labelvol.header.get_zooms()
    voxel_vol = dx*dy*dz

    volumes = np.zeros(nreg+1, dtype=float)
    for i in range(0, nreg+1):
        volumes[i] = voxel_vol * np.sum(labels == i)

    return volumes



def extract_vector(string: str, name: str) -> Optional[List[float]]:
    r"""
    Extract numerical vector from a block of text. The vector has to be on a single line with the format:
    <name> : (x0, x1, x2 [,...])
    If the vector in the correct format is missing, return None.

    >>> extract_vector("EXAMPLE\na: (1.0, 2.0, 3.0)\nb: (0.0, 0.0, 0.0)", "a")
    [1.0, 2.0, 3.0]

    >>> extract_vector("EMPTY", "a") is None
    True
    """

    for line in string.split("\n"):
        match = re.match(r"""^\s*
                         (.+?)              # name
                         \s*:\s*            # separator
                         \(([0-9.,\s-]+)\)   # vector: (x0, x1, ....)
                         \s*$""",
                         line, re.X)

        if match and match.group(1) == name:
            try:
                vector = [float(x) for x in match.group(2).split(",")]
                return vector
            except ValueError:
                pass

    return None


def pial_to_verts_and_triangs(pial_surf) -> (np.ndarray, np.ndarray):

    tmpdir = tempfile.TemporaryDirectory()
    pial_asc = os.path.join(tmpdir.name, os.path.basename(pial_surf + ".asc"))
    subprocess.run(['mris_convert', pial_surf, pial_asc])

    with open(pial_asc, 'r') as f:
        f.readline()
        nverts, ntriangs = [int(n) for n in f.readline().strip().split(' ')]

    vertices = np.genfromtxt(pial_asc, dtype=float, skip_header=2, skip_footer=ntriangs, usecols=(0, 1, 2))
    triangles = np.genfromtxt(pial_asc, dtype=int, skip_header=2+nverts, usecols=(0, 1, 2))
    assert vertices.shape == (nverts, 3)
    assert triangles.shape == (ntriangs, 3)

    completed_process = subprocess.run(["mris_info", pial_surf], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    mris_info = completed_process.stdout.decode('ascii')
    c_ras_list = extract_vector(mris_info, "c_(ras)")
    assert c_ras_list is not None
    vertices[:, 0:3] += np.array(c_ras_list)

    return vertices, triangles


def read_cortical_region_mapping(label_direc: os.PathLike, hemisphere: Hemisphere, fs_to_conn: RegionIndexMapping,
                                 parcellation) -> np.ndarray:
    """
    Read the region mapping of a cortical surface.

    When constructing the connectome, one needs to convert the FreeSurfer labels to MRtrix labels. For a volume, this
    can be done by `labelconvert`. For the surfaces, however, we have to read the FreeSurfer annotation, and then
    convert them to the MRtrix labels ourselves.

    To make things more complicated, the annotation file does not store the region indices, so some adjustment
    (hemisphere and parcellation dependent) is necessary. Ugh.
    """

    lut_lh_shift, lut_rh_shift = PARC_SHIFTS[parcellation]

    parc_file = "%s.aparc.%s.annot" % (hemisphere.value, parcellation)
    filename = os.path.join(label_direc, parc_file)
    region_mapping, _, _ = nibabel.freesurfer.io.read_annot(filename)

    region_mapping[region_mapping == -1] = 0   # Unknown regions in hemispheres

    # $FREESURFER_HOME/FreeSurferColorLUT.txt describes the shift
    if hemisphere == Hemisphere.lh:
        region_mapping += lut_lh_shift
    else:
        region_mapping += lut_rh_shift

    fs_to_conn_fun = np.vectorize(lambda n: fs_to_conn.source_to_target(n))
    region_mapping = fs_to_conn_fun(region_mapping)

    return region_mapping


def get_cortical_surfaces(cort_surf_direc: os.PathLike, label_direc: os.PathLike,
                          region_index_mapping: RegionIndexMapping, parcellation: str) -> Surface:
    verts_l, triangs_l = pial_to_verts_and_triangs(os.path.join(cort_surf_direc, Hemisphere.lh.value + ".pial"))
    verts_r, triangs_r = pial_to_verts_and_triangs(os.path.join(cort_surf_direc, Hemisphere.rh.value + ".pial"))

    region_mapping_l = read_cortical_region_mapping(label_direc, Hemisphere.lh, region_index_mapping, parcellation)
    region_mapping_r = read_cortical_region_mapping(label_direc, Hemisphere.rh, region_index_mapping, parcellation)

    surface = merge_surfaces([Surface(verts_l, triangs_l, region_mapping_l),
                              Surface(verts_r, triangs_r, region_mapping_r)])

    return surface


def get_subcortical_surfaces(subcort_surf_direc: os.PathLike, region_index_mapping: RegionIndexMapping,
                             parcellation: str) -> Surface:
    surfaces = []

    subcort_surf_files = glob.glob(os.path.join(subcort_surf_direc, "aseg_*.srf"))

    for filename in subcort_surf_files:
        match = re.match("^aseg_(\d+).srf$", os.path.split(filename)[1])
        fs_idx = int(match.group(1))
        conn_idx = region_index_mapping.source_to_target(fs_idx)
        with open(filename, 'r') as f:
            f.readline()
            nverts, ntriangs = [int(n) for n in f.readline().strip().split(' ')]

        vertices = np.genfromtxt(filename, dtype=float, skip_header=2, skip_footer=ntriangs, usecols=(0, 1, 2))
        triangles = np.genfromtxt(filename, dtype=int, skip_header=2 + nverts, usecols=(0, 1, 2))
        region_mapping = conn_idx * np.ones(nverts, dtype=int)
        surfaces.append(Surface(vertices, triangles, region_mapping))

    surface = merge_surfaces(surfaces)
    return surface


class MinimalSurfaceTest(unittest.TestCase):
    FLOAT_TOL = 1e-16

    def setUp(self):
        super().setUp()

        self.surf1 = Surface(np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]),
                             np.array([[0, 1, 2], [1, 3, 2]]),
                             np.array([1, 1, 1, 1]))
        self.surf2 = Surface(np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1]]),
                             np.array([[0, 1, 2], [1, 3, 2]]),
                             np.array([1, 1, 1, 1]))

        self.nverts1 = self.surf1.vertices.shape[0]
        self.nverts2 = self.surf2.vertices.shape[0]
        self.ntri1 = self.surf1.triangles.shape[0]
        self.ntri2 = self.surf2.triangles.shape[0]

    def test_merge_surfaces(self):
        surf_merged = merge_surfaces([self.surf1, self.surf2])

        self.assertEqual(surf_merged.vertices.shape, (self.nverts1 + self.nverts2, 3))
        self.assertEqual(surf_merged.triangles.shape, (self.ntri1 + self.ntri2, 3))
        self.assertEqual(surf_merged.region_mapping.shape, (self.nverts1 + self.nverts2,))

    def test_compute_triangle_normals(self):
        normals = compute_triangle_normals(self.surf1.triangles, self.surf1.vertices)
        self.assertTrue(np.allclose(normals,
                                    np.array([[0., 0., 1.], [0., 0., 1.]]),
                                    atol=self.FLOAT_TOL))

    def test_compute_vertex_normals(self):
        vert_triangles = compute_vertex_triangles(self.nverts1, self.ntri1, self.surf1.triangles)
        tri_angles = compute_triangle_angles(self.surf1.vertices, self.ntri1, self.surf1.triangles)
        tri_normals = compute_triangle_normals(self.surf1.triangles, self.surf1.vertices)
        vert_normals = compute_vertex_normals(self.nverts1, vert_triangles, self.surf1.triangles,
                                                   tri_angles, tri_normals, self.surf1.vertices)
        self.assertTrue(np.allclose(vert_normals,
                                    np.array([[0., 0., 1.], [0., 0., 1.], [0., 0., 1.], [0., 0., 1.]]),
                                    atol=self.FLOAT_TOL))

    def test_compute_triangle_areas(self):
        areas = compute_triangle_areas(self.surf1.vertices, self.surf1.triangles)
        self.assertTrue(np.allclose(areas,
                                    np.array([0.5, 0.5]),
                                    atol=self.FLOAT_TOL))

    def test_compute_region_params(self):
        surf = merge_surfaces([self.surf1, self.surf2])
        regions, areas, orientations, centers = compute_region_params(surf, False)

        self.assertTrue(np.equal(regions, [1]))
        self.assertTrue(np.allclose(areas, [2.0], atol=self.FLOAT_TOL))
        self.assertTrue(np.allclose(orientations, np.array([1, 0, 1])*np.sqrt(2)/2, atol=self.FLOAT_TOL))
        self.assertTrue(np.allclose(centers, [0.25, 0.5, 0.25], atol=self.FLOAT_TOL))


def create_tvb_dataset(cort_surf_direc: os.PathLike,
                       label_direc: os.PathLike,
                       subcort_surf_direc: os.PathLike,
                       source_lut: os.PathLike,
                       target_lut: os.PathLike,
                       parcellation: str,
                       weights_file: os.PathLike,
                       tract_lengths_file: os.PathLike,
                       label_volume_file: os.PathLike,
                       struct_zip_file: os.PathLike,
                       out_surfaces_dir: os.PathLike=None,
                       include_unknown: bool=False):
    """
    Parameters
    ----------
    cort_surf_direc: Directory that should contain:
                       rh.pial
                       lh.pial
    label_direc: Directory that should contain:
                   rh.aparc.annot
                   lh.aparc.annot
    subcort_surf_direc: Directory that should contain:
                          aseg_<NUM>.srf
                       for each <NUM> in SUBCORTICAL_REG_INDS
    source_lut: File with the color look-up table used for the original parcellation
    target_lut: File with the color look-up table used for the connectome generation
    parcellation: Parcellation type, currently either 'dk', 'destrieux', or 'vep'
    weights_file: text file with weights matrix (which should be upper triangular)
    tract_lengths_file: text file with tract length matrix (which should be upper triangular)
    label_volume_file: 3D volume with region labels
    struct_zip_file: zip file containing TVB structural dataset to be created
    out_surfaces_dir: directory where to put the surfaces and region mappings in TVB format
    include_unknown: include also unknown regions in the connectome
    """

    log = logging.getLogger('create_tvb_dataset')
    log.info('start')
    tic = time.time()

    region_index_mapping = RegionIndexMapping(source_lut, target_lut)

    surf_subcort = get_subcortical_surfaces(subcort_surf_direc, region_index_mapping, parcellation)
    surf_cort = get_cortical_surfaces(cort_surf_direc, label_direc, region_index_mapping, parcellation)

    region_params_subcort = compute_region_params(surf_subcort, True)
    region_params_cort = compute_region_params(surf_cort, False)

    nregions = max(region_index_mapping.trg_table.inds) + 1
    orientations = np.zeros((nregions, 3))
    areas = np.zeros(nregions)
    centers = np.zeros((nregions, 3))
    cortical = np.zeros(nregions, dtype=bool)
    volumes = get_region_volumes(label_volume_file)

    for region_params, is_cortical in [(region_params_subcort, False), (region_params_cort, True)]:
        regions, reg_areas, reg_orientations, reg_centers = region_params

        orientations[regions, :] = reg_orientations
        areas[regions] = reg_areas
        centers[regions, :] = reg_centers
        cortical[regions] = is_cortical

    weights = np.genfromtxt(os.fspath(weights_file))
    tract_lengths = np.genfromtxt(os.fspath(tract_lengths_file))

    if not include_unknown:
        # Remove the region from orientations, areas and centers
        indices = list(range(0, region_index_mapping.unknown_ind)) \
                  + list(range(region_index_mapping.unknown_ind + 1, nregions))

        names = region_index_mapping.trg_table.names[indices]
        orientations = orientations[indices]
        areas = areas[indices]
        centers = centers[indices]
        cortical = cortical[indices]
        volumes = volumes[indices]

        remap_dict = {ind: ind if ind < region_index_mapping.unknown_ind else ind - 1 for ind in range(nregions)}
        remap_dict[region_index_mapping.unknown_ind] = -1
        surf_subcort.remap(remap_dict)
        surf_cort.remap(remap_dict)

    else:
        # Add the region to weights and tract lengths
        names = region_index_mapping.trg_table.names

        np.insert(weights, region_index_mapping.unknown_ind, 0.0, axis=0)
        np.insert(weights, region_index_mapping.unknown_ind, 0.0, axis=1)
        np.insert(tract_lengths, region_index_mapping.unknown_ind, 0.0, axis=0)
        np.insert(tract_lengths, region_index_mapping.unknown_ind, 0.0, axis=1)

    dataset = StructuralDataset(orientations, areas, centers, cortical, weights, tract_lengths, names, volumes)
    dataset.save_to_txt_zip(struct_zip_file)

    if out_surfaces_dir:
        surf_subcort.save_region_mapping_txt(os.path.join(out_surfaces_dir, "region_mapping_subcort.%s.txt" % parcellation))
        surf_subcort.save_surf_zip(os.path.join(out_surfaces_dir, "surface_subcort.%s.zip" % parcellation))
        surf_cort.save_region_mapping_txt(os.path.join(out_surfaces_dir, "region_mapping_cort.%s.txt" % parcellation))
        surf_cort.save_surf_zip(os.path.join(out_surfaces_dir, "surface_cort.%s.zip" % parcellation))

    log.info('complete in %0.2fs', time.time() - tic)


def main():
    subject_dir, source_lut, target_lut, parcellation, weights_file, tract_lengths_file, out_file, out_surf = sys.argv[1:]

    logging.basicConfig(level=logging.INFO)

    create_tvb_dataset(
        os.path.join(subject_dir, "surf"),
        os.path.join(subject_dir, "label"),
        os.path.join(subject_dir, "aseg2srf", parcellation),
        source_lut,
        target_lut,
        parcellation,
        weights_file,
        tract_lengths_file,
        os.path.join(subject_dir, "dwi/label_in_T1.%s.nii.gz" % parcellation),
        out_file,
        out_surf
    )

if __name__ == "__main__":
    main()
