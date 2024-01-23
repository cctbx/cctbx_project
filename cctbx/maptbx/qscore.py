"""
This code provides methods to calculate the qscore metric for map-model validation,
as developed by Pintile et al.

As with the original implementation, multiple calculation options are provided.
The fastest (default) is to pass an mmtbx map-model-manager (mmm)
to the the function qscore_np, qscore_np(mmm)
"""

from __future__ import division
import math
import sys
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from itertools import chain

import numpy as np
import numpy.ma as ma

from scipy.spatial import KDTree

import cctbx
from cctbx.array_family import flex
from scitbx_array_family_flex_ext import bool as flex_bool

# from .qscore_utils import (
#     trilinear_interpolation,
#     rowwise_corrcoef,
#     sphere_points_np,
#     sphere_points_flex,
#     cdist_flex,
#     query_ball_point_flex,
#     query_atom_neighbors,
# )
# from .flex_utils import optimized_nd_to_1d_indices, nd_to_1d_indices, flex_std

try:
    from tqdm import tqdm
except ImportError:
    from .qscore_utils import DummyTQDM as tqdm


def radial_shell_worker_v1_np(args):
    """
    Calulate qscore for a single radial shell using version 1 (serial probe allocation) and numpy
    """
    (
        i,
        atoms_xyz,
        n_probes,
        radius_shell,
        tree,
        rtol,
        selection,
        n_probes_target,
    ) = args

    #
    # manage selection input
    if selection is None:
        selection = np.arange(len(atoms_xyz))
    else:
        assert selection.dtype == bool
    # do selection
    atoms_xyz_sel = atoms_xyz[selection]
    # print("sel_shape",atoms_xyz_sel.shape)
    n_atoms = atoms_xyz_sel.shape[0]

    if radius_shell == 0:
        radius_shell = 1e-9  # zero causes crash
    numPts = n_probes_target
    RAD = radius_shell
    outRAD = rtol
    all_pts = []  # list of probe arrays for each atom
    probe_xyz_r = np.full((n_atoms, n_probes_target, 3), -1.0)
    # print(atoms_xyz_sel)
    # print("n_atoms",n_atoms)
    for atom_i in range(n_atoms):
        coord = atoms_xyz_sel[atom_i]
        # print("atom_i",atom_i)
        # print("coord:",coord)
        pts = []

        # try to get at least [numPts] points at [RAD] distance
        # from the atom, that are not closer to other atoms
        for i in range(0, 50):
            # if we find the necessary number of probes in the first iteration, then i will never go to 1
            # points on a sphere at radius RAD...
            n_pts_to_grab = (
                numPts + i * 2
            )  # progressively more points are grabbed  with each failed iter
            # print("n_to_grab:",n_pts_to_grab)

            outPts = sphere_points_np(
                coord[None, :], RAD, n_pts_to_grab
            )  # get the points
            outPts = outPts.reshape(-1, 3)
            at_pts, at_pts_i = [None] * len(outPts), 0
            # print("probe candidates")

            for pt_i, pt in enumerate(
                outPts
            ):  # identify which ones to keep, progressively grow pts list
                # print(f"\t{pt[0]},{pt[1]},{pt[2]}")
                # query kdtree to find probe-atom interactions
                counts = tree.query_ball_point(
                    pt[None, :], RAD * outRAD, return_length=True
                )

                # each value in counts is the number of atoms within radius+tol of each probe
                count = counts.flatten()[0]
                ptsNear = count

                if ptsNear == 0:
                    at_pts[at_pts_i] = pt
                    at_pts_i += 1
                # if at_pts_i >= numPts:
                #     break

            if (
                at_pts_i >= numPts
            ):  # if we have enough points, take all the "good" points from this iter
                pts.extend(at_pts[0:at_pts_i])
                break
        # assert len(pts)>0, "Zero probes were found "
        pts = np.array(pts)  # should be shape (n_probes,3)
        all_pts.append(pts)

    # prepare output
    n_atoms = len(atoms_xyz)
    for i, r in enumerate(all_pts):
        if r.ndim == 2 and len(r) > 0:
            probe_xyz_r[i, :n_probes, :] = r[:n_probes_target, :]

    keep_sel = probe_xyz_r != -1.0
    keep_sel = np.mean(keep_sel, axis=-1, keepdims=True)
    keep_sel = np.squeeze(keep_sel, axis=-1)

    return probe_xyz_r, keep_sel.astype(bool)


def radial_shell_worker_v2_np(args):
    """
    Calulate qscore for a single radial shell using version 2 (parallel probe allocation) and numpy
    """
    # unpack args
    i, atoms_xyz, n_probes, radius_shell, tree, rtol, selection = args

    # manage selection input
    if selection is None:
        selection = np.arange(len(atoms_xyz))
    else:
        assert selection.dtype == bool

    # do selection
    atoms_xyz_sel = atoms_xyz[selection]
    n_atoms = atoms_xyz_sel.shape[0]

    # get probe coordinates
    probe_xyz = sphere_points_np(atoms_xyz_sel, radius_shell, n_probes)

    counts = tree.query_ball_point(
        probe_xyz, radius_shell * rtol, return_length=True
    )  # atom counts for each probe, for probes in shape (n_atoms,n_probes)
    probe_mask = (
        counts == 0
    )  # keep probes with 0 nearby atoms. The rtol ensures self is not counted

    return probe_xyz, probe_mask


def radial_shell_mp_np(

    model,
    n_probes=32,
    radii=np.linspace(0.1, 2, 12),
    rtol=0.9,
    num_processes=cpu_count(),
    selection=None,
    version=2,
    log=sys.stdout,
):
    """
    Generate probes for a model file using numpy
    """
    assert version in [1, 2], "Version must be 1 or 2"

    atoms_xyz = model.get_sites_cart().as_numpy_array()
    tree = KDTree(atoms_xyz)

    if version == 1:
        worker_func = radial_shell_worker_v1_np
        n_probes_target = n_probes
        # Create argument tuples for each chunk
        args = [
            (
                i,
                atoms_xyz,
                n_probes,
                radius_shell,
                tree,
                rtol,
                selection,
                n_probes_target,
            )
            for i, radius_shell in enumerate(radii)
        ]
    else:
        worker_func = radial_shell_worker_v2_np
        # Create argument tuples for each chunk
        args = [
            (i, atoms_xyz, n_probes, radius_shell, tree, rtol, selection)
            for i, radius_shell in enumerate(radii)
        ]

    # Create a pool of worker processes
    if num_processes > 1:
        with Pool(num_processes) as p:
            # Use the pool to run the trilinear_interpolation_worker function in parallel
            results = list(
                tqdm(p.imap(worker_func, args), total=len(args), file=log)
            )
    else:
        results = []
        for arg in tqdm(args, file=log):
            # for arg in args:
            result = worker_func(arg)
            results.append(result)

    probe_xyz_all = [result[0] for result in results]
    probe_mask_all = [result[1] for result in results]

    # debug
    # return probe_xyz_all, probe_mask_all

    # stack numpy
    probe_xyz = np.stack(probe_xyz_all)
    probe_mask = np.stack(probe_mask_all)

    return probe_xyz, probe_mask


def ndarray_to_nested_list(arr):
    """
    Convert a NumPy array of arbitrary dimensions into a nested list.
    :param arr: A NumPy array.
    :return: A nested list representing the array.
    """
    if arr.ndim == 1:
        return arr.tolist()
    return [ndarray_to_nested_list(sub_arr) for sub_arr in arr]

# Example usage


def qscore_np(
    mmm,
    selection=None,
    n_probes=32,
    shells=np.array(
        [
            0.1,
            0.27272727,
            0.44545455,
            0.61818182,
            0.79090909,
            0.96363636,
            1.13636364,
            1.30909091,
            1.48181818,
            1.65454545,
            1.82727273,
            2.0,
        ]
    ),
    version=2,
    nproc=cpu_count(),
    log=sys.stdout,
):
    """
    Calculate the qscore metric per-atom from an mmtbx map-model-manager, using numpy
    """
    model = mmm.model()
    mm = mmm.map_manager()
    volume = mm.map_data().as_numpy_array()
    radii = shells
    voxel_size = mm.pixel_sizes()

    probe_xyz, probe_mask = radial_shell_mp_np(
        model,
        n_probes=n_probes,
        num_processes=nproc,
        selection=selection,
        version=version,
        radii=radii,
        log=log,
    )

    # after the probe generation, versions 1 and 2 are the same

    # infer params from shape
    n_shells, n_atoms, n_probes, _ = probe_xyz.shape

    # flatten
    probe_xyz_flat = probe_xyz.reshape((n_atoms * n_shells * n_probes, 3))
    probe_mask_flat = probe_mask.reshape(-1)  # (n_shells*n_atoms*n_probes,)

    # select mask=True probes
    masked_probe_xyz_flat = probe_xyz_flat[probe_mask_flat]

    # interpolate
    masked_density = trilinear_interpolation(
        volume, masked_probe_xyz_flat, voxel_size=voxel_size
    )

    # reshape interpolated values to (n_shells,n_atoms, n_probes)

    d_vals = np.zeros((n_shells, n_atoms, n_probes))
    d_vals[probe_mask] = masked_density

    # reshape to (M,N*L) for rowwise correlation

    d_vals_2d = d_vals.transpose(1, 0, 2).reshape(d_vals.shape[1], -1)

    # create the reference data

    M = volume
    maxD = min(M.mean() + M.std() * 10, M.max())
    minD = max(M.mean() - M.std() * 1, M.min())
    A = maxD - minD
    B = minD
    u = 0
    sigma = 0.6
    x = np.array(radii)
    y = A * np.exp(-0.5 * ((x - u) / sigma) ** 2) + B

    # Stack and reshape data for correlation calc

    # stack the reference to shape (n_shells,n_atoms,n_probes)
    g_vals = np.repeat(y[:, None], n_probes, axis=1)
    g_vals = np.expand_dims(g_vals, 1)
    g_vals = np.tile(g_vals, (n_atoms, 1))

    # reshape
    g_vals_2d = g_vals.transpose(1, 0, 2).reshape(g_vals.shape[1], -1)
    d_vals_2d = d_vals.transpose(1, 0, 2).reshape(d_vals.shape[1], -1)
    mask_2d = probe_mask.transpose(1, 0, 2).reshape(probe_mask.shape[1], -1)

    # # CALCULATE Q

    # # numpy
    q = rowwise_corrcoef(g_vals_2d, d_vals_2d, mask=mask_2d)

    # # Log
    # import os

    # print("FINISHING...")
    # print(os.getcwd())
    # print("saving atoms xyz: atom_xyz_np.npy")
    # np.save("atoms_xyz_np.npy",model.get_sites_cart().as_numpy_array())
    # print("saving probe xyz: probe_xyz_np.npy")
    # np.save("probe_xyz_np.npy",probe_xyz)
    # print("saving probe mask: probe_mask_np.npy")
    # np.save("probe_mask_np.npy",probe_mask)
    # print("saving qscore: qscore_np.npy")
    # np.save("qscore_np.npy",q)

    return q


##############################################################################
# Functions that use only flex, no numpy
##############################################################################


def radial_shell_worker_v2_flex(args):
    """
    Calulate qscore for a single radial shell using version 2 (parallel probe allocation) and flex only
    """
    # unpack args
    i, atoms_xyz, n_probes, radius_shell, rtol, tree, selection = args

    # manage selection input
    if selection is None:
        selection = flex.size_t_range(len(atoms_xyz))
    else:
        assert isinstance(selection, flex_bool)

    # do selection
    n_atoms = selection.count(True)
    atoms_xyz_sel = atoms_xyz.select(selection)

    # get probe coordinates
    probe_xyz = sphere_points_flex(atoms_xyz_sel, radius_shell, n_probes)

    # query to find the number of atoms within the clash range of each probe
    counts = query_ball_point_flex(
        tree, atoms_xyz, probe_xyz, r=radius_shell * rtol
    )
    probe_mask = counts == 0
    return probe_xyz, probe_mask


def radial_shell_v2_mp_flex(
    model,
    n_probes=32,
    radii=np.linspace(0.1, 2, 12),
    rtol=0.9,
    num_processes=cpu_count(),
    selection=None,
    version=2,
    log=sys.stdout,
):
    """
    Generate probes for a model file using flex only
    """
    assert version in [1, 2], "Version must be 1 or 2"
    if version == 1:
        assert False
    else:
        worker_func = radial_shell_worker_v2_flex

    # get a "tree", which is just a dictionary of index:local neighbor indices
    tree, _ = query_atom_neighbors(model, radius=3.5)
    atoms_xyz = model.get_sites_cart()

    # i,atoms_xyz,n_probes,radius_shell,rtol, selection= args
    # Create argument tuples for each chunk

    args = [
        (i, atoms_xyz, n_probes, radius_shell, rtol, tree, selection)
        for i, radius_shell in enumerate(radii)
    ]

    # Create a pool of worker processes
    if num_processes > 1:
        with Pool(num_processes) as p:
            # Use the pool to run the trilinear_interpolation_worker function in parallel
            results = list(
                tqdm(p.imap(worker_func, args), total=len(args), file=log)
            )
    else:
        results = []
        for arg in tqdm(args, file=log):
            # for arg in args:
            result = worker_func(arg)
            results.append(result)

    # stack the results from each shell into single arrays
    probe_xyz_all = [result[0] for result in results]
    probe_mask_all = [result[1] for result in results]

    # # debug
    # return probe_xyz_all, probe_mask_all,tree

    n_atoms = probe_xyz_all[0].focus()[0]
    n_shells = len(probe_mask_all)
    out_shape = (n_shells, n_atoms, n_probes, 3)
    out_size = math.prod(out_shape)
    shell_size = math.prod(out_shape[1:])
    out_probes = flex.double(out_size, -1.0)
    out_mask = flex.bool(n_atoms * n_shells * n_probes, False)

    for i, p in enumerate(probe_xyz_all):
        start = i * shell_size
        stop = start + shell_size
        out_probes = out_probes.set_selected(
            flex.size_t_range(start, stop), p.as_1d()
        )
    out_probes.reshape(flex.grid(*out_shape))

    for i, k in enumerate(probe_mask_all):
        start = i * (n_atoms * n_probes)
        stop = start + (n_atoms * n_probes)
        out_mask = out_mask.set_selected(
            flex.size_t_range(start, stop), k.as_1d()
        )
    out_mask.reshape(flex.grid(n_shells, n_atoms, n_probes))

    return out_probes, out_mask


def qscore_flex(
    mmm,
    selection=None,
    n_probes=32,
    shells=[
        0.1,
        0.27272727,
        0.44545455,
        0.61818182,
        0.79090909,
        0.96363636,
        1.13636364,
        1.30909091,
        1.48181818,
        1.65454545,
        1.82727273,
        2.0,
    ],
    version=2,
    nproc=cpu_count(),
):
    """
    Calculate the qscore metric per-atom from an mmtbx map-model-manager, using flex only
    """
    model = mmm.model()
    mm = mmm.map_manager()
    radii = shells
    volume = mm.map_data()
    voxel_size = mm.pixel_sizes()

    probe_xyz, probe_mask = radial_shell_v2_mp_flex(
        model,
        n_probes=n_probes,
        num_processes=nproc,
        selection=selection,
        version=version,
        radii=radii,
    )

    # aliases
    probe_xyz_cctbx = probe_xyz
    probe_mask_cctbx = probe_mask

    # infer params from shape
    n_shells, n_atoms, n_probes, _ = probe_xyz.focus()

    # APPLY MASK BEFORE INTERPOLATION

    probe_mask_cctbx_fullflat = []

    for val in probe_mask_cctbx:
        for _ in range(3):  # since A has an additional dimension of size 3
            probe_mask_cctbx_fullflat.append(val)

    mask = flex.bool(probe_mask_cctbx_fullflat)
    # indices = flex.int([i for i in range(1, keep_mask_cctbx.size() + 1) for _ in range(3)])
    sel = probe_xyz_cctbx.select(mask)
    # sel_indices = indices.select(mask)
    masked_probe_xyz_flat_cctbx = flex.vec3_double(sel)

    # INTERPOLATE

    masked_density_cctbx = mm.density_at_sites_cart(
        masked_probe_xyz_flat_cctbx
    )

    # reshape interpolated values to (n_shells,n_atoms, n_probes)

    probe_mask_cctbx.reshape(flex.grid(n_shells * n_atoms * n_probes))
    d_vals_cctbx = flex.double(probe_mask_cctbx.size(), 0.0)
    d_vals_cctbx = d_vals_cctbx.set_selected(
        probe_mask_cctbx, masked_density_cctbx
    )
    d_vals_cctbx.reshape(flex.grid(n_shells, n_atoms, n_probes))

    # reshape to (M,N*L) for rowwise correlation

    def custom_reshape_indices(flex_array):
        N, M, L = flex_array.focus()
        result = flex.double(flex.grid(M, N * L))

        for i in range(N):
            for j in range(M):
                for k in range(L):
                    # Calculate the original flat index
                    old_index = i * M * L + j * L + k
                    # Calculate the new flat index after transpose and reshape
                    new_index = j * N * L + i * L + k
                    result[new_index] = flex_array[old_index]

        return result

    d_vals_2d_cctbx = custom_reshape_indices(d_vals_cctbx)

    # create the reference data
    M = volume
    M_std = flex_std(M)
    M_mean = flex.mean(M)
    maxD_cctbx = min(M_mean + M_std * 10, flex.max(M))
    minD_cctbx = max(M_mean - M_std * 1, flex.min(M))
    A_cctbx = maxD_cctbx - minD_cctbx
    B_cctbx = minD_cctbx
    u = 0
    sigma = 0.6
    x = flex.double(radii)
    y_cctbx = (
        A_cctbx * flex.exp(-0.5 * ((flex.double(x) - u) / sigma) ** 2)
        + B_cctbx
    )

    # Stack and reshape data for correlation calc

    # 1. Repeat y for n_probes (equivalent to np.repeat)
    g_vals_cctbx = [[val] * n_probes for val in y_cctbx]

    # 2. Add a new dimension (equivalent to np.expand_dims)
    g_vals_expanded = [[item] for item in g_vals_cctbx]

    # 3. Tile for each atom (equivalent to np.tile)
    g_vals_tiled = []
    for item in g_vals_expanded:
        g_vals_tiled.append(item * n_atoms)

    g_vals_cctbx = flex.double(np.array(g_vals_tiled))

    # # CALCULATE Q

    d_vals_cctbx = d_vals_cctbx.as_1d()
    g_vals_cctbx = g_vals_cctbx.as_1d()
    probe_mask_cctbx_double = probe_mask_cctbx.as_1d().as_double()
    q_cctbx = []
    for atomi in range(n_atoms):
        inds = nd_to_1d_indices(
            (None, atomi, None), (n_shells, n_atoms, n_probes)
        )
        # inds = optimized_nd_to_1d_indices(atomi,(n_shells,n_atoms,n_probes))
        inds = flex.uint32(inds)
        d_row = d_vals_cctbx.select(inds)
        g_row = g_vals_cctbx.select(inds)
        mask = probe_mask_cctbx.select(inds)

        d = d_row.select(mask)
        g = g_row.select(mask)
        qval = flex.linear_correlation(d, g).coefficient()
        q_cctbx.append(qval)

    q = flex.double(q_cctbx)
    return q


# qscore utils


class DummyTQDM:
    """A 'dummy' object that can be used for compatibility if tqdm is not installed"""

    def __init__(self, iterable=None, *args, **kwargs):
        self.iterable = iterable
        self.args = args
        self.kwargs = kwargs

    def __iter__(self):
        return iter(self.iterable)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


def trilinear_interpolation(voxel_grid, coords, voxel_size=None, offset=None):
    """Numpy trilinear interpolation"""
    assert voxel_size is not None, "Provide voxel size as an array or single value"

    # Apply offset if provided
    if offset is not None:
        coords = coords - offset

    # Transform coordinates to voxel grid index space
    index_coords = coords / voxel_size

    # Split the index_coords array into three arrays: x, y, and z
    x, y, z = index_coords.T

    # Truncate to integer values
    x0, y0, z0 = np.floor([x, y, z]).astype(int)
    x1, y1, z1 = np.ceil([x, y, z]).astype(int)

    # Ensure indices are within grid boundaries
    x0, y0, z0 = np.clip([x0, y0, z0], 0, voxel_grid.shape[0]-1)
    x1, y1, z1 = np.clip([x1, y1, z1], 0, voxel_grid.shape[0]-1)

    # Compute weights
    xd, yd, zd = [arr - arr.astype(int) for arr in [x, y, z]]

    # Interpolate along x
    c00 = voxel_grid[x0, y0, z0]*(1-xd) + voxel_grid[x1, y0, z0]*xd
    c01 = voxel_grid[x0, y0, z1]*(1-xd) + voxel_grid[x1, y0, z1]*xd
    c10 = voxel_grid[x0, y1, z0]*(1-xd) + voxel_grid[x1, y1, z0]*xd
    c11 = voxel_grid[x0, y1, z1]*(1-xd) + voxel_grid[x1, y1, z1]*xd

    # Interpolate along y
    c0 = c00*(1-yd) + c10*yd
    c1 = c01*(1-yd) + c11*yd

    # Interpolate along z
    c = c0*(1-zd) + c1*zd

    return c


def rowwise_corrcoef(A, B, mask=None):
    """Numpy masked array rowwise correlation coefficient"""
    assert A.shape == B.shape, f"A and B must have the same shape, got: {A.shape} and {B.shape}"

    if mask is not None:
        assert mask.shape == A.shape, "mask must have the same shape as A and B"
        A = ma.masked_array(A, mask=np.logical_not(mask))
        B = ma.masked_array(B, mask=np.logical_not(mask))

    # Calculate means
    A_mean = ma.mean(A, axis=1, keepdims=True)
    B_mean = ma.mean(B, axis=1, keepdims=True)

    # Subtract means
    A_centered = A - A_mean
    B_centered = B - B_mean

    # Calculate sum of products
    sumprod = ma.sum(A_centered * B_centered, axis=1)

    # Calculate square roots of the sum of squares
    sqrt_sos_A = ma.sqrt(ma.sum(A_centered**2, axis=1))
    sqrt_sos_B = ma.sqrt(ma.sum(B_centered**2, axis=1))

    # Return correlation coefficients
    cc = sumprod / (sqrt_sos_A * sqrt_sos_B)
    return cc.data


def sphere_points_np(ctr, rad, N, mode='SpherePts'):
    """
    Points on a sphere given centers, radius, number N

    TODO: Mode is confusing, it is not clear why there is a difference.
    mode='SpherePts' an attempt to literally copy the SpherePts function from mapq
    mode='original' a version that gives the same results as the QscorePt3 function from mapq
    """
    h = -1.0 + (2.0 * np.arange(N) / float(N-1))
    phis = np.arccos(h)

    if mode == 'original':
        thetas = np.zeros(N)
        a = (3.6 / np.sqrt(N * (1.0 - h[1:-1]**2)))
        thetas[1:-1] = a
        thetas = np.cumsum(thetas)
    elif mode == 'SpherePts':
        thetas = np.zeros(N)
        for k in range(1, N):
            if k == 1 or k == N - 1:
                thetas[k] = 0
            else:
                thetas[k] = (thetas[k-1] + 3.6 /
                             np.sqrt(N * (1 - h[k]**2))) % (2 * np.pi)
        thetas = thetas.cumsum()

    x = np.sin(phis) * np.cos(thetas)
    y = np.sin(phis) * np.sin(thetas)
    z = np.cos(phis)

    points = rad * np.stack([x, y, z], axis=-1)

    # Adjusting for multiple centers
    if ctr.ndim == 1:
        # Single center case
        points = points + ctr
    else:
        # Multiple centers case
        points = points.reshape(-1, 1, 3) + ctr.reshape(1, -1, 3)

    return points


def cdist_flex(A, B):
    """A flex implementation of the cdist function"""

    def indices_2d_flex(dimensions):
        N = len(dimensions)
        if N != 2:
            raise ValueError("Only 2D is supported for this implementation.")

        # Create the row indices
        row_idx = flex.size_t(chain.from_iterable(
            [[i] * dimensions[1] for i in range(dimensions[0])]))

        # Create the column indices
        col_idx = flex.size_t(chain.from_iterable(
            [list(range(dimensions[1])) for _ in range(dimensions[0])]))

        return row_idx, col_idx

    i_idxs, j_idxs = indices_2d_flex((A.focus()[0], B.focus()[0]))

    r = i_idxs
    xi = i_idxs*3
    yi = i_idxs*3 + 1
    zi = i_idxs*3 + 2

    xa = A.select(xi)
    ya = A.select(yi)
    za = A.select(zi)

    xj = j_idxs*3
    yj = j_idxs*3 + 1
    zj = j_idxs*3 + 2

    xb = B.select(xj)
    yb = B.select(yj)
    zb = B.select(zj)

    d = ((xb - xa)**2 + (yb - ya)**2 + (zb - za)**2)**0.5
    d.reshape(flex.grid((A.focus()[0], B.focus()[0])))

    return d


def query_atom_neighbors(model, radius=3.5, include_self=True, only_unit=True):
    """Perform radial nearest neighbor searches using cctbx tools, for atom coordinates in a model"""
    crystal_symmetry = model.crystal_symmetry()
    hierarchy = model.get_hierarchy()
    sites_cart = hierarchy.atoms().extract_xyz()
    sst = crystal_symmetry.special_position_settings().site_symmetry_table(
        sites_cart=sites_cart)
    conn_asu_mappings = crystal_symmetry.special_position_settings().\
        asu_mappings(buffer_thickness=5)
    conn_asu_mappings.process_sites_cart(
        original_sites=sites_cart,
        site_symmetry_table=sst)
    conn_pair_asu_table = cctbx.crystal.pair_asu_table(
        asu_mappings=conn_asu_mappings)
    conn_pair_asu_table.add_all_pairs(distance_cutoff=radius)
    pair_generator = cctbx.crystal.neighbors_fast_pair_generator(
        conn_asu_mappings,
        distance_cutoff=radius)
    fm = crystal_symmetry.unit_cell().fractionalization_matrix()
    om = crystal_symmetry.unit_cell().orthogonalization_matrix()

    pairs = list(pair_generator)
    inds = defaultdict(list)
    dists = defaultdict(list)

    for pair in pairs:
        i, j = pair.i_seq, pair.j_seq
        rt_mx_i = conn_asu_mappings.get_rt_mx_i(pair)
        rt_mx_j = conn_asu_mappings.get_rt_mx_j(pair)
        rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)

        if (only_unit and rt_mx_ji.is_unit_mx()) or (not only_unit):
            d = round(math.sqrt(pair.dist_sq), 6)
            inds[i].append(j)
            dists[i].append(d)

            # add reverse
            inds[j].append(i)
            dists[j].append(d)
            # print(pair.i_seq, pair.j_seq, rt_mx_ji, math.sqrt(pair.dist_sq), de)

    # add self
    if include_self:
        for key, value in list(inds.items()):
            dval = dists[key]
            dists[key] = dval+[0.0]
            inds[key] = value+[key]

    # sort
    for key, value in list(inds.items()):
        dval = dists[key]
        # sort
        sorted_pairs = sorted(set(list(zip(value, dval))))
        value_sorted, dval_sorted = zip(*sorted_pairs)
        inds[key] = value_sorted
        dists[key] = dval_sorted

    return inds, dists


def query_ball_point_flex(tree, tree_xyz, query_xyz, r=None):
    """
    Imitate the api of the scipy.spatial query_ball_point function, but using only flex arrays.
    Note: This just copies the api, it does not actually use a tree structure, so is much slower.
    """
    assert r is not None, "provide radius"
    n_atoms, n_probes, _ = query_xyz.focus()
    counts = []

    for atom_i in range(n_atoms):
        probe_range = (n_probes * atom_i * 3, n_probes * (atom_i+1) * 3)
        atom_probes_xyz = query_xyz.select(flex.size_t_range(*probe_range))
        atom_probes_xyz.reshape(flex.grid(n_probes, 3))
        nbrs = tree[atom_i]
        n_nbrs = len(nbrs)
        nbrs_xyz = tree_xyz.select(flex.size_t(nbrs)).as_1d().as_double()
        nbrs_xyz.reshape(flex.grid(len(nbrs), 3))
        d = cdist_flex(nbrs_xyz, atom_probes_xyz)
        sel = d < r
        count = []
        for nbr_i in range(n_probes):
            nbr_range = (slice(0, n_nbrs), slice(nbr_i, nbr_i+1))
            count_nbr = sel[nbr_range].count(True)
            count.append(count_nbr)

        counts.append(count)

    counts = flex_from_list(counts)
    return counts


# flex utils


def flex_from_list(lst, signed_int=False):
    """Generate a flex array from a list, try to infer type"""
    flat_list, shape = flatten_and_shape(lst)
    dtype = get_dtype_of_list(flat_list)
    type_mapper = {int: flex.size_t,
                   float: flex.double,
                   bool: flex.bool}
    if signed_int:
        type_mapper[int] = flex.int16

    # make flex array
    assert dtype in type_mapper, f"Unrecognized type: {dtype}"
    flex_func = type_mapper[dtype]
    flex_array = flex_func(flat_list)
    if len(shape) > 1:
        flex_array.reshape(flex.grid(*shape))
    return flex_array


def flatten_and_shape(lst):
    """Flatten a nested list and return its shape."""
    def helper(l):
        if not isinstance(l, list):
            return [l], ()
        flat = []
        shapes = []
        for item in l:
            f, s = helper(item)
            flat.extend(f)
            shapes.append(s)
        if len(set(shapes)) != 1:
            raise ValueError("Ragged nested list detected.")
        return flat, (len(l),) + shapes[0]

    flattened, shape = helper(lst)
    return flattened, shape


def get_dtype_of_list(lst):
    dtypes = {type(item) for item in lst}

    if len(dtypes) > 1:
        raise ValueError("Multiple data types detected.")
    elif len(dtypes) == 0:
        raise ValueError("Empty list provided.")
    else:
        return dtypes.pop()


def nd_to_1d_indices(indices, shape):
    """Generate the 1d indices given nd indices and an array shape"""
    # Normalize indices to always use slice objects
    normalized_indices = []
    for dim, idx in enumerate(indices):
        if idx is None:
            normalized_indices.append(slice(0, shape[dim]))
        else:
            normalized_indices.append(idx)

    # If any index is a slice, recursively call function for each value in slice
    for dim, (i, s) in enumerate(zip(normalized_indices, shape)):
        if isinstance(i, slice):
            result_indices = []
            start, stop, step = i.indices(s)
            for j in range(start, stop, step):
                new_indices = list(normalized_indices)
                new_indices[dim] = j
                result_indices.extend(nd_to_1d_indices(new_indices, shape))
            return result_indices

    # If no slices, calculate single 1D index
    index = 0
    stride = 1
    for i, dim in reversed(list(zip(normalized_indices, shape))):
        index += i * stride
        stride *= dim
    return [index]


def optimized_nd_to_1d_indices(i, shape):
    """Similar to above, but hardcoded to select a single index on dimension 1"""
    # For fixed input of (None, i, None), we directly compute based on given structure
    result_indices = []

    # Pre-compute for 1st dimension which is always a slice
    start1, stop1 = 0, shape[0]

    # Pre-compute for 3rd dimension which is always a slice
    start3, stop3 = 0, shape[2]
    stride3 = 1

    # Directly compute for 2nd dimension which is variable
    stride2 = shape[2]
    index2 = i * stride2 * shape[0]

    for val1 in range(start1, stop1):
        for val3 in range(start3, stop3):
            result_indices.append(val1 * stride2 + index2 + val3 * stride3)

    return result_indices


def flex_std(flex_array):
    """Standard deviation"""
    n = flex_array.size()
    if n <= 1:
        raise ValueError("Sample size must be greater than 1")

    # Compute the mean
    mean_value = flex.mean(flex_array)

    # Compute the sum of squared deviations
    squared_deviations = (flex_array - mean_value) ** 2
    sum_squared_deviations = flex.sum(squared_deviations)

    # Compute the standard deviation
    std_dev = (sum_squared_deviations / (n - 1)) ** 0.5
    return std_dev
