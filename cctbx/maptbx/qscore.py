"""
This code provides methods to calculate the qscore metric for map-model validation,
as developed by Pintile et al.

As with the original implementation, multiple calculation options are provided.
The fastest (default) is to pass an mmtbx map-model-manager (mmm)
to the the function qscore_np, qscore_np(mmm)
"""

from __future__ import division
<<<<<<< Updated upstream
import math
import sys
from collections import defaultdict
from multiprocessing import Pool, cpu_count
=======
from libtbx.utils import null_out
from multiprocessing import Pool
>>>>>>> Stashed changes
from itertools import chain
from collections import defaultdict

import numpy as np
import numpy.ma as ma
from scipy.spatial import KDTree
import pandas as pd

from cctbx.array_family import flex

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

<<<<<<< Updated upstream
try:
    from tqdm import tqdm
except ImportError:
    from .qscore_utils import DummyTQDM as tqdm


def radial_shell_worker_v1_np(args):
=======
master_phil_str = """
  qscore
  {

    nproc = 1
        .type = int
        .help = Number of processors to use
        .short_caption = Number of processors to use
        .expert_level = 1
    n_probes_target = 8
        .type = int
        .help = Number of radial probes to use
        .short_caption = Number of radial probes to use
        .expert_level = 1
    n_probes_max = 16
        .type = int
        .help = Max number of radial probes to use
        .short_caption = Number of radial probes to use
        .expert_level = 1
    n_probes_min = 4
        .type = int
        .help = Min number of radial probes to use
        .short_caption = Number of radial probes to use
        .expert_level = 1
    selection_str = None
      .type = str
      .help = Only calculate q for atoms within this selection
      .short_caption = Only calculate q for atoms within this selection
      .expert_level = 1

    shell_radius_start = 0.1
      .type = float
      .help = Start testing density at this radius from atom
      .short_caption = Start testing density at this radius from atom
      .expert_level = 1

    shell_radius_stop = 2
      .type = float
      .help = Stop testing density at this radius from atom
      .short_caption = Stop testing density at this radius from atom
      .expert_level = 1

    shell_radius_num = 20
      .type = int
      .help = The number of radial shells
      .short_caption = The number of radial shells\
                         (includes start/stop, so minimum 2)
      .expert_level = 1

    shells = None
      .type = float
      .multiple = True
      .help = Explicitly provide radial shells

    rtol = 0.9
      .type = float
      .help = Mapq rtol value, the "real" shell radii are r*rtol

    probe_allocation_method = precalculate
      .type = str
      .help = The method used to allocate radial probes
      .short_caption = Either 'progressive' or 'precalculate'. \
                       Progressive is the original method \
                       where probes are proposed and rejected iteratively. \
                       Precalculate is where probes are pre-allocated and \
                       rejected once. Parallelization is done by radial shell.\
                       Precalculate is faster but with slightly varied results.

    backend = numpy
      .type = str
      .help = Either 'numpy' or 'flex'
      .short_caption = Determines which backend mdoe is used
      .expert_level = 1

    debug = False
      .type = bool
      .help = Return much more debug information
      .short_caption = Returns a dictionary with additional debug information
      .expert_level = 1
  }

  """



################################################################################
#### Probe generation functions
################################################################################

# Original mapq version
def SpherePtsVectorized ( ctr, rad, N ) :
  """
  Function for generating points on a sphere. For testing only.
  It retains the original mapq code pattern
  """
  thetas, phis = [], []
  from math import acos, sin, cos, sqrt, pi
  for k in range ( 1, N+1 ) :
    h = -1.0 + ( 2.0*float(k-1)/float(N-1) )
    phis.append ( acos(h) )
    thetas.append ( 0 if k == 1 or k == N else
                    (thetas[k-2] + 3.6/sqrt(N*(1.0-h**2.0))) % (2*pi) )

  pts = [None] * N
  for i, theta, phi in zip ( range(N), thetas, phis ):
    v = np.array([ sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)])

    pt = ctr + v * rad
    pts[i] = pt
  pts = np.array(pts)
  pts = pts.swapaxes(0,1)
  return pts


# Fast numpy version
def generate_probes_np(atoms_xyz, rad, n_probes):
  """
  Generate probes using the same methodology as Pintile mapq, but vectorized

  atoms_xyz: np array of shape (n_atoms,3)
  rad: the radius at which to place the probes
  N: the number of probes per atom

  Returns:
    probes (np.ndarray): shape (n_atoms,n_probes,3)
  """
  assert atoms_xyz.ndim == 2 and atoms_xyz.shape[-1]==3, (
    "Provide coordinates in shape (n_atoms,3)")

  N = n_probes
  h = -1.0 + (2.0 * np.arange(N) / float(N-1))
  phis = np.arccos(h)

  thetas = np.zeros(N)
  a = (3.6 / np.sqrt(N * (1.0 - h[1:-1]**2)))
  thetas[1:-1] = a
  thetas = np.cumsum(thetas)


  x = np.sin(phis) * np.cos(thetas)
  y = np.sin(phis) * np.sin(thetas)
  z = np.cos(phis)

  probes = rad * np.stack([x, y, z], axis=-1)

  # Adjusting location of generated points relative to point ctr
  probes = probes.reshape(-1, 1, 3) + atoms_xyz.reshape(1, -1, 3)

  # reshape (n_atoms,n_probes,3)
  probes = probes.swapaxes(0,1)
  return probes


# Generate Points with flex

def cumsum_flex(arr):
  """
  Return an array that is the cumulative sum of arr
  Analogous to np.cumsum
  """
  result = flex.double(len(arr))
  running_sum = 0.0
  for i, x in enumerate(arr):
      running_sum += x
      result[i] = running_sum
  return result

def broadcast_add_vec3(ctr, points):
  """
  Broadcast add two flex.vec3_double arrays.

  Params:
    ctr (flex.vec3_double): the 'center' coordinates
    points (flex.vec3_double): the points that will be added to each ctr

  Returns:
    result (flex.vec3_double): array of shape (N*M,3), 1 point for each center
  """
  N = points.size()
  M = ctr.size()
  extended_ctr = flex.vec3_double()
  extended_points = flex.vec3_double()

  # Extend ctr and points
  for point in ctr:
    extended_ctr.extend(flex.vec3_double([point] * N))
  for _ in range(M):
    extended_points.extend(points)

  # Perform addition
  result = flex.vec3_double(M * N)
  space =flex.size_t_range(M*N)
  for i in space:
    pt = extended_ctr[i:i+1] + extended_points[i:i+1]
    result=result.set_selected(space[i:i+1],pt)
  return result


def generate_probes_flex(ctr, rad, N):
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
    #
    # manage selection input
    if selection is None:
        selection = np.arange(len(atoms_xyz))
=======
    # Generating spherical coordinates
    x = flex.sin(phis) * flex.cos(thetas)
    y = flex.sin(phis) * flex.sin(thetas)
    z = flex.cos(phis)
    points = rad * flex.vec3_double(x, y, z)
    # Adding the generated points to center points
    probes = broadcast_add_vec3(ctr, points)
    return probes


################################################################################
#### Progressive mode functions
################################################################################

def get_probe_mask(
      atom_tree,
      probes_xyz,
      r=None,
      expected=None,
      log=null_out(),
      ):
  """
  atoms_xyz shape  (n_atoms,3)
  probes_xyz shape (n_atoms,n_probes,3)

  If expected is None, infer atom indices from probes_xyz
  Else expected should be a single value, or have shape  (n_atoms,n_probes)

  """

  assert r is not None, "Provide a radius"
  assert probes_xyz.ndim ==3 and probes_xyz.shape[-1] == 3,(
    "Provide probes_xyz as shape: (n_atoms,n_probes,3)")

  n_atoms_probe,n_probes,_ = probes_xyz.shape
  dim = probes_xyz.shape[-1] # 3 for cartesian coords


  # reshaped_probes shape (n_atoms*n_probes,3)
  reshaped_probes = probes_xyz.reshape(-1, 3)
  atom_indices = np.tile(np.arange(n_atoms_probe), (probes_xyz.shape[1], 1)).T

  if not expected:
    atom_indices = np.tile(np.arange(n_atoms_probe), (probes_xyz.shape[1], 1)).T
  else:
    atom_indices = np.full(probes_xyz.shape,expected)

  associated_indices = atom_indices.reshape(-1)


  # query
  # Check if any other tree points are within r of each query point
  query_points = reshaped_probes
  other_points_within_r = []
  for i, (query_point,idx) in enumerate(zip(query_points,associated_indices)):

    indices_within_r = atom_tree.query_ball_point(query_point, r)

    # Exclude the associated point
    associated_index = associated_indices[i]
    other_indices = []
    for  idx in indices_within_r:
      if idx != associated_index:
        other_indices.append(idx)
      if len(indices_within_r)==0:
        other_indices.append(-1)


    print(other_indices,file=log)

    other_points_within_r.append(other_indices)

  # true are points that don't get rejected
  num_nbrs_other = np.array(
     [len(inds) for i,inds in enumerate(other_points_within_r)])

  num_nbrs_other = num_nbrs_other.reshape((n_atoms_probe,n_probes))
  mask = num_nbrs_other==0

  return mask


def shell_probes_progressive(
      atoms_xyz=None,   # A numpy array of shape (N,3)
      atoms_tree=None,  # An atom_xyz scipy kdtree
      selection_bool=None,# Boolean atom selection
      n_probes_target=8,# The desired number of probes per shell
      n_probes_max=16,  # The maximum number of probes allowed
      n_probes_min=4,
      RAD=1.5,          # The nominal radius of this shell
      rtol=0.9,         # Multiplied with RAD to get actual radius
      log = null_out(),
      ):
  """
  Generate probes progressively for a single shell (radius)
  """

  # Do input validation
  if not atoms_tree:
    assert atoms_tree is None, (
      "If not providing an atom tree, \
        provide a 2d atom coordinate array to build tree")

    atoms_tree = KDTree(atoms_xyz)

  # Manage log
  if log is None:
    log = null_out()

  # manage selection input
  if selection_bool is None:
    selection_bool = np.full(len(atoms_xyz),True)

  # do selection
  atoms_xyz_sel = atoms_xyz[selection_bool]
  n_atoms = atoms_xyz_sel.shape[0]

  all_pts = []  # list of probe arrays for each atom
  for atom_i in range(n_atoms):
    coord = atoms_xyz_sel[atom_i:atom_i+1]
    outRAD = RAD * rtol


    print(coord,file=log)
    pts = []
    i_log = []
    # try to get at least numPts] points at [RAD] distance
    # from the atom, that are not closer to other atoms
    N_i = 50

    # If we find the necessary number of probes in the first iteration,
    #   then i will never go to 1
    for i in range(0, N_i):
      rejections = 0



      # progressively more points are grabbed  with each failed iter
      n_pts_to_grab = (n_probes_target + i * 2)

      # get the points in shape (n_atoms,n_pts_to_grab,3)
      outPts = generate_probes_np(coord, RAD, n_pts_to_grab)

      # initialize points to keep
      at_pts, at_pts_i = [None] * outPts.shape[1], 0

      # mask for outPts, are they are closest to the expected atom
      # mask shape (n_atoms,n_pts_to_grab)
      # NOTE: n_atoms != len(outPts)

      # will get mask of shape (n_atoms,n_probes)
      mask = get_probe_mask(atoms_tree,outPts,r=outRAD,expected=atom_i,log=log)

      # identify which ones to keep, progressively grow pts list
      for pt_i, pt in enumerate(outPts[0]):
        keep = mask[0,pt_i] # only one atom TODO: vectorize atoms
        if keep:
          at_pts[at_pts_i] = pt
          at_pts_i += 1
        else:
          #print("REJECTING...",pt,file=log)
          rejections+=1
          pass

      # check if we have enough points to break the search loop
      if ( at_pts_i >= n_probes_target):
        pts.extend(at_pts[0:at_pts_i])
        pts = pts + [np.array([np.nan,np.nan,np.nan])]*(n_probes_max-len(pts))
        #print(pts)
        break

      i_log.append(i)
      if i>=N_i:
          assert False, "Too many iterations to get probes"
      if i>0:
          print("Going another round..",file=log)
      # End sampling iteration


    #Finish working on a single atom
    pts = np.array(pts)
    if pts.shape == (0,): # all probes clashed
      pts = np.full((n_probes_max,3),np.nan)
>>>>>>> Stashed changes
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


<<<<<<< Updated upstream
def radial_shell_worker_v2_np(args):
=======
################################################################################
#### numpy/scipy-based functions (precalculate mode)
################################################################################

def shell_probes_precalculate(
      atoms_xyz=None,   # A numpy array of shape (N,3)
      atoms_tree=None,  # An atom_xyz scipy kdtree
      selection_bool=None,# Boolean atom selection
      n_probes_target=8,# The desired number of probes per shell
      n_probes_max=16,  # The maximum number of probes allowed
      n_probes_min=4,   # The min number of probes allowed without error
      RAD=1.5,          # The nominal radius of this shell
      rtol=0.9,         # Multiplied with RAD to get actual radius
      log = null_out(),
      strict = False,
      ):
  """
  Generate probes by precalculating for a single shell (radius)
  """

  # Do input validation
  if not atoms_tree:
    assert atoms_tree is None, ("If not providing an atom tree,\
      provide a 2d atom coordinate array to build tree")

    # make atom kdtree
    atoms_tree = KDTree(atoms_xyz)

  # Manage log
  if log is None:
    log = null_out()

  # manage selection input
  if selection_bool is None:
    selection_bool = np.full(len(atoms_xyz),True)


  # do selection
  atoms_xyz_sel = atoms_xyz[selection_bool]

  # get probe coordinates
  probe_xyz = generate_probes_np(atoms_xyz_sel, RAD, n_probes_max)
  n_atoms, n_probes, _ = probe_xyz.shape
  probe_xyz_flat = probe_xyz.reshape(-1,3)

  # modify "real" radius as in mapq
  outRAD = RAD*rtol

  # query kdtree to get neighbors and their distances
  dists, atom_indices = atoms_tree.query(probe_xyz_flat, k=2)
  dists = dists.reshape((n_atoms,n_probes,2))
  atom_indices = atom_indices.reshape((n_atoms,n_probes,2))

  # Build an index array that would be expected if each probe is near "its" atom
  row_indices = np.arange(n_atoms)[:, np.newaxis]

  # Mask for whether each probe's nearest atom is the one expected
  expected_atom_mask = atom_indices[:,:,0]==row_indices

  # A second mask to determine if the second nearest neighbor should be rejected
  #  (whether the second nearest neighbor is within the rejection radius)
  within_r_mask = dists[:,:,1]<outRAD #

  # Combine masks
  probe_mask = expected_atom_mask & ~within_r_mask

  # Debug/Validation on number of probes per atom
  n_probes_per_atom = probe_mask.sum(axis=1)
  insufficient_probes = np.where(n_probes_per_atom<n_probes_target)[0]
  problematic_probes = np.where(n_probes_per_atom<n_probes_min)[0]
  if strict:
    assert n_probes_per_atom.min() >= n_probes_min, (
      f"Some atoms have less than {n_probes_min} probes. \
          ({len(problematic_probes)}). Consider raising n_probes_max")

  return probe_xyz, probe_mask


def get_probes(
    atoms_xyz=None,
    sites_cart = None,
    atoms_tree = None,
    params=None,
    selection_bool_np = None,
    selection_bool_flex=None,
    worker_func=None,
    log=None):
>>>>>>> Stashed changes
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


<<<<<<< Updated upstream
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
=======
def calc_qscore(mmm,params,log=null_out(),debug=False):
  """
  Calculate qscore from map model manager
  """
  model = mmm.model()
  # never do hydrogen
  model = model.select(model.selection("not element H"))
  mm = mmm.map_manager()
>>>>>>> Stashed changes


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

<<<<<<< Updated upstream
    # select mask=True probes
    masked_probe_xyz_flat = probe_xyz_flat[probe_mask_flat]

    # interpolate
    masked_density = trilinear_interpolation(
        volume, masked_probe_xyz_flat, voxel_size=voxel_size
=======
  # interpolate
  volume = mm.map_data().as_numpy_array()
  voxel_size = mm.pixel_sizes()
  # masked_density = trilinear_interpolation(
  #    volume, masked_probe_xyz_flat, voxel_size=voxel_size)
  masked_density = mm.density_at_sites_cart(
     flex.vec3_double(masked_probe_xyz_flat)).as_numpy_array()

  d_vals = np.full((n_shells, n_atoms, n_probes),np.nan)
  d_vals[probe_mask] = masked_density

  # g vals
  # create the reference data
  radii = params.shells
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

  # set masked area to nan
  g_vals[~probe_mask] = np.nan

  # reshape
  g_vals_2d = g_vals.transpose(1, 0, 2).reshape(g_vals.shape[1], -1)
  d_vals_2d = d_vals.transpose(1, 0, 2).reshape(d_vals.shape[1], -1)
  mask_2d = probe_mask.transpose(1, 0, 2).reshape(probe_mask.shape[1], -1)

  # CALCULATE Q
  q = rowwise_corrcoef(g_vals_2d, d_vals_2d, mask=mask_2d)

  # round sensibly
  q = np.around(q,4)

  # aggregate per residue
  qscore_df = aggregate_qscore_per_residue(model,q,window=3)
  q = flex.double(q)
  qscore_per_residue = flex.double(qscore_df["Q-scorePerResidue"].values)

  # Output
  if debug or params.debug:
    # Collect debug data
    result = {
      "atom_xyz":atom_xyz,
      "probe_xyz":probe_xyz,
      "probe_mask":probe_mask,
      "d_vals":d_vals,
      "g_vals":g_vals,
      "qscore_per_atom":q,
      "qscore_per_residue":qscore_per_residue,
      "qscore_dataframe":qscore_df,
    }
  else:
    result = {
    "qscore_per_atom":q,
    "qscore_per_residue":qscore_per_residue,
    "qscore_dataframe":qscore_df
    }
  return result

################################################################################
#### CCTBX flex-based functions (precalculate mode)
################################################################################
def shell_probes_precalculate_flex(
      atoms_xyz=None,   # sites_cart. Flex vec3_double
      atoms_tree=None,  # A KDTree
      selection_bool=None, # Boolean atom selection
      n_probes_target=8,# The desired number of probes per shell
      n_probes_max=16,  # The maximum number of probes allowed
      n_probes_min=4,   # The min number of probes allowed without error
      RAD=1.5,          # The nominal radius of this shell
      rtol=0.9,         # Multiplied with RAD to get actual radius
      log = null_out(),
      strict = False,
      ):
  """
  Generate probes by precalculating for a single shell (radius)
  """

  # Do input validation
  if not atoms_tree:
    assert atoms_tree is None, ("If not providing an atom tree,\
      provide a 2d atom coordinate array to build tree")

    # make atom kdtree
    atoms_tree = KDTreeFlex(atoms_xyz)

  # Manage log
  if log is None:
      log = null_out()

  # manage selection input

  if selection_bool is None:
      selection_bool = flex.bool(len(atoms_xyz),True)


  # do selection
  sites_sel = atoms_xyz.select(selection_bool)
  n_atoms = len(sites_sel)

  # get probe coordinates
  probe_sites = generate_probes_flex(sites_sel, RAD, n_probes_max)


  # modify "real" radius as in mapq
  outRAD = RAD*rtol

  # query kdtree to get neighbors and their distances
  dists,atom_indices = atoms_tree.query_nearest_neighbors(probe_sites,k=2)
  atom_indices_flat = flex_from_list(atom_indices)

  # Perform equivalent to atom_indices[:,:,0] if (n_atoms,n_probes,k)
  dim0_indices = flex.size_t_range(0, n_atoms*n_probes_max*2, 2)
  atom_indices_flat = atom_indices_flat.select(dim0_indices)

  # Build an index array that would be expected if each probe is near "its" atom
  row_indices_flat = flex.size_t([
     i for i in range(n_atoms) for _ in range(n_probes_max)])


  # Mask for whether each probe's nearest atom is the one expected
  expected_mask = row_indices_flat == atom_indices_flat


  # A second mask to determine if the second nearest neighbor should be rejected
  #  (whether the second nearest neighbor is within the rejection radius)
  dists = flex_from_list(dists)
  # perform equivalent selcetion to dists[:,:,1] if (n_atoms,n_probes,k)
  dist_dim1_sel = flex.size_t([i * 2 + 1 for i in range(n_atoms * n_probes_max)])
  dists_dim1 = dists.select(dist_dim1_sel)
  within_r_mask = dists_dim1<outRAD

  # Combine masks
  probe_mask_flat = expected_mask & (~within_r_mask)

  # Reshape (this only stores the anticipated shape)
  probe_xyz = probe_sites
  probe_xyz.reshape(flex.grid((n_atoms,n_probes_max))) # vec3
  probe_mask = probe_mask_flat
  probe_mask.reshape(flex.grid((n_atoms,n_probes_max)))

  return probe_xyz, probe_mask

def calc_qscore_flex(mmm,params,log=null_out(),debug=False):
  """
  Calculate qscore from map model manager
  """
  model = mmm.model()
  # never do hydrogen
  model = model.select(model.selection("not element H"))
  mm = mmm.map_manager()


  # Get atoms
  atoms_xyz = model.get_sites_cart()

  # do selection
  if params.selection_str != None:
    selection_bool = mmm.model().selection(params.selection_str) # boolean
    if selection_bool.count(True) ==0:
      print("Finished... nothing selected")
      return {"qscore_per_atom":None}
  else:
    selection_bool = flex.bool(len(atoms_xyz),True)


  # Get probes and probe mask (probes to reject)
  probe_xyzs,probe_masks = get_probes(
    atoms_xyz=atoms_xyz,
    atoms_tree = None,
    params=params,
    selection_bool_flex = selection_bool,
    worker_func=shell_probes_precalculate_flex,
    log = log,
>>>>>>> Stashed changes
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


<<<<<<< Updated upstream
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
=======
  # for each "shell row"
  for shell_idx,(probe_xyz,probe_mask) in enumerate(zip(probe_xyzs,probe_masks)):
    probe_mask.reshape(flex.grid(n_atoms*params.n_probes_max))
    probe_xyz_sel = probe_xyz.select(probe_mask)

    # get d_vals for a "shell row"
    d_vals = mm.density_at_sites_cart(probe_xyz_sel)

    # get reference for the "shell row"
    g_vals = flex.double(n_atoms*params.n_probes_max,y_cctbx[shell_idx])
    g_vals = g_vals.select(probe_mask)

    # calc flat indices for the shell
    start_idx = shell_idx * n_atoms * n_probes
    stop_idx = start_idx + n_atoms * n_probes
    shell_space = flex.size_t_range(start_idx,stop_idx)

    # apply shell mask
    masked_shell_space = shell_space.select(probe_mask)

    # put d,g, mask into results
    full_flat_d.set_selected(masked_shell_space,d_vals)
    full_flat_g.set_selected(masked_shell_space,g_vals)
    full_flat_mask  = full_flat_mask.set_selected(shell_space,probe_mask)

  # calculate q
  def calculate_1d_indices_for_atom(n_shells, n_atoms, n_probes, atom_idx):
    indices_1d = []
    for shell_idx in range(n_shells):
      for probe_idx in range(n_probes):
        index_1d = shell_idx * (n_atoms * n_probes) + atom_idx * n_probes + probe_idx
        indices_1d.append(index_1d)
    return indices_1d

  qscore_per_atom = []
  mask_check_1d = []
  for atomi in range(n_atoms):
    inds = calculate_1d_indices_for_atom(n_shells,n_atoms,n_probes,atomi)
    inds = flex.size_t(inds)

    # select an "atom row" for d,g,mask values
    d_row = full_flat_d.select(inds)
    g_row = full_flat_g.select(inds)
    mask = full_flat_mask.select(inds)

    # subset the d,g rows using the mask
    d = d_row.select(mask)
    g = g_row.select(mask)

    # calculate correlation between masked d and g
    qval = flex.linear_correlation(d, g).coefficient()
    qscore_per_atom.append(qval)

  qscore_per_atom = flex.double(qscore_per_atom)

  qscore_df = aggregate_qscore_per_residue(model,np.array(qscore_per_atom),window=3)
  qscore_per_residue = flex.double(qscore_df["Q-scorePerResidue"].values)

  # Output
  if debug or params.debug:
    # Collect debug data
    result = {
      "atom_xyz":atoms_xyz,
      "probe_xyz":probe_xyzs,
      "probe_mask":probe_masks,
      "d_vals":full_flat_d,
      "g_vals":full_flat_g,
      "qscore_per_atom":qscore_per_atom,
      "mask_check_1d":mask_check_1d,
      "qscore_per_residue":qscore_per_residue,
      "qscore_dataframe":qscore_df
    }
  else:
    result = {
    "qscore_per_atom":qscore_per_atom,
    "qscore_per_residue":qscore_per_residue,
    "qscore_dataframe":qscore_df

    }
  return result
>>>>>>> Stashed changes


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

<<<<<<< Updated upstream
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
=======
################################################################################
#### KDTree implementation using flex arrays (no numpy/scipy)
################################################################################
class KDTreeFlexNode:
  def __init__(self, index, point, left=None, right=None):
    self.index = index
    self.point = point
    self.left = left
    self.right = right
>>>>>>> Stashed changes


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

<<<<<<< Updated upstream
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
=======

  def pre_sort_indices(self, points):
    # Sort indices for each axis and return the sorted indices
    x,y,z = points.parts()
    sorted_indices = [
       flex.sort_permutation(x),
       flex.sort_permutation(y),
       flex.sort_permutation(z)]
    return sorted_indices

  def build_tree(self, indices, points, depth):
    if not indices:
      return None

    axis = depth % self.dims


    sorted_indices = self.axis_sorted_indices[axis]


    # Step 4: Create an empty boolean mask of length N, initially set to False
    mask = flex.bool(len(sorted_indices),False)

    # Directly set mask to True for positions in your query array
    mask.set_selected(indices,True)

    # Step 5: Apply the mask to the sorted indices,
    #   then use this to create a sorted mask
    # This step seems to be where you're looking to optimize.
    # To directly use the sorted_indices to index into 'mask' and maintain sorting:
    sorted_mask = mask.select(sorted_indices)

    # Now, apply this sorted_mask to select from the sorted_indices
    sorted_indices_this_axis = sorted_indices.select(sorted_mask)


    if len(sorted_indices_this_axis) == 0:
      return None

    median_idx = len(sorted_indices_this_axis) // 2
    median_index = sorted_indices_this_axis[median_idx]

    left_indices = sorted_indices_this_axis[:median_idx]
    right_indices = sorted_indices_this_axis[median_idx + 1:]

    return KDTreeFlexNode(
        median_index,
        points[median_index:median_index+1],
        left=self.build_tree(left_indices, points, depth + 1),
        right=self.build_tree(right_indices, points, depth + 1)
    )
  def _nearest_neighbor(self, root, point, depth=0, best=None, k=1):
    if root is None:
      return best

    if best is None:
      best = []
>>>>>>> Stashed changes

    # infer params from shape
    n_shells, n_atoms, n_probes, _ = probe_xyz.focus()

    # APPLY MASK BEFORE INTERPOLATION

<<<<<<< Updated upstream
    probe_mask_cctbx_fullflat = []

    for val in probe_mask_cctbx:
        for _ in range(3):  # since A has an additional dimension of size 3
            probe_mask_cctbx_fullflat.append(val)
=======
    # Check the current root distance
    current_dist = root.point.max_distance(point)
    if len(best) < k or current_dist < best[-1]['dist']:
      best.append({'index': root.index, 'dist': current_dist})
      best.sort(key=lambda x: x['dist'])
      best = best[:k]  # Keep only k nearest

    # Check if we need to search the opposite branch
    if len(best) < k or (
      abs(point[0][axis] - root.point[0][axis]) < best[-1]['dist']):
      best = self._nearest_neighbor(opposite_branch, point, depth + 1, best, k)
>>>>>>> Stashed changes

    mask = flex.bool(probe_mask_cctbx_fullflat)
    # indices = flex.int([i for i in range(1, keep_mask_cctbx.size() + 1) for _ in range(3)])
    sel = probe_xyz_cctbx.select(mask)
    # sel_indices = indices.select(mask)
    masked_probe_xyz_flat_cctbx = flex.vec3_double(sel)

<<<<<<< Updated upstream
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
=======
  def query_nearest_neighbors(self, query_points, k=1):
    dists, inds = [], []
    for i,point in enumerate(query_points):
      point = query_points[i:i+1]
      nearest = self._nearest_neighbor(self.root, point, k=k)
      dists.append([n['dist'] for n in nearest])
      inds.append([n['index'] for n in nearest])
    return dists, inds
>>>>>>> Stashed changes


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


<<<<<<< Updated upstream
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
=======
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

def trilinear_interpolation(voxel_grid, coords, voxel_size=None, offset=None):
  """Numpy trilinear interpolation"""
  assert voxel_size is not None,(
      "Provide voxel size as an array or single value")

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
  assert A.shape == B.shape, (
      f"A and B must have the same shape, got: {A.shape} and {B.shape}")

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

def starmap_wrapper(task):
  """
  A generic wrapper function to call any worker function with specified arguments.
  Params:
    task (dict): A dictionary containing the 'func', 'args', and 'kwargs' keys.

  Returns:
    return: The result of the worker function call.
  """
  # Extract the function to call, its args, and kwargs from the task
  func = task['func']
  args = task.get('args', ())
  kwargs = task.get('kwargs', {})

  # Dynamically call the function with args and kwargs
  return func(*args, **kwargs)


def aggregate_qscore_per_residue(model,qscore_per_atom,window=3):
  # assign residue indices to each atom

  model = model.select(model.selection("not element H"))
  atoms = model.get_atoms()
  res_seqs = [atom.parent().parent().resseq_as_int() for atom in atoms]
  chain_ids = [atom.parent().parent().parent().id for atom in atoms]
  res_names = ["".join(atom.parent().parent().unique_resnames()) for atom in atoms]
  names = [atom.name.strip() for atom in atoms]
  df = pd.DataFrame({"resseq":res_seqs,
                    "chain_id":chain_ids,
                    "resname":res_names,
                    "name":names,
                    "Q-score":qscore_per_atom,

                    })

  # group atoms into residues
  df["rg_index"] = df.groupby(["chain_id","resseq","resname"]).ngroup()

  # average qscore by residue
  grouped_means = df.groupby(['chain_id',"resseq","resname","rg_index"],as_index=False)['Q-score'].mean()

  # roll over residue means
  for chain_id,group in grouped_means.groupby("chain_id"):
    grouped_means.loc[group.index, "Q-scorePerResidue"] = variable_neighbors_rolling_mean(group["Q-score"],window).values


  # now grouped means is a df with each row being a "residue" "QscoreRollingMean" is the per-residue value to match mapq

  # place back in atom df
  df = df.merge(grouped_means[['rg_index', 'Q-scorePerResidue']], on='rg_index', how='left')
  df.drop("rg_index",axis=1,inplace=True)
  return df

def variable_neighbors_rolling_mean(series, window=3):
  """
  Aggregate per-atom qscore to per-residue in the same
    manner as the original mapq program
  """
  # Container for the rolling means
  rolling_means = []

  # Calculate rolling mean for each index
  for i in range(len(series)):
    # Determine the start and end indices of the window
    start_idx = max(0, i - window)
    end_idx = min(len(series), i + window + 1)

    # Calculate mean for the current window
    window_mean = series.iloc[start_idx:end_idx].mean()
    rolling_means.append(window_mean)

  return pd.Series(rolling_means)


def write_bild_spheres(xyz,filename="sphere.bild",r=0.5):
  """
  Write a chimerax .bild file with spheres

  Params:
    xyz (np.array): Cartesian coordinates (N,3)
    filename (str): The filename for the .bild file
    r (float): sphere radius
  """
  out = ""
  for x,y,z in xyz:
    s = f".sphere {x} {y} {z} {r}\n"
    out+=s

  with open(filename,"w") as fh:
    fh.write(out)



def cctbx_atoms_to_df(atoms):
  """
  Build a pandas dataframe from a cctbx shared atoms object

  Params:
    atoms (iotbx_pdb_hierarchy_ext.af_shared_atom): The atom array

  Returns:
    df_atoms (pd.DataFrame): A pandas dataframe with the core data present
  """
  # Define values composition
  func_mapper = {
                        #"model_id", # model
                        "id":lambda atom: atom.i_seq,
                        "Chain":lambda atom: atom.parent().parent().parent().id,
                        "Resseq":lambda atom: atom.parent().parent().resseq_as_int(),
                        "Resname":lambda atom: atom.parent().parent().unique_resnames()[0],
                        "Name":lambda atom: atom.name.strip(),
                        "Element":lambda atom: atom.element,
                        "Alt Loc": lambda atom: atom.parent().altloc
  }

  # Build as dictionaries
  data = defaultdict(list)
  for atom in atoms:
    for key,func in func_mapper.items():
      data[key].append(func(atom))


  # Include values non-composition
  xyz = atoms.extract_xyz().as_numpy_array()
  data["x"] = xyz[:,0]
  data["y"] = xyz[:,1]
  data["z"] = xyz[:,2]

  # Build dataframe from dictionaries
  df_atoms = pd.DataFrame(data,index=list(range(len(atoms))))

  return df_atoms
>>>>>>> Stashed changes
