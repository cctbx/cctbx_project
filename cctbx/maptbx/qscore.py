import sys
from multiprocessing import Pool
import numpy as np
import numpy.ma as ma
from scipy.spatial import KDTree
import pandas as pd
from libtbx.utils import null_out
from cctbx.array_family import flex


################################################################################
#### Probe generation functions
################################################################################

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
    pts = np.array(pts)  # should be shape (n_probes,3)
    assert pts.shape == (n_probes_max,3),(

    f"Pts shape must be ({n_probes_max},3), not {pts.shape})")

    if pts.shape == (0,): # all probes clashed
      pts = np.full((n_probes_max,3),np.nan)

    all_pts.append(pts)


  # prepare output
  probe_xyz = np.stack(all_pts)
  probe_mask = ~(np.isnan(probe_xyz))[:,:,0]

  return probe_xyz, probe_mask

################################################################################
#### Run shell functions for multiple shells(possibly in parallel)
################################################################################
def get_probes(
    atoms_xyz=None,
    sites_cart = None,
    atoms_tree = None,
    params=None,
    selection_bool_np = None,
    selection_bool_flex=None,
    worker_func=None,
    log=None):
  """
  Generate probes for multiple radial shells (params.shells)
  """
  if params.backend == "numpy":
    assert worker_func in [shell_probes_progressive,shell_probes_precalculate]
    if atoms_tree is None:
      atoms_tree = KDTree(atoms_xyz)
    selection_bool = selection_bool_np
  else:
    assert False, f"Unrecognized backend: {params.backend}"


  kwargs_list = [
    {
        "atoms_xyz":atoms_xyz,   # A numpy array of shape (N,3)
        "atoms_tree":atoms_tree,  # An atom_xyz scipy kdtree
        "selection_bool":selection_bool,# Boolean atom selection
        "n_probes_target":params.n_probes_target,# The desired number of probes per shell
        "n_probes_max":params.n_probes_max,  # The maximum number of probes allowed
        "n_probes_min":params.n_probes_min,
        "RAD":RAD,          # The nominal radius of this shell
        "rtol":params.rtol,         # Multiplied with RAD to get actual radius
        "log":log,
      }
     for RAD in params.shells]
  #DEBUG
  if params.nproc=="DEBUG":
    with Pool() as pool:
      results = pool.map(starmap_wrapper, kwargs_list)
  else:
    results = []
    for kwarg in kwargs_list:
        result = worker_func(**kwarg)
        results.append(result)


  probe_xyz = [result[0] for result in results]
  probe_mask = [result[1] for result in results]

  if params.backend == "numpy":
    return np.stack(probe_xyz), np.array(probe_mask)
  else:
    return probe_xyz, probe_mask
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

def calc_qscore(mmm,params,log=null_out(),debug=False):
  """
  Calculate qscore from map model manager
  """
  model = mmm.model()
  # never do hydrogen
  model = model.select(model.selection("not element H"))
  mm = mmm.map_manager()


  # Get atoms
  atoms_xyz = model.get_sites_cart().as_numpy_array()

  # do selection
  if params.selection_str != None:
    selection_bool = mmm.model().selection(params.selection_str).as_numpy_array() # boolean
    if selection_bool.sum() ==0:
      print("Finished... nothing selected")
      return {"qscore_per_atom":None}
  else:
    selection_bool = np.full(len(atoms_xyz),True)

  # determine worker func
  if params.probe_allocation_method == "precalculate":
    worker_func=shell_probes_precalculate
  else:
    worker_func=shell_probes_progressive


  # Get probes and probe mask (probes to reject)
  probe_xyz,probe_mask = get_probes(
    atoms_xyz=atoms_xyz,
    atoms_tree = None,
    params=params,
    selection_bool_np = selection_bool,
    worker_func=worker_func,
    log = log,
    )

  # after the probe generation, versions 1 and 2 are the same

  # infer params from shape
  n_shells, n_atoms, n_probes, _ = probe_xyz.shape

  # flatten
  probe_xyz_flat = probe_xyz.reshape((n_atoms * n_shells * n_probes, 3))
  probe_mask_flat = probe_mask.reshape(-1)  # (n_shells*n_atoms*n_probes,)
  masked_probe_xyz_flat = probe_xyz_flat[probe_mask_flat]

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
      "atom_xyz":atoms_xyz,
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
                    "Q-score":np.array(qscore_per_atom),

                    })

  # group atoms into residues
  df["rg_index"] = df.groupby(["chain_id","resseq","resname"]).ngroup()

  # average qscore by residue
  grouped_means = df.groupby(['chain_id',"resseq","resname","rg_index"],as_index=False)['Q-score'].mean()
  # DEBUG:

  if 'Q-score' not in grouped_means.columns:
    import pdb
    pdb.set_trace()


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


