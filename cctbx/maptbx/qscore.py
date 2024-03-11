from __future__ import division
from collections import defaultdict
import warnings

from libtbx.utils import null_out
from cctbx.array_family import flex
from libtbx import easy_mp
import numpy as np
import numpy.ma as ma
from scipy.spatial import KDTree
import pandas as pd


master_phil_str = """
  qscore
  {

    nproc = 1
        .type = int
        .help = Number of processors to use
        .short_caption = Number of processors to use
        .expert_level = 1
    n_probes = 32
        .type = int
        .help = Number of radial probes to use
        .short_caption = Number of radial probes to use
        .expert_level = 1

    selection = None
      .type = str
      .help = Only calculate atoms within this selection
      .short_caption = Only test atoms within this selection
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
      .short_caption = The number of radial shells (includes start/stop, so minimum 2)
      .expert_level = 1

    rtol = 0.9
      .type = float
      .help = Mapq rtol value, the "real" shell radii are r*rtol

    write_probes = False
      .type = bool
      .help = Write the qscore probes as a .bild file to visualize in Chimera

    write_to_bfactor_pdb = False
      .type = bool
      .help = Write out a pdb file with the Q-score per atom in the B-factor field
  }

  """

################################################################################
#### Probe generation functions
################################################################################



# Fast numpy version
def generate_probes_np(sites_cart, rad, n_probes):
  """
  Generate probes using the same methodology as Pintile mapq, but vectorized

  sites_cart: np array of shape (n_atoms,3)
  rad: the radius at which to place the probes
  N: the number of probes per atom

  Returns:
    probes (np.ndarray): shape (n_atoms,n_probes,3)
  """
  assert sites_cart.ndim == 2 and sites_cart.shape[-1]==3, (
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
  probes = probes.reshape(-1, 1, 3) + sites_cart.reshape(1, -1, 3)

  # reshape (n_atoms,n_probes,3)
  probes = probes.swapaxes(0,1)
  return probes


################################################################################
#### Run shell functions for multiple shells(possibly in parallel)
################################################################################

class GatherProbes:
  def __init__(self,
               func,
               fixed_kwargs,
               ):
    self.func = func
    self.fixed_kwargs = fixed_kwargs

  def __call__(self,RAD):

    return self.func(RAD=RAD,**self.fixed_kwargs)


def get_probes(
    sites_cart=None,
    atoms_tree = None,
    shells = None,
    n_probes = None,
    rtol=None,
    nproc=1,
    selection_bool = None,
    worker_func=None,
    log = null_out()):

  """
  Generate probes for multiple radial shells (shells)
  """
  # Create before multiprocessing
  atoms_tree = KDTree(sites_cart)



  assert shells is not None, "Must provide explicit radial shells"
  fixed_kwargs = {
      "sites_cart": sites_cart,  # A numpy array of shape (N,3)
      "atoms_tree": atoms_tree,  # An atom_xyz scipy kdtree
      "selection_bool": selection_bool,  # Boolean atom selection
      "n_probes": n_probes,  # The desired number of probes per shell
      "rtol": rtol,  # Multiplied with RAD to get actual radius
      "log": log,
      }

  gather = GatherProbes(worker_func,fixed_kwargs)

  results = easy_mp.pool_map(
      processes=nproc,
      fixed_func=gather,
      args=shells)

  probe_xyz = [result[0] for result in results]
  probe_mask = [result[1] for result in results]

  n_shells = len(shells)
  probe_xyz_stacked = np.empty((n_shells,*probe_xyz[0].shape))
  for i,array in enumerate(probe_xyz):
    probe_xyz_stacked[i] = array
  return probe_xyz_stacked, np.array(probe_mask)


def shell_probes_precalculate(
      sites_cart=None,   # A numpy array of shape (N,3)
      atoms_tree=None,  # An atom_xyz scipy kdtree
      selection_bool=None,# Boolean atom selection
      n_probes=8,# The desired number of probes per shell
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
    atoms_tree = KDTree(np.array(sites_cart))

  # Manage log
  if log is None:
    log = null_out()

  # manage selection input
  if selection_bool is None:
    selection_bool = np.full(len(sites_cart),True)


  # do selection
  sites_cart_sel = sites_cart[selection_bool]

  # get probe coordinates
  probe_xyz = generate_probes_np(sites_cart_sel, RAD, n_probes)
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
  n_probes_min = 4
  strict=False
  n_probes_per_atom = probe_mask.sum(axis=1)
  insufficient_probes = np.where(n_probes_per_atom<n_probes)[0]
  problematic_probes = np.where(n_probes_per_atom<n_probes_min)[0]
  if strict:
    if n_probes_per_atom.min() >= n_probes_min:
      print(
      f"Some atoms have less than {n_probes_min} probes. \
          ({len(problematic_probes)}). Consider raising n_probes",file=log)

  return probe_xyz, probe_mask


def calc_qscore(mmm,
                selection=None,
                shells=None,
                n_probes=8,
                rtol=0.9,
                nproc=1,
                log=null_out(),
                debug=False):
  """
  Calculate qscore from map model manager
  """
  model = mmm.model()
  # never do hydrogen
  model = model.remove_hydrogens()
  mmm.set_model(model)
  mm = mmm.map_manager()


  # Get atoms
  sites_cart = model.get_sites_cart().as_numpy_array()

  # do selection
  if selection != None:
    selection_bool = mmm.model().selection(selection) # boolean
    if selection_bool.sum() ==0:
      print("Finished... nothing selected",file=log)
      return {"qscore_per_atom":None}
  else:
    selection_bool = flex.bool(model.get_number_of_atoms(),True)


  # determine worker func
  worker_func=shell_probes_precalculate


  # Get probes and probe mask (probes to reject)
  probe_xyz,probe_mask = get_probes(
    sites_cart=sites_cart,
    atoms_tree = None,
    shells=shells,
    n_probes=n_probes,
    rtol=rtol,
    nproc=nproc,
    selection_bool = selection_bool,
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
  masked_density = mm.density_at_sites_cart(
    flex.vec3_double(masked_probe_xyz_flat)).as_numpy_array()

  d_vals = np.full((n_shells, n_atoms, n_probes),np.nan)
  d_vals[probe_mask] = masked_density

  # g vals
  # create the reference data
  M = volume
  maxD = min(M.mean() + M.std() * 10, M.max())
  minD = max(M.mean() - M.std() * 1, M.min())
  A = maxD - minD
  B = minD
  u = 0
  sigma = 0.6
  x = np.array(shells)
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
  model = model.select(flex.bool(selection_bool))
  qscore_df = aggregate_qscore_per_residue(model,q,window=3)
  q = flex.double(q)
  qscore_per_residue = flex.double(qscore_df["Q-Residue"].values)

  # Output
  result = {
    "qscore_per_atom":q,
    "qscore_per_residue":qscore_per_residue,
    "qscore_dataframe":qscore_df
    }

  if debug:
    # Collect debug data
    result.update({
      "atom_xyz":sites_cart,
      "probe_xyz":probe_xyz,
      "probe_mask":probe_mask,
      "d_vals":d_vals,
      "g_vals":g_vals,
    })

  return result


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



def aggregate_qscore_per_residue(model,qscore_per_atom,window=3):
  # assign residue indices to each atom

  atoms = model.get_atoms()
  df = cctbx_atoms_to_df(atoms)
  df["Q-score"] = qscore_per_atom
  df["rg_index"] = df.groupby(["chain_id", "resseq", "resname"]).ngroup()
  grouped_means = df.groupby(['chain_id', "resseq", "resname", "rg_index"],
                             as_index=False)['Q-score'].mean().rename(
                               columns={'Q-score': 'Q-Residue'})

  #grouped_means['RollingMean'] = None  # Initialize column to avoid KeyError

  # Until pandas is updated, need to suppress warning
  warnings.filterwarnings("ignore", category=FutureWarning)
  # for chain_id, group in grouped_means.groupby("chain_id"):
  #   rolling_means = variable_neighbors_rolling_mean(group['Q-Residue'], window)
  #   grouped_means.loc[group.index, 'RollingMean'] = rolling_means.values



  # Merge the updated 'Q-Residue' and 'RollingMean' back into the original DataFrame
  df = df.merge(grouped_means[['rg_index', 'Q-Residue']], on='rg_index', how='left')
  df.drop("rg_index", axis=1, inplace=True)
  # df["Q-ResidueRolling"] = df["RollingMean"].astype(float)
  # df.drop(columns=["RollingMean"],inplace=True)
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
                        "model_id":lambda atom: atom.parent().parent().parent().parent().id,
                        "chain_id":lambda atom: atom.parent().parent().parent().id,
                        "resseq":lambda atom: atom.parent().parent().resseq_as_int(),
                        "resname":lambda atom: atom.parent().parent().unique_resnames()[0],
                        "name":lambda atom: atom.name.strip(),
                        "element":lambda atom: atom.element,
                        "altloc": lambda atom: atom.parent().altloc
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
