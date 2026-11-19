"""
Deviations of bond lengths and angles from ideal values versus map values at
atom centers, for macromolecule (protein and nucleotide) non-H atoms.
"""

from cctbx.array_family import flex
from libtbx import group_args

def rank_0_1(values):
  """
  Rank of each value among all values, scaled to be between 0 (smallest) and
  1 (largest). Equal values get the same (average) rank.

  Parameters
  ----------
  values : flex.double

  Returns
  -------
  flex.double
  """
  n = values.size()
  if n == 0: return flex.double()
  if n == 1: return flex.double([1.])
  perm = flex.sort_permutation(values)
  ranks = flex.double(n, 0)
  i = 0
  while i < n:
    j = i
    while j+1 < n and values[perm[j+1]] == values[perm[i]]: j += 1
    rank = (i + j) / 2.
    for k in range(i, j+1): ranks[perm[k]] = rank
    i = j + 1
  return ranks / (n - 1)

def map_values(fmodel, sites_frac, map_type="2mFo-DFmodel",
               sigma_scale=False):
  """
  Values of a map at given positions (eight-point interpolation).

  Parameters
  ----------
  fmodel : mmtbx.f_model.manager
  sites_frac : flex.vec3_double
      Fractional coordinates of the positions.
  map_type : str
      Any map type known to mmtbx.map_tools.electron_density_map.
  sigma_scale : bool
      Scale the map to have zero mean and unit standard deviation.

  Returns
  -------
  flex.double
  """
  from cctbx import maptbx
  from mmtbx import map_tools
  import cctbx.miller
  crystal_gridding = fmodel.f_obs().crystal_gridding(
    d_min             = fmodel.f_obs().d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 1./4)
  map_coefficients = map_tools.electron_density_map(
    fmodel = fmodel).map_coefficients(
      map_type     = map_type,
      isotropize   = True,
      fill_missing = False)
  fft_map = cctbx.miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = map_coefficients)
  if sigma_scale: fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  result = flex.double()
  for site_frac in sites_frac:
    result.append(map_data.eight_point_interpolation(site_frac))
  return result

def _per_restraint(proxies, restraint_type, value_name, sites_cart,
                   i_seqs_all, maps):
  """
  Table of restraints: one entry per proxy in each field (see run).
  """
  result = group_args(
    i_seqs       = [],
    ideal        = flex.double(),
    model        = flex.double(),
    delta        = flex.double(),
    weight       = flex.double(),
    origin_id    = flex.int(),
    map_2mfo_dfc = [],
    map_mfo_dfc  = [],
    rank         = [])
  for proxy in proxies:
    restraint = restraint_type(sites_cart = sites_cart, proxy = proxy)
    result.i_seqs.append(tuple([i_seqs_all[i] for i in proxy.i_seqs]))
    result.ideal.append(getattr(restraint, value_name+"_ideal"))
    result.model.append(getattr(restraint, value_name+"_model"))
    result.delta.append(-restraint.delta)
    result.weight.append(proxy.weight)
    result.origin_id.append(proxy.origin_id)
    for key in ["map_2mfo_dfc", "map_mfo_dfc", "rank"]:
      values = getattr(maps, key)
      getattr(result, key).append(tuple([values[i] for i in proxy.i_seqs]))
  return result

def _per_atom(proxies, restraints, rank, bins, d_min):
  """
  Mean absolute deviation from ideal over restraints of each atom (deltas),
  with density rank of the atom (maps), and mean of deltas in bins of rank.
  Atoms are in order of their first appearance in proxies.
  """
  atom_deltas = {}
  for proxy, delta in zip(proxies, restraints.delta):
    for i_seq in proxy.i_seqs:
      if i_seq not in atom_deltas: atom_deltas[i_seq] = flex.double()
      atom_deltas[i_seq].append(abs(delta))
  deltas = flex.double()
  maps   = flex.double()
  for i_seq, values in atom_deltas.items():
    deltas.append(flex.mean(values))
    maps  .append(rank[i_seq])
  x, y = flex.double(), flex.double()
  for low, high in bins:
    in_bin = deltas.select((maps >= low) & (maps <= high))
    if in_bin.size() == 0: continue
    x.append(high)
    y.append(flex.mean(in_bin))
  return group_args(
    bins       = bins,
    x          = x,
    y          = y,
    deltas     = deltas,
    maps       = maps,
    d_min      = "%3.1f"%d_min,
    restraints = restraints)

def run(model, fmodel, n_bins=5):
  """
  Deviations of bonds and angles from ideal values and map values at atom
  centers, for protein and nucleotide non-H atoms. All restraints with all
  atoms among these are used (not symmetry bonds).

  Parameters
  ----------
  model : mmtbx.model.manager
      Model with restraints; it is not changed.
  fmodel : mmtbx.f_model.manager
      Source of maps and R-factors.
  n_bins : int
      Number of equal bins of density rank between 0 and 1.

  Returns
  -------
  libtbx.group_args or None
      None if there are no protein or nucleotide atoms. Otherwise:

      r_work, r_free : R-factors of fmodel.

      g : geometry statistics of the atoms used.

      atoms : one entry per atom used, in the order of the model:
        i_seqs (indices in the given model), id_strs, map_2mfo_dfc and
        map_mfo_dfc (2mFo-DFc and mFo-DFc map values in sigma units),
        rank (rank of 2mFo-DFc value among the atoms used, from 0: lowest
        to 1: highest).

      bonds, angles : group_args with

        deltas, maps : one entry per atom that has restraints: mean absolute
          deviation from ideal over the restraints of the atom and rank of
          the atom;

        bins, x, y : bins of rank, upper limits of non-empty bins and mean
          of deltas in them;

        d_min : resolution, as a string;

        restraints : one entry per restraint: i_seqs (tuple of indices in
          the given model), ideal, model, delta (model - ideal), weight,
          origin_id, map_2mfo_dfc, map_mfo_dfc, rank (tuples of values at
          the atoms of the restraint).
  """
  from cctbx import geometry_restraints
  selection = model.selection(
    "(protein or nucleotide) and not (element H or element D)")
  selection = selection & ~model.solvent_selection()
  if selection.count(True) == 0: return None
  i_seqs_all = selection.iselection()
  model = model.select(selection)
  sites_cart = model.get_sites_cart()
  sites_frac = model.get_sites_frac()
  maps = group_args(
    map_2mfo_dfc = map_values(fmodel = fmodel, sites_frac = sites_frac,
      map_type = "2mFo-DFmodel", sigma_scale = True),
    map_mfo_dfc = map_values(fmodel = fmodel, sites_frac = sites_frac,
      map_type = "mFo-DFmodel", sigma_scale = True))
  maps.rank = rank_0_1(maps.map_2mfo_dfc)
  geometry = model.get_restraints_manager().geometry
  bond_proxies, asu = geometry.get_all_bond_proxies(sites_cart = sites_cart)
  bins = [[i/n_bins, (i+1)/n_bins] for i in range(n_bins)]
  d_min = fmodel.f_obs().d_min()
  results = []
  for proxies, restraint_type, value_name in [
      (bond_proxies,           geometry_restraints.bond,  "distance"),
      (geometry.angle_proxies, geometry_restraints.angle, "angle")]:
    restraints = _per_restraint(
      proxies        = proxies,
      restraint_type = restraint_type,
      value_name     = value_name,
      sites_cart     = sites_cart,
      i_seqs_all     = i_seqs_all,
      maps           = maps)
    results.append(_per_atom(
      proxies    = proxies,
      restraints = restraints,
      rank       = maps.rank,
      bins       = bins,
      d_min      = d_min))
  return group_args(
    r_work = fmodel.r_work(),
    r_free = fmodel.r_free(),
    g      = model.geometry_statistics().result(slim = True),
    atoms  = group_args(
      i_seqs       = i_seqs_all,
      id_strs      = [a.id_str() for a in model.get_hierarchy().atoms()],
      map_2mfo_dfc = maps.map_2mfo_dfc,
      map_mfo_dfc  = maps.map_mfo_dfc,
      rank         = maps.rank),
    bonds  = results[0],
    angles = results[1])
