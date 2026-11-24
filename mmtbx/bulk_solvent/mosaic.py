"""
Mosaic bulk-solvent model.

Reference: Afonine, Adams, Sobolev & Urzhumtsev (2024). Protein Science,
33, e4909.

The standard flat bulk-solvent model fills every part of the unit cell that
is not occupied by the atomic model with the same solvent. The solvent mask
consists of one large connected region (main) and, usually, a number of
isolated regions (secondary): pockets and cavities of the macromolecule.
Their content differs from one to another: empty, partially or fully filled
with disordered solvent, or filled with something denser than bulk solvent.

Here each secondary region n gets its own occupancy q_n relative to the main
bulk solvent:

  F_bulk = k_mask(s) * [ F_whole + sum_n (q_n - 1) * F_n ],

F_whole: structure factors of the whole mask (as used by the flat model),
F_n: structure factors of region n. Regions smaller than volume_min are
emptied (q_n=0); all others are fitted to the data, 0 <= q_n <= 1, with
work reflections only.

Entry point: run(fmodel, params, log). Used by
mmtbx.f_model.manager.update_all_scales(mosaic=True).
"""

from __future__ import absolute_import, division, print_function

from cctbx.array_family import flex
from cctbx import maptbx
from scitbx import matrix
from scitbx import minimizers
from scitbx import fftpack
from libtbx import adopt_init_args
from libtbx import group_args
from mmtbx.bulk_solvent import kbu_refinery
import mmtbx.masks
import boost_adaptbx.boost.python as bp
import math
mosaic_ext = bp.import_ext("mmtbx_mosaic_ext")

# A bound that is not a bound.
NO_LIMIT = 1.e30

# Upper limit of a fitted occupancy: a region has at most as much solvent as
# the main bulk solvent. What a region has in excess of that (for example a
# ligand that is not in the model yet) is not bulk solvent: it is left in
# the difference map, as it is with the flat model. Without the limit such
# density is absorbed by the solvent model, and at low resolution R-free
# gets worse.
OCCUPANCY_MAX = 1.0

# Number of calls of run() that did the region analysis (used by tests to
# make sure the analysis is not repeated when it should not be).
number_of_runs = 0

# Cost of one (grid point, reflection) pair of direct summation and of one
# grid point of FFT, in arbitrary but equal units; used only to choose the
# faster of the two ways to get the structure factors of a region.
_COST_DIRECT = 1.0
_COST_FFT    = 80.0

# ------------------------------------------------------------------------------
# Regions of the bulk-solvent mask
# ------------------------------------------------------------------------------

class regions(object):
  """
  Connected regions of the bulk-solvent mask of a model.

  The mask is computed in the whole unit cell (P1) on the same grid and with
  the same radii as mmtbx.masks.manager does for structure factors.
  Symmetry-related regions are one region.

  Parameters
  ----------
  xray_structure : cctbx.xray.structure
      Atoms that define the mask, as mmtbx.masks.manager holds them (hydrogens
      and zero-occupancy atoms already removed if the mask ignores them).
  mask_params : libtbx.phil.scope_extract
      mmtbx.masks.mask_master_params.
  d_min : float, optional
      Resolution of the data; needed only if
      mask_params.use_resolution_based_gridding is True.

  Attributes
  ----------
  n_real : tuple of int
  n_grid : int
  volume_per_point : float
  solvent_fraction : float
  main : libtbx.group_args
      The largest region: label, n_points, volume (its points are not kept,
      see solvent_mask).
  secondary : list of libtbx.group_args
      All other regions, largest first. Fields: label (int), n_points (int),
      volume (A^3, all symmetry copies in the unit cell together), points
      (flex.size_t, 1d grid indices of all copies, sorted), centroid_frac
      (fractional coordinates of the centre of one connected copy).
  """

  def __init__(self, xray_structure, mask_params, d_min=None):
    mp = mask_params
    self.unit_cell = xray_structure.unit_cell()
    self.space_group = xray_structure.space_group()
    self._xray_structure = xray_structure
    self._mask_params = mask_params
    if(not mp.use_resolution_based_gridding): d_min = None
    self.n_real = tuple(mmtbx.masks.gridding_n_real(
      mask_params = mp, xray_structure = xray_structure, d_min = d_min))
    mask = self.solvent_mask()
    self.n_grid = self.n_real[0]*self.n_real[1]*self.n_real[2]
    self.volume_per_point = self.unit_cell.volume()/self.n_grid
    co = maptbx.connectivity(
      map_data                   = mask,
      threshold                  = 0.01,
      preprocess_against_shallow = False,
      wrapping                   = True)
    del mask
    # Labels of connected pieces; then the same with all symmetry copies of
    # a region under one label
    pieces = co.result().as_1d().deep_copy()
    co.merge_symmetry_related_regions(space_group = self.space_group)
    labels = co.result().as_1d()
    counts = co.regions() # label 0: not solvent
    del co
    self.main = None
    self.secondary = []
    solvent = [(c, i) for i, c in enumerate(counts) if i != 0 and c > 0]
    self.solvent_fraction = sum([c for c, i in solvent])/self.n_grid
    if(len(solvent) == 0): return
    main_count, main_label = max(solvent)
    self.main = group_args(
      label    = main_label,
      n_points = main_count,
      volume   = main_count*self.volume_per_point)
    if(len(solvent) == 1): return
    # Points of all secondary regions in one pass, grouped by label and
    # sorted within a label.
    points = ((labels != 0) & (labels != main_label)).iselection()
    labels_sel = labels.select(points)
    key = labels_sel.as_double()*float(self.n_grid) + points.as_double()
    perm = flex.sort_permutation(key)
    points = points.select(perm)
    labels_sel = labels_sel.select(perm)
    pieces_sel = pieces.select(points)
    start = 0
    for label, count in enumerate(counts):
      if(label == 0 or label == main_label or count == 0): continue
      end = start + count
      assert labels_sel[start] == label and labels_sel[end-1] == label
      p = points[start:end]
      self.secondary.append(group_args(
        label         = label,
        n_points      = count,
        volume        = count*self.volume_per_point,
        points        = p,
        centroid_frac = self._centroid(
          points = p.select(pieces_sel[start:end] == pieces_sel[start]))))
      start = end
    assert start == points.size()
    self.secondary.sort(key = lambda r: (-r.n_points, r.label))

  def solvent_mask(self, padded=False):
    """
    The whole bulk-solvent mask (main and secondary regions) in the unit
    cell: flex.double on n_real, one in the solvent, zero elsewhere.
    Computed at each call. padded: with the padding that an in-place Fourier
    transform needs (see structure_factors.of_map).
    """
    mp = self._mask_params
    mask = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = self._xray_structure,
      p1                       = True,
      for_structure_factors    = True,
      solvent_radius           = mp.solvent_radius,
      shrink_truncation_radius = mp.shrink_truncation_radius,
      n_real                   = self.n_real,
      in_asu                   = False).mask_data
    assert mask.is_padded() and mask.focus() == self.n_real
    if(not padded): maptbx.unpad_in_place(map = mask)
    return mask

  def _centroid(self, points):
    """
    Centre (fractional) of a connected set of grid points; the set may cross
    the unit cell boundaries.
    """
    nx,ny,nz = self.n_real
    p = points.as_double()
    jz = flex.fmod(p, nz)
    jy = flex.fmod(flex.floor(p/nz), ny)
    jx = flex.floor(p/(ny*nz))
    result = []
    for j, n in [(jx,nx), (jy,ny), (jz,nz)]:
      # shortest shifts from the first point, then their mean
      d = (j - j[0])/n
      d = d - flex.floor(d + 0.5)
      c = j[0]/n + flex.mean(d)
      result.append(c - math.floor(c))
    return tuple(result)

  def mask_of(self, regions):
    """
    Map (flex.double on n_real) that is one in the given regions and zero
    elsewhere.
    """
    return self.mask_of_points(points = [r.points for r in regions])

  def mask_of_points(self, points, values=None, padded=False):
    """
    Map that has values[i] (one if values is None) at the grid points
    points[i] (flex.size_t), zero elsewhere. padded: with the padding that
    an in-place Fourier transform needs (see structure_factors.of_map).
    """
    if(values is None): values = [1]*len(points)
    nz = self.n_real[2]
    grid = flex.grid(self.n_real)
    extra = 0
    if(padded):
      m_real = fftpack.real_to_complex_3d(self.n_real).m_real()
      grid = flex.grid(m_real).set_focus(self.n_real)
      extra = m_real[2] - nz
    m = flex.double(grid.size_1d(), 0)
    for p, v in zip(points, values):
      if(v == 0 or p.size() == 0): continue
      # position in the padded array: each row along z is longer by extra
      if(extra > 0): p = p + (p/nz)*extra # integer division (flex.size_t)
      m.set_selected(p, v)
    m.reshape(grid)
    return m

class structure_factors(object):
  """
  Structure factors of regions of the mask for the reflections of
  miller_array with d >= d_min (zero for the others, and for those not in
  selection if that is given), by direct summation
  over the region grid points (cctbx.maptbx.loft) or by FFT of the unit
  cell, whichever is faster for the region.

  Parameters
  ----------
  regions : mmtbx.bulk_solvent.mosaic.regions
  miller_array : cctbx.miller.array
      Defines the reflections; data are not used.
  d_min : float
      High-resolution limit of the bulk-solvent contribution
      (mask_params.Fmask_res_high).
  selection : flex.bool, optional
      Reflections of miller_array to compute structure factors for.

  Attributes
  ----------
  selection : flex.bool
      Reflections that structure factors are computed for.
  n_selected : int
  """

  def __init__(self, regions, miller_array, d_min, selection=None):
    self.regions = regions
    self.size = miller_array.indices().size()
    self.selection = miller_array.d_spacings().data() >= d_min
    if(selection is not None): self.selection = self.selection & selection
    self.n_selected = self.selection.count(True)
    ms = miller_array.select(self.selection)
    self.low = ms.customized_copy(
      data = flex.complex_double(ms.indices().size(), 0), sigmas = None)

  def of_points(self, points, method=None, compact=False):
    """
    Parameters
    ----------
    points : flex.size_t
        1d grid indices of all symmetry copies of the region(s).
    method : str, optional
        "direct" or "fft"; default: the faster one.
    compact : bool
        Return values for the selected reflections only.

    Returns
    -------
    flex.complex_double
        One value per reflection of miller_array, or per selected
        reflection (compact=True).
    """
    assert method in [None, "direct", "fft"]
    if(self.n_selected == 0 or points.size() == 0):
      if(compact): return flex.complex_double(self.n_selected, 0)
      return flex.complex_double(self.size, 0)
    # Direct summation is over one of the points related by lattice
    # centring (same result, see maptbx.loft)
    space_group = self.regions.space_group
    if(method is None):
      method = "fft"
      if(_COST_DIRECT*points.size()*self.n_selected/space_group.n_ltr() <
         _COST_FFT*self.regions.n_grid): method = "direct"
    if(method == "direct"):
      data = maptbx.loft(
        points           = points,
        n_real           = self.regions.n_real,
        miller_indices   = self.low.indices(),
        volume_per_point = self.regions.volume_per_point,
        space_group      = space_group)
    else:
      data = self.of_map(
        map_data = self.regions.mask_of_points(points = [points],
          padded = True),
        compact  = True)
    if(compact): return data
    return flex.complex_double(self.size, 0).set_selected(
      self.selection, data)

  def of_map(self, map_data, compact=False):
    """
    FFT of a map given in the whole unit cell (it must obey the space group
    symmetry). Returns one value per reflection of miller_array, or per
    selected reflection (compact=True).

    A padded map (see regions.mask_of_points, regions.solvent_mask) is
    transformed in place: no copy of the map is made, and the map cannot be
    used after this call. Whole-cell maps are what limits memory for large
    unit cells.
    """
    if(self.n_selected == 0):
      if(compact): return flex.complex_double(0)
      return flex.complex_double(self.size, 0)
    data = self.low.structure_factors_from_map(
      map          = map_data,
      in_place_fft = map_data.is_padded(),
      use_scale    = True,
      # the map is real: one half of its transform is computed, the rest
      # follows from Friedel's law (also for anomalous sets of indices)
      anomalous_flag = False).data()
    if(compact): return data
    return flex.complex_double(self.size, 0).set_selected(
      self.selection, data)

  def mask_mismatch(self, f_mask):
    """
    Relative difference between f_mask (cctbx.miller.array, the mask
    structure factors that fmodel uses) and the structure factors of all
    regions found here taken together. Zero if the regions are the regions
    of that mask.
    """
    r = self.regions
    if(r.main is None): return 1.
    f = self.of_map(map_data = r.solvent_mask(padded = True))
    den = flex.sum(flex.abs(f_mask.data()))
    if(den == 0): return 1.
    return flex.sum(flex.abs(f - f_mask.data()))/den

# ------------------------------------------------------------------------------
# Search for region occupancies
# ------------------------------------------------------------------------------

def alg4(f_obs, base, G, lower=None, upper=None, x_start=None, n_cycles=5,
         normal_matrix=None):
  """
  Phased simultaneous search (Alg4 of Afonine et al., 2024) for the
  coefficients x_n of

    F_model = base + sum_n x_n * G_n

  minimizing sum |F_obs*exp(i*phi) - F_model|^2, phi being the phases of the
  current F_model. All x_n are found at once from a system of linear
  equations, solved with bounds lower <= x_n <= upper; phases are then
  updated and the search repeated, n_cycles times in total.

  Parameters
  ----------
  f_obs : flex.double
  base : flex.complex_double
      Part of the model that is not fitted.
  G : list of flex.complex_double
  lower, upper : flex.double, optional
      Bounds for x; no bound if None.
  x_start : flex.double, optional
      Starting values, zero if None; they define the first set of phases.
  n_cycles : int
  normal_matrix : flex.double (n x n), optional
      mosaic_ext.normal_matrix(G), if it has been computed already (it is
      the expensive part and does not depend on f_obs and base).

  Returns
  -------
  flex.double
  """
  n = len(G)
  if(normal_matrix is None): normal_matrix = flex.double(flex.grid(0, 0))
  if(lower is None):   lower   = flex.double(n, -NO_LIMIT)
  if(upper is None):   upper   = flex.double(n,  NO_LIMIT)
  if(x_start is None): x_start = flex.double(n, 0)
  return mosaic_ext.alg4(
    G        = G,
    f_obs    = f_obs,
    base     = base,
    lower    = lower,
    upper    = upper,
    x_start  = x_start,
    n_cycles = n_cycles,
    normal_matrix = normal_matrix)

class minimizer2(object):
  """
  L-BFGS minimizer (optionally with curvatures) used by algorithm_2.
  """

  def __init__(self, calculator, min_iterations=0, max_iterations=2000):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.n = self.x.size()
    self.cntr=0

  def run(self, use_curvatures):
    self.minimizer = kbu_refinery.lbfgs_run(
      target_evaluator=self,
      min_iterations=self.min_iterations,
      max_iterations=self.max_iterations,
      use_curvatures=use_curvatures)
    self(requests_f_and_g=True, requests_diag=use_curvatures)
    return self

  def __call__(self, requests_f_and_g, requests_diag):
    self.cntr+=1
    self.calculator.update_target_and_grads(x=self.x)
    if (not requests_f_and_g and not requests_diag):
      requests_f_and_g = True
      requests_diag = True
    if (requests_f_and_g):
      self.f = self.calculator.target()
      self.g = self.calculator.gradients()
      self.d = None
    if (requests_diag):
      self.d = self.calculator.curvatures()
      if(self.d.all_eq(0)): self.d=None
      else:
        self.d = 1 / self.d
    return self.x, self.f, self.g, self.d

class tg(object):
  """
  Target, gradients and curvatures of the unphased search (algorithm_2):

    T(x) = sum_i (I_model_i - I_obs_i)^2 / (4 * sum_i I_obs_i),
    I_model = | sum_n x_n F_n |^2.

  i_obs and F are arrays with data() and select() (cctbx.miller.array).
  """

  def __init__(self, x, i_obs, F, use_curvatures,
               bound_flags=None, lower_bound=None, upper_bound=None):
    self.x = x
    self.i_obs = i_obs
    self.F = F
    self.t = None
    self.g = None
    self.d = None
    self.bound_flags=bound_flags
    self.lower_bound=lower_bound
    self.upper_bound=upper_bound
    # Sum from small to large values to reduce loss of precision
    s = flex.sort_permutation(self.i_obs.data())
    self.i_obs = self.i_obs.select(s)
    self.F = [f.select(s) for f in self.F]
    #
    self.sum_i_obs = flex.sum(self.i_obs.data())
    self.use_curvatures=use_curvatures
    self.tgo = mosaic_ext.alg2_tg(
      F     = [f.data() for f in self.F],
      i_obs = self.i_obs.data())
    self.update_target_and_grads(x=x)

  def update(self, x):
    self.update_target_and_grads(x = x)

  def update_target_and_grads(self, x):
    self.x = x
    self.tgo.update(self.x)
    self.t = self.tgo.target()
    self.g = self.tgo.gradient()
    if self.use_curvatures:
      d = flex.double()
      for j in range(len(self.F)):
        tmp1 = flex.double(self.i_obs.data().size(),0)
        tmp2 = flex.double(self.i_obs.data().size(),0)
        for m, km in enumerate(self.x):
          zz = flex.real( self.F[j].data()*flex.conj(self.F[m].data()) )
          tmp1 += km * zz
          tmp2 += zz
        d.append(flex.sum(tmp1*tmp1 + tmp2))
      self.d=d

  def target(self): return self.t

  def gradients(self): return self.g

  def gradient(self): return self.gradients()

  def curvatures(self): return self.d/self.sum_i_obs

def algorithm_2(i_obs, F, x, use_curvatures=False, use_lbfgsb=False,
                macro_cycles=10, max_iterations=100, lower=None, upper=None):
  """
  Unphased one-step search (Alg2 of Afonine et al., 2024): coefficients x_n
  of F_model = sum_n x_n F_n that minimize the residual between I_obs and
  |F_model|^2 (see class tg), by L-BFGS minimization from the start x.

  Parameters
  ----------
  i_obs : cctbx.miller.array
      Observed intensities.
  F : list of cctbx.miller.array
      Complex arrays F_n; the first one is usually the atomic model.
  x : flex.double
      Starting values, one per array of F.
  use_curvatures : bool or None
      True: L-BFGS with curvatures; False: without; None: with curvatures
      first, then without.
  use_lbfgsb : bool
      Use L-BFGS-B in the steps that do not use curvatures.
  macro_cycles : int
  max_iterations : int
  lower, upper : flex.double, optional
      Bounds for x. If given, minimization is done by L-BFGS-B with these
      bounds, without curvatures. If not given, x is kept positive by
      resetting negative values between macro-cycles, and with
      use_curvatures=None the last step of each macro-cycle uses the bounds
      0.9 <= x[0] <= 1.1, 0 <= x[n] <= 12.

  Returns
  -------
  flex.double
  """
  assert use_curvatures in [True, False, None]
  assert use_lbfgsb in [True, False, None]
  n = x.size()
  assert n == len(F)
  if(lower is not None or upper is not None):
    if(lower is None): lower = flex.double(n, -NO_LIMIT)
    if(upper is None): upper = flex.double(n,  NO_LIMIT)
    assert lower.size() == n and upper.size() == n
    # 0: unbounded, 1: lower only, 2: both, 3: upper only
    flags = flex.int()
    for lo, hi in zip(lower, upper):
      assert lo <= hi
      has_lo, has_hi = lo > -NO_LIMIT, hi < NO_LIMIT
      if(has_lo and has_hi): flags.append(2)
      elif(has_lo):          flags.append(1)
      elif(has_hi):          flags.append(3)
      else:                  flags.append(0)
    def clip(x):
      return flex.double([min(hi, max(lo, v))
        for v, lo, hi in zip(x, lower, upper)])
    x = clip(x)
    for it in range(macro_cycles):
      # the minimizer works on the array it is given
      calculator = tg(i_obs = i_obs, F=F, x = x.deep_copy(),
        use_curvatures=False,
        bound_flags=flags, lower_bound=lower, upper_bound=upper)
      m = minimizers.lbfgs(
        mode='lbfgsb', max_iterations=max_iterations, calculator=calculator)
      x_new = clip(m.x)
      converged = flex.max(flex.abs(x_new-x)) < 1.e-6
      x = x_new
      if(converged): break
    return x
  upper_default = flex.double([1.1] + [12]*(n-1))
  lower_default = flex.double([0.9] + [0]*(n-1))
  mode = 'lbfgs'
  if(use_lbfgsb): mode = 'lbfgsb'
  for it in range(macro_cycles):
    x = x.set_selected(x<0,1.e-6)
    if(use_curvatures is True): # Curvatures only
      calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=True)
      m = minimizer2(
        max_iterations=max_iterations, calculator=calculator).run(
          use_curvatures=True)
      x = m.x
      x = x.set_selected(x<0,1.e-6)
    elif(use_curvatures is False): # No curvatures at all
      calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=False)
      m = minimizers.lbfgs(
        mode=mode, max_iterations=max_iterations, calculator=calculator)
      x = m.x
      x = x.set_selected(x<0,1.e-6)
    elif(use_curvatures is None):
      calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=True)
      m = minimizer2(
        max_iterations=max_iterations, calculator=calculator).run(
          use_curvatures=True)
      x = m.x
      x = x.set_selected(x<0,1.e-6)
      calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=False,
        bound_flags=flex.int(x.size(),2), lower_bound=lower_default,
        upper_bound=upper_default)
      m = minimizers.lbfgs(
        mode=mode, max_iterations=max_iterations, calculator=calculator)
      x = m.x
      x = x.set_selected(x<0,1.e-6)
  return m.x

def algorithm_3(i_obs, fc, f_masks):
  """
  Unphased two-step search: first the products X_nm = x_n*x_m are found
  from a system of linear equations, then x_n are derived from them.

  Parameters
  ----------
  i_obs : cctbx.miller.array
      Observed intensities.
  fc : cctbx.miller.array
      Atomic model structure factors.
  f_masks : list of cctbx.miller.array

  Returns
  -------
  list of float
      Coefficients of [fc]+f_masks.
  """
  F = [fc]+f_masks
  Gnm = []
  cs = {}
  cntr=0
  nm=[]
  # Compute and store Gnm
  for n, Fn in enumerate(F):
    for m, Fm in enumerate(F):
      if m < n:
        continue
      Gnm.append( flex.real( Fn.data()*flex.conj(Fm.data()) ) )
      cs[(n,m)] = cntr
      cntr+=1
      nm.append((n,m))
  # Keep track of indices for "upper triangular matrix vs full"
  for k,v in zip(list(cs.keys()), list(cs.values())):
    i,j=k
    if i==j: continue
    else: cs[(j,i)]=v
  # Generate and solve system Ax=b, x = A_1*b
  A = []
  b = []
  for u, Gnm_u in enumerate(Gnm):
    for v, Gnm_v in enumerate(Gnm):
      scale = 2
      n,m=nm[v]
      if n==m: scale=1
      A.append( flex.sum(Gnm_u*Gnm_v)*scale )
    b.append( flex.sum(Gnm_u * i_obs.data()) )
  A = matrix.sqr(A)
  A_1 = A.inverse()
  b = matrix.col(b)
  x = A_1 * b
  # Expand Xmn from solution x
  Xmn = []
  for n, Fn in enumerate(F):
    rows = []
    for m, Fm in enumerate(F):
      x_ = x[cs[(n,m)]]
      rows.append(x_)
    Xmn.append(rows)
  # x_n from X_nm
  lnK = []
  for j, Fj in enumerate(F):
    t1 = flex.sum( flex.log( flex.double(Xmn[j]) ) )
    t2 = 0
    for n, Fn in enumerate(F):
      for m, Fm in enumerate(F):
        t2 += math.log(Xmn[n][m])
    t2 = t2 / (2*len(F))
    lnK.append( 1/len(F)*(t1-t2) )
  return [math.exp(x) for x in lnK]

# ------------------------------------------------------------------------------
# Mosaic model for an fmodel
# ------------------------------------------------------------------------------

def default_params():
  """
  Default settings: the mosaic scope of
  mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params.
  """
  from mmtbx.bulk_solvent import bulk_solvent_and_scaling as bss
  return bss.master_params.extract().mosaic

def mask_xray_structure(fmodel):
  """
  Atoms that define the bulk-solvent mask of fmodel.
  """
  mm = fmodel.mask_manager
  if(mm is not None and mm.xray_structure is not None):
    return mm.xray_structure
  xrs = fmodel.xray_structure
  if(xrs is None): return None
  mp = fmodel.mask_params
  selection = flex.bool(xrs.scatterers().size(), True)
  if(mp.ignore_zero_occupancy_atoms):
    selection &= xrs.scatterers().extract_occupancies() > 0
  if(mp.ignore_hydrogens):
    selection &= ~xrs.hd_selection()
  return xrs.select(selection)

def work_selection(fmodel):
  """
  Reflections used to fit occupancies: work reflections; all reflections if
  no reflection is flagged as free.
  """
  flags = fmodel.r_free_flags()
  n = fmodel.f_obs().data().size()
  if(flags is None or flags.data().count(True) == 0):
    return flex.bool(n, True)
  return ~flags.data()

def not_applicable(fmodel):
  """
  Reason (str) why the mosaic model cannot be used with this fmodel, or None
  if it can.
  """
  if(fmodel.twinned()): return "twinned data"
  xrs = fmodel.xray_structure
  if(xrs is None): return "no atomic model"
  if(xrs.scattering_type_registry().last_table() == "neutron"):
    return "neutron data"
  if(len(fmodel.f_masks()) != 1): return "more than one bulk-solvent mask"
  d = fmodel.f_obs().d_spacings().data()
  low = d >= fmodel.mask_params.Fmask_res_high
  n_low_work = (low & work_selection(fmodel)).count(True)
  if(n_low_work < 50):
    return "%d work reflections at d >= %s A, at least 50 needed"%(
      n_low_work, str(fmodel.mask_params.Fmask_res_high))
  if(flex.abs(fmodel.f_masks()[0].data()).all_eq(0)):
    return "no bulk-solvent contribution"
  return None

class result(object):
  """
  Result of run().

  Attributes
  ----------
  correction : flex.complex_double
      sum_n (q_n - 1) F_n over all secondary regions, one value per
      reflection of fmodel.f_obs(). F_whole + correction is the mask term of
      the mosaic model.
  regions : list of libtbx.group_args
      Regions with volume >= volume_min, largest first: volume, occupancy,
      status ("fitted"), centroid_frac, centroid_cart, near (str, set by
      annotate()).
  n_small, volume_small : int, float
      Number and total volume of the emptied small regions.
  n_not_fitted, volume_not_fitted : int, float
      Regions left as in the flat model because there are more regions than
      can be fitted in reasonable time (see _MAX_FITTED_REGIONS).
  n_fitted : int
  volume_min : float
  algorithm : str
  """

  def __init__(self, correction, regions, n_small,
               volume_small, volume_min, algorithm, solvent_fraction,
               n_not_fitted=0, volume_not_fitted=0):
    adopt_init_args(self, locals())
    self.n_fitted = len(regions)

  def n_regions(self):
    return self.n_small + self.n_fitted + self.n_not_fitted

  def mean_occupancy(self):
    if(self.n_fitted == 0): return None
    return sum([r.occupancy for r in self.regions])/self.n_fitted

  def summary(self):
    s = "mosaic: regions %d, small %d, fitted %d"%(self.n_regions(),
      self.n_small, self.n_fitted)
    if(self.n_not_fitted > 0): s += ", not fitted %d"%self.n_not_fitted
    if(self.n_fitted > 0): s += ", <q> %.2f"%self.mean_occupancy()
    return s

  def annotate(self, model):
    """
    Set the near field of each region: the residues of model
    (mmtbx.model.manager) nearest to the region centre.
    """
    sites_frac = model.get_xray_structure().sites_frac()
    if(sites_frac.size() == 0): return self
    uc = model.crystal_symmetry().unit_cell()
    ops = model.crystal_symmetry().space_group().all_ops()
    atoms = model.get_hierarchy().atoms()
    for r in self.regions:
      # symmetry copy of the region centre that is nearest to the model
      best = None
      for op in ops:
        c = op * r.centroid_frac
        d = uc.min_mod_short_distance(sites_frac, c)
        if(best is None or d < best[0]): best = (d, c)
      c = best[1]
      x,y,z = (sites_frac - c).parts()
      x = x - flex.floor(x + 0.5)
      y = y - flex.floor(y + 0.5)
      z = z - flex.floor(z + 0.5)
      dist = uc.orthogonalize(flex.vec3_double(x,y,z)).norms()
      near = []
      for i in flex.sort_permutation(dist):
        ag = atoms[i].parent()
        rg = ag.parent()
        label = "%s %s%s%s"%(rg.parent().id.strip(), ag.resname.strip(),
          rg.resseq.strip(), rg.icode.strip())
        if(not label in [n[0] for n in near]):
          near.append((label, dist[i]))
        if(len(near) == 3): break
      r.near = ", ".join(["%s (%.1f)"%n for n in near])
      r.centroid_cart = uc.orthogonalize(c)
    return self

  def show(self, log, prefix="", max_regions=None):
    """
    Print the summary and the table of regions (the max_regions largest ones
    if given).
    """
    if(log is None): return
    p = prefix
    print(p+"Mosaic bulk-solvent model (%s)"%self.algorithm, file=log)
    print(p+"  isolated regions of the solvent mask: %d"%self.n_regions(),
      file=log)
    print(p+"  smaller than %s A^3, emptied: %d (%.0f A^3 in total)"%(
      ("%.1f"%self.volume_min).rstrip("0").rstrip("."), self.n_small,
      self.volume_small), file=log)
    print(p+"  occupancy fitted: %d"%self.n_fitted, file=log)
    if(self.n_not_fitted > 0):
      print(p+"  too many regions, the smaller ones are not fitted "
        "(occupancy 1): %d (%.0f A^3 in total)"%(self.n_not_fitted,
        self.volume_not_fitted), file=log)
    if(self.n_fitted == 0): return
    print(p+"  Occupancy is relative to the main bulk solvent: 0 - empty, "
"1 - same\n%s  (upper limit)."%p, file=log)
    has_near = len([r for r in self.regions if r.near is not None]) > 0
    head = p+"    #  volume(A^3)  occupancy"
    if(has_near): head += "          nearest residues (distance, A)"
    print(head, file=log)
    shown = self.regions
    if(max_regions is not None): shown = self.regions[:max_regions]
    for i, r in enumerate(shown):
      line = p+"  %3d  %11.0f  %9.2f      "%(i+1, r.volume, r.occupancy)
      if(r.near is not None): line += "    %s"%r.near
      print(line, file=log)
    if(len(shown) < self.n_fitted):
      print(p+"  ... and %d smaller regions, <q> %.2f"%(
        self.n_fitted-len(shown),
        sum([r.occupancy for r in self.regions[len(shown):]])/
          (self.n_fitted-len(shown))), file=log)

# Smallest number of work reflections (with bulk-solvent contribution) per
# fitted occupancy.
_REFLECTIONS_PER_REGION = 20

# Limits on the number of regions with fitted occupancy. Time grows with the
# square of this number. Memory: one array of structure factors per region,
# for the reflections of the fit; _MAX_FIT_ELEMENTS is the largest total
# number of elements of these arrays (16 bytes each, about 1 GB).
# Structures that reach these limits are rare (virus capsids). The largest
# regions are fitted, the others stay as in the flat bulk-solvent model
# (occupancy 1).
_MAX_FITTED_REGIONS = 300
_MAX_FIT_ELEMENTS = 60000000

# Limits for the two overall scales of fit_occupancies; outside of them the
# scales of the flat model are used as they are.
_SCALE_MIN, _SCALE_MAX = 0.5, 2.0

def _fit_with_scales(f_obs, atoms, rest, G, algorithm, normal_matrix=None):
  """
  a, b and q_n >= 0 of F_model = a*atoms + b*(rest + sum_n q_n G_n).
  Returns (a, b, q), or None if a or b comes out unreasonable.
  normal_matrix: of [atoms, rest] + G (alg4 only).
  """
  n = len(G)
  if(algorithm == "alg4"):
    x = alg4(
      f_obs   = f_obs,
      base    = flex.complex_double(f_obs.size(), 0),
      G       = [atoms, rest] + G,
      lower   = flex.double([-NO_LIMIT, -NO_LIMIT] + [0]*n),
      x_start = flex.double(n+2, 1),
      normal_matrix = normal_matrix)
  else:
    x = algorithm_2(
      i_obs = _array(f_obs*f_obs),
      F     = [_array(atoms), _array(rest)] + [_array(g) for g in G],
      x     = flex.double(n+2, 1),
      lower = flex.double([_SCALE_MIN, _SCALE_MIN] + [0]*n),
      upper = flex.double([_SCALE_MAX, _SCALE_MAX] + [NO_LIMIT]*n))
  if(_SCALE_MIN < x[0] < _SCALE_MAX and _SCALE_MIN < x[1] < _SCALE_MAX):
    return x[0], x[1], x[2:]/x[1]
  return None

def _fit_bounded(f_obs, base, G, q_max, q_start, algorithm,
                 normal_matrix=None):
  """
  0 <= q_n <= q_max of F_model = base + sum_n q_n G_n.
  normal_matrix: of G (alg4 only).
  """
  n = len(G)
  if(algorithm == "alg4"):
    return alg4(f_obs = f_obs, base = base, G = G,
      lower = flex.double(n, 0), upper = flex.double(n, q_max),
      x_start = q_start, normal_matrix = normal_matrix)
  x = flex.double([1])
  x.extend(q_start)
  x = algorithm_2(
    i_obs = _array(f_obs*f_obs),
    F     = [_array(base)] + [_array(g) for g in G],
    x     = x,
    lower = flex.double([1] + [0]*n),
    upper = flex.double([1] + [q_max]*n))
  return x[1:]

def fit_occupancies(f_obs, atoms, rest, G, algorithm="alg4",
                    refine_scales=True, q_max=OCCUPANCY_MAX):
  """
  Occupancies 0 <= q_n <= q_max of regions in

    F_model = a*atoms + b*(rest + sum_n q_n G_n),

  fitted to the amplitudes f_obs.

  atoms, rest and G_n carry the overall and bulk-solvent scales of the flat
  model. Those scales were found with all regions at q=1, so they are
  biased when the regions are far from that; a and b (two numbers for all
  reflections) take this up, which makes q_n relative to the main bulk
  solvent as it is with the regions adjusted. a and b are not used after
  that: all scales are determined again with the fitted q_n.

  If a region comes out above q_max, the occupancies are fitted once more
  within the limits, with a and b kept: these are the scales of the model
  that describes the data best, and they should not change to make up for
  what the limit leaves unexplained in one region.

  Parameters
  ----------
  f_obs : flex.double
  atoms : flex.complex_double
      Atomic model contribution with its overall scale.
  rest : flex.complex_double
      Scaled contribution of the part of the mask that is not fitted (main
      region; small regions excluded).
  G : list of flex.complex_double
      Scaled contributions of the regions.
  algorithm : str
      "alg4" (phased) or "alg2" (unphased).
  refine_scales : bool
      If False a=b=1.
  q_max : float, optional
      None: no upper limit.

  Returns
  -------
  flex.double
  """
  assert algorithm in ["alg4", "alg2"]
  n = len(G)
  if(n == 0): return flex.double()
  if(q_max is None): q_max = NO_LIMIT
  assert q_max > 0
  # The normal matrix of the search (the expensive part) is computed once
  # for all terms; its lower right block is the matrix of the regions alone.
  A, A_regions = None, None
  if(algorithm == "alg4"):
    A = mosaic_ext.normal_matrix(G = [atoms, rest] + G)
    A_regions = A.matrix_copy_block(i_row = 2, i_column = 2, n_rows = n,
      n_columns = n)
  q = None
  if(refine_scales):
    fit = _fit_with_scales(f_obs = f_obs, atoms = atoms, rest = rest, G = G,
      algorithm = algorithm, normal_matrix = A)
    if(fit is not None):
      a, b, q = fit
      if((q > q_max).count(True) > 0):
        q.set_selected(q > q_max, q_max)
        if(algorithm == "alg4"):
          # same as fitting q with the regions scaled by b: b*q is fitted
          # with the regions as they are
          x_max = q_max
          if(q_max < NO_LIMIT): x_max = q_max*b
          x = _fit_bounded(f_obs = f_obs, base = atoms*a + rest*b, G = G,
            q_max = x_max, q_start = q*b, algorithm = algorithm,
            normal_matrix = A_regions)
          at_limit = x >= x_max
          q = x/b
          q.set_selected(at_limit, q_max) # (q_max*b)/b may not be q_max
        else:
          q = _fit_bounded(f_obs = f_obs, base = atoms*a + rest*b,
            G = [g*b for g in G], q_max = q_max, q_start = q,
            algorithm = algorithm)
  if(q is None):
    q = _fit_bounded(f_obs = f_obs, base = atoms + rest, G = G,
      q_max = q_max, q_start = flex.double(n, 1), algorithm = algorithm,
      normal_matrix = A_regions)
  q.set_selected(q < 0, 0)
  q.set_selected(q > q_max, q_max)
  return q

def _same_mask(fmodel, regions, sf):
  """
  True if the regions are the regions of the mask that fmodel uses. The
  solvent fraction is compared if the mask manager of fmodel has it for the
  current mask (no cost), else the structure factors (one FFT).
  """
  mm = fmodel.mask_manager
  fraction = getattr(mm, "solvent_content_via_mask", None)
  if(fraction is not None and mm.xray_structure is not None and
     mm.xray_structure.scatterers().size() ==
     regions._xray_structure.scatterers().size()):
    return abs(fraction - regions.solvent_fraction) < 1.e-6
  return sf.mask_mismatch(f_mask = fmodel.f_masks()[0]) < 1.e-3

def _log(log, message):
  if(log is not None): print(message, file=log)

def run(fmodel, params=None, log=None):
  """
  Mosaic bulk-solvent model for an fmodel that has the scales of the flat
  model (mmtbx.f_model.manager after update_all_scales): isolated regions of
  the mask, and their occupancies fitted to the work reflections.

  fmodel is not changed.

  Parameters
  ----------
  fmodel : mmtbx.f_model.manager
  params : libtbx.phil.scope_extract, optional
      The mosaic scope of
      mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params; defaults if
      None.
  log : file-like, optional

  Returns
  -------
  mmtbx.bulk_solvent.mosaic.result or None
      None if the model is not applicable (see not_applicable) or the mask
      has no isolated regions; the reason is printed to log.
  """
  if(params is None): params = default_params()
  if(params.volume_min is None or params.volume_min < 0):
    from libtbx.utils import Sorry
    raise Sorry("mosaic.volume_min must be zero or positive.")
  reason = not_applicable(fmodel)
  if(reason is None and fmodel.k_masks()[0].all_eq(0)):
    # bulk solvent switched off, or no contribution found by scaling
    reason = "no bulk-solvent contribution"
  if(reason is not None):
    _log(log, "mosaic: skipped (%s)"%reason)
    return None
  global number_of_runs
  number_of_runs += 1
  mp = fmodel.mask_params
  f_obs = fmodel.f_obs()
  regs = regions(
    xray_structure = mask_xray_structure(fmodel),
    mask_params    = mp,
    d_min          = f_obs.d_min())
  if(len(regs.secondary) == 0):
    _log(log, "mosaic: skipped (no isolated regions in the solvent mask)")
    return None
  sf = structure_factors(regions = regs, miller_array = f_obs,
    d_min = mp.Fmask_res_high)
  f_whole = fmodel.f_masks()[0]
  if(not _same_mask(fmodel = fmodel, regions = regs, sf = sf)):
    _log(log, "mosaic: skipped (cannot reproduce the bulk-solvent mask)")
    return None
  small = [r for r in regs.secondary if r.volume <  params.volume_min]
  fit   = [r for r in regs.secondary if r.volume >= params.volume_min]
  # Reflections of the fit
  sf_fit = structure_factors(regions = regs, miller_array = f_obs,
    d_min = mp.Fmask_res_high, selection = work_selection(fmodel))
  sel, n_work = sf_fit.selection, sf_fit.n_selected
  # regs.secondary is sorted by size, largest first
  n_max = min(_MAX_FITTED_REGIONS, _MAX_FIT_ELEMENTS//max(1, n_work))
  not_fitted = fit[n_max:]
  fit = fit[:n_max]
  if(n_work < _REFLECTIONS_PER_REGION*len(fit)):
    _log(log, "mosaic: skipped (%d work reflections for %d regions, at "
      "least %d per region needed)"%(n_work, len(fit),
      _REFLECTIONS_PER_REGION))
    return None
  points_small = flex.size_t()
  for r in small: points_small.extend(r.points)
  q = flex.double(len(fit), 1)
  if(len(fit) > 0):
    # Scales of the flat model. Structure factors of the regions are kept
    # for the reflections of the fit only, and only until the fit is done.
    k_total = (fmodel.k_isotropic()*fmodel.k_anisotropic()).select(sel)
    scale = k_total*fmodel.k_masks()[0].select(sel)
    G = [sf_fit.of_points(points = r.points, compact = True)*scale
         for r in fit]
    # Current F_model = atoms + all of the mask (every region at q=1)
    f_model = fmodel.f_model().data().select(sel)
    whole = f_whole.data().select(sel)*scale
    atoms = f_model - whole
    # Mask without the small regions (emptied) and the regions being fitted
    rest = whole - sf_fit.of_points(points = points_small,
      compact = True)*scale
    for g in G: rest = rest - g
    q = fit_occupancies(
      f_obs     = f_obs.data().select(sel),
      atoms     = atoms,
      rest      = rest,
      G         = G,
      algorithm = params.algorithm)
    del G, rest, atoms, whole, f_model
  # sum_n (q_n - 1) F_n for all reflections with a bulk-solvent
  # contribution, in one Fourier transform
  correction = sf.of_map(map_data = regs.mask_of_points(
    points = [points_small] + [r.points for r in fit],
    values = [-1.] + [qn-1. for qn in q],
    padded = True))
  uc = regs.unit_cell
  return result(
    correction       = correction,
    regions          = [group_args(
      volume        = r.volume,
      occupancy     = qn,
      status        = "fitted",
      centroid_frac = r.centroid_frac,
      centroid_cart = uc.orthogonalize(r.centroid_frac),
      near          = None) for r, qn in zip(fit, q)],
    n_small          = len(small),
    volume_small     = sum([r.volume for r in small]),
    n_not_fitted      = len(not_fitted),
    volume_not_fitted = sum([r.volume for r in not_fitted]),
    volume_min       = params.volume_min,
    algorithm        = params.algorithm,
    solvent_fraction = regs.solvent_fraction)

class _array(object):
  """
  Minimal array with data() and select(), as class tg needs.
  """
  def __init__(self, data): self._data = data
  def data(self): return self._data
  def select(self, selection): return _array(self._data.select(selection))
