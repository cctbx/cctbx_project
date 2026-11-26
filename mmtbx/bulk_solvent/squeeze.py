"""
Squeeze: unmodelled density in the solvent region as a contribution to
F_model.

After all interpretable density has been modelled, the mFo-DFc map of a
high-resolution structure still has positive peaks in the solvent region
(partially ordered water, buffer, cryoprotectant, smeared ligands) that
cannot be built as atoms. Here this density is Fourier transformed and added
to the model structure factors:

  F_model = k_total * (F_calc + F_bulk + F_squeeze),
  F_squeeze = k * exp(-B*s^2/4) * FT[ rho ],

where rho is the mFo-DFc map (work reflections only, absolute scale) inside
the solvent mask and above a threshold, zero elsewhere. k and B are the only
fitted parameters (least squares against work reflections). k is larger than
one because only the tops of the peaks are used.

The solvent mask is computed without water: such peaks are often next to
modelled water molecules, and a mask computed with water does not contain
them. The whole model (water included) is used for F_calc and the map.

This is what PLATON SQUEEZE (Spek, A. L. (2015). Acta Cryst. C71, 9-18) does
for small molecules, with a thresholded map and a fitted scale and B.

The contribution describes density that the atomic model does not explain.
It is meant for a complete model (final refinement before deposition) at
high resolution: used earlier it removes from the map what still has to be
built.
"""
from cctbx.array_family import flex
from cctbx import maptbx
import mmtbx.masks

def default_params():
  """
  Default settings: the squeeze scope of
  mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params.
  """
  from mmtbx.bulk_solvent import bulk_solvent_and_scaling as bss
  return bss.master_params.extract().squeeze

def validate_params(params):
  """
  Raise Sorry if the settings (squeeze scope) make no sense.
  """
  from libtbx.utils import Sorry
  p = params
  if(p.threshold is None or p.threshold <= 0):
    raise Sorry("squeeze.threshold must be positive.")
  if(p.number_of_cycles is None or p.number_of_cycles < 1):
    raise Sorry("squeeze.number_of_cycles must be 1 or more.")
  if(p.grid_step_factor is None or p.grid_step_factor <= 0 or
     p.grid_step_max is None or p.grid_step_max <= 0):
    raise Sorry("squeeze.grid_step_factor and squeeze.grid_step_max must be "
      "positive.")
  if(p.k_max is None or p.k_max <= 0):
    raise Sorry("squeeze.k_max must be positive.")
  if(p.b_min is None or p.b_max is None or p.b_min >= p.b_max or
     p.b_start is None or not (p.b_min <= p.b_start <= p.b_max)):
    raise Sorry("squeeze: b_min <= b_start <= b_max is required.")

def work_selection(fmodel):
  """
  Reflections used for the map and the fit: work reflections; all
  reflections if no reflection is flagged as free.
  """
  flags = fmodel.r_free_flags()
  n = fmodel.f_obs().data().size()
  if(flags is None or flags.data().count(True) == 0):
    return flex.bool(n, True)
  return ~flags.data()

def not_applicable(fmodel, params=None):
  """
  Reason (str) why squeeze cannot be used with this fmodel, or None if it
  can.
  """
  if(params is None): params = default_params()
  if(fmodel.twinned()): return "twinned data"
  xrs = fmodel.xray_structure
  if(xrs is None): return "no atomic model"
  if(xrs.scattering_type_registry().last_table() == "neutron"):
    return "neutron data"
  d_min = fmodel.f_obs().d_min()
  if(params.d_min_max is not None and d_min > params.d_min_max+1.e-6):
    return "data resolution %.2f A, %s A or better needed"%(d_min,
      ("%.2f"%params.d_min_max))
  return None

def crystal_gridding(fmodel, params):
  """
  Grid of the map: step is d_min/grid_step_factor, at most grid_step_max.
  """
  f_obs = fmodel.f_obs()
  step = min(params.grid_step_max, f_obs.d_min()/params.grid_step_factor)
  return maptbx.crystal_gridding(
    unit_cell        = f_obs.unit_cell(),
    space_group_info = f_obs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = step)

def difference_map(fmodel, gridding):
  """
  mFo-DFc map of fmodel in the whole unit cell: work reflections only,
  absolute scale (e/A^3, no F000), anisotropy removed.

  Returns
  -------
  flex.double (3d, not padded)
  """
  mc = fmodel.electron_density_map().map_coefficients(
    map_type                   = "mFo-DFc",
    isotropize                 = True,
    exclude_free_r_reflections = False)
  sel = work_selection(fmodel)
  if(mc.indices().size() == sel.size()):
    mc = mc.customized_copy(
      data = mc.data().deep_copy().set_selected(~sel, 0j))
  else: # Bijvoet mates were merged
    mc = mc.common_set(fmodel.f_obs().select(sel).average_bijvoet_mates())
  fft_map = mc.fft_map(
    symmetry_flags   = maptbx.use_space_group_symmetry,
    crystal_gridding = gridding)
  fft_map.apply_volume_scaling()
  return fft_map.real_map_unpadded()

def solvent_mask(fmodel, params, n_real, water_selection=None):
  """
  Solvent mask in the whole unit cell on the grid n_real: flex.double, one in
  the solvent, zero elsewhere. Atoms of water_selection (flex.bool) do not
  define the mask if params.exclude_water_from_mask is True.
  """
  xrs = fmodel.xray_structure
  n = xrs.scatterers().size()
  selection = flex.bool(n, True)
  if(params.exclude_water_from_mask and water_selection is not None):
    assert water_selection.size() == n, [water_selection.size(), n]
    selection &= ~water_selection
  if(params.ignore_hydrogens):
    selection &= ~xrs.hd_selection()
  if(params.ignore_zero_occupancy_atoms):
    selection &= xrs.scatterers().extract_occupancies() > 0
  mp = fmodel.mask_params
  r_solvent, r_shrink = params.solvent_radius, params.shrink_truncation_radius
  if(r_solvent is None): r_solvent = mp.solvent_radius
  if(r_shrink  is None): r_shrink  = mp.shrink_truncation_radius
  mask = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xrs.select(selection),
    p1                       = True,
    for_structure_factors    = False,
    solvent_radius           = r_solvent,
    shrink_truncation_radius = r_shrink,
    n_real                   = n_real,
    in_asu                   = False).mask_data
  maptbx.unpad_in_place(map = mask)
  return mask

class _kb_target(object):
  """
  Least-squares target for scitbx.minimizers.lbfgs (mode lbfgsb):

    T(k,B) = scale * sum (|k_total*(F0 + k*exp(-B*ss)*F)| - F_obs)^2
                   / sum F_obs^2,   ss = s^2/4,

  with variables x = (k, B/b_scale) so that both are of order one.
  """

  b_scale = 10.
  scale = 1.e4

  def __init__(self, f_obs, k_total, f0, f, ss, k, b, k_max, b_min, b_max):
    self.f_obs, self.k_total, self.f0, self.f, self.ss = \
      f_obs, k_total, f0, f, ss
    self.norm = self.scale/flex.sum(f_obs*f_obs)
    self.x = flex.double([k, b/self.b_scale])
    self.bound_flags = flex.int(2, 2)
    self.lower_bound = flex.double([0, b_min/self.b_scale])
    self.upper_bound = flex.double([k_max, b_max/self.b_scale])
    self.update(x = self.x)

  def kb(self):
    return self.x[0], self.x[1]*self.b_scale

  def update(self, x):
    self.x = x
    k, b = self.kb()
    e = flex.exp(-b*self.ss)
    fe = self.f*self.k_total*e                 # d f_model / d k
    f_model = self.f0*self.k_total + fe*k
    a = flex.abs(f_model)
    a.set_selected(a < 1.e-9, 1.e-9)
    d = a - self.f_obs
    self.t = flex.sum(d*d)*self.norm
    # d|f_model|/dp = Re(conj(f_model) * d f_model/dp) / |f_model|
    q = flex.real(flex.conj(f_model)*fe)*(2*self.norm*d/a)
    self.g = flex.double([
      flex.sum(q),
      flex.sum(q*self.ss)*(-k)*self.b_scale])

  def target(self): return self.t

  def gradients(self): return self.g

def fit_kb(f_obs, k_total, f0, f, ss, params, k=1., b=None):
  """
  Scale k and B of the contribution k*exp(-B*ss)*f added to f0, so that
  |k_total*(f0 + k*exp(-B*ss)*f)| matches f_obs (least squares).

  Parameters
  ----------
  f_obs : flex.double
  k_total : flex.double
      Overall scale (isotropic and anisotropic).
  f0 : flex.complex_double
      Model structure factors without the overall scale and without the
      contribution.
  f : flex.complex_double
      The contribution without scale.
  ss : flex.double
      s^2/4 = 1/(4 d^2).
  params : libtbx.phil.scope_extract
      squeeze scope: k_max, b_min, b_max, b_start are used.
  k, b : float, optional
      Start values.

  Returns
  -------
  tuple of float
      k, B
  """
  from scitbx import minimizers
  if(b is None): b = params.b_start
  calculator = _kb_target(f_obs = f_obs, k_total = k_total, f0 = f0, f = f,
    ss = ss, k = k, b = b, k_max = params.k_max, b_min = params.b_min,
    b_max = params.b_max)
  minimizers.lbfgs(mode = "lbfgsb", calculator = calculator,
    max_iterations = 500)
  return calculator.kb()

class result(object):
  """
  Result of run().

  Attributes
  ----------
  f : flex.complex_double
      Fourier transform of the selected density, one value per reflection of
      f_obs of the fmodel, without scale.
  k, b : float
      Scale and B-factor of the contribution.
  threshold : float
      Threshold in sigma of the map.
  sigma : float
      Standard deviation of the mFo-DFc map (e/A^3) the threshold refers to.
  density_max : float
      Largest value of the map in the solvent region (e/A^3).
  volume : float
      Volume of the selected density (A^3, whole unit cell).
  n_electrons : float
      Integral of the selected density (whole unit cell, without scale).
  solvent_fraction : float
      Fraction of the unit cell that is the solvent mask used here.
  water_excluded : bool
      True if water did not define the mask.
  """

  def __init__(self, f, k, b, threshold, sigma, density_max, volume,
               n_electrons, solvent_fraction, water_excluded):
    self.f = f
    self.k, self.b = k, b
    self.threshold = threshold
    self.sigma = sigma
    self.density_max = density_max
    self.volume = volume
    self.n_electrons = n_electrons
    self.solvent_fraction = solvent_fraction
    self.water_excluded = water_excluded

  def scaled(self, ss):
    """
    The contribution k*exp(-B*ss)*f (flex.complex_double), ss = s^2/4 of the
    reflections of f.
    """
    return self.f*(self.k*flex.exp(-self.b*ss))

  def select(self, selection):
    """
    The same for a selection (flex.bool) of the reflections.
    """
    return self._copy(f = self.f.select(selection))

  def customized_copy(self, f=None, k=None, b=None):
    return self._copy(f = f, k = k, b = b)

  def _copy(self, f=None, k=None, b=None):
    if(f is None): f = self.f
    if(k is None): k = self.k
    if(b is None): b = self.b
    return result(f = f, k = k, b = b, threshold = self.threshold,
      sigma = self.sigma, density_max = self.density_max,
      volume = self.volume, n_electrons = self.n_electrons,
      solvent_fraction = self.solvent_fraction,
      water_excluded = self.water_excluded)

  def summary(self):
    return ("squeeze: density above %s sigma (%.3f e/A^3) in the solvent "
      "mask: %.1f A^3, k %.2f, B %.1f"%(("%.2f"%self.threshold).rstrip(
      "0").rstrip("."), self.threshold*self.sigma, self.volume, self.k,
      self.b))

  def show(self, log, prefix=""):
    if(log is None): return
    mask = "solvent mask"
    if(self.water_excluded): mask += " computed without water"
    lines = [
      "Squeeze: unmodelled density in the solvent region added to F_model",
      "  %s: %.1f %% of the unit cell"%(mask, self.solvent_fraction*100),
      "  mFo-DFc map: sigma %.3f e/A^3, largest value in the mask %.3f "
        "e/A^3 (%.1f sigma)"%(self.sigma, self.density_max,
        self.density_max/max(self.sigma, 1.e-9)),
      "  density used: above %.2f sigma, %.1f A^3 in the unit cell, "
        "%.1f e"%(self.threshold, self.volume, self.n_electrons),
      "  scale k: %.3f  B: %.2f A^2"%(self.k, self.b)]
    for line in lines: print(prefix+line, file=log)

def _log(log, message):
  if(log is not None): print(message, file=log)

def refit(fmodel, result, params=None):
  """
  Scale and B of an existing squeeze contribution with the current scales
  of fmodel, which has that contribution in F_model already (see
  mmtbx.f_model.manager.f_squeeze).

  Returns
  -------
  mmtbx.bulk_solvent.squeeze.result
      Copy of result with new k and b.
  """
  if(params is None): params = default_params()
  sel = work_selection(fmodel)
  k_total = (fmodel.k_isotropic()*fmodel.k_anisotropic()).select(sel)
  f0 = fmodel.f_model().data().select(sel)*(1./k_total) - \
       fmodel.f_squeeze().data().select(sel)
  k, b = fit_kb(
    f_obs   = fmodel.f_obs().data().select(sel),
    k_total = k_total,
    f0      = f0,
    f       = result.f.select(sel),
    ss      = fmodel.ss.select(sel),
    params  = params,
    k       = result.k,
    b       = result.b)
  return result.customized_copy(k = k, b = b)

def run(fmodel, params=None, water_selection=None, log=None):
  """
  Squeeze contribution for an fmodel that has its scales and no squeeze
  contribution (mmtbx.f_model.manager after update_all_scales).

  fmodel is not changed.

  Parameters
  ----------
  fmodel : mmtbx.f_model.manager
  params : libtbx.phil.scope_extract, optional
      The squeeze scope of
      mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params; defaults if
      None.
  water_selection : flex.bool, optional
      Water atoms of fmodel.xray_structure. They do not define the solvent
      mask (params.exclude_water_from_mask). If None the mask is computed
      from all atoms.
  log : file-like, optional

  Returns
  -------
  mmtbx.bulk_solvent.squeeze.result or None
      None if not applicable (see not_applicable) or there is no density
      above the threshold in the solvent mask; the reason is printed to log.
  """
  if(params is None): params = default_params()
  validate_params(params = params)
  reason = not_applicable(fmodel = fmodel, params = params)
  if(reason is not None):
    _log(log, "squeeze: skipped (%s)"%reason)
    return None
  gridding = crystal_gridding(fmodel = fmodel, params = params)
  n_real = gridding.n_real()
  map_data = difference_map(fmodel = fmodel, gridding = gridding)
  sigma = map_data.sample_standard_deviation()
  mask = solvent_mask(fmodel = fmodel, params = params, n_real = n_real,
    water_selection = water_selection)
  assert mask.all() == map_data.all()
  solvent_fraction = flex.mean(mask.as_1d())
  # density in the mask, zero elsewhere and below the threshold
  rho = map_data.as_1d()*mask.as_1d()
  del mask
  density_max = flex.max(rho)
  rho.set_selected(rho < params.threshold*sigma, 0)
  n_points = (rho > 0).count(True)
  if(n_points == 0 or sigma <= 0):
    _log(log, "squeeze: skipped (no density above %.2f sigma in the solvent "
      "mask)"%params.threshold)
    return None
  f_obs = fmodel.f_obs()
  volume_per_point = f_obs.unit_cell().volume()/rho.size()
  n_electrons = flex.sum(rho)*volume_per_point
  rho.reshape(map_data.accessor())
  f = f_obs.structure_factors_from_map(
    map            = rho,
    use_scale      = True,
    anomalous_flag = False,
    use_sg         = False).data()
  del rho
  sel = work_selection(fmodel)
  k_total = (fmodel.k_isotropic()*fmodel.k_anisotropic()).select(sel)
  k, b = fit_kb(
    f_obs   = f_obs.data().select(sel),
    k_total = k_total,
    f0      = fmodel.f_model().data().select(sel)*(1./k_total),
    f       = f.select(sel),
    ss      = fmodel.ss.select(sel),
    params  = params)
  if(k <= 0):
    _log(log, "squeeze: skipped (fitted scale is zero)")
    return None
  return result(
    f                = f,
    k                = k,
    b                = b,
    threshold        = params.threshold,
    sigma            = sigma,
    density_max      = density_max,
    volume           = n_points*volume_per_point,
    n_electrons      = n_electrons,
    solvent_fraction = solvent_fraction,
    water_excluded   = bool(params.exclude_water_from_mask and
                            water_selection is not None))
