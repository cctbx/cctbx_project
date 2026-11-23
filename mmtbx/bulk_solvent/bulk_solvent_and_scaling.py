from __future__ import absolute_import, division, print_function
from mmtbx import bulk_solvent
import iotbx.phil
from cctbx.array_family import flex
from cctbx import adptbx
from libtbx import group_args

import boost_adaptbx.boost.python as bp
from six.moves import range
ext = bp.import_ext("mmtbx_f_model_ext")

# Settings of the mosaic bulk-solvent model (mmtbx.bulk_solvent.mosaic), also
# used by programs that have their own mosaic scope.
mosaic_params_str = """\
    volume_min = 50
      .type = float
      .help = Isolated regions smaller than this volume (A^3, all symmetry \
              copies in the unit cell counted together) get no bulk \
              solvent at all. Occupancies between 0 (empty) and 1 (as the \
              main bulk solvent) are fitted for the larger ones. \
              50 A^3 is about the volume of two water molecules, smaller \
              regions cannot hold disordered solvent. Increase to fit \
              fewer regions.
    algorithm = *alg4 alg2
      .type = choice(multi=False)
      .expert_level = 3
      .help = Search for region occupancies. alg4 (default) uses phases of \
              the current model and solves for all regions at once. alg2 \
              uses no phases (minimization against intensities) and is \
              slower. Both use work reflections only and keep occupancies \
              between 0 and 1.
"""

master_params = iotbx.phil.parse("""\
  mode = slow *fast
    .type = choice(multi=False)
    .help = fast: analytical, resolution-bin-wise bulk-solvent and overall \
            anisotropic scales (no parameters other than bulk_solvent and \
            anisotropic_scaling). slow: flat bulk-solvent model k_sol, b_sol \
            and overall anisotropic scale b_cart found by grid search and \
            minimization, controlled by the parameters in the slow scope. \
            Twinned data use a separate protocol that reads neither.
  bulk_solvent = True
    .type = bool
    .help = Account for bulk solvent (both modes).
  anisotropic_scaling = True
    .type = bool
    .help = Determine the overall anisotropic scale (both modes).
  apply_back_trace = True
    .type = bool
    .help = Move the isotropic part of the overall anisotropic scale into \
            atomic B-factors (both modes). Only done if the calling program \
            asks for it, one mask is used and a model is available.
  slow
    .help = Parameters of mode=slow. Not used by mode=fast.
  {
    k_sol_b_sol_grid_search = True
      .type = bool
      .expert_level=2
      .help = Start with a grid search over k_sol and b_sol (one mask only).
    minimization_k_sol_b_sol = True
      .type = bool
      .expert_level=2
      .help = Refine k_sol and b_sol by minimization (together with b_cart \
              if minimization_b_cart=True).
    minimization_b_cart = True
      .type = bool
      .expert_level=2
      .help = Refine b_cart by minimization (together with k_sol and b_sol \
              if minimization_k_sol_b_sol=True).
    symmetry_constraints_on_b_cart = True
      .type = bool
      .expert_level=2
      .help = Constrain b_cart to the crystal symmetry. If False all six \
              components are independent.
    k_sol_max = 0.6
      .type = float
      .expert_level=2
      .help = Upper limit of k_sol in minimization.
    k_sol_min = 0.0
      .type = float
      .expert_level=2
      .help = Lower limit of k_sol in minimization.
    b_sol_max = 300.0
      .type = float
      .expert_level=2
      .help = Upper limit of b_sol in minimization.
    b_sol_min = 0.0
      .type = float
      .expert_level=2
      .help = Lower limit of b_sol in minimization.
    k_sol_grid_search_max = 0.6
      .type = float
      .expert_level=2
      .help = Largest k_sol of the grid search.
    k_sol_grid_search_min = 0.0
      .type = float
      .expert_level=2
      .help = Smallest k_sol of the grid search.
    b_sol_grid_search_max = 80.0
      .type = float
      .expert_level=2
      .help = Largest b_sol of the grid search.
    b_sol_grid_search_min = 10.0
      .type = float
      .expert_level=2
      .help = Smallest b_sol of the grid search.
    k_sol_step = 0.1
      .type = float
      .expert_level=2
      .help = Step of the grid search in k_sol.
    b_sol_step = 5.0
      .type = float
      .expert_level=2
      .help = Step of the grid search in b_sol.
    max_iterations = 100
      .type = int
      .expert_level=3
      .help = Maximal number of iterations of the minimization (k_sol, \
              b_sol and b_cart are refined together, within the limits \
              above).
    fix_k_sol = None
      .type = float
      .expert_level=2
      .help = Use this k_sol and do not determine it. Requires grid search \
              and k_sol, b_sol minimization to be off.
    fix_b_sol = None
      .type = float
      .expert_level=2
      .help = Use this b_sol and do not determine it. Requires grid search \
              and k_sol, b_sol minimization to be off.
    fix_b_cart
      .expert_level=2
      .style = box
      .help = Use this b_cart and do not determine it (all six values are \
              required). Requires b_cart minimization to be off.
    {
      b11 = None
        .type = float
      b22 = None
        .type = float
      b33 = None
        .type = float
      b12 = None
        .type = float
      b13 = None
        .type = float
      b23 = None
        .type = float
    }
  }
  mosaic
    .help = Settings of the mosaic bulk-solvent model (Afonine et al., \
            Protein Sci. 2024, 33, e4909): isolated regions of the solvent \
            mask (pockets, cavities) get individual occupancies relative to \
            the main bulk solvent. These are settings only. Whether the \
            model is used is decided by the program, for example \
            main.mosaic in phenix.refine. X-ray data, one mask, not twinned.
  {
""" + mosaic_params_str + """  }
""")

_fix_b_cart_names = ["b11", "b22", "b33", "b12", "b13", "b23"]

def non_default_slow_parameters(params):
  """
  List parameters of the slow scope that differ from their defaults.

  Parameters
  ----------
  params : libtbx.phil.scope_extract
      Extracted master_params.

  Returns
  -------
  list of str
      Names such as "slow.k_sol_max" or "slow.fix_b_cart.b11", in the order
      of definition.
  """
  defaults = master_params.extract().slow
  result = []
  for name in defaults.__dict__.keys():
    if(name.startswith("_") or name == "fix_b_cart"): continue
    if(getattr(params.slow, name) != getattr(defaults, name)):
      result.append("slow.%s"%name)
  for name in _fix_b_cart_names:
    if(getattr(params.slow.fix_b_cart, name) is not None):
      result.append("slow.fix_b_cart.%s"%name)
  return result

def _sorry(problem, options):
  from libtbx.utils import Sorry
  lines = ["bulk_solvent_and_scale: %s"%problem]
  for i, option in enumerate(options):
    if(i == 0): lines.append("  Either %s"%option)
    else:       lines.append("  or %s"%option)
  raise Sorry("\n".join(lines))

def validate_slow_params(params):
  """
  Stop with Sorry if parameters of mode=slow contradict each other.

  Parameters
  ----------
  params : libtbx.phil.scope_extract
      Extracted master_params.
  """
  s = params.slow
  kb_on = []
  if(s.k_sol_b_sol_grid_search):
    kb_on.append("slow.k_sol_b_sol_grid_search=True")
  if(s.minimization_k_sol_b_sol):
    kb_on.append("slow.minimization_k_sol_b_sol=True")
  for name in ["fix_k_sol", "fix_b_sol"]:
    value = getattr(s, name)
    if(value is None): continue
    what = name[4:]
    if(len(kb_on) > 0):
      _sorry("slow.%s=%s cannot be combined with %s."%(name, str(value),
        " and ".join(kb_on)), [
        "set slow.k_sol_b_sol_grid_search=False and "
        "slow.minimization_k_sol_b_sol=False to keep %s fixed,"%what,
        "remove slow.%s to let %s be determined."%(name, what)])
    if(not params.bulk_solvent):
      _sorry("slow.%s=%s cannot be combined with bulk_solvent=False."%(
        name, str(value)), [
        "set bulk_solvent=True to use the fixed %s,"%what,
        "remove slow.%s to have no bulk-solvent contribution."%name])
  fbc = [getattr(s.fix_b_cart, name) for name in _fix_b_cart_names]
  n_set = len(fbc) - fbc.count(None)
  if(n_set > 0 and n_set < 6):
    missing = ["slow.fix_b_cart.%s"%n for n, v in zip(_fix_b_cart_names, fbc)
      if v is None]
    _sorry("slow.fix_b_cart needs all six values, missing: %s."%(
      ", ".join(missing)), [
      "give all six of b11, b22, b33, b12, b13, b23,",
      "give none of them to let b_cart be determined."])
  if(n_set == 6):
    if(not params.anisotropic_scaling):
      _sorry("slow.fix_b_cart cannot be combined with "
        "anisotropic_scaling=False.", [
        "set anisotropic_scaling=True to use the fixed b_cart,",
        "remove slow.fix_b_cart to have no anisotropic scale."])
    if(s.minimization_b_cart):
      _sorry("slow.fix_b_cart cannot be combined with "
        "slow.minimization_b_cart=True.", [
        "set slow.minimization_b_cart=False to keep b_cart fixed,",
        "remove slow.fix_b_cart to let b_cart be determined."])
  for lo, hi in [("k_sol_min", "k_sol_max"), ("b_sol_min", "b_sol_max"),
                 ("k_sol_grid_search_min", "k_sol_grid_search_max"),
                 ("b_sol_grid_search_min", "b_sol_grid_search_max")]:
    v_lo, v_hi = getattr(s, lo), getattr(s, hi)
    if(v_lo is None or v_hi is None):
      _sorry("slow.%s and slow.%s must both be defined."%(lo, hi), [
        "give both values,", "leave both at their defaults."])
    if(v_lo > v_hi):
      _sorry("slow.%s=%s is larger than slow.%s=%s."%(lo, str(v_lo), hi,
        str(v_hi)), [
        "decrease slow.%s,"%lo, "increase slow.%s."%hi])
  for name in ["k_sol_step", "b_sol_step"]:
    value = getattr(s, name)
    if(value is None or value < 0.001):
      _sorry("slow.%s=%s must be 0.001 or larger."%(name, str(value)), [
        "give a step of 0.001 or more,",
        "set slow.k_sol_b_sol_grid_search=False and leave slow.%s at its "
        "default to skip the grid search."%name])
  if(s.max_iterations is None or s.max_iterations < 1):
    _sorry("slow.max_iterations=%s must be at least 1."%str(
      s.max_iterations), [
      "give a number of 1 or more,",
      "set slow.minimization_k_sol_b_sol=False and "
      "slow.minimization_b_cart=False to skip the minimization."])

def validate_params(params, fast, twinned=False, n_masks=1):
  """
  Stop with Sorry if bulk-solvent and scaling parameters contradict each
  other or the data they are used with.

  Parameters
  ----------
  params : libtbx.phil.scope_extract
      Extracted master_params.
  fast : bool
      Protocol that is going to be used: True for mode=fast, False for
      mode=slow.
  twinned : bool
      Data are treated as twinned.
  n_masks : int
      Number of bulk-solvent masks.
  """
  mosaic = getattr(params, "mosaic", None)
  if(mosaic is not None and
     (mosaic.volume_min is None or mosaic.volume_min < 0)):
    from libtbx.utils import Sorry
    raise Sorry("mosaic.volume_min must be zero or positive (it is %s)."%
      str(mosaic.volume_min))
  changed = non_default_slow_parameters(params)
  if(twinned):
    if(params.mode == "slow" or len(changed) > 0):
      given = changed[:]
      if(params.mode == "slow"): given.insert(0, "mode=slow")
      _sorry("twinned data use the twin scaling protocol, which does not "
        "use: %s."%", ".join(given), [
        "remove these parameters (leave mode=fast) to refine against "
        "twinned data,",
        "remove the twin law to use them."])
    return
  if(fast):
    if(len(changed) > 0):
      _sorry("mode=fast does not use the parameters of the slow scope, but "
        "these were changed: %s."%", ".join(changed), [
        "set mode=slow to use them,",
        "remove them to use mode=fast."])
    if(n_masks > 1):
      _sorry("mode=fast works with one bulk-solvent mask, but %d masks are "
        "used (mask.n_radial_shells)."%n_masks, [
        "set mode=slow to use several masks,",
        "set mask.n_radial_shells=1 to use mode=fast."])
    return
  validate_slow_params(params)

def _extract_fix_b_cart(fix_b_cart_scope):
  fbs = fix_b_cart_scope
  b_cart = [fbs.b11,fbs.b22,fbs.b33,fbs.b12,fbs.b13,fbs.b23]
  if(b_cart.count(None) > 0): return None
  else: return b_cart

class bulk_solvent_and_scales(object):

  def __init__(self,
               fmodel_kbu = None,
               params = None,
               log    = None):
    self.fmodel_kbu = fmodel_kbu
    self.params = params
    self.log = log
    if(self.params is None): self.params = master_params.extract()
    validate_slow_params(params = self.params)
    p = self.params.slow
    # What to do. Decided here and kept local: params may be shared by many
    # calls and must not be changed.
    mask_ok = not self.fmodel_kbu.check_f_mask_all_zero()
    n_masks = len(self.fmodel_kbu.f_masks)
    bulk_solvent = self.params.bulk_solvent and mask_ok
    self._anisotropic_scaling = self.params.anisotropic_scaling
    fix_b_cart = _extract_fix_b_cart(fix_b_cart_scope = p.fix_b_cart)
    self._b_cart_is_free = self._anisotropic_scaling and fix_b_cart is None
    grid_search = bulk_solvent and p.k_sol_b_sol_grid_search and n_masks==1
    minimization_k_sol_b_sol = bulk_solvent and p.minimization_k_sol_b_sol
    minimization_b_cart = self._anisotropic_scaling and p.minimization_b_cart
    if(not (self.params.bulk_solvent or self.params.anisotropic_scaling)):
      return
    if(bulk_solvent and p.fix_k_sol is not None):
      self.fmodel_kbu.update(k_sols = [p.fix_k_sol]*n_masks)
    if(bulk_solvent and p.fix_b_sol is not None):
      self.fmodel_kbu.update(b_sols = [p.fix_b_sol]*n_masks)
    if(fix_b_cart is not None):
      self.fmodel_kbu.update(b_cart = fix_b_cart)
    if(grid_search):
      ksol,bsol,b_cart,target = self._ksol_bsol_grid_search()
      self.fmodel_kbu.update(k_sols = ksol, b_sols = bsol, b_cart = b_cart)
      assert abs(target-self.fmodel_kbu.r_factor()) < 1.e-6
    if(minimization_k_sol_b_sol or minimization_b_cart):
      self._minimization(
        refine_kb = minimization_k_sol_b_sol,
        refine_u  = minimization_b_cart)
    # Tiny k_sol means no bulk solvent
    ksols = list(self.k_sols())
    bsols = list(self.b_sols())
    do_update = False
    for ik in range(len(ksols)):
      if(abs(ksols[ik]) < 0.01 and ksols[ik] != 0):
        ksols[ik] = 0.
        do_update = True
    if(ksols.count(0) == len(ksols) and bsols.count(0) != len(bsols)):
      bsols = [0.]*len(ksols)
      do_update = True
    if(do_update):
      self.fmodel_kbu.update(k_sols = ksols, b_sols = bsols)

  def k_sols(self):
    return self.fmodel_kbu.k_sols()

  def b_sols(self):
    return self.fmodel_kbu.b_sols()

  def b_cart(self):
    return self.fmodel_kbu.b_cart()

  def format_scale_matrix(self, m=None, log=None):
    sm = m
    if(sm is None): sm = self.b_cart()
    out = log
    if(sm is None):
      print("  k_anisotropic=1", file=log)
      return
    if(len(sm)<=6):
      print("      b_cart(11,22,33,12,13,23):",\
        ",".join([str("%8.4f"%i).strip() for i in sm]), file=out)

  def u_star(self):
    return self.fmodel_kbu.u_star()

  def _symmetrized_b_cart(self, b_cart):
    """
    b_cart obeying crystal symmetry, if constraints are requested.
    """
    if(not self.params.slow.symmetry_constraints_on_b_cart): return b_cart
    uc = self.fmodel_kbu.f_obs.unit_cell()
    u_star = adptbx.u_cart_as_u_star(uc, adptbx.b_as_u(b_cart))
    u_star = self.fmodel_kbu.f_obs.space_group().average_u_star(
      u_star = u_star)
    return adptbx.u_as_b(adptbx.u_star_as_u_cart(uc, u_star))

  def _within_bounds(self, k_sols, b_sols):
    p = self.params.slow
    k_sols = [max(p.k_sol_min, min(p.k_sol_max, v)) for v in k_sols]
    b_sols = [max(p.b_sol_min, min(p.b_sol_max, v)) for v in b_sols]
    return k_sols, b_sols

  def _is_within_bounds(self, k_sols, b_sols):
    k, b = self._within_bounds(k_sols = k_sols, b_sols = b_sols)
    return list(k) == list(k_sols) and list(b) == list(b_sols)

  def _grid_range(self, x_max, x_min, step, bound_max, bound_min):
    """
    Grid points that are within bounds; the nearest bound if there are none.
    """
    result = [x for x in kb_range(x_max, x_min, step)
      if x >= bound_min and x <= bound_max]
    if(len(result) == 0):
      if(x_max < bound_min): result = [bound_min]
      else:                  result = [bound_max]
    return result

  def _ksol_bsol_grid_search(self):
    kbu = self.fmodel_kbu
    p = self.params.slow
    assert len(kbu.f_masks)==1
    start_r_factor = kbu.r_factor()
    start_ksol     = list(kbu.k_sols())
    start_bsol     = list(kbu.b_sols())
    start_b_cart   = kbu.b_cart()
    k_sols = self._grid_range(p.k_sol_grid_search_max, p.k_sol_grid_search_min,
      p.k_sol_step, p.k_sol_max, p.k_sol_min)
    b_sols = self._grid_range(p.b_sol_grid_search_max, p.b_sol_grid_search_min,
      p.b_sol_step, p.b_sol_max, p.b_sol_min)
    if(self._b_cart_is_free):
      # Overall anisotropic scale is determined for each k_sol, b_sol
      o = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
        f_obs          = kbu.f_obs.data(),
        f_calc         = kbu.f_calc.data(),
        f_mask         = kbu.f_masks[0].data(),
        ss             = kbu.ss,
        k_sol_range    = flex.double(k_sols),
        b_sol_range    = flex.double(b_sols),
        miller_indices = kbu.f_obs.indices(),
        r_ref          = start_r_factor)
    else:
      # Overall anisotropic scale is kept as is
      o = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
        f_obs       = kbu.f_obs.data(),
        f_calc      = kbu.f_calc.data(),
        f_mask      = kbu.f_masks[0].data(),
        k_total     = kbu.k_anisotropic(),
        ss          = kbu.ss,
        k_sol_range = flex.double(k_sols),
        b_sol_range = flex.double(b_sols),
        r_ref       = start_r_factor)
    def restore_start():
      kbu.update(k_sols=start_ksol, b_sols=start_bsol, b_cart=start_b_cart)
      return start_ksol, start_bsol, start_b_cart, kbu.r_factor()
    if(not o.updated()): return restore_start()
    ksol   = [o.k_sol()]
    bsol   = [o.b_sol()]
    b_cart = start_b_cart
    if(self._b_cart_is_free):
      b_cart = self._symmetrized_b_cart(b_cart = adptbx.u_as_b(
        adptbx.u_star_as_u_cart(kbu.f_obs.unit_cell(), o.u_star())))
    kbu.update(k_sols=ksol, b_sols=bsol, b_cart=b_cart)
    r_factor = kbu.r_factor()
    # Possible only if symmetry constraints made the grid result worse
    if(r_factor >= start_r_factor): return restore_start()
    if(ksol[0] < 0.01 and p.k_sol_min <= 0):
      ksol = [0]
      kbu.update(k_sols=ksol)
      r_factor = kbu.r_factor()
    assert abs(kbu.r_factor()-r_factor) < 1.e-6
    assert r_factor <= start_r_factor
    return ksol, bsol, b_cart, r_factor

  def r_factor(self):
    return self.fmodel_kbu.r_factor()

  def r_all(self):
    return self.r_factor()

  def _minimization(self, refine_kb, refine_u):
    """
    Refine k_sols, b_sols (refine_kb) and b_cart (refine_u) together, within
    slow.k_sol_min/max and slow.b_sol_min/max. The result replaces the start
    if it has a lower R-factor, or if the start is outside the bounds.
    """
    import mmtbx.bulk_solvent.kbu_refinery as kbu_refinery
    kbu = self.fmodel_kbu
    p = self.params.slow
    k_start = list(kbu.k_sols())
    b_start = list(kbu.b_sols())
    b_cart_start = kbu.b_cart()
    if(refine_u):
      # Start from b_cart that obeys the constraints, so that the result is
      # compared with (and may be replaced by) nothing else.
      b_cart_start = self._symmetrized_b_cart(b_cart = b_cart_start)
      kbu.update(b_cart = b_cart_start)
    start_is_ok = True
    if(refine_kb):
      start_is_ok = self._is_within_bounds(k_sols=k_start, b_sols=b_start)
    r_start = kbu.r_factor()
    # k_sol=0 is a stationary point for b_sol (no gradient), and b_sol=0
    # next to it a poor place to start from: if nothing is known yet (no
    # grid search), minimization starts from typical values.
    k_first, b_first = k_start, b_start
    if(refine_kb and k_start.count(0) == len(k_start)):
      k_first, b_first = self._within_bounds(
        k_sols = [typical_k_sol]*len(k_start),
        b_sols = [typical_b_sol]*len(b_start))
    obj = kbu_refinery.tgc(
      f_obs   = kbu.f_obs,
      f_calc  = kbu.f_calc,
      f_masks = kbu.f_masks,
      ss      = kbu.ss,
      k_sols  = k_first,
      b_sols  = b_first,
      u_star  = kbu.u_star(),
      symmetry_constraints = p.symmetry_constraints_on_b_cart)
    obj.minimize_joint(
      k_sol_min      = p.k_sol_min,
      k_sol_max      = p.k_sol_max,
      b_sol_min      = p.b_sol_min,
      b_sol_max      = p.b_sol_max,
      refine_kb      = refine_kb,
      refine_u       = refine_u,
      max_iterations = p.max_iterations)
    k_end, b_end, b_cart_end = k_start, b_start, b_cart_start
    if(refine_kb):
      k_end, b_end = list(obj.kbu.k_sols()), list(obj.kbu.b_sols())
    if(refine_u):
      b_cart_end = adptbx.u_as_b(adptbx.u_star_as_u_cart(
        kbu.f_obs.unit_cell(), obj.kbu.u_star()))
    kbu.update(k_sols = k_end, b_sols = b_end, b_cart = b_cart_end)
    r_end = kbu.r_factor()
    if(r_end >= r_start and start_is_ok):
      if(not start_is_ok):
        k_start, b_start = self._within_bounds(k_sols=k_start, b_sols=b_start)
      kbu.update(k_sols = k_start, b_sols = b_start, b_cart = b_cart_start)

  def k_masks(self):
    return self.fmodel_kbu.k_masks()

  def k_anisotropic(self):
    return self.fmodel_kbu.k_anisotropic()

  def k_isotropic(self):
    return self.fmodel_kbu.k_isotropic()

  def apply_back_trace_of_overall_exp_scale_matrix(self, xray_structure=None):
    if(xray_structure is None): return None
    k_sol, b_sol, b_cart = self.k_sols(), self.b_sols(), self.b_cart()
    assert len(k_sol)==1 # XXX Only one mask!
    k_sol = k_sol[0]
    b_sol = b_sol[0]
    #
    xrs = xray_structure
    if(xrs is None): return
    b_min = min(b_sol, xrs.min_u_cart_eigenvalue()*adptbx.u_as_b(1.))
    if(b_min < 0): xrs.tidy_us()
    b_iso = (b_cart[0]+b_cart[1]+b_cart[2])/3.0
    b_test = b_min+b_iso
    if(b_test < 0.0): b_adj = b_iso + abs(b_test) + 0.001
    else: b_adj = b_iso
    b_cart_new = [b_cart[0]-b_adj,b_cart[1]-b_adj,b_cart[2]-b_adj,
                  b_cart[3],      b_cart[4],      b_cart[5]]
    b_sol_new = b_sol + b_adj
    xrs.shift_us(b_shift = b_adj)
    b_min = min(b_sol_new, xrs.min_u_cart_eigenvalue()*adptbx.u_as_b(1.))
    assert b_min >= 0.0
    xrs.tidy_us()
    #
    assert self.fmodel_kbu
    k_masks = [ext.k_mask(self.fmodel_kbu.ss, k_sol, b_sol_new)]
    u_star=adptbx.u_cart_as_u_star(
      self.fmodel_kbu.f_obs.unit_cell(), adptbx.b_as_u(b_cart_new))
    k_anisotropic = ext.k_anisotropic(self.fmodel_kbu.f_obs.indices(), u_star)
    self.fmodel_kbu = self.fmodel_kbu.update(
      b_sols = [b_sol_new],
      b_cart = b_cart_new)
    return group_args(
      xray_structure = xrs,
      b_adj          = b_adj,
      b_sol          = b_sol_new,
      b_cart         = b_cart_new)

# Mean values over the PDB (Fokine & Urzhumtsev, 2002, Acta Cryst. D58,
# 1387-1392): where minimization of k_sol and b_sol starts from if nothing
# better is known.
typical_k_sol = 0.35
typical_b_sol = 46.

def kb_range(x_max, x_min, step):
  sc = 1000.
  return [i/sc for i in range(int(x_min*sc), int(x_max*sc)+1, int(step*sc))]

class ksol_bsol_result(object):
  def __init__(self, k_sol, b_sol, b_cart, r_factor):
    self.k_sol = k_sol
    self.b_sol = b_sol
    self.b_cart = b_cart
    self.r_factor = r_factor
