from __future__ import division, print_function
from libtbx import adopt_init_args
from scitbx import minimizers
import mmtbx.refinement.data
from mmtbx.refinement import calculators
from cctbx import adptbx
from libtbx.utils import user_plus_sys_time
from libtbx.str_utils import format_value
import sys
import libtbx.log
from scitbx.array_family import flex
import mmtbx.refinement.restraints
import mmtbx.refinement.weights
from libtbx import group_args
import mmtbx.refinement.refinement_flags

class occupancy_and_u_iso_calculator(object):
  """
  Target and gradients for simultaneous refinement of occupancies and
  isotropic ADP of selected atoms; for scitbx.minimizers.lbfgs (lbfgsb).

  Structure factors are a sum over atoms: those of atoms that are not
  selected are calculated once, each update calculates only selected atoms.
  For this sum to be equal to structure factors calculated from all atoms
  at once, sampling of atoms in FFT calculations is fixed by the structure
  the calculator is made with (as in mmtbx.refinement.occupancies), with
  u_iso of one refined atom set to u_min.

  fmodel and its xray_structure are changed in-place. x: occupancies of
  selected atoms followed by u_iso of atoms selected by selection_u_iso
  (must be a subset of selection). Target is scaled by its starting value
  (L-BFGS-B stopping rules are absolute).
  """
  def __init__(self, fmodel, selection, selection_u_iso, q_min, q_max,
               u_min, u_max, refine_q=True, refine_u_iso=True):
    adopt_init_args(self, locals())
    import cctbx.xray.structure_factors
    assert fmodel.twin_law is None
    assert (selection_u_iso & ~selection).count(True) == 0
    self.i_selection = selection.iselection()
    self.selection_u_iso_part = selection_u_iso.select(selection)
    self.n_q = self.i_selection.size()
    self.n_u = self.selection_u_iso_part.count(True)
    xrs = fmodel.xray_structure
    p = fmodel.sfg_params
    # Sampling is fixed for a structure in which all refined atoms scatter
    # (atoms with zero occupancy are not looked at) and one of them has the
    # smallest u_iso the refinement can get to: an atom with u_iso smaller
    # than the smallest one of that structure cannot be sampled.
    xrs_sampling = xrs.deep_copy_scatterers()
    scatterers_sampling = xrs_sampling.scatterers()
    for i in self.i_selection:
      if(scatterers_sampling[i].occupancy <= 0):
        scatterers_sampling[i].occupancy = 1
    if(self.n_u > 0):
      u_iso = xrs_sampling.extract_u_iso_or_u_equiv()
      u_sel = u_iso.select(selection_u_iso)
      smallest = flex.min(u_sel)
      if(refine_u_iso): smallest = min(u_min, smallest)
      i_smallest = selection_u_iso.iselection()[flex.min_index(u_sel)]
      scatterers_sampling[i_smallest].u_iso = max(smallest, 0)
    self.sf_managers = []
    for manager_type, miller_set in [
          (cctbx.xray.structure_factors.from_scatterers, fmodel.f_obs()),
          (cctbx.xray.structure_factors.gradients, fmodel.f_obs_work())]:
      manager = manager_type(
        miller_set                   = miller_set,
        cos_sin_table                = p.cos_sin_table,
        grid_resolution_factor       = p.grid_resolution_factor,
        quality_factor               = p.quality_factor,
        u_base                       = p.u_base,
        b_base                       = p.b_base,
        wing_cutoff                  = p.wing_cutoff,
        exp_table_one_over_step_size = p.exp_table_one_over_step_size)
      manager.fix_sampling(xray_structure = xrs_sampling)
      self.sf_managers.append(manager)
    xrs_fixed = xrs.deep_copy_scatterers()
    xrs_fixed.set_occupancies(0, self.i_selection)
    self.f_calc_fixed = self._f_calc(xray_structure = xrs_fixed)
    self.xrs_part = xrs.select(self.i_selection)
    scatterers = self.xrs_part.scatterers()
    scatterers.flags_set_grads(state = False)
    scatterers.flags_set_grad_occupancy(
      iselection = flex.size_t_range(self.n_q))
    scatterers.flags_set_grad_u_iso(
      iselection = self.selection_u_iso_part.iselection())
    self.target_functor = fmodel.target_functor()
    self.x = self.xrs_part.scatterers().extract_occupancies()
    u_iso = self.xrs_part.extract_u_iso_or_u_equiv()
    self.x.extend(u_iso.select(self.selection_u_iso_part))
    q_lower, q_upper = flex.double(self.n_q, q_min), flex.double(self.n_q, q_max)
    u_lower, u_upper = flex.double(self.n_u, u_min), flex.double(self.n_u, u_max)
    if(not refine_q):
      q_lower, q_upper = self.x[:self.n_q], self.x[:self.n_q]
    if(not refine_u_iso):
      u_lower, u_upper = self.x[self.n_q:], self.x[self.n_q:]
    self.lower_bound = q_lower.deep_copy()
    self.lower_bound.extend(u_lower)
    self.upper_bound = q_upper.deep_copy()
    self.upper_bound.extend(u_upper)
    for i in range(self.x.size()):
      self.x[i] = min(self.upper_bound[i], max(self.lower_bound[i], self.x[i]))
    self.bound_flags = flex.int(self.x.size(), 2)
    # FFT gradient wrt occupancy vanishes as occupancy goes to zero
    self.occupancy_small = 0.01
    self.n_calls = 0
    self.target_scale = None
    self.t, self.g = None, None
    self.update(x = self.x)

  def _f_calc(self, xray_structure):
    return self.sf_managers[0](
      xray_structure = xray_structure,
      miller_set     = self.fmodel.f_obs(),
      algorithm      = self.fmodel.sfg_params.algorithm).f_calc()

  def _set(self, x):
    q = x[:self.n_q]
    self.xrs_part.scatterers().set_occupancies(q)
    if(self.n_u > 0):
      u_iso = self.xrs_part.extract_u_iso_or_u_equiv()
      u_iso.set_selected(self.selection_u_iso_part, x[self.n_q:])
      self.xrs_part.set_u_iso(values = u_iso,
        selection = self.selection_u_iso_part)

  def update(self, x):
    self.x = x
    fmodel = self.fmodel
    self._set(x = x)
    f_calc = self._f_calc(xray_structure = self.xrs_part)
    fmodel.update(f_calc = f_calc.customized_copy(
      data = f_calc.data() + self.f_calc_fixed.data()))
    tg = self.target_functor(compute_gradients = True)
    t = tg.target_work()
    d_target_d_f_calc = tg.d_target_d_f_calc_work()
    q = x[:self.n_q]
    small = q < self.occupancy_small
    if(small.count(True) > 0):
      self.xrs_part.scatterers().set_occupancies(
        q.deep_copy().set_selected(small, self.occupancy_small))
    g = self.sf_managers[1](
      xray_structure         = self.xrs_part,
      u_iso_refinable_params = None,
      miller_set             = d_target_d_f_calc,
      d_target_d_f_calc      = d_target_d_f_calc.data(),
      n_parameters           = 0,
      algorithm              = fmodel.sfg_params.algorithm)
    self.g = g.d_target_d_occupancy()
    self.g.extend(g.d_target_d_u_iso().select(self.selection_u_iso_part))
    if(small.count(True) > 0):
      self.xrs_part.scatterers().set_occupancies(q)
    if(self.target_scale is None):
      self.target_scale = 1
      if(t > 0): self.target_scale = 1/t
    self.t = t*self.target_scale
    self.g = self.g*self.target_scale
    self.n_calls += 1

  def target(self): return self.t

  def gradients(self): return self.g

  def finalize(self):
    """
    Put refined parameters into fmodel.xray_structure and make fmodel
    consistent with it the regular way.
    """
    self._set(x = self.x)
    xrs = self.fmodel.xray_structure
    occ = xrs.scatterers().extract_occupancies()
    occ.set_selected(self.i_selection,
      self.xrs_part.scatterers().extract_occupancies())
    xrs.scatterers().set_occupancies(occ)
    if(self.n_u > 0):
      u_iso = xrs.extract_u_iso_or_u_equiv()
      u_iso.set_selected(self.i_selection,
        self.xrs_part.extract_u_iso_or_u_equiv())
      xrs.set_u_iso(values = u_iso, selection = self.selection_u_iso)
    self.fmodel.update_xray_structure(update_f_calc = True)

class unrestrained_qbr_fsr(object):
  """
  Unrestrained sequential occupancy (q) and isotropic ADP (b) refinement.
  Use example: refinement of ocupancy and B-factors of newly placed water only.
  Input fmodel is changed in-place.

  selection: atoms to refine. selection_b: atoms to refine B of (default:
  selection); anisotropic atoms must not be in it.

  freeze_unselected=True: structure factors of atoms that are not refined
  are computed once and stay fixed, each iteration computes only the refined
  atoms, and occupancies and B are refined simultaneously in one
  minimization (max_iterations_joint iterations at most), see
  occupancy_and_u_iso_calculator. Same result as macro_cycles alternations
  of occupancy and B refinement with several times fewer structure factor
  calculations, each of them faster if a small part of the model is
  refined. Not used (regular calculation is done) if coordinates are
  refined, data are twinned or structure factor algorithm is not fft or
  direct.
  """
  def __init__(self,
               fmodel,
               model,
               selection,
               refine_q = True,
               refine_b = True,
               refine_xyz = True,
               macro_cycles=3,
               max_iterations=50,
               q_min = 0.,
               q_max = 1,
               b_min = 5,
               b_max = 45,
               max_xyz_shift = 0.5,
               selection_b = None,
               freeze_unselected = False,
               max_iterations_joint = 150,
               log   = None):
    adopt_init_args(self, locals())
    if(self.selection_b is None): self.selection_b = self.selection
    if(freeze_unselected and not refine_xyz and fmodel.twin_law is None and
       fmodel.sfg_params.algorithm in ["fft", "direct"] and
       selection.count(True) > 0 and (refine_q or refine_b)):
      self._run_frozen()
    else:
      self._run_all()

  def _run_all(self):
    model = self.model
    rt_old = model.get_refinement_flags()
    rf = mmtbx.refinement.refinement_flags.manager( # This is ugly!!
      individual_sites   = True,
      individual_adp     = True,
      occupancies        = True,
      sites_individual   = flex.bool(model.size(), True),
      adp_individual_iso = flex.bool(model.size(), True),
      s_occupancies      = flex.bool(model.size(), True))
    model.set_refinement_flags(flags = rf)
    self._refine(
      fmodel      = self.fmodel,
      model       = model,
      selection   = self.selection,
      selection_b = self.selection_b)
    # UGLY. Make sure to re-store flags
    model.set_refinement_flags(flags = rt_old)

  def _run_frozen(self):
    fmodel = self.fmodel
    calculator = occupancy_and_u_iso_calculator(
      fmodel          = fmodel,
      selection       = self.selection,
      selection_u_iso = self.selection_b,
      q_min           = self.q_min,
      q_max           = self.q_max,
      u_min           = adptbx.b_as_u(self.b_min),
      u_max           = adptbx.b_as_u(self.b_max),
      refine_q        = self.refine_q,
      refine_u_iso    = self.refine_b)
    minimizers.lbfgs(
      calculator     = calculator,
      max_iterations = self.max_iterations_joint,
      mode           = 'lbfgsb')
    calculator.finalize()
    self.n_calls = calculator.n_calls
    if self.log is not None:
      print("occ+adp: r_work=%6.4f r_free=%6.4f"%(
        fmodel.r_work(), fmodel.r_free()), file = self.log)

  def _refine(self, fmodel, model, selection, selection_b):
    log = self.log
    data = mmtbx.refinement.data.fs(fmodel = fmodel)
    restraints = mmtbx.refinement.restraints.manager(
      model = model, use_target=False)
    refine_b = self.refine_b and selection_b.count(True) > 0
    for micro_cycle in range(self.macro_cycles):
      # Occupancy
      if self.refine_q:
        calculator = calculators.occ(
          data       = data,
          restraints = restraints,
          selection  = selection,
          q_min      = self.q_min,
          q_max      = self.q_max).calculator()
        minimized = minimizers.lbfgs(
          calculator     = calculator,
          max_iterations = self.max_iterations,
          mode           = 'lbfgsb')
        if log is not None:
          print("occ: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)
      # ADP
      if refine_b:
        calculator = calculators.adp(
          data       = data,
          restraints = restraints,
          selection  = selection_b,
          u_min      = adptbx.b_as_u(self.b_min),
          u_max      = adptbx.b_as_u(self.b_max)).calculator()
        minimized = minimizers.lbfgs(
          calculator     = calculator,
          max_iterations = self.max_iterations,
          mode           = 'lbfgsb')
        if log is not None:
          print("adp: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)
      # XYZ
      if self.refine_xyz:
        calculator = calculators.xyz(
          data       = data,
          restraints = restraints,
          selection  = selection,
          max_shift  = self.max_xyz_shift).calculator()
        minimized = minimizers.lbfgs(
          calculator     = calculator,
          max_iterations = self.max_iterations,
          mode           = 'lbfgsb')
        if log is not None:
          print("xyz: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)

class simple_fsr(object):
  def __init__(self,
               model,
               fmodel,
               sites_cart_list,
               refine_xyz        = True,
               refine_adp        = True,
               twin_laws         = [None,],
               macro_cycles      = 3,
               max_iterations    = 25,
               update_all_scales = True,
               update_mask       = True,
               max_xyz_shift     = 5.,
               b_min             = 1,
               b_max             = 200,
               log               = sys.stdout):
    adopt_init_args(self, locals())
    self.ma             = libtbx.log.manager(log = self.log)
    self.total_time     = 0
    self.r_work         = None
    self.r_free         = None
    self.data           = None
    self.u_iso_start    = None
    self.xyz_restraints = False
    self.weights        = None
    self.results        = []
    self._call(msg="Check and setup", func=self._check_and_setup)

  def _check_and_setup(self):
    assert self.model.get_xray_structure() == self.fmodel.xray_structure
    self.data = mmtbx.refinement.data.fs(fmodel = self.fmodel)
    assert len(self.sites_cart_list)>0 and isinstance(
      self.sites_cart_list[0], flex.vec3_double)
    assert len(set([_.size() for _ in self.sites_cart_list])) == 1
    assert self.model.size() == self.sites_cart_list[0].size()
    assert self.model.get_xray_structure().is_similar(self.fmodel.xray_structure)
    self.u_iso_start = \
      self.fmodel.xray_structure.extract_u_iso_or_u_equiv().deep_copy()
    self.restraints = mmtbx.refinement.restraints.manager(model = self.model)
    self._update_r_factors()

  def run(self):
    n = len(self.sites_cart_list)
    for i, sites_cart in enumerate(self.sites_cart_list):
      self._call(
        msg="Refining sites %d out of %d"%(i,n), func=None)
      self._call(msg="set sites ", func=self._set_sites_cart, args=sites_cart)
      self._macrocycle()
      self.results.append(
        group_args(
          sites_cart_start = sites_cart,
          sites_cart_final = self.fmodel.xray_structure.sites_cart(),
          r_work           = self.fmodel.r_work(),
          r_free           = self.fmodel.r_free()
        )
      )
      self.model.set_sites_cart(
        sites_cart=self.fmodel.xray_structure.sites_cart())
      self.model.set_b_iso(values=
        self.fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))

  def _macrocycle(self):
    for mc in range(self.macro_cycles):
      self._call(msg="weights   ", func=self._compute_weights)
      self._call(msg="refine xyz", func=self._refine_xyz)
      self._call(msg="refine adp", func=self._refine_adp)
      self.fmodel.update_all_scales(remove_outliers=False)

  def _call(self, msg, func = None, args=None):
    timer = user_plus_sys_time()
    if(func is not None):
      if args is None: func()
      else:            func(args)
    #
    self.model.set_sites_cart(
      sites_cart=self.fmodel.xray_structure.sites_cart())
    self.model.set_b_iso(values=
      self.fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
    #
    assert self.model.get_xray_structure() == self.fmodel.xray_structure
    self._update_r_factors()
    t = timer.elapsed()
    self.total_time += t
    self.ma.add_and_show( self._format_msg(m=msg, t=t) )

  def _format_msg(self, m, t):
    m1 = "r_work: %6s r_free: %6s"%(format_value("%6.4f",self.r_work),
      format_value("%6.4f",self.r_free))
    m2 = "time: %6.3f total time: %6.3f"%(t, self.total_time)
    l = "%s | %s | %s"%(m, m1, m2)
    return l

  def _compute_weights(self):
    params = mmtbx.refinement.weights.master_params.extract()
    self.weights = mmtbx.refinement.weights.weight(
      fmodel                             = self.fmodel,
      model                              = self.model,
      correct_special_position_tolerance = 1.0,
      target_weights_params              = params,
      macro_cycle                        = 0,
      show_summary                       = False,
      log                                = self.log)

  def _update_r_factors(self):
    self.r_work, self.r_free = self.fmodel.r_work(), self.fmodel.r_free()

  def _set_sites_cart(self, sites_cart):
    self.model.set_sites_cart(sites_cart = sites_cart)
    self.model.set_b_iso(values = self.u_iso_start * adptbx.u_as_b(1))
    self.fmodel.update_xray_structure(
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True,
      update_f_mask  = self.update_mask)
    if(self.update_all_scales):
      self.fmodel.update_all_scales(remove_outliers=False)

  def _refine_xyz(self):
    if not self.refine_xyz: return
    weights = self.weights.xyz_weights_result
    wx = weights.wx * weights.wx_scale
    calculator = calculators.xyz(
      data              = self.data,
      data_weight       = wx,
      restraints_weight = 1, # XXX
      restraints        = self.restraints,
      max_shift         = self.max_xyz_shift).calculator()
    minimized = minimizers.lbfgs(
      calculator     = calculator,
      max_iterations = self.max_iterations,
      mode           = "lbfgsb")

  def _refine_adp(self):
    if not self.refine_adp: return
    weights = self.weights.adp_weights_result
    wx = weights.wx * weights.wx_scale
    calculator = calculators.adp(
      data              = self.data,
      data_weight       = wx,
      restraints_weight = 1, # XXX
      restraints        = self.restraints,
      u_min             = adptbx.b_as_u(self.b_min),
      u_max             = adptbx.b_as_u(self.b_max)).calculator()
    minimized = minimizers.lbfgs(
      calculator     = calculator,
      max_iterations = self.max_iterations,
      mode           = "lbfgsb")
