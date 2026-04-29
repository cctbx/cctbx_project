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

class unrestrained_qbr_fsr(object):
  """
  Unrestrained sequential occupancy (q) and isotropic ADP (b) refinement.
  Use example: refinement of ocupancy and B-factors of newly placed water only.
  Input fmodel is changed in-place.
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
               log   = None):
    adopt_init_args(self, locals())
    data = mmtbx.refinement.data.fs(fmodel = fmodel)
    restraints = mmtbx.refinement.restraints.manager(
      model = model, use_target=False)
    #
    rt_old = model.get_refinement_flags()
    rf = mmtbx.refinement.refinement_flags.manager( # This is ugly!!
      individual_sites   = True,
      individual_adp     = True,
      occupancies        = True,
      sites_individual   = flex.bool(model.size(), True),
      adp_individual_iso = flex.bool(model.size(), True),
      s_occupancies      = flex.bool(model.size(), True))
    model.set_refinement_flags(flags = rf)
    #
    for micro_cycle in range(macro_cycles):
      # Occupancy
      if refine_q:
        calculator = calculators.occ(
          data       = data,
          restraints = restraints,
          selection  = selection,
          q_min      = q_min,
          q_max      = q_max).calculator()
        minimized = minimizers.lbfgs(
          calculator     = calculator,
          max_iterations = max_iterations,
          mode           = 'lbfgsb')
        if log is not None:
          print("occ: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)
      # ADP
      if refine_b:
        calculator = calculators.adp(
          data       = data,
          restraints = restraints,
          selection  = selection,
          u_min      = adptbx.b_as_u(b_min),
          u_max      = adptbx.b_as_u(b_max)).calculator()
        minimized = minimizers.lbfgs(
          calculator     = calculator,
          max_iterations = max_iterations,
          mode           = 'lbfgsb')
        if log is not None:
          print("adp: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)
      # XYZ
      if refine_xyz:
        calculator = calculators.xyz(
          data       = data,
          restraints = restraints,
          selection  = selection,
          max_shift  = max_xyz_shift).calculator()
        minimized = minimizers.lbfgs(
          calculator     = calculator,
          max_iterations = max_iterations,
          mode           = 'lbfgsb')
        if log is not None:
          print("xyz: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)
    # UGLY. Make sure to re-store flags
    model.set_refinement_flags(flags = rt_old)

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
