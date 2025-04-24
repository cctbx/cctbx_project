from __future__ import absolute_import, division, print_function

import mmtbx.f_model.f_model_info
import libtbx.load_env
from cctbx.array_family import flex
import math, sys, os
from cctbx import miller
from cctbx import adptbx
from libtbx import adopt_init_args
from mmtbx import masks
from cctbx import xray
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from mmtbx.scaling.sigmaa_estimation import sigmaa_estimator
from cctbx import miller
import cctbx.xray.structure_factors
from cctbx.array_family import flex
import math
from cctbx import xray
from cctbx import adptbx
import boost_adaptbx.boost.python as bp
import mmtbx
from libtbx.math_utils import iround
from libtbx.utils import user_plus_sys_time, date_and_time, Sorry
from libtbx.str_utils import format_value, show_string
import libtbx.path
from six.moves import cStringIO as StringIO
import iotbx.phil
from mmtbx.scaling import outlier_rejection
from mmtbx.scaling import absolute_scaling
from cctbx import sgtbx
from mmtbx import map_tools
from copy import deepcopy
from libtbx import group_args
import mmtbx.scaling.ta_alpha_beta_calc
import mmtbx.refinement.targets
from libtbx import Auto
import mmtbx.arrays
import scitbx.math
from cctbx import maptbx
from libtbx.test_utils import approx_equal
import libtbx
import mmtbx.bulk_solvent
import six
from six.moves import zip, range

ext = bp.import_ext("mmtbx_f_model_ext")

time_bulk_solvent_and_scale         = 0.0
time_mask                           = 0.0
time_f_calc                         = 0.0
time_alpha_beta                     = 0.0
time_target                         = 0.0
time_gradients_wrt_atomic_parameters = 0.0
time_fmodel_core_data               = 0.0
time_r_factors                      = 0.0
time_phase_errors                   = 0.0
time_foms                           = 0.0
time_show                           = 0.0

def show_times(out = None):
  if(out is None): out = sys.stdout
  total = time_mask                           +\
          time_f_calc                         +\
          time_alpha_beta                     +\
          time_target                         +\
          time_gradients_wrt_atomic_parameters +\
          time_fmodel_core_data               +\
          time_r_factors                      +\
          time_phase_errors                   +\
          time_foms
  print("  Micro-tasks:", file=out)
  print("    mask                            = %-7.2f" % (time_mask), file=out)
  print("    f_calc                          = %-7.2f" % time_f_calc, file=out)
  print("    alpha_beta                      = %-7.2f" % time_alpha_beta, file=out)
  print("    target                          = %-7.2f" % time_target, file=out)
  print("    gradients_wrt_atomic_parameters = %-7.2f" % \
    time_gradients_wrt_atomic_parameters, file=out)
  print("    fmodel                         = %-7.2f" % time_fmodel_core_data, file=out)
  print("    r_factors                      = %-7.2f" % time_r_factors, file=out)
  print("    phase_errors                   = %-7.2f" % time_phase_errors, file=out)
  print("    foms                           = %-7.2f" % time_foms, file=out)
  print("    TOTAL for micro-tasks          = %-7.2f" % total, file=out)
  return total

class arrays(object):

  def __init__(self,
               core,
               f_obs,
               r_free_flags,
               hl_coeffs=None,
               core_twin=None,
               twin_fraction=None):
    adopt_init_args(self, locals())
    self.work_sel = ~self.r_free_flags.data()
    self.free_sel = self.r_free_flags.data()
    if(self.free_sel.all_ne(True)): self.free_sel = ~self.free_sel
    self.d_spacings = self.f_obs.d_spacings().data()
    self.d_spacings_work = self.d_spacings.select(self.work_sel)
    self.d_spacings_free = self.d_spacings.select(self.free_sel)
    self._update_derived_arrays()

  def _update_derived_arrays(self):
    if(self.core_twin is None):
      self.f_model = self.core.f_model
    else:
      self.f_model = self.f_obs.array(data =
        mmtbx.utils.apply_twin_fraction(
          amplitude_data_part_one = flex.abs(self.core.f_model.data()),
          amplitude_data_part_two = flex.abs(self.core_twin.f_model.data()),
          twin_fraction           = self.twin_fraction))
    self.f_model_work = self.f_model.select(self.work_sel)
    self.f_model_free = self.f_model.select(self.free_sel)
    self.k_anisotropic_work = self.core.data.k_anisotropic.select(self.work_sel) #XXX no-twin
    self.k_isotropic_work   = self.core.k_isotropic.select(self.work_sel)
    self.f_obs_work = self.f_obs.select(self.work_sel)
    self.f_obs_free = self.f_obs.select(self.free_sel)

  def update_core(self, core=None, core_twin=None, twin_fraction=None):
    if(core is not None): self.core = core
    if(core_twin is not None): self.core_twin = core_twin
    if(twin_fraction is not None): self.twin_fraction = twin_fraction
    if([core, core_twin].count(None)<2): self._update_derived_arrays()

class manager_mixin(object):

  def target_w(self, alpha_beta=None):
    global time_target
    timer = user_plus_sys_time()
    result = self.target_functor(
      alpha_beta=alpha_beta)(compute_gradients=False).target_work()
    time_target += timer.elapsed()
    return result

  def target_t(self, alpha_beta=None):
    global time_target
    timer = user_plus_sys_time()
    result = self.target_functor(
      alpha_beta=alpha_beta)(compute_gradients=False).target_test()
    time_target += timer.elapsed()
    return result

  def one_time_gradients_wrt_atomic_parameters(self, **keyword_args):
    return self.target_functor()(compute_gradients=True) \
      .gradients_wrt_atomic_parameters(**keyword_args)

sf_and_grads_accuracy_master_params = iotbx.phil.parse("""\
  algorithm = *fft direct taam
    .type = choice
  cos_sin_table = False
    .type = bool
  grid_resolution_factor=1/3.
    .type = float
  quality_factor = None
    .type = float
  u_base = None
    .type = float
  b_base = None
    .type = float
  wing_cutoff = None
    .type = float
  exp_table_one_over_step_size = None
    .type = float
  extra
    .help = Extra parameters for additional algorithms, e.g. TAAM using DiSCaMB
  {
  }
""")

if cctbx.xray.structure_factors.pydiscamb_is_installed:
  sf_and_grads_accuracy_master_params.adopt_scope(
    iotbx.phil.parse(
      "extra  { discamb { include scope pydiscamb.cctbx_interface.pydiscamb_master_params } }",
      process_includes=True,
    )
  )

alpha_beta_master_params = iotbx.phil.parse("""\
  include scope mmtbx.max_lik.maxlik.alpha_beta_params
  sigmaa_estimator
    .expert_level=2
  {
    include scope mmtbx.scaling.sigmaa_estimation.sigmaa_estimator_params
  }
""", process_includes=True)

def _scale_helper(num, den, selection=None, num_num=False):
  if (selection is not None):
    num = num.select(selection)
    den = den.select(selection)
  if (den.size() == 0):
    raise RuntimeError("No data for scale calculation.")
  denom = flex.sum(den*den)
  if (denom == 0):
    raise RuntimeError("Zero denominator in scale calculation.")
  if (num_num): return flex.sum(num*num) / denom
  return flex.sum(num*den) / denom


class manager_kbu(object):
  def __init__(self,
               f_obs,
               f_calc,
               f_masks,
               ss,
               f_part1=None,
               f_part2=None,
               k_sols = [0.],
               b_sols = [0.],
               u_star = [0,0,0,0,0,0]):
    self.f_obs   = f_obs
    self.f_calc  = f_calc
    self.f_masks = f_masks
    self.f_part1 = f_part1
    self.f_part2 = f_part2
    self.ss      = ss
    if(u_star is None): u_star = [0,0,0,0,0,0]
    if not ((type(k_sols) is list) or (k_sols is None)):
      assert (type(k_sols) is float) or \
        (type(k_sols) is int), type(k_sols)
      k_sols = [k_sols]
    assert len(k_sols) >= 1
    #
    if(len(f_masks)>1):
      if(len(k_sols)==1 and abs(k_sols[0])<1.e-6):
        k_sols = [0]*len(f_masks)
        b_sols = [0]*len(f_masks)
    #
    typfm = type(f_masks)
    if( not((typfm is list) or (typfm is None)) ):
      self.f_masks = [self.f_masks]
    assert len(k_sols) == len(self.f_masks), \
        "VALS: %d %d"%(len(k_sols),len(self.f_masks))
    self.fmdata = []
    for fm in self.f_masks:
      self.fmdata.append(fm.data())
      assert self.f_calc.indices().all_eq(fm.indices())
    if(self.f_part1 is None):
      self.f_part1 = self.f_calc.customized_copy(data =
        flex.complex_double(self.f_calc.data().size(), 0))
    if(self.f_part2 is None):
      self.f_part2 = self.f_calc.customized_copy(data =
        flex.complex_double(self.f_calc.data().size(), 0))
    assert len(k_sols)==len(b_sols), [k_sols, b_sols]
    self.data = ext.core(
      f_calc        = self.f_calc.data(),
      shell_f_masks = self.fmdata,
      k_sols        = k_sols,
      b_sols        = b_sols,
      f_part1       = self.f_part1.data(),
      f_part2       = self.f_part2.data(),
      u_star        = u_star,
      hkl           = self.f_calc.indices(),
      ss            = self.ss)
    self.f_model = miller.array(miller_set=self.f_calc, data=self.data.f_model)
    self.uc = self.data.uc

  def update(self, k_sols=None, b_sols=None, b_cart=None, u_star=None,
                   f_calc=None):
    if(k_sols is None): k_sols = self.data.k_sols()
    if(b_sols is None): b_sols = self.data.b_sols()
    if(u_star is None): u_star = self.data.u_star
    if(b_cart is not None):
      u_star = u_star = adptbx.u_cart_as_u_star(
        self.f_obs.unit_cell(),adptbx.b_as_u(b_cart))
    if(type(k_sols) is int or type(k_sols) is float): k_sols=[k_sols]
    if(type(b_sols) is int or type(b_sols) is float): b_sols=[b_sols]
    assert len(k_sols)==len(b_sols), [k_sols, b_sols]
    if(f_calc is not None): self.f_calc = f_calc
    self.__init__(
      f_obs   = self.f_obs,
      f_calc  = self.f_calc,
      f_masks = self.f_masks,
      f_part1 = self.f_part1,
      f_part2 = self.f_part2,
      ss      = self.ss,
      k_sols  = list(k_sols),
      b_sols  = list(b_sols),
      u_star  = u_star)
    return self

  def select(self, selection):
    return manager_kbu(
      f_obs   = self.f_obs.select(selection=selection),
      f_calc  = self.f_calc.select(selection=selection),
      f_masks = [fm.select(selection=selection) for fm in self.f_masks],
      f_part1 = self.f_part1.select(selection=selection),
      f_part2 = self.f_part2.select(selection=selection),
      ss      = self.ss.select(selection),
      k_sols  = list(self.k_sols()),
      b_sols  = list(self.b_sols()),
      u_star  = self.u_star())

  def deep_copy(self):
    return self.select(selection=flex.bool(self.f_obs.data().size(),True))

  def k_sols(self):
    return self.data.k_sols()

  def b_sols(self):
    return self.data.b_sols()

  def u_star(self): return self.data.u_star

  def r_factor(self):
    return mmtbx.bulk_solvent.r_factor(self.f_obs.data(), self.data.f_model)

  def check_f_mask_all_zero(self):
    for fm in self.f_masks:
      if(not flex.abs(fm.data()).all_eq(0)):
        return False
    return True

  def b_cart(self):
    uc = self.f_obs.unit_cell()
    b_cart = adptbx.u_as_b(
      adptbx.u_star_as_u_cart(uc, self.data.u_star))
    return b_cart

  def core_data_work(self):
    return self

  def k_mask(self):
    assert len(self.k_sols()) == 1
    return ext.k_mask(self.ss, self.data.k_sol(0), self.b_sol(0))

  def k_masks(self):
    result = []
    for k_sol, b_sol in zip(self.data.k_sols(), self.data.b_sols()):
      result.append( ext.k_mask(self.ss, k_sol, b_sol) )
    return result

  def k_anisotropic(self):
    return ext.k_anisotropic(self.f_calc.indices(), self.u_star())

  def k_isotropic(self):
    return flex.double(self.f_obs.data().size(),
      mmtbx.bulk_solvent.scale(self.f_obs.data(), self.f_model.data()))

class manager(manager_mixin):

  def __init__(self,
         f_obs                        = None,
         i_obs                        = None,
         r_free_flags                 = None,
         f_mask                       = None,
         f_part1                      = None,
         f_part2                      = None,
         f_calc                       = None,
         abcd                         = None,
         epsilons                     = None,
         sf_and_grads_accuracy_params = None,
         target_name                  = "ml",
         alpha_beta_params            = None,
         xray_structure               = None,
         mask_params                  = None,
         use_f_model_scaled           = False,
         mask_manager                 = None,
         twin_law                     = None,
         twin_fraction                = 0,
         max_number_of_bins           = 30,
         bin_selections = None,
         k_mask=None,
         k_isotropic=None,
         k_anisotropic=None,
         k_sol = None,
         b_sol = None,
         b_cart =None,
         _target_memory               = None,
         n_resolution_bins_output     = None,
         scale_method="combo",
         origin=None,
         data_type=None):
    self._origin = origin
    self._data_type = data_type
    self.russ = None
    if(twin_law is not None): target_name = "twin_lsq_f"
    self.k_sol, self.b_sol, self.b_cart = k_sol, b_sol, b_cart
    self.bin_selections = bin_selections
    self.scale_method = scale_method
    self.arrays = None
    self.alpha_beta_cache = None
    self.twin_law = twin_law
    self.twin_law_str = twin_law
    self.twin_fraction = twin_fraction
    self.twin_set = self._set_twin_set(miller_array = f_obs)
    self.ss = 1./flex.pow2(f_obs.d_spacings().data()) / 4.
    assert [k_sol, b_sol].count(None) in [0,2]
    #XXX check consistency instead XXX assert [k_sol, k_mask].count(None) in [2,1]
    #XXX check consistency instead XXX assert [b_cart, k_anisotropic].count(None) in [2,1]
    if([k_sol, b_sol].count(None)==0):
      k_mask = ext.k_mask(self.ss, k_sol, b_sol)
    if(b_cart is not None):
      u_star = adptbx.u_cart_as_u_star(f_obs.unit_cell(), adptbx.b_as_u(b_cart))
      k_anisotropic = ext.k_anisotropic(f_obs.indices(), u_star)
    if(f_part1 is None): # XXX ???
      f_part1 = f_obs.array(
        data = flex.complex_double(f_obs.data().size(),0+0j))
    if(f_part2 is None): # XXX ???
      f_part2 = f_obs.array(
        data = flex.complex_double(f_obs.data().size(),0+0j))
    if(r_free_flags is None):
      r_free_flags = f_obs.array(
        data = flex.bool(f_obs.data().size(),False))
    if abcd is not None and not abcd.indices().all_eq(f_obs.indices()):
        abcd = abcd.complete_array(
          new_data_value=(0,0,0,0), d_min=f_obs.d_min()).common_set(f_obs)
    self._f_obs = f_obs
    self._i_obs = i_obs
    self._r_free_flags = r_free_flags
    assert type(f_obs) == type(r_free_flags)
    self._hl_coeffs = abcd
    if(sf_and_grads_accuracy_params is None):
      sf_and_grads_accuracy_params = sf_and_grads_accuracy_master_params.extract()
    self.sfg_params = sf_and_grads_accuracy_params
    if(alpha_beta_params is None):
      alpha_beta_params = alpha_beta_master_params.extract()
    self.twin = False
    if(twin_law is not None): self.twin=True
    assert f_obs is not None
    assert f_obs.is_real_array()
    assert f_obs.is_unique_set_under_symmetry()
    assert r_free_flags.is_unique_set_under_symmetry()
    if(twin_law is None):
      # twin mate of mapped-to-asu does not have to obey this!
      assert f_obs.is_in_asu()
      assert r_free_flags.is_in_asu()
    self.alpha_beta_params = alpha_beta_params
    self.xray_structure    = xray_structure
    self.use_f_model_scaled= use_f_model_scaled
    self.max_number_of_bins = max_number_of_bins
    self.n_resolution_bins_output = n_resolution_bins_output
    if(mask_params is not None):
      self.mask_params = mask_params
    else:
      self.mask_params = mmtbx.masks.mask_master_params.extract()
    if(self._data_type == "neutron"):
      self.mask_params.ignore_hydrogens=False
    self.mask_manager = mask_manager
    if(self.mask_manager is None):
      self.mask_manager = masks.manager(
        miller_array      = f_obs,
        miller_array_twin = self.twin_set,
        mask_params       = self.mask_params)
    if(r_free_flags is not None):
      assert r_free_flags.indices().all_eq(f_obs.indices())
    self.uc = f_obs.unit_cell()
    if(self.xray_structure is None):
      assert [f_calc, f_mask].count(None) == 0
      if not type(f_mask) is list:
        f_mask = [f_mask]
      assert f_calc.is_complex_array()
      assert f_calc.data().size() <= f_obs.data().size()
      for fm in f_mask:
        assert fm.is_complex_array()
        assert fm.data().size() <= f_obs.data().size(), \
          [fm.data().size(),f_obs.data().size()]
        assert fm.data().size() == f_calc.data().size()
      f_calc_twin = None
      f_mask_twin = None
    else:
      if(f_calc is None): f_calc = self.compute_f_calc(miller_array=f_obs)
      f_calc_twin = None
      if(self.twin_set is not None):
        f_calc_twin = self.compute_f_calc(miller_array = self.twin_set)
      if(f_mask is None):
        f_mask = self.mask_manager.shell_f_masks(
          xray_structure = self.xray_structure,
          force_update   = True)
      f_mask_twin = self.mask_manager.shell_f_masks_twin()
    self.epsilons = epsilons
    if(self.epsilons is None): self.epsilons = f_obs.epsilons().data().as_double()
    self.update_core(
      f_calc        = f_calc,
      f_mask        = f_mask,
      f_calc_twin   = f_calc_twin,
      f_mask_twin   = f_mask_twin,
      k_mask        = k_mask,
      k_isotropic   = k_isotropic,
      k_anisotropic = k_anisotropic,
      f_part1       = f_part1.common_set(f_calc),  # XXX WHY COMMON_SET ?
      f_part2       = f_part2.common_set(f_calc))  # XXX WHY COMMON_SET ?
    self.set_target_name(target_name=target_name)
    self._target_memory = _target_memory
    self._structure_factor_gradients_w = None
    self._wilson_b = None
    self.k_h = None
    self.b_h = None
    if(self.bin_selections is None):
      self.bin_selections = self.f_obs().log_binning()

  def __getstate__(self):
    self._structure_factor_gradients_w = None
    return (self.__dict__,)

  def __setstate__(self, state):
    assert len(state) == 1
    self.__dict__.update(state[0])

  def origin(self):
    return self._origin

  def data_type(self):
    return self._data_type

  def is_twin_fmodel_manager(self):
    return False

  def _get_structure_factor_gradients_w(self):
    if (self._structure_factor_gradients_w is None):
      self._structure_factor_gradients_w \
        = cctbx.xray.structure_factors.gradients(
         miller_set                   = self.f_obs_work(),
         cos_sin_table                = self.sfg_params.cos_sin_table,
         grid_resolution_factor       = self.sfg_params.grid_resolution_factor,
         quality_factor               = self.sfg_params.quality_factor,
         u_base                       = self.sfg_params.u_base,
         b_base                       = self.sfg_params.b_base,
         wing_cutoff                  = self.sfg_params.wing_cutoff,
         exp_table_one_over_step_size =
                                  self.sfg_params.exp_table_one_over_step_size)
    return self._structure_factor_gradients_w

  structure_factor_gradients_w = property(_get_structure_factor_gradients_w)

  def _set_twin_set(self, miller_array):
    result = None
    if(self.twin_law is not None):
      twin_law_xyz = sgtbx.rt_mx(symbol=self.twin_law, r_den=12, t_den=144)
      twin_law_matrix = twin_law_xyz.as_double_array()[0:9]
      twin_mi = mmtbx.utils.create_twin_mate(
        miller_indices  = miller_array.indices(),
        twin_law_matrix = twin_law_matrix)
      result = miller_array.customized_copy(
        indices = twin_mi,
        crystal_symmetry = miller_array.crystal_symmetry())
    return result

  def compute_f_calc(self, miller_array = None, xray_structure=None):
    xrs = xray_structure
    if(xrs is None): xrs = self.xray_structure
    if(miller_array is None): miller_array = self.f_obs()
    p = self.sfg_params
    if(miller_array.indices().size()==0):
      raise RuntimeError("Empty miller_array.")
    manager = miller_array.structure_factors_from_scatterers(
      xray_structure               = xrs,
      algorithm                    = p.algorithm,
      cos_sin_table                = p.cos_sin_table,
      grid_resolution_factor       = p.grid_resolution_factor,
      quality_factor               = p.quality_factor,
      u_base                       = p.u_base,
      b_base                       = p.b_base,
      wing_cutoff                  = p.wing_cutoff,
      exp_table_one_over_step_size = p.exp_table_one_over_step_size,
      extra_params = p.extra
      )
    m = manager.manager()
    return manager.f_calc()

  def update_twin_fraction(self):
    if(self.twin_set is None):
      assert self.arrays.core_twin is None
      return
    tfb = self.twin_fraction
    r_work = self.r_work()
    twin_fraction = 0.0
    tf_best = tfb
    while twin_fraction <= 1.0:
      r_work_= abs(mmtbx.bulk_solvent.r_factor(
        self.f_obs_work().data(),
        self.arrays.core.data.f_model.select(self.arrays.work_sel),
        self.arrays.core_twin.data.f_model.select(self.arrays.work_sel),
        twin_fraction))
      if(r_work_ < r_work):
        r_work = r_work_
        tf_best = twin_fraction
      twin_fraction += 0.01
    self.update(twin_fraction = tf_best)

  def core_data_work(self): # XXX WHICH one : twin? non-twin?
    return core(fmodel = self.select(self.arrays.work_sel))

  def outlier_selection(self, log = None, use_model=True):
    n_free = self.r_free_flags().data().count(True)
    fo = self.f_obs()
    fo.set_observation_type_xray_amplitude() # XXX WHY ?
    result = outlier_rejection.outlier_manager(
      miller_obs   = fo,
      r_free_flags = self._r_free_flags,
      out          = "silent")
    s1 = result.basic_wilson_outliers().data()
    s2 = result.extreme_wilson_outliers().data()
    s3 = result.beamstop_shadow_outliers().data()
    s4 = None
    if(n_free > 0 and use_model):
      s4 = result.model_based_outliers(f_model = self.f_model()).data()
      result = s1 & s2 & s3 & s4
    else: result = s1 & s2 & s3
    if(log is not None):
      print(file=log)
      print("basic_wilson_outliers    =", s1.count(False), file=log)
      print("extreme_wilson_outliers  =", s2.count(False), file=log)
      print("beamstop_shadow_outliers =", s3.count(False), file=log)
      if(n_free > 0 and s4 is not None):
        print("model_based_outliers     =", s4.count(False), file=log)
      print("total                    =", result.count(False), file=log)
    return result

  def remove_outliers(self, log = None, use_model=True):
    if(log is not None):
      print("Distribution of F-obs values:", file=log)
      show_histogram(data = self.f_obs().data(), n_slots = 10, log = log)
    o_sel = self.outlier_selection(log = log, use_model=use_model)
    if(o_sel.count(False) > 0):
      #
      indices = self.f_obs().indices()
      data = self.f_obs().data()
      sigmas = self.f_obs().sigmas()
      o_sel_neg = (~o_sel).iselection()
      d_spacings = self.f_obs().d_spacings().data()
      if(log is not None):
        print("Discarded reflections:", file=log)
        if(sigmas is not None):
          print("    h    k    l           f_obs       sigma   d_min", file=log)
          fmt = "%5d%5d%5d %15.4f %10.4f %7.4f"
          for si in o_sel_neg:
            i = indices[si]
            print(fmt%(i[0],i[1],i[2], data[si],
              sigmas[si], d_spacings[si]), file=log)
        else:
          print("    h    k    l           f_obs    d_min", file=log)
          fmt = "%5d%5d%5d %15.4f %7.4f"
          for si in o_sel_neg:
            i = indices[si]
            print(fmt%(i[0],i[1],i[2], data[si],
              d_spacings[si]), file=log)
      #
      new = self.select(selection = o_sel, in_place=True,
        deep_copy_xray_structure=False)
      if(log is not None):
        print("\nDistribution of active F-obs values after outliers rejection:", file=log)
        show_histogram(data=new.arrays.f_obs.data(),n_slots=10,log=log)

  def wilson_b(self, force_update = False):
    if(self.xray_structure is not None and (self._wilson_b is None or
       force_update)):
      result = absolute_scaling.ml_iso_absolute_scaling(
        miller_array = self.f_obs(),
        n_residues   = self.xray_structure.scatterers().size()/8).b_wilson
      self._wilson_b = result
    else: result = self._wilson_b
    return result

  def top_largest_f_obs_f_model_differences(self, threshold_percent,
        n_slots=10000):
    # XXX make a method of scitbx
    def build_histogram(data, n_slots, threshold):
      hm = flex.histogram(data = data, n_slots = n_slots)
      lc_1 = hm.data_min()
      size = data.size()
      values = flex.double()
      counts = flex.double()
      s_1 = enumerate(hm.slots())
      for (i_1,n_1) in s_1:
        hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
        count = n_1*100./size
        values.append((lc_1+hc_1)/2)
        counts.append(count)
        lc_1 = hc_1
      new_counts = flex.double()
      for i, c in enumerate(counts):
        new_counts.append(flex.sum(counts[i:]))
      result = None
      for v, c in zip(values, new_counts):
        if(c<threshold and result is None): result = v
      return result
    d = flex.abs(flex.abs(self.f_obs().data()) -
      flex.abs(self.f_model_scaled_with_k1().data()))
    r = build_histogram(data=d, n_slots=n_slots, threshold=threshold_percent)
    return r, d

  def deep_copy(self):
    return self.select()

  def select(self, selection=None, in_place=False, deep_copy_xray_structure=True):
    if(selection is None): selection = flex.bool(self.f_obs().data().size(),True)
    if(self.hl_coeffs() is not None):
      new_hl_coeffs = self.hl_coeffs().select(selection=selection)
    else:
      new_hl_coeffs = None
    if(self.mask_manager is not None):
      new_mask_manager = self.mask_manager.select(selection = selection) # XXX
    else:
      new_mask_manager = None
    xrs = self.xray_structure
    if(xrs is not None and deep_copy_xray_structure):
      xrs = self.xray_structure.deep_copy_scatterers()
    new_f_masks = [fm.select(selection=selection) for fm in self.f_masks()]
    if(self.arrays.core.k_mask is None):
      new_k_masks = None
    else:
      new_k_masks = [km.select(selection) for km in self.arrays.core.k_masks]
    if(self.arrays.core.k_isotropic is None):
      k_isotropic = None
    else:
      k_isotropic = self.arrays.core.k_isotropic.select(selection)
    if(self.arrays.core.k_anisotropic is None):
      k_anisotropic = None
    else:
      k_anisotropic = self.arrays.core.k_anisotropic.select(selection)
    #
    bin_selections = self.bin_selections
    if(bin_selections is not None):
      tmp = []
      for bs in bin_selections:
        bs_ = bs.select(selection)
        tmp.append(bs_)
        if(bs_.count(True)==0):
          tmp = None
          break
      bin_selections = tmp
    #
    result = manager(
      f_obs                        = self.f_obs().select(selection=selection),
      r_free_flags                 = self._r_free_flags.select(selection=selection),
      f_part1                      = self.f_part1().select(selection=selection),
      f_part2                      = self.f_part2().select(selection=selection),
      sf_and_grads_accuracy_params = deepcopy(self.sfg_params),
      target_name                  = self.target_name,
      abcd                         = new_hl_coeffs,
      epsilons                     = self.epsilons.select(selection),
      alpha_beta_params            = deepcopy(self.alpha_beta_params),
      xray_structure               = xrs,
      f_calc                       = self.f_calc().select(selection=selection),
      f_mask                       = new_f_masks,
      mask_params                  = deepcopy(self.mask_params),
      mask_manager                 = new_mask_manager,
      twin_law                     = self.twin_law,
      twin_fraction                = self.twin_fraction,
      max_number_of_bins           = self.max_number_of_bins,
      k_mask                       = new_k_masks,
      k_isotropic                  = k_isotropic,
      k_anisotropic                = k_anisotropic,
      bin_selections               = bin_selections,
      _target_memory               = self._target_memory,
      n_resolution_bins_output     = self.n_resolution_bins_output,
      k_sol                        = self.k_sol,
      b_sol                        = self.b_sol,
      b_cart                       = self.b_cart,
      origin                       = self.origin(),
      data_type                    = self.data_type())
    result.twin = self.twin
    result.twin_law_str = self.twin_law_str
    result.k_h = self.k_h
    result.b_h = self.b_h
    if(in_place): # XXX USE THIS INSTEAD OF ABOVE
      # TODO: six.moves.zip this file
      for k, v in six.iteritems(self.__dict__):
        try: self.__dict__[k] = result.__dict__[k]
        except KeyError:
          raise RuntimeError("Corrupt refinement: ad hoc key", k)
    return result

  def resolution_filter(self,
        d_max=0,
        d_min=0):
    if(d_min is None): d_min = 0
    if(d_max is None): d_max = 0
    return self.select(
      selection=self.f_obs().resolution_filter_selection(d_max=d_max,d_min=d_min))

  def update_xray_structure(self,
                            xray_structure      = None,
                            update_f_calc       = False,
                            update_f_mask       = False,
                            force_update_f_mask = False):
    self.alpha_beta_cache = None
    if(xray_structure is not None):
      self.xray_structure = xray_structure
    if(self.arrays.core is None or
       (self.twin_set is not None)):
      force_update_f_mask = True
    f_calc = None
    f_calc_twin = None
    if(update_f_calc):
      f_calc = self.compute_f_calc()
      if(self.twin_law is not None):
        f_calc_twin = self.compute_f_calc(miller_array = self.twin_set)
    f_masks = None
    f_mask_twin = None
    if(update_f_mask):
      mngmsks = self.mask_manager.shell_f_masks(
        xray_structure = self.xray_structure,
        force_update   = force_update_f_mask)
      curfmsks = self.f_masks()
      if(mngmsks is not None):
        f_masks = mngmsks[:] # copy
      if(mngmsks is not None and curfmsks is not None):
        assert len(curfmsks) == len(mngmsks), [len(curfmsks), len(mngmsks)]
        for i in range(len(curfmsks)):
          if( mngmsks[i].data().size() != curfmsks[i].data().size() ):
            f_masks[i] = mngmsks[i].common_set(curfmsks[i])
      f_mask_twin = self.mask_manager.shell_f_masks_twin()
    if([f_calc, f_masks, f_calc_twin, f_mask_twin].count(None) == 4):
      set_core_flag = False
    else: set_core_flag = True
    if(f_calc is None and self.arrays.core is not None):
      f_calc = self.f_calc()
    if(f_masks is None and self.arrays.core is not None):
      f_masks = self.f_masks()
    if(set_core_flag):
      self.update_core(
        f_calc      = f_calc,
        f_mask      = f_masks,
        f_calc_twin = f_calc_twin,
        f_mask_twin = f_mask_twin)

  def fmodel_kbu(self):
    return manager_kbu(
      f_obs   = self.f_obs(),
      f_calc  = self.f_calc(),
      f_masks = self.f_masks(),
      ss      = self.ss,
      f_part1 = self.f_part1(),
      f_part2 = self.f_part2())

  def fmodel_kbu_twin(self):
    if(self.arrays.core_twin is None):
      assert self.twin_fraction == 0 or self.twin_fraction is None
      return None
    return manager_kbu(
      f_obs   = self.f_obs(),
      f_calc  = self.f_calc_twin(),
      f_masks = self.f_masks_twin(),
      ss      = self.ss,
      f_part1 = self.f_part1_twin(),
      f_part2 = self.f_part2_twin())

  def update_core(self,
                  f_calc        = None,
                  f_mask        = None,
                  f_calc_twin   = None,
                  f_mask_twin   = None,
                  f_part1       = None,
                  f_part2       = None,
                  f_part1_twin  = None,
                  f_part2_twin  = None,
                  k_anisotropic = None,
                  k_isotropic   = None,
                  k_mask        = None):
    core_ = None
    self.alpha_beta_cache = None
    core_twin_ = None
    if(self.arrays is None):
      core_ = mmtbx.arrays.init(
        f_calc        = f_calc,
        f_masks       = f_mask,
        f_part1       = f_part1,
        f_part2       = f_part2,
        k_masks       = k_mask,
        k_isotropic   = k_isotropic,
        k_anisotropic = k_anisotropic)
      if(self.twin_set is not None):
        core_twin_ = mmtbx.arrays.init(
          f_calc        = self.compute_f_calc(miller_array = self.twin_set),
          f_masks       = self.mask_manager.shell_f_masks_twin(),
          k_masks       = k_mask,
          k_isotropic   = k_isotropic,
          k_anisotropic = k_anisotropic,
          f_part1       = f_part1_twin,
          f_part2       = f_part2_twin)
      f_obs = self.f_obs()
      r_free_flags = self._r_free_flags
      hl_coeffs = self._hl_coeffs
      if(f_obs.data().size() != core_.data.f_model.size()):
        f_obs = f_obs.common_set(core_.f_model)
        r_free_flags = r_free_flags.common_set(core_.f_model)
        if(hl_coeffs is not None):
          hl_coeffs = hl_coeffs.common_set(core_.f_model)
      self.arrays = arrays(
        core          = core_,
        core_twin     = core_twin_,
        f_obs         = f_obs,
        r_free_flags  = r_free_flags,
        hl_coeffs     = hl_coeffs,
        twin_fraction = self.twin_fraction)
    else:
      core_ = self.arrays.core.update(
        f_calc        = f_calc,
        f_masks       = f_mask,
        k_masks       = k_mask,
        k_isotropic   = k_isotropic,
        k_anisotropic = k_anisotropic,
        f_part1       = f_part1,
        f_part2       = f_part2)
      if(self.twin_set is not None):
        core_twin_ = self.arrays.core_twin.update(
          f_calc        = f_calc_twin,
          f_masks       = f_mask_twin,
          k_masks       = k_mask,
          k_isotropic   = k_isotropic,
          k_anisotropic = k_anisotropic,
          f_part1       = f_part1_twin,
          f_part2       = f_part2_twin)
      self.arrays.update_core(core = core_, core_twin = core_twin_,
        twin_fraction = self.twin_fraction)

  def show_mask_optimization_statistics(self, prefix="", out=None):
    if(out is None): return
    line_len = len("|-"+"|"+prefix)
    fill_len = 80-line_len-1
    print("|-"+prefix+"-"*(fill_len)+"|", file=out)
    print("| Solvent (probe) radius= %4.2f Shrink truncation radius= %4.2f%s|"%(
      self.mask_params.solvent_radius,
      self.mask_params.shrink_truncation_radius," "*17), file=out)
    print("| all data:                         500 lowest resolution reflections:        |", file=out)
    fmtl = "| r_work= %s r_free= %s     r_work= %s  "%(
      format_value(format="%6.4f", value=self.r_work()).strip(),
      format_value(format="%6.4f", value=self.r_free()).strip(),
      format_value(format="%6.4f", value=self.r_work_low()).strip())
    pad = " "*(78-len(fmtl))
    print(fmtl + pad + "|", file=out)
    print("|"+"-"*77+"|", file=out)
    print(file=out)

  def optimize_mask(self, params = None, thorough=False, out = None,
      nproc=1):
    if(len(self.f_masks()) != 1):
      return False
    if (nproc is None) : nproc = 1
    rw_ = self.r_work()
    rf_ = self.r_free()
    rwl_= self.r_work_low()
    r_solv_   = self.mask_params.solvent_radius
    r_shrink_ = self.mask_params.shrink_truncation_radius
    self.show_mask_optimization_statistics(prefix="Mask optimization: start",
      out = out)
    trial_range = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4]
    r_shrinks = []
    r_solvs = []
    for tr1 in trial_range:
      for tr2 in trial_range:
        r_shrinks.append(tr1)
        r_solvs.append(tr2)
    r_solv = trial_range[:]
    r_shrink = trial_range[:]
    trial_params = list(zip(r_solvs, r_shrinks))
    parallel = False
    mask_results = []
    if(nproc is Auto) or (nproc > 1):
      parallel = True
      from libtbx import easy_mp
      stdout_and_results = easy_mp.pool_map(
        processes=nproc,
        fixed_func=self.try_mask_params,
        args=trial_params,
        func_wrapper="buffer_stdout_stderr") # XXX safer for phenix GUI
      mask_results = [ r for so, r in stdout_and_results ]
    else :
      for r_solv, r_shrink in trial_params :
        result = self.try_mask_params((r_solv, r_shrink))
        result.show(out=out)
        mask_results.append(result)
    for result in mask_results :
      if (parallel):
        result.show(out=out)
      if(result.r_work < rw_):
        rw_       = result.r_work
        rf_       = result.r_free
        r_solv_   = result.r_solv
        r_shrink_ = result.r_shrink
    if(out is not None): print(file=out)
    self.mask_params.solvent_radius = r_solv_
    self.mask_params.shrink_truncation_radius = r_shrink_
    self.mask_manager.mask_params = self.mask_params
    self.update_xray_structure(
      xray_structure      = self.xray_structure,
      update_f_calc       = False,
      update_f_mask       = True,
      force_update_f_mask = True)
    self.show_mask_optimization_statistics(prefix="Mask optimization: final",
      out = out)
    return True

  # XXX parallel routine
  def try_mask_params(self, args):
    r_solv, r_shrink = args
    self.mask_params.solvent_radius = r_solv
    self.mask_params.shrink_truncation_radius = r_shrink
    self.mask_manager.mask_params = self.mask_params
    self.update_xray_structure(
      xray_structure      = self.xray_structure,
      update_f_calc       = False,
      update_f_mask       = True,
      force_update_f_mask = True)
    return mask_result(
      r_solv=r_solv,
      r_shrink=r_shrink,
      r_work=self.r_work(),
      r_free=self.r_free(),
      r_work_low=self.r_work_low())

  def check_f_mask_all_zero(self):
    for fm in self.f_masks():
      if(not flex.abs(fm.data()).all_eq(0)):
        return False
    return True

  def f_obs_f_model_abs_differences(self):
    return flex.abs(flex.abs(self.f_obs().data()) -
      flex.abs(self.f_model_scaled_with_k1().data()))

  def f_model_full_set(self, n_real, d_min=None, method="sphere",
                       include_f000=True):
    assert method in ["sphere","box"]
    if(method == "sphere"): assert d_min is not None
    assert [self.k_sol, self.b_sol, self.b_cart].count(None) == 0
    assert self.twin_law is None
    u_star = adptbx.u_cart_as_u_star(
      self.f_calc().unit_cell(), adptbx.b_as_u(self.b_cart))
    core = ext.core(
      f_calc        = self.f_calc().data(),
      shell_f_masks = [self.f_masks()[0].data()],
      k_sols        = [self.k_sol],
      b_sol         = self.b_sol,
      f_part1       = self.f_part1().data(),
      f_part2       = self.f_part2().data(),
      u_star        = u_star,
      hkl           = self.f_calc().indices(),
      uc            = self.f_calc().unit_cell(),
      ss            = self.ss)
    k1 = _scale_helper(num=self.f_obs().data(), den=flex.abs(core.f_model))
    # consistency check BEGIN
    fmodel_data = core.f_model*k1
    fm1 = abs(self.f_model_scaled_with_k1()).data()
    fm2 = flex.abs(fmodel_data)
    d = flex.abs(fm1-fm2)
    assert approx_equal(d.min_max_mean().as_tuple(), [0,0,0])
    # consistency check END
    if(method == "box"):
      full_set = miller.structure_factor_box_from_map(
        crystal_symmetry = self.f_calc().crystal_symmetry(),
        n_real           = n_real,
        include_000      = include_f000)
    else:
      full_set = miller.build_set(
        crystal_symmetry = self.f_calc().crystal_symmetry(),
        anomalous_flag   = self.f_obs().anomalous_flag(),
        d_min            = d_min)
    ##
    zero = flex.complex_double(full_set.indices().size(), 0)
    f_calc_full_set = self.compute_f_calc(miller_array = full_set)
    f_mask_full_set = masks.manager(
      miller_array = f_calc_full_set,
      mask_params  = self.mask_params).shell_f_masks(
        xray_structure = self.xray_structure,
        force_update   = True)[0]
    ss_full_set = 1./flex.pow2(f_calc_full_set.d_spacings().data()) / 4.
    core = ext.core(
      f_calc        = f_calc_full_set.data(),
      shell_f_masks = [f_mask_full_set.data()],
      k_sols        = [self.k_sol],
      b_sol         = self.b_sol,
      f_part1       = zero,
      f_part2       = zero,
      u_star        = u_star,
      hkl           = full_set.indices(),
      uc            = full_set.unit_cell(),
      ss            = ss_full_set)
    fmodel_data_full_set = core.f_model*k1
    return miller.array(
      miller_set = full_set,
      data       = fmodel_data_full_set)

  def update_f_part1_all(self, purpose,
                     refinement_neg_cutoff=-2.5):
    zero = flex.complex_double(self.f_calc().data().size(),0)
    zero_a = self.f_calc().customized_copy(data = zero)
    self.update_core(f_part1 = zero_a)
    mp = mmtbx.masks.mask_master_params.extract()
    mp.solvent_radius = mp.solvent_radius/2
    mp.shrink_truncation_radius = mp.shrink_truncation_radius/2
    mmtbx_masks_asu_mask_obj = mmtbx.masks.asu_mask(
      xray_structure = self.xray_structure.expand_to_p1(sites_mod_positive=True),
      d_min          = self.f_obs().d_min(),
      mask_params    = mp)
    bulk_solvent_mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = self.xray_structure.unit_cell(),
      space_group_info      = self.f_obs().space_group_info(),
      pre_determined_n_real = bulk_solvent_mask.focus())
    mc = self.electron_density_map().map_coefficients(
      map_type   = "mFo-DFc",
      isotropize = True,
      exclude_free_r_reflections = False)
    fft_map = mc.fft_map(crystal_gridding = crystal_gridding)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded() # important: not scaled!
    maptbx.truncate_between_min_max(
      map_data = map_data,
      min      = 0,
      max      = 2)
    map_data = map_data*bulk_solvent_mask
    min_map = flex.min(map_data)
    if(min_map == 0): return
    assert min_map < 0
    f_diff = mc.structure_factors_from_map(
      map            = map_data,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    fm = manager(
      f_obs        = self.f_obs(),
      r_free_flags = self.r_free_flags(),
      f_calc       = self.f_model(), # important: not f_model_scaled_with_k1() !
      f_mask       = f_diff)
    fm.update_all_scales(remove_outliers=False, update_f_part1=False)
    rw, rf = fm.r_work(), fm.r_free()
    # put back into self: tricky!
    kt  = self.k_isotropic()*self.k_anisotropic()
    ktp = fm.k_isotropic()*fm.k_anisotropic()
    assert kt.all_ne(0)
    self.update_core(
      k_isotropic   = flex.double(kt.size(),1),
      k_anisotropic = kt*ktp,
      f_part1       = f_diff.customized_copy(
        data = f_diff.data()*(1/kt)*fm.k_masks()[0]))
    # normally it holds up to 1.e-6, so if assertion below breakes then
    # there must be something terribly wrong!
    assert approx_equal(rw, self.r_work(), 1.e-4)
    assert approx_equal(rf, self.r_free(), 1.e-4)

  def update_f_part1(self):
    """
    Identify negative blobs in mFo-DFc synthesis in solvent region only,
    then leave only those blobs (set everything else to zero),
    then FT modified synthesis into f_diff s.f. and add them to Fmodel as scaled
    Fmask contribution. Finally, extract updated scales and scaled f_diff and
    update self with them. Store f_diff in F_part_1 partial contribution array
    of Fmodel, which is an arbitrary choice (could be F_part_2, for instance).
    Rationale: there should not be negative peaks in bulk-solvent region unless
               they are noise or errors of bulk-solvent model, so they are
               always good to remove.
    """
    zero = flex.complex_double(self.f_calc().data().size(),0)
    zero_a = self.f_calc().customized_copy(data = zero)
    self.update_core(f_part1 = zero_a)
    if(self.xray_structure is None): return # need mask
    sgt = self.f_obs().space_group().type()
    d_spacings = self.f_calc().d_spacings().data()
    if(flex.max(d_spacings)>3.0 and (d_spacings>3.0).count(True)>500):
      S_E_L_F = self.resolution_filter(d_min=3.0)
    else:
      S_E_L_F = self.deep_copy()
    crystal_gridding = S_E_L_F.f_obs().crystal_gridding(
      d_min              = S_E_L_F.f_obs().d_min(),
      symmetry_flags     = maptbx.use_space_group_symmetry,
      resolution_factor  = 1./4)
    # bulk-solvent mask filter
    mmtbx_masks_asu_mask_obj = mmtbx.masks.asu_mask(
      xray_structure = S_E_L_F.xray_structure.expand_to_p1(sites_mod_positive=True),
      n_real         = crystal_gridding.n_real())
    asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
    bulk_solvent_mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
    # residual map
    mc = S_E_L_F.electron_density_map().map_coefficients(
      map_type   = "mFo-DFc",
      isotropize = True,
      exclude_free_r_reflections = False)
    #
    exclude_free_r_reflections = True
    r_free_flags_data = S_E_L_F.r_free_flags().data()
    if(exclude_free_r_reflections and
       0 not in [r_free_flags_data.count(True),r_free_flags_data.count(False)]):
      mc = mc.customized_copy(
        data = mc.data().set_selected(r_free_flags_data, 0+0j))
    fft_map = mc.fft_map(
      symmetry_flags   = maptbx.use_space_group_symmetry,
      crystal_gridding = crystal_gridding)
    map_data = fft_map.real_map_unpadded() # important: not scaled!
    # in P1
    map_filtered_inverted = map_data * bulk_solvent_mask *(-1.)
    # in ASU
    asu_map = asu_map_ext.asymmetric_map(sgt, map_filtered_inverted)
    map_data_asu = asu_map.data()
    map_data_asu = map_data_asu.shift_origin()
    # connectivity analysis; work at 0.2 e/A**3 (lower end of solvent density)
    co = maptbx.connectivity(map_data=map_data_asu/mc.unit_cell().volume(),
      threshold=0.2)
    conn = co.result()
    z = zip(co.regions(),range(0,co.regions().size()))
    #z = zip(co.maximum_values(),range(0,co.regions().size()))
    sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
    good = []
    r_free = S_E_L_F.r_free()
    #
    cntr = 0
    for p in sorted_by_volume:
      v, i = p
    #
    #for i, v in enumerate(co.regions()):
    #
      cntr += 1
      if(cntr>50): break
      if(i==0): continue
      if(v>40):
        map_data_asu_ = maptbx.update_f_part1_helper(
          connectivity_map = conn,
          map_data         = map_data_asu,
          region_id        = i)
        asu_map = asu_map_ext.asymmetric_map(sgt, map_data_asu_,
          crystal_gridding.n_real())
        f_diff = mc.customized_copy(
          data = asu_map.structure_factors(mc.indices()))
        fm = manager(
          f_obs        = S_E_L_F.f_obs(),
          r_free_flags = S_E_L_F.r_free_flags(),
          f_calc       = S_E_L_F.f_model(), # not f_model_scaled_with_k1() !
          f_mask       = f_diff)
        fm.update_all_scales(remove_outliers=False, update_f_part1=False)
        if(fm.r_free() < r_free): good.append(i)
    # may be inefficient?
    F = conn.deep_copy()
    if(not 1 in good): F = F.set_selected(F==1, 0)
    for i in good:
      F = F.set_selected(F==i, 1)
    F = F.set_selected(F!=1, 0)
    # filter and invert back to negative
    map_data_asu_ = map_data_asu * (F.as_double() * (-1.))
    asu_map = asu_map_ext.asymmetric_map(sgt, map_data_asu_,
      crystal_gridding.n_real())
    f_diff = mc.customized_copy(
      indices = mc.indices(),
      data    = asu_map.structure_factors(mc.indices()))
    fm = manager(
      f_obs        = S_E_L_F.f_obs(),
      r_free_flags = S_E_L_F.r_free_flags(),
      f_calc       = S_E_L_F.f_model(), # not f_model_scaled_with_k1() !
      f_mask       = f_diff)
    fm.update_all_scales(remove_outliers=False, update_f_part1=False)
    kt  = S_E_L_F.k_isotropic()*S_E_L_F.k_anisotropic()
    ktp = fm.k_isotropic()*fm.k_anisotropic()
    f_part1 = f_diff.customized_copy(data = f_diff.data()*(1/kt)*fm.k_masks()[0])
    # concatenate
    fc = self.f_calc()
    fc_lp = fc.lone_set(f_part1)
    fc_lp = fc_lp.customized_copy(data = fc_lp.data()*0)
    f_part1 = f_part1.concatenate(fc_lp)
    f_part1, fc = f_part1.common_sets(fc)
    # add to self (fmodel)
    self.update(f_part1 = f_part1)

  def bins(self):
    k_masks       = self.k_masks()
    k_isotropic   = self.k_isotropic()
    k_anisotropic = self.k_anisotropic()
    f_model       = self.f_model_scaled_with_k1()
    f_obs         = self.f_obs()
    work_flags    =~self.r_free_flags().data()
    free_flags    = self.r_free_flags().data()
    result = []
    for sel in self.bin_selections:
      d        = self.arrays.d_spacings.select(sel)
      d_min    = flex.min(d)
      d_max    = flex.max(d)
      sel_work = sel & work_flags
      nw       = sel_work.count(True)
      nf       = (sel & free_flags).count(True)
      fo       = f_obs.select(sel)
      fm       = f_model.select(sel)
      fo_mean  = flex.mean(fo.data())
      fm_mean  = flex.mean(flex.abs(fm.data()))
      cmpl     = fo.completeness(d_max=d_max)*100.
      ki       = flex.mean(k_isotropic.select(sel))
      ka       = flex.mean(k_anisotropic.select(sel))
      r        = mmtbx.bulk_solvent.r_factor(
        f_obs.select(sel_work).data(),
        f_model.select(sel_work).data())
      km = " ".join(["%5.3f"%flex.mean(km_.select(sel)) for km_ in k_masks])
      result.append(group_args(
        d_min   = d_min,
        d_max   = d_max,
        nw      = nw,
        nf      = nf,
        fo_mean = fo_mean,
        fm_mean = fm_mean,
        cmpl    = cmpl,
        ki      = ki,
        ka      = ka,
        r       = r,
        km      = km))
    return result

  def show_short(self, show_k_mask=True, log=None, prefix=""):
    if(log is None): log = sys.stdout
    if(show_k_mask):
      fmt="%s%7.3f-%-7.3f %6.2f %5d %5d %6.4f %9.3f %9.3f %5.3f   %s"
      print("%s   Resolution    Compl Nwork Nfree R_work    <Fobs>  <Fmodel> ktotal  kmask"%prefix, file=log)
      for b in self.bins():
        print(fmt % (prefix,b.d_max,b.d_min,b.cmpl,b.nw,b.nf,b.r,b.fo_mean,b.fm_mean,b.ki*b.ka,b.km), file=log)
    else:
      fmt="%s%7.3f-%-7.3f %6.2f %5d %5d %6.4f %9.3f %9.3f %5.3f"
      print("%s   Resolution    Compl Nwork Nfree R_work    <Fobs>  <Fmodel> ktotal"%prefix, file=log)
      for b in self.bins():
        print(fmt % (prefix,b.d_max,b.d_min,b.cmpl,b.nw,b.nf,b.r,b.fo_mean,b.fm_mean,b.ki*b.ka), file=log)

  def show(self, log=None, suffix=None, show_header=True, show_approx=True):
    if(log is None): log = sys.stdout
    l="Statistics in resolution bins"
    if(suffix is not None): l += " %s"%suffix
    f1 = "Total model structure factor:"
    f2 = "  F_model = k_total * (F_calc + k_mask * F_mask)\n"
    f3 = "    k_total = k_isotropic * k_anisotropic"
    m=(77-len(l))//2
    if(show_header):
      print("\n","="*m,l,"="*m,"\n", file=log)
      print(f1, file=log)
      print(f2, file=log)
      print(f3, file=log)
    fmt="%7.3f-%-7.3f %6.2f %5d %5d %6.4f %9.3f %9.3f %5.3f %5.3f %s"
    print("   Resolution    Compl Nwork Nfree R_work    <Fobs>  <Fmodel> kiso   kani kmask", file=log)
    for b in self.bins():
      print(fmt % (b.d_max,b.d_min,b.cmpl,b.nw,b.nf,b.r,b.fo_mean,b.fm_mean,b.ki,b.ka,b.km), file=log)
    def overall_isotropic_kb_estimate(self):
      k_total = self.k_isotropic() * self.k_anisotropic()
      assert self.ss.size() == self.k_isotropic().size()
      r = scitbx.math.gaussian_fit_1d_analytical(x=flex.sqrt(self.ss), y=k_total)
      return r.a, r.b
    k_overall, b_overall = overall_isotropic_kb_estimate(self)
    print(file=log)
    f3="  Approximation of k_total with k_overall*exp(-b_overall*s**2/4)"
    f4="    k_overall=%-8.4f b_overall=%-8.4f"%(k_overall, b_overall)
    if(show_approx):
      print(f3, file=log)
      print(f4, file=log)
    #
    if(self.russ is not None):
      if(self.russ.b_cart is not None and self.twin_law is None):
        b_cart = " ".join(["%6.4f"%b for b in self.russ.b_cart])
        print("  k_sol: %6.4f"%self.russ.k_sol[0], file=log)
        print("  b_sol: %6.4f"%self.russ.b_sol[0], file=log)
        print("  b_cart:", b_cart, file=log)

  def remove_unreliable_atoms_and_update(self, min_map_value=0.5, min_cc=0.7):
    """
    Identify and remove 'unreliably' pleaced atoms, and create a new fmodel
    with updated xray_structure.
    """
    # XXX see map_tools.py: duplication! Consolidate.
    coeffs = map_tools.electron_density_map(
      fmodel=self).map_coefficients(
        map_type         = "2mFo-DFc",
        isotropize       = True,
        fill_missing     = False)
    crystal_gridding = self.f_obs().crystal_gridding(
      d_min              = self.f_obs().d_min(),
      resolution_factor  = 1./3)
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = coeffs)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    rho_atoms = flex.double()
    for site_frac in self.xray_structure.sites_frac():
      rho_atoms.append(map_data.eight_point_interpolation(site_frac))
    #rho_mean = flex.mean_default(
    #  rho_atoms.select(rho_atoms>min_map_value), min_map_value)
    rho_mean = flex.mean_default(rho_atoms, min_map_value)
    sel_exclude = rho_atoms > rho_mean/6.
    ##
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = self.f_model())
    fft_map.apply_sigma_scaling()
    map_data2 = fft_map.real_map_unpadded()
    #
    sites_cart = self.xray_structure.sites_cart()
    #sel_exclude = flex.bool(sites_cart.size(), True)
    for i_seq, site_cart in enumerate(sites_cart):
      selection = maptbx.grid_indices_around_sites(
        unit_cell  = coeffs.unit_cell(),
        fft_n_real = map_data.focus(),
        fft_m_real = map_data.all(),
        sites_cart = flex.vec3_double([site_cart]),
        site_radii = flex.double([1.5]))
      cc = flex.linear_correlation(x=map_data.select(selection),
        y=map_data2.select(selection)).coefficient()
      if(cc<min_cc): sel_exclude[i_seq] = False
    #
    xray_structure_truncated = self.xray_structure.select(sel_exclude)
    fmodel_result = self.deep_copy()
    fmodel_result.update_xray_structure(
      xray_structure = xray_structure_truncated,
      update_f_calc  = True,
      update_f_mask  = True)
    fmodel_result.update_all_scales(update_f_part1=False)
    return fmodel_result

  def classical_dm_map_coefficients(self):
    import mmtbx.density_modification
    import iotbx.phil
    from libtbx.utils import null_out
    fmodel_truncated = self.remove_unreliable_atoms_and_update()
    coeffs = map_tools.electron_density_map(
      fmodel=fmodel_truncated).map_coefficients(
        map_type         = "2mFo-DFc",
        isotropize       = True,
        fill_missing     = False)
    import mmtbx.utils
    solvent_content = mmtbx.utils.f_000(xray_structure =
      self.xray_structure).solvent_fraction # here we use the original xrs
    params = iotbx.phil.parse(mmtbx.density_modification.master_params_str).extract()
    params.solvent_fraction = solvent_content + 0.15
    params.grid_resolution_factor = 1./3
    params.initial_steps = 3
    params.shrink_steps = 5
    params.final_steps = 3
    hl_model = miller.set(crystal_symmetry=self.f_obs().crystal_symmetry(),
        indices = self.f_obs().indices(),
        anomalous_flag=False).array(
          data=self.f_model_phases_as_hl_coefficients(
            map_calculation_helper=None))
    dm = mmtbx.density_modification.density_modification(
      params           = params,
      f_obs            = self.f_obs(),
      hl_coeffs_start  = hl_model,
      map_coeffs       = coeffs,
      log              = null_out(),
      as_gui_program   = False)
    return dm.map_coeffs_in_original_setting

  def resolve_dm_map_coefficients(self):
    if(not libtbx.env.has_module("solve_resolve")):
      raise Sorry("solve_resolve not available.")
    fmodel_truncated = self.remove_unreliable_atoms_and_update()
    coeffs = map_tools.electron_density_map(
      fmodel=fmodel_truncated).map_coefficients(
        map_type         = "2mFo-DFc",
        isotropize       = True,
        fill_missing     = False)
    import mmtbx.utils
    solvent_content = mmtbx.utils.f_000(xray_structure =
      self.xray_structure).solvent_fraction # here we use the original xrs
    return mmtbx.map_tools.resolve_dm_map(
        fmodel       = fmodel_truncated,
        map_coeffs   = coeffs,
        xrs          = self.xray_structure,
        add_mask     = True,  # use add_mask instead of model_mask
        mask_cycles  = 2,
        minor_cycles = 2,
        use_model_hl = True,
        fill         = True)

  def k_sol_b_sol_from_k_mask(self):
    sel = self.f_obs().d_spacings().data()>=4.0
    k_mask = self.k_masks()[0].select(sel)
    ss = self.ss.select(sel)
    r = scitbx.math.gaussian_fit_1d_analytical(x=flex.sqrt(ss), y=k_mask)
    k_sol_0, b_sol_0 = flex.max(k_mask), r.b
    if(k_sol_0<0):   k_sol_0=0
    if(k_sol_0>0.6): k_sol_0=0.6
    if(b_sol_0<0):   b_sol_0=0
    if(b_sol_0>150): b_sol_0=150
    k_range = [k/100. for k in range(-10,11)]
    b_range = [b for b in range(-20,21)]
    t = 1.e+9
    k_best, b_best = None, None
    for k_sh in k_range:
      for b_sh in b_range:
        k_sol = k_sol_0 + k_sh
        b_sol = b_sol_0 + b_sh
        if(k_sol<0 or k_sol>0.6 or b_sol<0 or b_sol>150): continue
        k_mask_ = k_sol * flex.exp(-b_sol * ss)
        delta = k_mask - k_mask_
        t_ = flex.sum(delta*delta)
        if(t_<t):
          t = t_
          k_best, b_best = k_sol, b_sol
    return k_best, b_best

  def need_to_refine_hd_scattering_contribution(self):
    refine_hd_scattering = True
    hd_selection = self.xray_structure.hd_selection()
    occ_h_all_zero = self.xray_structure.select(
      hd_selection).scatterers().extract_occupancies().all_eq(0.0) # riding H
    if(self.xray_structure is None or
       self.xray_structure.guess_scattering_type_neutron() or
       hd_selection.count(True)==0 or
       not occ_h_all_zero):
      refine_hd_scattering = False
    return refine_hd_scattering

  def update_all_scales(
        self,
        update_f_part1 = True,
        apply_back_trace=False,
        params = None, # XXX DUMMY
        nproc = None,  # XXX DUMMY
        cycles = None,
        fast = True,
        optimize_mask = False,
        refine_hd_scattering = True,
        refine_hd_scattering_method = "fast",
        bulk_solvent_and_scaling = True,
        remove_outliers = True,
        apply_scale_k1_to_f_obs=True,
        show = False,
        verbose=None,
        log = None):
    self.alpha_beta_cache = None
    if(apply_scale_k1_to_f_obs):
      self.apply_scale_k1_to_f_obs()
    from mmtbx.bulk_solvent import f_model_all_scales
    o = f_model_all_scales.run(
      fmodel               = self,
      apply_back_trace     = apply_back_trace,
      remove_outliers      = remove_outliers,
      fast                 = fast,
      params               = params,
      log                  = log,
      refine_hd_scattering = refine_hd_scattering)
    if(optimize_mask): self.optimize_mask(out=log)
    if(remove_outliers and self.twin_law is None):
      for i in [1,2,3,4,5]:
        if(self.r_work_low()>0.7 and len(self.bin_selections)>1):
          ds_low = flex.min(
            self.f_obs().d_spacings().data().select(self.bin_selections[0]))
          f_obs = self.f_obs().resolution_filter(d_max=ds_low)
          abcd = self.hl_coeffs()
          if(abcd is not None):
            abcd = self.hl_coeffs().common_set(f_obs)
          self.__init__(
             f_obs                        = f_obs,
             r_free_flags                 = self.r_free_flags().common_set(f_obs),
             f_part1                      = self.f_part1()     .common_set(f_obs),
             f_part2                      = self.f_part2()     .common_set(f_obs),
             f_calc                       = self.f_calc()      .common_set(f_obs),
             f_mask = [f.common_set(f_obs) for f in self.f_masks()],
             abcd                         = abcd,
             epsilons                     = None,
             sf_and_grads_accuracy_params = self.sfg_params,
             target_name                  = "ml",
             alpha_beta_params            = self.alpha_beta_params ,
             xray_structure               = self.xray_structure    ,
             mask_params                  = self.mask_params       ,
             use_f_model_scaled           = self.use_f_model_scaled,
             twin_law                     = self.twin_law          ,
             twin_fraction                = self.twin_fraction     ,
             max_number_of_bins           = self.max_number_of_bins,
             bin_selections               = None    ,
             n_resolution_bins_output     = self.n_resolution_bins_output,
             scale_method                 = self.scale_method)
          o = f_model_all_scales.run(
            fmodel=self, apply_back_trace = apply_back_trace,
            remove_outliers = remove_outliers, fast = fast,
            params = params, log = log, refine_hd_scattering=refine_hd_scattering)
    self.russ = o.russ
    return o.russ

  def apply_scale_k1_to_f_obs(self):
    r_start = self.r_work()
    fo = self.f_obs().data()
    fc = abs(self.f_model()).data()
    sc = flex.sum(fo*fc)/flex.sum(fc*fc)
    if sc < 1.1 and sc > 0.9: return
    if sc == 0: return
    if self.twin_law is not None: return
    sigmas = self.f_obs().sigmas()
    if(sigmas is not None):
      sigmas = sigmas/sc
    f_obs_new = self.f_obs().customized_copy(
      data   = self.f_obs().data()/sc,
      sigmas = sigmas)
    self.update(f_obs = f_obs_new)
    r_final = self.r_work()
    assert approx_equal(r_start, r_final), [r_start, r_final]

  def r_work4(self):
    s4 = self.resolution_filter(d_min=4)
    if(s4.f_obs().data().size()==0): return 0
    return s4.r_work()

  def _get_target_name(self): return self._target_name
  target_name = property(_get_target_name)

  def target_attributes(self):
    if (self.target_name is None): return None
    result = mmtbx.refinement.targets.target_names.get(self.target_name)
    if (result is None):
      raise RuntimeError(
        "Unknown target name: %s" % show_string(self.target_name))
    return result

  def target_functor(self, alpha_beta=None):
    return mmtbx.refinement.targets.target_functor(manager=self,
      alpha_beta=alpha_beta)

  def set_target_name(self, target_name):
    if (target_name == "ls"): target_name = "ls_wunit_k1"
    if(target_name=="mlhl" and self._hl_coeffs is None):
      raise Sorry("Using MLHL target requires HL present: no HL provided.")
    self._target_name = target_name

  def determine_n_bins(self,
        free_reflections_per_bin,
        max_n_bins=None,
        min_n_bins=1,
        min_refl_per_bin=100):
    assert free_reflections_per_bin > 0
    n_refl = self.arrays.free_sel.size()
    n_free = self.arrays.free_sel.count(True)
    n_refl_per_bin = free_reflections_per_bin
    if (n_free != 0):
      n_refl_per_bin *= n_refl / n_free
    n_refl_per_bin = min(n_refl, iround(n_refl_per_bin))
    result = max(1, iround(n_refl / max(1, n_refl_per_bin)))
    if (min_n_bins is not None):
      result = max(result, min(min_n_bins, iround(n_refl / min_refl_per_bin)))
    if (max_n_bins is not None):
      result = min(max_n_bins, result)
    return result

  def target_unscaled_w(self, alpha_beta=None):
    return self.target_w(alpha_beta=alpha_beta)*self.f_calc_w().data().size()

  def target_unscaled_t(self, alpha_beta=None):
    return self.target_t(alpha_beta=alpha_beta)*self.f_calc_t().data().size()

  def _update_mask_params(self, mask_params):
    self.mask_params = mask_params
    self.mask_manager = masks.manager(
      miller_array      = self.f_obs(),
      miller_array_twin = self.twin_set,
      mask_params       = self.mask_params)

  def update(self, f_calc                       = None,
                   f_obs                        = None,
                   f_mask                       = None,
                   r_free_flags                 = None,
                   f_part1                      = None,
                   k_mask                       = None,
                   k_anisotropic                = None,
                   sf_and_grads_accuracy_params = None,
                   target_name                  = None,
                   abcd                         = None,
                   alpha_beta_params            = None,
                   twin_fraction                = None,
                   xray_structure               = None,
                   epsilons                     = None,
                   mask_params                  = None):
    self.update_core(
      f_calc = f_calc,
      f_mask = f_mask,
      f_part1 = f_part1,
      k_mask = k_mask,
      k_anisotropic=k_anisotropic)
    if(epsilons is not None):
      self.epsilons = epsilons
    if(mask_params is not None): # Reguires updating mask manager!
      self._update_mask_params(mask_params=mask_params)
    if(f_obs is not None):
      self.arrays.f_obs = f_obs
      self.update_core()
    if(r_free_flags is not None):
      self._r_free_flags = r_free_flags
    if(sf_and_grads_accuracy_params is not None):
      self.sfg_params = sf_and_grads_accuracy_params
      self.update_xray_structure(update_f_calc  = True)
    if(target_name is not None):
      self.set_target_name(target_name=target_name)
    if abcd is not None:
      if not abcd.indices().all_eq(self.f_obs().indices()):
        abcd = abcd.complete_array(
          new_data_value=(0,0,0,0), d_min=self.f_obs().d_min())\
             .common_set(self.f_obs())
      self.abcd = abcd
      self._hl_coeffs = abcd
      self.arrays.hl_coeffs = abcd.common_set(self.f_obs())
    if(twin_fraction is not None):
      self.twin_fraction = twin_fraction
      self.update_core()
    if(alpha_beta_params is not None):
      self.alpha_beta_params = alpha_beta_params
    if(xray_structure is not None):
       self.update_xray_structure(xray_structure = xray_structure,
                                  update_f_mask  = True,
                                  update_f_calc  = True)
    return self

  def checksum(self,
               k_isotropic=True,
               k_anisotropic=True,
               f_calc=True,
               f_masks=True,
               k_masks=True,
               f_obs=True):
    result = 0
    sz = self.k_isotropic().size()
    if(k_isotropic):   result += flex.sum(self.k_isotropic())/sz
    if(k_anisotropic): result += flex.sum(self.k_anisotropic())/sz
    if(f_calc):        result += flex.sum(abs(self.f_calc()).data())/sz
    if(f_masks):       result += sum([flex.sum(abs(m).data())/sz
                                     for m in self.f_masks()])
    if(k_masks):       result += sum([flex.sum(m)/sz for m in self.k_masks()])
    if(f_obs):         result += flex.sum(self.f_obs().data())/sz
    return result

  def hl_coeffs(self):
    if(self.arrays is None):
      return self._hl_coeffs
    else:
      return self.arrays.hl_coeffs

  def f_obs_scaled(self, include_fom=False):
    scale = 1.0 / self.k_total()
    if include_fom:
      scale = scale * self.figures_of_merit()
    return miller.array(
      miller_set = self.f_obs(),
      data       = self.f_obs().data()*scale)

  def f_obs(self):
    if(self.arrays is not None):
      return self.arrays.f_obs
    else:
      return self._f_obs

  def i_obs(self):
    return self._i_obs

  def r_free_flags(self):
      return self._r_free_flags

  def f_bulk(self):
    return miller.array(miller_set = self.f_obs(),
      data = self.arrays.core.data.f_bulk)

  def f_bulk_w(self):
    if(self.r_free_flags().data().count(True) > 0):
      return self.f_bulk().select(self.arrays.work_sel)
    else:
      return self.f_bulk()

  def f_bulk_t(self):
    if(self.r_free_flags().data().count(True) > 0):
      return self.f_bulk().select(self.arrays.free_sel)
    else:
      return self.f_bulk()

  def k_anisotropic(self):
    return self.arrays.core.data.k_anisotropic

  def k_anisotropic_work(self):
    return self.arrays.k_anisotropic_work

  def k_anisotropic_test(self):
    return self.k_anisotropic().select(self.arrays.free_sel)

  def k_masks(self):
    return self.arrays.core.k_masks

  def k_masks_twin(self):
    return self.arrays.core_twin.k_masks

  def k_isotropic(self):
    return self.arrays.core.k_isotropic

  def k_isotropic_work(self):
    return self.arrays.k_isotropic_work

  def k_total(self):
    return self.k_isotropic()*self.k_anisotropic()*self.scale_k1()*\
      self.arrays.core.k_isotropic_exp

  def k_total_work(self):
    return self.k_isotropic_work()*self.k_anisotropic_work()*self.scale_k1_w()*\
      self.arrays.core.k_isotropic_exp.select(self.arrays.work_sel)

  def f_obs_work(self):
    return self.arrays.f_obs_work

  def f_obs_free(self):
    return self.arrays.f_obs_free

  def f_model(self):
    return self.arrays.f_model

  def f_model_no_scales(self):
    # or better name: f_calc + f_bulk
    f_masks = self.f_masks()
    k_masks = self.k_masks()
    assert len(f_masks)==1
    assert len(k_masks)==1
    d = self.f_calc().data() + f_masks[0].data()*k_masks[0]
    return f_masks[0].customized_copy(data = d)

  def f_model_work(self):
    return self.arrays.f_model_work

  def f_masks(self):
    return self.arrays.core.f_masks

  def f_masks_twin(self):
    return self.arrays.core_twin.f_masks

  def f_model_free(self):
    return self.arrays.f_model_free

  def f_model_scaled_with_k1(self):
    return miller.array(
      miller_set = self.f_obs(),
      data       = self.scale_k1()*self.f_model().data())

  def f_model_scaled_with_k1_composite_work_free(self):
    ma_w = self.f_model_scaled_with_k1_w()
    ma_f = self.f_model_scaled_with_k1_t()
    if(ma_w.indices().size() == ma_f.indices().size()): return ma_w
    return ma_w.concatenate(ma_f)

  def f_model_scaled_with_k1_t(self):
    return miller.array(
      miller_set = self.f_obs_free(),
      data       = self.scale_k1_t()*self.f_model_free().data())

  def f_model_scaled_with_k1_w(self):
    return miller.array(
      miller_set = self.f_obs_work(),
      data       = self.scale_k1_w()*self.f_model_work().data())

  def f_star_w_star_obj(self):
    alpha, beta = self.alpha_beta()
    obj = max_lik.f_star_w_star_mu_nu(
      f_obs          = self.f_obs().data(),
      f_model        = flex.abs(self.f_model().data()),
      alpha          = alpha.data(),
      beta           = beta.data(),
      space_group    = self.f_obs().space_group(),
      miller_indices = self.f_obs().indices())
    return obj

  def f_star_w_star(self):
    obj = self.f_star_w_star_obj()
    f_star = miller.array(miller_set = self.f_obs(), data = obj.f_star())
    w_star = miller.array(miller_set = self.f_obs(), data = obj.w_star())
    return f_star, w_star

  def f_part1(self):
    return self.arrays.core.f_part1

  def f_part1_twin(self):
    if(self.arrays.core_twin is not None):
      return self.arrays.core_twin.f_part1

  def f_part2(self):
    return self.arrays.core.f_part2

  def f_part2_twin(self):
    if(self.arrays.core_twin is not None):
      return self.arrays.core_twin.f_part2

  def f_calc(self):
    return self.arrays.core.f_calc

  def f_calc_twin(self):
    if(self.arrays.core_twin is not None):
      return self.arrays.core_twin.f_calc

  def f_calc_w(self):
    if(self.r_free_flags().data().count(True) > 0):
      return self.f_calc().select(~self.r_free_flags().data())
    else:
      return self.f_calc()

  def f_calc_t(self):
    if(self.r_free_flags().data().count(True) > 0):
      return self.f_calc().select(self.r_free_flags().data())
    else:
      return self.f_calc()

  def f_obs_vs_f_model(self, log=None):
    if(log is None): log = sys.stdout
    for fo, fm, d in zip(self.f_obs().data(),
                         abs(self.f_model_scaled_with_k1()),
                         self.f_obs().d_spacings().data()):
      print("Fobe   Fmodel   resolution (A)", file=log)
      print("%12.3f %12.3f %6.3f"%(fo, fm, d), file=log)

  #Ensemble refinement alpha beta parameters
  set_sigmaa = None
  n_obs      = None
  n_calc     = None
  def n_obs_n_calc(self, update_nobs_ncalc = True):
    p = self.alpha_beta_params.sigmaa_estimator
    fmodel = self.f_model()
    n_obs, n_calc = mmtbx.scaling.ta_alpha_beta_calc.ta_alpha_beta_calc(
             miller_obs = self.f_obs(),
             miller_calc = fmodel,
             r_free_flags = self.r_free_flags(),
             ta_d = self.set_sigmaa,
             kernel_width_free_reflections=p.kernel_width_free_reflections,
             kernel_on_chebyshev_nodes=p.kernel_on_chebyshev_nodes,
             n_sampling_points=p.number_of_sampling_points,
             n_chebyshev_terms=p.number_of_chebyshev_terms,
             use_sampling_sum_weights=p.use_sampling_sum_weights).eobs_and_ecalc_miller_array_normalizers()
    if update_nobs_ncalc:
      self.n_obs  = n_obs
      self.n_calc = n_calc
    return n_obs, n_calc

  def alpha_beta_with_restrained_n_calc(self):
    assert self.n_obs is not None
    assert self.n_calc is not None
    alpha = self.set_sigmaa * flex.sqrt(self.n_obs / self.n_calc)
    beta  = (1.0 - self.set_sigmaa*self.set_sigmaa) * self.n_obs
    alpha = self.f_obs().array(data=alpha)
    beta  = self.f_obs().array(data=beta)
    return alpha, beta

  def alpha_beta(self, f_obs = None, f_model = None):
    global time_alpha_beta
    timer = user_plus_sys_time()
    if self.alpha_beta_cache is not None:
      return self.alpha_beta_cache
    if(f_obs is None): f_obs = self.f_obs()
    if(f_model is None): f_model = self.f_model()
    # Use work reflections if test set is too small
    flags = self.r_free_flags().data()
    d_max = self.f_obs().d_max_min()[0]
    if(flags.count(True)<300 and d_max>15 and
       flags.count(True)!=flags.size()): # XXX ad hoc. Make a parameter?
      alpha, beta = maxlik.alpha_beta_est_manager(
           f_obs                    = f_obs,
           f_calc                   = f_model,
           free_reflections_per_bin = 140,
           flags                    = ~flags,
           interpolation            = False,
           epsilons                 = self.epsilons).alpha_beta()
      self.alpha_beta_cache = (alpha, beta)
      return self.alpha_beta_cache
    alpha, beta = None, None
    ab_params = self.alpha_beta_params
    #Nobs and Ncalc restrained for ensemble refinement
    if self.set_sigmaa != None:
      alpha, beta = self.alpha_beta_with_restrained_n_calc()
      self.alpha_beta_cache = (alpha, beta)
      return self.alpha_beta_cache
    if(self.alpha_beta_params is not None):
       assert self.alpha_beta_params.method in ("est", "calc")
       assert self.alpha_beta_params.estimation_algorithm in [
         "analytical", "iterative"]
       if (self.alpha_beta_params.method == "est"):
         if (self.alpha_beta_params.estimation_algorithm == "analytical"):
           alpha, beta = maxlik.alpha_beta_est_manager(
             f_obs           = f_obs,
             f_calc          = f_model,
             free_reflections_per_bin
               = self.alpha_beta_params.free_reflections_per_bin,
             flags           = self.r_free_flags().data(),
             interpolation   = self.alpha_beta_params.interpolation,
             epsilons        = self.epsilons).alpha_beta()
         else:
           p = self.alpha_beta_params.sigmaa_estimator
           alpha, beta = sigmaa_estimator(
             miller_obs  = f_obs,
             miller_calc = f_model,
             r_free_flags=self.r_free_flags(),
             kernel_width_free_reflections=p.kernel_width_free_reflections,
             kernel_on_chebyshev_nodes=p.kernel_on_chebyshev_nodes,
             n_sampling_points=p.number_of_sampling_points,
             n_chebyshev_terms=p.number_of_chebyshev_terms,
             use_sampling_sum_weights=p.use_sampling_sum_weights).alpha_beta()
       else:
         n_atoms_missed = ab_params.number_of_macromolecule_atoms_absent + \
                          ab_params.number_of_waters_absent
         alpha, beta = maxlik.alpha_beta_calc(
           f                = f_obs,
           n_atoms_absent   = n_atoms_missed,
           n_atoms_included = ab_params.n_atoms_included,
           bf_atoms_absent  = ab_params.bf_atoms_absent,
           final_error      = ab_params.final_error,
           absent_atom_type = ab_params.absent_atom_type).alpha_beta()
    else:
       alpha, beta = maxlik.alpha_beta_est_manager(
         f_obs                    = f_obs,
         f_calc                   = f_model,
         free_reflections_per_bin = 140,
         flags                    = self.r_free_flags().data(),
         interpolation            = False,
         epsilons                 = self.epsilons).alpha_beta()
    time_alpha_beta += timer.elapsed()
    self.alpha_beta_cache = (alpha, beta)
    return self.alpha_beta_cache

  def alpha_beta_w(self, only_if_required_by_target=False):
    if(only_if_required_by_target):
      if(self.target_name not in ["ml", "mlhl"]): return None, None
    alpha, beta = self.alpha_beta()
    if(self.r_free_flags().data().count(True) > 0):
      return alpha.select(self.arrays.work_sel), \
             beta.select(self.arrays.work_sel)
    else:
      return alpha, beta

  def alpha_beta_t(self):
    alpha, beta = self.alpha_beta()
    if(self.r_free_flags().data().count(True) > 0):
      return alpha.select(self.arrays.free_sel), \
             beta.select(self.arrays.free_sel)
    else:
      return alpha, beta

  def sigmaa(self, f_obs = None, f_model = None):
    p = self.alpha_beta_params.sigmaa_estimator
    sigmaa_obj = sigmaa_estimator(
      miller_obs                    = self.f_obs(),
      miller_calc                   = self.f_model(),
      r_free_flags                  = self.r_free_flags(),
      kernel_width_free_reflections = p.kernel_width_free_reflections,
      kernel_on_chebyshev_nodes     = p.kernel_on_chebyshev_nodes,
      n_sampling_points             = p.number_of_sampling_points,
      n_chebyshev_terms             = p.number_of_chebyshev_terms,
      use_sampling_sum_weights      = p.use_sampling_sum_weights)
    return sigmaa_obj

  def model_error_ml(self):
    #XXX needs clean solution / one more unfinished project
    if (self.alpha_beta_params is None):
      free_reflections_per_bin = 140
      estimation_algorithm = "analytical"
    else:
      free_reflections_per_bin=self.alpha_beta_params.free_reflections_per_bin
      estimation_algorithm = self.alpha_beta_params.estimation_algorithm
    assert estimation_algorithm in ["analytical", "iterative"]
    est_exceptions = []
    if(estimation_algorithm == "analytical"):
      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = self.f_obs(),
        f_calc          = self.f_model_scaled_with_k1(),
        free_reflections_per_bin = free_reflections_per_bin,
        flags           = self.r_free_flags().data(),
        interpolation   = True,
        epsilons        = self.epsilons).alpha_beta()
    else:
      p = self.alpha_beta_params.sigmaa_estimator
      alpha, beta = sigmaa_estimator(
        miller_obs=self.f_obs(),
        miller_calc=self.f_model_scaled_with_k1(),
        r_free_flags=self.r_free_flags(),
        kernel_width_free_reflections=p.kernel_width_free_reflections,
        kernel_on_chebyshev_nodes=p.kernel_on_chebyshev_nodes,
        n_sampling_points=p.number_of_sampling_points,
        n_chebyshev_terms=p.number_of_chebyshev_terms,
        use_sampling_sum_weights=p.use_sampling_sum_weights).alpha_beta()
    sel = alpha.data() > 0
    alpha_data = alpha.data().select(sel)
    sj = self.ss.select(sel)
    sel = alpha_data > 1.
    alpha_data = alpha_data.set_selected(sel, 1.)
    aj = -math.pi**3 * sj
    bj = flex.log(alpha_data)
    den = flex.sum(aj*aj)
    if(den != 0):
      omega = math.sqrt( flex.sum(aj*bj) / den )
    else: omega = None
    return omega

  def _r_factor(self,
                type="work",
                d_min=None,
                d_max=None,
                d_spacings=None,
                selection=None):
    global time_r_factors
    if(type == "work"):
      f_obs = self.f_obs_work().data()
      f_model = self.f_model_scaled_with_k1_w().data()
    elif(type == "free"):
      f_obs = self.f_obs_free().data()
      f_model = self.f_model_scaled_with_k1_t().data()
    elif(type == "all"):
      f_obs = self.f_obs().data()
      f_model = self.f_model_scaled_with_k1().data()
    else: raise RuntimeError
    if(selection is not None): assert [d_min, d_max].count(None) == 2
    if([d_min, d_max].count(None) < 2):
      assert selection is None and d_spacings is not None
    timer = user_plus_sys_time()
    if(d_min is not None or d_max is not None):
      keep = flex.bool(d_spacings.size(), True)
      if (d_max is not None): keep &= d_spacings <= d_max
      if (d_min is not None): keep &= d_spacings >= d_min
      f_obs   = f_obs.select(keep)
      f_model = f_model.select(keep)
    if(selection is not None):
      f_obs   = f_obs.select(selection)
      f_model = f_model.select(selection)
    result = abs(mmtbx.bulk_solvent.r_factor(f_obs, f_model, 1.0))
    time_r_factors += timer.elapsed()
    if(result >= 1.e+9): result = None
    return result

  def r_work(self, d_min = None, d_max = None, selection = None):
    return self._r_factor(
      type       = "work",
      d_min      = d_min,
      d_max      = d_max,
      d_spacings = self.arrays.d_spacings_work,
      selection  = selection)

  def r_free(self, d_min = None, d_max = None, selection = None):
    return self._r_factor(
      type       = "free",
      d_min      = d_min,
      d_max      = d_max,
      d_spacings = self.arrays.d_spacings_free,
      selection  = selection)

  def r_all(self):
    return self._r_factor(type="all")

  def scale_k1_w_for_twin_targets(self):
    s = ~self.r_free_flags().data()
    fo = self.f_obs_work().data()
    fmask = self.f_masks()[0].data().select(s)
    fm = self.k_anisotropic_work() * (self.f_calc_w().data() +
      self.k_masks()[0].select(s) * fmask)
    return _scale_helper(num=fo, den=flex.abs(fm), selection=None)

  #XXX Fix k1 option for TA
  set_scale_switch = 0
  def scale_k1(self, selection = None):
    if self.set_scale_switch is not None:
      if (self.set_scale_switch == 0):
        return _scale_helper(
        num=self.f_obs().data(),
        den=flex.abs(self.f_model().data()),
        selection=selection)
      if (self.set_scale_switch > 0):
        return self.set_scale_switch
    else:
      return _scale_helper(
        num=self.f_obs().data(),
        den=flex.abs(self.f_model().data()),
        selection=selection)

  def scale_k1_w(self, selection = None):
    if self.set_scale_switch is not None:
      if (self.set_scale_switch == 0):
        return _scale_helper(
          num=self.f_obs_work().data(),
          den=flex.abs(self.f_model_work().data()),
          selection=selection)
      if (self.set_scale_switch > 0):
        return self.set_scale_switch

  def scale_k1_t(self, selection = None):
    if self.set_scale_switch is not None:
      if (self.set_scale_switch == 0):
        return _scale_helper(
          num=self.f_obs_free().data(),
          den=flex.abs(self.f_model_free().data()),
          selection=selection)
      if (self.set_scale_switch > 0):
        return self.set_scale_switch

  def scale_k2(self, selection = None):
    return _scale_helper(
      num=flex.abs(self.arrays.f_model.data()),
      den=self.f_obs().data(),
      selection=selection)

  def scale_k2_w(self, selection = None):
    return _scale_helper(
      num=flex.abs(self.f_model_work().data()),
      den=self.f_obs_work().data(),
      selection=selection)

  def scale_k2_t(self, selection = None):
    return _scale_helper(
      num=flex.abs(self.f_model_free().data()),
      den=self.f_obs_free().data(),
      selection=selection)

  def scale_k3_w(self, selection = None):
    mul_eps_sqrt = flex.sqrt(
        self.f_obs_work().multiplicities().data().as_double()
      / self.f_obs_work().epsilons().data().as_double())
    result = _scale_helper(
      num=self.f_obs_work().data() * mul_eps_sqrt,
      den=flex.abs(self.f_model_work().data()) * mul_eps_sqrt,
      selection=selection,
      num_num=True)
    if (result is None): return None
    return result**0.5

  def scale_k3_t(self, selection = None):
    mul_eps_sqrt = flex.sqrt(
        self.f_obs_free().multiplicities().data().as_double()
      / self.f_obs_free().epsilons().data().as_double())
    result = _scale_helper(
      num=self.f_obs_free().data() * mul_eps_sqrt,
      den=flex.abs(self.f_model_free().data()) * mul_eps_sqrt,
      selection=selection,
      num_num=True)
    if (result is None): return None
    return result**0.5

  def k_mask_low(self):
    k_masks = self.k_masks()
    assert len(k_masks) == 1
    k_mask_l = k_masks[0].select(self.bin_selections[0])
    assert approx_equal(flex.min(k_mask_l), flex.max(k_mask_l), 1.e-3)
    return k_mask_l[0]

  def r_n_lowest(self, n=500):
    ds = self.f_obs().d_spacings().data()
    sel = flex.sort_permutation(ds,reverse=True)
    ds = ds.select(sel)
    d_min = ds[n]
    return self.resolution_filter(d_min=d_min).r_work()

  def r_two_bins_lowest(self):
    f = self.select(self.bin_selections[0] | self.bin_selections[1])
    return mmtbx.bulk_solvent.r_factor(
      f.f_obs_work().data(), f.f_model_work().data())

  def r_work_low(self):
    f = self.select(self.bin_selections[0])


    # Scan for best k
    #fo, fc = f.f_obs_work().data(), f.f_model_work().data()
    #scale = mmtbx.bulk_solvent.scale(fo, fc)
    #r1 = mmtbx.bulk_solvent.r_factor(fo, fc, scale)
    #r1_ = mmtbx.bulk_solvent.r_factor(fo, fc)
    #assert approx_equal(r1, r1_)
    #
    #r2=r1
    #k_best = scale
    #k=scale-0.5
    #while k<=scale+0.5:
    #  r2_ = mmtbx.bulk_solvent.r_factor(fo, fc,k)
    #  if(r2_<r2):
    #    k_best = k
    #    r2 = r2_
    #  k+=0.01
    #return mmtbx.bulk_solvent.r_factor(
    #  f.f_obs_work().data(), f.f_model_work().data(), k_best)

    return mmtbx.bulk_solvent.r_factor(
      f.f_obs_work().data(), f.f_model_work().data())

  def r_work_high(self):
    f = self.select(self.bin_selections[len(self.bin_selections)-1])
    return mmtbx.bulk_solvent.r_factor(
      f.f_obs_work().data(), f.f_model_work().data())

  def completeness_high(self):
    return self.f_obs().select(
      self.bin_selections[len(self.bin_selections)-1]).completeness()

  def completeness_low(self):
    return self.f_obs().select(self.bin_selections[0]).completeness()

  def r_factors(self, prefix="", as_string=True):
    rw,rf,rh,rl,r4=self.r_work(),self.r_free(),self.r_work_high(),\
      self.r_work_low(), self.r_work4()
    f = "%s r_work=%6.4f r_free=%6.4f r_high=%6.4f r_low=%6.4f r_4=%6.4f"
    if(as_string): return f%(prefix, rw,rf,rh,rl,r4)
    else:          return group_args(rw=rw,rf=rf,rh=rh,rl=rl,r4=r4)

  def r_overall_low_high(self, d = 6.0):
    r_work = self.r_work()
    d_max, d_min = self.f_obs_work().d_max_min()
    if(d_max < d): d = d_max
    if(d_min > d): d = d_min
    n_low = self.f_obs_work().resolution_filter(d_min = d, d_max = 999.9).data().size()
    if(n_low > 0):
       r_work_l = self.r_work(d_min = d, d_max = 999.9)
    else:
       r_work_l = None
    n_high = self.f_obs_work().resolution_filter(d_min = 0.0, d_max = d).data().size()
    if(n_high > 0):
       r_work_h = self.r_work(d_min = 0.0, d_max = d)
    else:
       r_work_h = None
    if(r_work_l is not None):
       r_work_l = r_work_l
    else:
       r_work_l = 0.0
    if(r_work_h is not None):
       r_work_h = r_work_h
    else:
       r_work_h = 0.0
    return r_work, r_work_l, r_work_h, n_low, n_high

  def scale_ml_wrapper(self):
    if (self.alpha_beta_params is None): return 1.0
    if (self.alpha_beta_params.method != "calc"): return 1.0

  def figures_of_merit(self):
    alpha, beta = self.alpha_beta()
    global time_foms
    timer = user_plus_sys_time()
    result = max_lik.fom_and_phase_error(
      f_obs          = self.f_obs().data(),
      f_model        = flex.abs(self.arrays.f_model.data()),
      alpha          = alpha.data(),
      beta           = beta.data(),
      epsilons       = self.epsilons,
      centric_flags  = self.f_obs().centric_flags().data()).fom()
    time_foms += timer.elapsed()
    return result

  def fom(self):
    r = self.figures_of_merit()
    return miller.array(
      miller_set = self.f_obs(),
      data       = r)

  def figures_of_merit_work(self):
    fom = self.figures_of_merit()
    if(self.r_free_flags().data().count(True) > 0):
      return fom.select(self.arrays.work_sel)
    else:
      return fom

  def phase_errors(self):
    alpha, beta = self.alpha_beta()
    global time_phase_errors
    timer = user_plus_sys_time()
    result = max_lik.fom_and_phase_error(
      f_obs          = self.f_obs().data(),
      f_model        = flex.abs(self.arrays.f_model.data()),
      alpha          = alpha.data(),
      beta           = beta.data(),
      epsilons       = self.epsilons,
      centric_flags  = self.f_obs().centric_flags().data()).phase_error()
    time_phase_errors += timer.elapsed()
    return result

  def phase_errors_test(self):
    pher = self.phase_errors()
    if(self.r_free_flags().data().count(True) > 0):
      return pher.select(self.r_free_flags().data())
    else:
      return pher

  def phase_errors_work(self):
    pher = self.phase_errors()
    if(self.r_free_flags().data().count(True) > 0):
      return pher.select(self.arrays.work_sel)
    else:
      return pher

  def filter_by_fom(self, fom_threshold = 0.2):
    fom = self.figures_of_merit()
    sel = fom > fom_threshold
    dsel = self.f_obs().d_spacings().data() < 5.
    sel = sel.set_selected(~dsel,True)
    self = self.select(sel)
    return self

  def filter_by_delta_fofc(self, n):
    deltas = flex.double()
    assert self.f_model_scaled_with_k1().data().size() == self.f_obs().data().size()
    for fc, fo in zip(self.f_model_scaled_with_k1().data(), self.f_obs().data()):
      deltas.append(abs(fo - abs(fc)))
    sel = flex.sort_permutation(deltas, reverse=True)
    self = self.select(sel)
    sel = flex.size_t(range(n, deltas.size()))
    sel = flex.bool(deltas.size(), sel)
    self = self.select(sel)
    return self

  def twinned(self):
    return self.arrays.core_twin is not None

  def map_calculation_helper(self,
                             free_reflections_per_bin = 100,
                             interpolation = True):
    class result(object):
      def __init__(self, fmodel, free_reflections_per_bin, interpolation):
        self.f_obs = fmodel.f_obs()
        self.f_model = fmodel.f_model()
        assert fmodel.epsilons.size()==self.f_obs.data().size()
        self.alpha, self.beta = maxlik.alpha_beta_est_manager(
          f_obs                    = self.f_obs,
          f_calc                   = self.f_model,
          free_reflections_per_bin = free_reflections_per_bin,
          flags                    = fmodel.r_free_flags().data(),
          interpolation            = interpolation,
          epsilons                 = fmodel.epsilons).alpha_beta()
        self.fom = max_lik.fom_and_phase_error(
          f_obs          = self.f_obs.data(),
          f_model        = flex.abs(self.f_model.data()),
          alpha          = self.alpha.data(),
          beta           = self.beta.data(),
          epsilons       = fmodel.epsilons,
          centric_flags  = self.f_obs.centric_flags().data()).fom()
    return result(
      fmodel                   = self,
      free_reflections_per_bin = free_reflections_per_bin,
      interpolation            = interpolation)

  def f_model_phases_as_hl_coefficients(self, map_calculation_helper,
        k_blur=None, b_blur=None):
    if(map_calculation_helper is not None):
      mch = map_calculation_helper
    else:
      mch = self.map_calculation_helper()
    f_model_phases = mch.f_model.phases().data()
    sin_f_model_phases = flex.sin(f_model_phases)
    cos_f_model_phases = flex.cos(f_model_phases)
    if([k_blur, b_blur].count(None) == 2):
      t = maxlik.fo_fc_alpha_over_eps_beta(
        f_obs   = mch.f_obs,
        f_model = mch.f_model,
        alpha   = mch.alpha,
        beta    = mch.beta)
    else:
      assert [k_blur, b_blur].count(None) == 0
      t = 2*k_blur * flex.exp(-b_blur*self.ss)
    hl_a_model = t * cos_f_model_phases
    hl_b_model = t * sin_f_model_phases
    return flex.hendrickson_lattman(a = hl_a_model, b = hl_b_model)

  def f_model_phases_as_hl_coeffs_array(self):
    hl_coeffs = miller.set(crystal_symmetry=self.f_obs().crystal_symmetry(),
        indices = self.f_obs().indices(),
        anomalous_flag=False).array(
          data=self.f_model_phases_as_hl_coefficients(
            map_calculation_helper=None))
    return hl_coeffs

  def combined_hl_coefficients(self, map_calculation_helper):
    result = None
    if(self.hl_coeffs() is not None):
      result = self.hl_coeffs().data() + self.f_model_phases_as_hl_coefficients(
        map_calculation_helper)
    return result

  def combine_phases(self, n_steps = 360, map_calculation_helper=None):
    result = None
    if(self.hl_coeffs() is not None):
      integrator = miller.phase_integrator(n_steps = n_steps)
      phase_source = integrator(
        space_group= self.f_obs().space_group(),
        miller_indices = self.f_obs().indices(),
        hendrickson_lattman_coefficients =
          self.combined_hl_coefficients(map_calculation_helper))
      class tmp:
        def __init__(self, phase_source):
          self.phase_source = phase_source
        def phases(self):
          return flex.arg(self.phase_source)
        def fom(self):
          return flex.abs(self.phase_source)
        def f_obs_phase_and_fom_source(self):
          return self.phase_source
      result = tmp(phase_source = phase_source)
    return result

  def electron_density_map(self):
    return map_tools.electron_density_map(fmodel = self)

  def map_coefficients(self, **kwds):
    emap = self.electron_density_map()
    return emap.map_coefficients(**kwds)

  def _get_real_map(self, **kwds):
    map_coeffs = self.map_coefficients(**kwds)
    return map_coeffs.fft_map(
      resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()

  def two_fofc_map(self, **kwds):
    kwds['map_type'] = "2mFo-DFc"
    return self._get_real_map(**kwds)

  def fofc_map(self, **kwds):
    kwds['map_type'] = "mFo-DFc"
    return self._get_real_map(**kwds)

  def anomalous_map(self, **kwds):
    if (not self.f_obs().anomalous_flag()) : return None
    kwds['map_type'] = "anom"
    return self._get_real_map(**kwds)

  def info(self, free_reflections_per_bin = None, max_number_of_bins = None,
      n_bins=None):
    if(free_reflections_per_bin is None):
      free_reflections_per_bin= self.alpha_beta_params.free_reflections_per_bin
    if(max_number_of_bins is None):
      max_number_of_bins = self.max_number_of_bins
    if (n_bins is None):
      n_bins = self.n_resolution_bins_output
    return mmtbx.f_model.f_model_info.info(
      fmodel                   = self,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins       = max_number_of_bins,
      n_bins                   = n_bins)

  def fft_vs_direct(self, reflections_per_bin = 250,
                          n_bins              = 0,
                          out                 = None):
    if(out is None): out = sys.stdout
    f_obs_w = self.f_obs_work()
    f_obs_w.setup_binner(reflections_per_bin = reflections_per_bin,
                         n_bins              = n_bins)
    fmodel_dc = self.deep_copy()
    fft_p = sf_and_grads_accuracy_master_params.extract()
    fft_p.algorithm = "fft"
    dir_p = sf_and_grads_accuracy_master_params.extract()
    dir_p.algorithm = "direct"
    if(self.sfg_params.algorithm == "fft"):
       fmodel_dc.update(sf_and_grads_accuracy_params = dir_p)
    elif(self.sfg_params.algorithm == "direct"):
       fmodel_dc.update(sf_and_grads_accuracy_params = fft_p)
    print("|"+"-"*77+"|", file=out)
    print("| Bin     Resolution   Compl.  No.       Scale_k1             R-work          |", file=out)
    print("|number     range              Refl.      direct       fft direct fft-direct,%|", file=out)
    deltas = flex.double()
    for i_bin in f_obs_w.binner().range_used():
        sel         = f_obs_w.binner().selection(i_bin)
        r_work_1    = self.r_work(selection = sel)
        scale_k1_1  = self.scale_k1_w(selection = sel)
        r_work_2    = fmodel_dc.r_work(selection = sel)
        scale_k1_2  = fmodel_dc.scale_k1_w(selection = sel)
        f_obs_sel   = f_obs_w.select(sel)
        d_max,d_min = f_obs_sel.d_max_min()
        compl       = f_obs_sel.completeness(d_max = d_max)
        n_ref       = sel.count(True)
        delta       = abs(r_work_1-r_work_2)*100.
        deltas.append(delta)
        d_range     = f_obs_w.binner().bin_legend(
                   i_bin = i_bin, show_bin_number = False, show_counts = False)
        format = "|%3d: %-17s %4.2f %6d %14.3f %6.4f %6.4f   %9.4f  |"
        if(self.sfg_params.algorithm == "fft"):
           print(format % (i_bin,d_range,compl,n_ref,
                                            scale_k1_2,r_work_1,r_work_2,delta), file=out)
        else:
           print(format % (i_bin,d_range,compl,n_ref,
                                            scale_k1_1,r_work_2,r_work_1,delta), file=out)
    print("|"+"-"*77+"|", file=out)
    out.flush()

  def explain_members(self, out=None, prefix="", suffix=""):
    if (out is None): out = sys.stdout
    def zero_if_almost_zero(v, eps=1.e-6):
      if (abs(v) < eps): return 0
      return v
    for line in [
       "Fmodel        = k_isotropic * k_anisotropic * (Fcalc + k_mask * Fmask)",
       "Fcalc         = structure factors calculated from atomic model",
       "Fmask         = structure factors calculated from bulk-solvent mask",
       "k_isotropic   = overall resolution-dependent scale factor",
       "k_anisotropic = overall Millerindex-dependent scale factor",
       "k_mask        = bulk-solvent scale factor"]:
      print(prefix + line + suffix, file=out)

  def export_f_obs_flags_as_mtz(self,
      file_name,
      merge_anomalous=False,
      include_hendrickson_lattman=True):
    """
    Dump all input data to an MTZ file using standard column labels.  This may
    be useful when running modules or programs that require an MTZ file as
    input (rather than taking f_model.manager or the Miller arrays directly).
    """
    f_obs = self.f_obs()
    flags = self.r_free_flags()
    hl_coeffs = self.hl_coeffs()
    if (merge_anomalous):
      f_obs = f_obs.average_bijvoet_mates()
      flags = flags.average_bijvoet_mates()
      if (hl_coeffs is not None):
        hl_coeffs = hl_coeffs.average_bijvoet_mates()
    mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F")
    if (hl_coeffs is not None) and (include_hendrickson_lattman):
      mtz_dataset.add_miller_array(hl_coeffs, column_root_label="HL")
    mtz_dataset.add_miller_array(flags, column_root_label="FreeR_flag")
    mtz_dataset.mtz_object().write(file_name)

  def export(self, out=None, format="mtz"):
    assert format in ["mtz"]
    file_name = None
    if (out is None):
      out = sys.stdout
    elif (hasattr(out, "name") and file_name is None):
      file_name = libtbx.path.canonical_path(file_name=out.name)
    if(format == "mtz"):
      assert file_name is not None
      mtz_dataset = self.f_obs().as_mtz_dataset(column_root_label="FOBS")
      mtz_dataset.add_miller_array(
        miller_array=self.r_free_flags(), column_root_label="R_FREE_FLAGS")
      mtz_dataset.add_miller_array(
        miller_array=self.f_model_scaled_with_k1(), column_root_label="FMODEL")
      mtz_dataset.add_miller_array(
        miller_array=self.f_calc(), column_root_label="FCALC")
      for ifm,fm in enumerate(self.f_masks()):
        if( len(self.f_masks())==1 ):
          label = "FMASK"
        else:
          label= "FMASK%d"%(ifm+1)
        mtz_dataset.add_miller_array(
          miller_array=fm, column_root_label=label)
      mtz_dataset.add_miller_array(
        miller_array=self.f_calc().customized_copy(data=self.k_isotropic()),
        column_root_label="K_ISOTROPIC")
      mtz_dataset.add_miller_array(
        miller_array=self.f_calc().customized_copy(data=self.k_anisotropic()),
        column_root_label="K_ANISOTROPIC")
      for ifm,fm in enumerate(self.k_masks()):
        if( len(self.f_masks())==1 ):
          label = "K_MASK"
        else:
          label= "K_MASK%d"%(ifm+1)
        mtz_dataset.add_miller_array(
          miller_array=self.f_obs().customized_copy(data=fm),
          column_root_label=label)
      mtz_dataset.add_miller_array(
        miller_array=self.f_obs().array(data=self.figures_of_merit()),
        column_root_label="FOM")
      alpha, beta = self.alpha_beta()
      mtz_dataset.add_miller_array(
        miller_array=alpha, column_root_label="ALPHA")
      mtz_dataset.add_miller_array(
        miller_array=beta, column_root_label="BETA")
      mtz_dataset.add_miller_array(
        miller_array=self.f_obs().centric_flags(),
        column_root_label="CENTRIC_FLAGS")
      if(self.hl_coeffs() is not None):
        mtz_dataset.add_miller_array(
          miller_array=self.hl_coeffs(), column_root_label="HL")
      hl_model = miller.set(crystal_symmetry=self.f_obs().crystal_symmetry(),
        indices = self.f_obs().indices(),
        anomalous_flag=False).array(
          data=self.f_model_phases_as_hl_coefficients(
            map_calculation_helper=None))
      mtz_dataset.add_miller_array(
        miller_array = hl_model,
        column_root_label="HLmodel")
      hl_comb = miller.set(crystal_symmetry=self.f_obs().crystal_symmetry(),
        indices = self.f_obs().indices(),
        anomalous_flag=False).array(
          data=self.combined_hl_coefficients(map_calculation_helper=None))
      mtz_dataset.add_miller_array(
        miller_array = hl_model,
        column_root_label="HLcomb")
      mtz_dataset.add_miller_array(
        miller_array = self.f_obs().d_spacings(),
        column_root_label="RESOLUTION", column_types="R")
      mtz_history_buffer = flex.std_string()
      ha = mtz_history_buffer.append
      ha(date_and_time())
      ha("file name: %s" % os.path.basename(file_name))
      ha("directory: %s" % os.path.dirname(file_name))
      s = StringIO()
      self.explain_members(out=s)
      for line in s.getvalue().splitlines():
        ha(line)
      mtz_object = mtz_dataset.mtz_object()
      mtz_object.add_history(lines=mtz_history_buffer)
      out.close()
      mtz_object.write(file_name=file_name)

def show_histogram(data, n_slots, log):
  hm = flex.histogram(data = data, n_slots = n_slots)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    print("%10.3f - %-10.3f : %d" % (lc_1, hc_1, n_1), file=log)
    lc_1 = hc_1

class mask_result(object):
  def __init__(self, r_solv, r_shrink, r_work, r_free, r_work_low):
    adopt_init_args(self, locals())

  def show(self, out):
    if (out is None) : return
    print("r_solv=%s r_shrink=%s r_work=%s r_free=%s r_work_low=%s"%(
          format_value("%6.2f", self.r_solv),
          format_value("%6.2f", self.r_shrink),
          format_value("%6.4f", self.r_work),
          format_value("%6.4f", self.r_free),
          format_value("%6.4f", self.r_work_low)), file=out)

# XXX backwards compatibility 2011-02-08
# bad hack... - relative import
# Change to absolute import
from . import f_model_info
class info(f_model_info.info):
  pass

class resolution_bin(f_model_info.resolution_bin):
  pass
