from __future__ import division
try: import phaser
except ImportError: phaser = None
else: import phaser.phenix_adaptors.sad_target
from cctbx.array_family import flex
import math, time, sys, os, random
from cctbx import miller
from cctbx import crystal
from cctbx import adptbx
from scitbx import lbfgs
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal, not_approx_equal
from mmtbx import bulk_solvent
from mmtbx import masks
from cctbx import xray
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from mmtbx.scaling.sigmaa_estimation import sigmaa_estimator
from mmtbx.refinement import print_statistics
from cctbx.eltbx.xray_scattering import wk1995
from mmtbx.max_lik import max_like_non_uniform
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from cctbx import miller
from iotbx.cns.miller_array import crystal_symmetry_as_cns_comments
import cctbx.xray.structure_factors
from cctbx.array_family import flex
from stdlib import math
from cctbx import xray
from cctbx import adptbx
import boost.python
import mmtbx
from libtbx.math_utils import iround
from libtbx.utils import Sorry, user_plus_sys_time, date_and_time
from libtbx.str_utils import format_value, show_string
import libtbx.path
from cStringIO import StringIO
import iotbx.phil
from mmtbx.scaling import outlier_rejection
from mmtbx.scaling import absolute_scaling
import mmtbx.scaling.twin_analyses
from cctbx import sgtbx


ext = boost.python.import_ext("mmtbx_f_model_ext")

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
  print >> out, "  Micro-tasks:"
  print >> out, "    mask                            = %-7.2f" % (time_mask)
  print >> out, "    f_calc                          = %-7.2f" % time_f_calc
  print >> out, "    alpha_beta                      = %-7.2f" % time_alpha_beta
  print >> out, "    target                          = %-7.2f" % time_target
  print >> out, "    gradients_wrt_atomic_parameters = %-7.2f" % \
    time_gradients_wrt_atomic_parameters
  print >> out, "    fmodel                         = %-7.2f" % time_fmodel_core_data
  print >> out, "    r_factors                      = %-7.2f" % time_r_factors
  print >> out, "    phase_errors                   = %-7.2f" % time_phase_errors
  print >> out, "    foms                           = %-7.2f" % time_foms
  print >> out, "    TOTAL for micro-tasks          = %-7.2f" % total
  return total

class set_core(object):
  def __init__(self, f_calc,
                     f_mask,
                     b_cart,
                     k_sol,
                     b_sol,
                     overall_scale,
                     uc,
                     ss,
                     work,
                     test):
    adopt_init_args(self, locals())
    global time_fmodel_core_data
    timer = user_plus_sys_time()
    self.f_mask.indices().all_eq(self.f_calc.indices())
    self.core = ext.core(f_calc        = self.f_calc.data(),
                         f_mask        = self.f_mask.data(),
                         b_cart        = self.b_cart,
                         k_sol         = self.k_sol,
                         b_sol         = self.b_sol,
                         overall_scale = self.overall_scale,
                         hkl           = self.f_calc.indices(),
                         uc            = self.uc,
                         ss            = self.ss)
    self.f_model   = miller.array(miller_set = self.f_calc,
                                  data       = self.core.f_model)
    self.f_model_w = self.f_model.select(self.work)
    self.f_model_t = self.f_model.select(self.test)
    self.fb_cart_w = self.core.fb_cart.select(self.work)
    self.fb_cart_t = self.core.fb_cart.select(self.test)
    time_fmodel_core_data += timer.elapsed()

  def __getstate__(self):
    return {"args": (
      self.f_calc,
      self.f_mask,
      self.b_cart,
      self.k_sol,
      self.b_sol,
      self.overall_scale,
      self.uc,
      self.ss,
      self.work,
      self.test)}

  def __setstate__(self, state):
    assert len(state) == 1
    self.__init__(*state["args"])

class target_attributes(object):

  def __init__(self,
        family,
        specialization=None,
        requires_external_scale=False):
    adopt_init_args(self, locals())
    assert self.validate()

  def validate(self):
    if (self.family == "lsm"):
      self.family = "ls"
      self.pseudo_ml = True
    else:
      self.pseudo_ml = False
    if (self.family == "ls"):
      return self.specialization in [None, "k1", "k2"]
    if (self.family == "ml"):
      return self.specialization in [None, "hl", "sad"]
    return False

  def requires_experimental_phases(self):
    return (self.family == "ml" and self.specialization == "hl")

  def ls_apply_scale_to_f_calc(self):
    if (self.family == "ls"):
      return (self.specialization == "k1")
    return None

class manager_mixin(object):

  def target_w(self):
    global time_target
    timer = user_plus_sys_time()
    result = self.target_functor()(compute_gradients=False).target_work()
    time_target += timer.elapsed()
    return result

  def target_t(self):
    global time_target
    timer = user_plus_sys_time()
    result = self.target_functor()(compute_gradients=False).target_work()
    time_target += timer.elapsed()
    return result

  def one_time_gradients_wrt_atomic_parameters(self, **keyword_args):
    return self.target_functor()(compute_gradients=True) \
      .gradients_wrt_atomic_parameters(**keyword_args)

sf_and_grads_accuracy_params = iotbx.phil.parse("""\
  algorithm = *fft direct
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
""")

alpha_beta_params = iotbx.phil.parse("""\
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

class manager(manager_mixin):

  target_names = {
    "ls_wunit_k1": target_attributes("ls", "k1"),
    "ls_wunit_k2": target_attributes("ls", "k2"),
    "ls_wunit_kunit": target_attributes("ls", "k1", True),
    "ls_wunit_k1_fixed": target_attributes("ls", "k1"),
    "ls_wunit_k1ask3_fixed": target_attributes("ls", "k1", True),
    "ls_wexp_k1": target_attributes("ls", "k1"),
    "ls_wexp_k2": target_attributes("ls", "k2"),
    "ls_wexp_kunit": target_attributes("ls", "k1", True),
    "ls_wff_k1": target_attributes("ls", "k1"),
    "ls_wff_k2": target_attributes("ls", "k2"),
    "ls_wff_kunit": target_attributes("ls", "k1", True),
    "ls_wff_k1_fixed": target_attributes("ls", "k1"),
    "ls_wff_k1ask3_fixed": target_attributes("ls", "k1", True),
    "lsm_k1": target_attributes("lsm", "k1"),
    "lsm_k2": target_attributes("lsm", "k2"),
    "lsm_kunit": target_attributes("lsm", "k1", True),
    "lsm_k1_fixed": target_attributes("lsm", "k1"),
    "lsm_k1ask3_fixed": target_attributes("lsm", "k1", True),
    "ml": target_attributes("ml"),
    "mlhl": target_attributes("ml", "hl"),
    "ml_sad": target_attributes("ml", "sad")}

  def __init__(self,
         f_obs                        = None,
         r_free_flags                 = None,
         b_cart                       = [0.,0.,0.,0.,0.,0.],
         k_sol                        = 0.0,
         b_sol                        = 0.0,
         sf_and_grads_accuracy_params = sf_and_grads_accuracy_params.extract(),
         target_name                  = None,
         abcd                         = None,
         alpha_beta_params            = alpha_beta_params.extract(),
         xray_structure               = None,
         f_mask                       = None,
         f_calc                       = None,
         mask_params                  = None,
         trust_xray_structure         = False,
         update_xray_structure        = True,
         use_f_model_scaled           = False,
         overall_scale                = 1.0,
         f_ordered_solvent            = None,
         f_ordered_solvent_dist       = None,
         mask_manager                 = None,
         n_ordered_water              = 0,
         b_ordered_water              = 0.0,
         max_number_of_bins           = 30):
    self.twin = False
    assert f_obs is not None
    assert f_obs.is_real_array()
    self.f_obs             = f_obs
    self.r_free_flags      = None
    self.sfg_params = sf_and_grads_accuracy_params
    self.abcd              = abcd
    self.alpha_beta_params = alpha_beta_params
    self.xray_structure    = xray_structure
    self.use_f_model_scaled= use_f_model_scaled
    self.overall_scale     = overall_scale
    self.max_number_of_bins = max_number_of_bins
    self.mask_manager = mask_manager
    if (f_ordered_solvent is None):
      self.f_ordered_solvent = self.f_obs.array(
        data=flex.complex_double(self.f_obs.data().size(), 0.0))
    else:
      self.f_ordered_solvent = f_ordered_solvent
    if (f_ordered_solvent_dist is None):
      self.f_ordered_solvent_dist = self.f_obs.array(
        data=flex.complex_double(self.f_obs.data().size(), 0.0))
    else:
      self.f_ordered_solvent_dist = f_ordered_solvent_dist
    self.n_ordered_water = n_ordered_water
    self.b_ordered_water = b_ordered_water
    if(mask_params is not None):
       self.mask_params = mask_params
    else:
       self.mask_params = mmtbx.masks.mask_master_params.extract()
    if(r_free_flags is not None):
       assert r_free_flags.indices().all_eq(self.f_obs.indices())
       self.update_r_free_flags(r_free_flags)
    self.f_obs_w = self.f_obs.select(self.work)
    self.f_obs_t = self.f_obs.select(self.test)
    self.uc = self.f_obs.unit_cell()
    self.d_spacings = self.f_obs.d_spacings().data()
    self.d_spacings_w = self.d_spacings.select(self.work)
    self.d_spacings_t = self.d_spacings.select(self.test)
    self.ss = 1./flex.pow2(self.d_spacings) / 4.
    if(self.xray_structure is None):
       assert [f_calc, f_mask].count(None) == 0
       assert f_mask.is_complex_array()
       assert f_calc.is_complex_array()
       assert f_mask.indices().all_eq(self.f_obs.indices())
       assert f_calc.indices().all_eq(self.f_obs.indices())
       self.update_core(f_calc       = f_calc,
                        f_mask       = f_mask,
                        b_cart       = b_cart,
                        k_sol        = k_sol,
                        b_sol        = b_sol)
    else:
       if(self.mask_manager is None):
         self.mask_manager = masks.manager(
           miller_array   = self.f_obs,
           xray_structure = self.xray_structure,
           mask_params    = self.mask_params)
       if(not trust_xray_structure):
          assert [f_calc, f_mask].count(None) == 2
       if(update_xray_structure):
          self.update_xray_structure(xray_structure       = self.xray_structure,
                                     update_f_calc        = True,
                                     update_f_mask        = True,
                                     k_sol                = k_sol,
                                     b_sol                = b_sol,
                                     b_cart               = b_cart)
       else:
          self.update_core(f_calc       = f_calc,
                           f_mask       = f_mask,
                           b_cart       = b_cart,
                           k_sol        = k_sol,
                           b_sol        = b_sol)
    assert len(b_cart) == 6
    if(self.abcd is not None):
       assert self.abcd.indices().all_eq(self.f_obs.indices()) == 1
    self.set_target_name(target_name=target_name)
    self._structure_factor_gradients_w = None
    self._wilson_b = None

  def __getstate__(self):
    self._structure_factor_gradients_w = None
    return (self.__dict__,)

  def __setstate__(self, state):
    assert len(state) == 1
    self.__dict__.update(state[0])

  def _get_structure_factor_gradients_w(self):
    if (self._structure_factor_gradients_w is None):
      self._structure_factor_gradients_w \
        = cctbx.xray.structure_factors.gradients(
         miller_set                   = self.f_obs_w,
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

  def update_core(self, f_calc        = None,
                        f_mask        = None,
                        b_cart        = None,
                        k_sol         = None,
                        b_sol         = None,
                        r_free_flags  = None):
    if(f_calc is not None): f_calc_ = f_calc
    else: f_calc_ = self.f_calc()
    if(f_mask is not None): f_mask_ = f_mask
    else: f_mask_ = self.f_mask()
    if(b_cart is not None): b_cart_ = b_cart
    else: b_cart_ = self.b_cart()
    if(k_sol is not None): k_sol_ = k_sol
    else: k_sol_ = self.k_sol()
    if(b_sol is not None): b_sol_ = b_sol
    else: b_sol_ = self.b_sol()
    if(r_free_flags is not None):
       work = ~r_free_flags.data()
       test =  r_free_flags.data()
    else:
       work = self.work
       test = self.test
    self.core = set_core(f_calc        = f_calc_,
                         f_mask        = f_mask_,
                         b_cart        = b_cart_,
                         k_sol         = k_sol_,
                         b_sol         = b_sol_,
                         overall_scale = self.overall_scale,
                         uc            = self.uc,
                         ss            = self.ss,
                         work          = work,
                         test          = test)

  def outlier_selection(self, show = False, log = None):
    if(log is None): log = sys.stdout
    result = outlier_rejection.outlier_manager(
      miller_obs   = self.f_obs,
      r_free_flags = self.r_free_flags,
      out          = "silent")
    s1 = result.basic_wilson_outliers().data()
    s2 = result.extreme_wilson_outliers().data()
    s3 = result.beamstop_shadow_outliers().data()
    result = s1 & s2 & s3
    if(show):
      print >> log
      print >> log, "basic_wilson_outliers    =", s1.count(False)
      print >> log, "extreme_wilson_outliers  =", s2.count(False)
      print >> log, "beamstop_shadow_outliers =", s3.count(False)
      print >> log, "total                    =", result.count(False)
    return result

  def remove_outliers(self, show = False, log = None):
    if(log is None): log = sys.stdout
    if(show):
      print >> log, "Distribution of F-obs values before outliers rejection:"
      show_histogram(data = self.f_obs.data(), n_slots = 10, log = log)
    result = self.select(selection = self.outlier_selection(show=show,log=log))
    if(show):
      print >> log, "\nDistribution of F-obs values after outliers rejection:"
      show_histogram(data = result.f_obs.data(), n_slots = 10, log = log)
    return result

  def wilson_b(self, force_update = False):
    if(self._wilson_b is None or force_update):
      result = absolute_scaling.ml_iso_absolute_scaling(
        miller_array = self.f_obs,
        n_residues   = self.xray_structure.scatterers().size()/8).b_wilson
      self._wilson_b = result
    else: result = self._wilson_b
    return result

  def twin_test(self, cut_off = 3.5):
    result = None
    if(self.f_obs.d_min() < 9.0):
      result = mmtbx.scaling.twin_analyses.twin_analyses_brief(
        miller_array = self.f_obs, cut_off = cut_off)
      if(result): result = "twinned"
      elif(result is None): result = "unknown"
      elif(not result): result = "not twinned"
      else: raise Sorry("Twin analysis failed.")
    return result

  def deep_copy(self):
    if(self.abcd is not None):
       abcd = self.abcd.deep_copy()
    else:
       abcd = None
    if(self.xray_structure is None):
       xrs = None
    else:
       xrs = self.xray_structure.deep_copy_scatterers()
    if(self.mask_manager is not None):
      new_mask_manager = self.mask_manager.deep_copy()
    else: new_mask_manager = None
    return manager(
              f_obs                        = self.f_obs.deep_copy(),
              r_free_flags                 = self.r_free_flags.deep_copy(),
              b_cart                       = self.b_cart(),
              f_mask                       = self.f_mask().deep_copy(),
              k_sol                        = self.k_sol(),
              b_sol                        = self.b_sol(),
              sf_and_grads_accuracy_params = self.sfg_params,  # XXX deep_copy ?
              target_name                  = self.target_name,
              abcd                         = abcd,
              alpha_beta_params            = self.alpha_beta_params,
              xray_structure               = xrs,
              mask_manager                 = new_mask_manager,
              f_calc                       = self.f_calc().deep_copy(),
              mask_params                  = self.mask_params,  # XXX deep_copy ?
              trust_xray_structure         = True,
              update_xray_structure        = False,
      overall_scale         = self.overall_scale,
      f_ordered_solvent      = self.f_ordered_solvent.deep_copy(),
      f_ordered_solvent_dist = self.f_ordered_solvent_dist.deep_copy(),
      n_ordered_water        = self.n_ordered_water,
      b_ordered_water        = self.b_ordered_water,
      max_number_of_bins = self.max_number_of_bins)

  def select(self, selection, update_xray_structure=False):
    if (self.abcd is not None):
      abcd = self.abcd.select(selection=selection)
    else:
      abcd = None
    if(self.mask_manager is not None):
      new_mask_manager = self.mask_manager.select(selection = selection)
    else:
      new_mask_manager = None
    result = manager(
      f_obs=self.f_obs.select(selection=selection),
      r_free_flags=self.r_free_flags.select(selection=selection),
      b_cart=tuple(self.b_cart()),
      k_sol=self.k_sol(),
      b_sol=self.b_sol(),
      sf_and_grads_accuracy_params = self.sfg_params,
      target_name=self.target_name,
      abcd=abcd,
      alpha_beta_params=self.alpha_beta_params,
      xray_structure=self.xray_structure,
      f_calc=self.f_calc().select(selection=selection),
      f_mask=self.f_mask().select(selection=selection),
      mask_params=self.mask_params,
      mask_manager = new_mask_manager,
      trust_xray_structure=True,
      update_xray_structure=update_xray_structure,
      overall_scale=self.overall_scale,
      f_ordered_solvent=self.f_ordered_solvent.select(
        selection=selection),
      f_ordered_solvent_dist=self.f_ordered_solvent_dist.select(
        selection=selection),
      n_ordered_water=self.n_ordered_water,
      b_ordered_water=self.b_ordered_water,
      max_number_of_bins=self.max_number_of_bins)
    return result

  def resolution_filter(self,
        d_max=0,
        d_min=0,
        update_xray_structure=False):
    return self.select(
      selection=self.f_obs.resolution_filter_selection(
        d_max=d_max, d_min=d_min),
      update_xray_structure=update_xray_structure)

  def apply_back_b_iso(self):
    eps = math.pi**2*8
    unit_cell = self.xray_structure.unit_cell()
    b_min = min(self.b_sol(), self.xray_structure.min_u_cart_eigenvalue())
    if(b_min < 0):
      self.xray_structure.tidy_us(u_min = 1.e-6)
    b_iso = self.b_iso()
    b_test = b_min+b_iso
    if(b_test < 0.0): b_adj = b_iso + abs(b_test) + 0.001
    else: b_adj = b_iso
    if(abs(b_adj) <= 300.0):
      b_cart = self.b_cart()
      b_cart_new = [b_cart[0]-b_adj,b_cart[1]-b_adj,b_cart[2]-b_adj,
                    b_cart[3],      b_cart[4],      b_cart[5]]
      self.update(b_cart = b_cart_new)
      self.update(b_sol = self.k_sol_b_sol()[1] + b_adj)
      self.xray_structure.shift_us(b_shift = b_adj)
      b_min = min(self.b_sol(), self.xray_structure.min_u_cart_eigenvalue())
      assert b_min >= 0.0
      self.xray_structure.tidy_us(u_min = 1.e-6)
      self.update_xray_structure(
        xray_structure           = self.xray_structure,
        update_f_calc            = True,
        update_f_mask            = False,
        update_f_ordered_solvent = False,
        out                      = None)

  def set_f_ordered_solvent(self, params):
    raise RuntimeError("Not implemented.")
    if(params.nu_fix_b_atoms is not None):
       self.n_ordered_water = params.nu_fix_n_atoms
       self.b_ordered_water = params.nu_fix_b_atoms
       self.f_ordered_solvent = max_like_non_uniform.f_ordered_solvent(
                            f                    = self.f_ordered_solvent_dist,
                            n_water_atoms_absent = self.n_ordered_water,
                            bf_atoms_absent      = self.b_ordered_water,
                            absent_atom_type     = "O")
    else:
       r = self.target_w()
       f_ordered_solvent = self.f_ordered_solvent
       n_ordered_water   = self.n_ordered_water
       b_ordered_water   = self.b_ordered_water
       n_atoms_prot = self.xray_structure.scatterers().size()
       n_residues = n_atoms_prot / 10
       n_solvent_max = n_residues * 2
       n_solvent_min = n_residues / 2
       u_isos = self.xray_structure.extract_u_iso_or_u_equiv()
       b_iso_mean = flex.mean(u_isos * math.pi**2*8)
       b_solvent_max = int(b_iso_mean + 35.0)
       b_solvent_min = int(b_iso_mean - 5.0)
       for n in range(n_solvent_min, n_solvent_max+1, 10):
           for b in range(b_solvent_min, b_solvent_max+1, 5):
               self.f_ordered_solvent = max_like_non_uniform.f_ordered_solvent(
                            f                    = self.f_ordered_solvent_dist,
                            n_water_atoms_absent = n,
                            bf_atoms_absent      = b,
                            absent_atom_type     = "O")
               r_i = self.target_w()
               if(r_i < r):
                  r = r_i
                  f_ordered_solvent = self.f_ordered_solvent
                  n_ordered_water = n
                  b_ordered_water = b
       self.n_ordered_water = n_ordered_water
       self.b_ordered_water = b_ordered_water
       self.f_ordered_solvent = f_ordered_solvent
       assert approx_equal(self.target_w(), r)
       ############## ????
       self.alpha_beta_params.n_water_atoms_absent = self.n_ordered_water
       self.alpha_beta_params.bf_atoms_absent = self.b_ordered_water

  def update_xray_structure(self,
                            xray_structure           = None,
                            update_f_calc            = False,
                            update_f_mask            = False,
                            update_f_ordered_solvent = False,
                            force_update_f_mask      = False,
                            out                      = None,
                            k_sol                    = None,
                            b_sol                    = None,
                            b_cart                   = None):
    if (xray_structure is not None):
      self.xray_structure = xray_structure
    f_calc = None
    if(update_f_calc):
       global time_f_calc
       timer = user_plus_sys_time()
       assert self.xray_structure is not None
       f_calc = self.f_obs.structure_factors_from_scatterers(
         xray_structure = self.xray_structure,
         algorithm                    = self.sfg_params.algorithm,
         cos_sin_table                = self.sfg_params.cos_sin_table,
         grid_resolution_factor       = self.sfg_params.grid_resolution_factor,
         quality_factor               = self.sfg_params.quality_factor,
         u_base                       = self.sfg_params.u_base,
         b_base                       = self.sfg_params.b_base,
         wing_cutoff                  = self.sfg_params.wing_cutoff,
         exp_table_one_over_step_size =
                         self.sfg_params.exp_table_one_over_step_size
                         ).f_calc()
       time_f_calc += timer.elapsed()
    if(update_f_ordered_solvent):
       nu_manager = max_like_non_uniform.ordered_solvent_distribution(
                                               structure = self.xray_structure,
                                               fo        = self.f_obs,
                                               grid_step = step,
                                               rad       = 0.0)
       nu_map = nu_manager.distribution_as_array()
       self.f_ordered_solvent_dist = nu_manager.fcalc_from_distribution()
    f_mask = None
    if(update_f_mask):
       global time_mask
       timer = user_plus_sys_time()
       f_mask = self.mask_manager.f_mask(
         xray_structure_new = self.xray_structure,
         force_update       = force_update_f_mask)
       time_mask += timer.elapsed()
    if([f_calc, f_mask].count(None) == 2): set_core_flag = False
    else: set_core_flag = True
    if(f_calc is None): f_calc = self.f_calc()
    if(f_mask is None): f_mask = self.f_mask()
    if(set_core_flag):
       self.update_core(f_calc = f_calc,
                        f_mask = f_mask,
                        b_cart = b_cart,
                        k_sol  = k_sol,
                        b_sol  = b_sol)

  def optimize_mask_and_update_solvent_and_scale(
                                self, params = None, out = None, verbose=-1):
    rw = self.r_work()
    r_solv_   = self.mask_params.solvent_radius
    r_shrink_ = self.mask_params.shrink_truncation_radius
    if(verbose > 0):
       self.show_mask_optimization_statistics(prefix="Mask optimization start",
                                              out   = out)
    for r_solv in [0.8, 1.0, 1.2, 1.4]:
        for r_shrink in [0.8, 1.0, 1.2, 1.4]:
          self.mask_params.solvent_radius = r_solv
          self.mask_params.shrink_truncation_radius = r_shrink
          self.update_xray_structure(
                                xray_structure           = self.xray_structure,
                                update_f_calc            = False,
                                update_f_mask            = True,
                                update_f_ordered_solvent = False,
                                force_update_f_mask      = True,
                                out                      = None)
          self.update_solvent_and_scale(params=params, out=None, verbose=-1)
          rw_ = self.r_free()
          if(rw_ is not None and (rw is None or rw_ < rw)):
             rw = rw_
             r_solv_ = r_solv
             r_shrink_ = r_shrink
    self.mask_params.solvent_radius = r_solv_
    self.mask_params.shrink_truncation_radius = r_shrink_
    self.update_xray_structure(xray_structure           = self.xray_structure,
                               update_f_calc            = False,
                               update_f_mask            = True,
                               update_f_ordered_solvent = False,
                               force_update_f_mask      = True,
                               out                      = None)
    self.update_solvent_and_scale(params = params, out = out, verbose = -1)
    if(verbose > 0):
       self.show_mask_optimization_statistics(prefix="Mask optimization final",
                                              out   = out)

  def show_mask_optimization_statistics(self, prefix="", out=None):
    if(out is None): out = sys.stdout
    line_len = len("|-"+"|"+prefix)
    fill_len = 80-line_len-1
    print >> out, "|-"+prefix+"-"*(fill_len)+"|"
    print >> out, "| r_work= %s     r_free= %s     Rad_solv= %4.2f     Rad_shrink= %4.2f   |"%\
     (format_value(format="%6.4f", value=self.r_work()),
      format_value(format="%6.4f", value=self.r_free()),
      self.mask_params.solvent_radius,
      self.mask_params.shrink_truncation_radius)
    print >> out, "|"+"-"*77+"|"
    print >> out

  def update_solvent_and_scale(self, params = None, out = None, verbose=None):
    global time_bulk_solvent_and_scale
    timer = user_plus_sys_time()
    if(out is None): out = sys.stdout
    if(params is None):
       params = bss.master_params.extract()
    if(verbose is not None): params.verbose=verbose
    bss.bulk_solvent_and_scales(fmodel = self, params = params, log = out)
    overall_scale = self.scale_k1_w()
    if(self.use_f_model_scaled):
       self.overall_scale = overall_scale
       self.update_core()
    elif(overall_scale < 0.01 or overall_scale > 100.0):
       self.overall_scale = overall_scale
       self.update_core()
    time_bulk_solvent_and_scale += timer.elapsed()

  def _get_target_name(self): return self._target_name
  target_name = property(_get_target_name)

  def target_attributes(self):
    if (self.target_name is None): return None
    try: return manager.target_names[self.target_name]
    except AttributeError:
      raise RuntimeError(
        "Unknown target name: %s" % show_string(self.target_name))

  def target_functor(self):
    return target_functor(manager=self)

  def set_target_name(self, target_name):
    if (target_name == "ls"): target_name = "ls_wunit_k1"
    self._target_name = target_name

  def update_r_free_flags(self, r_free_flags):
    assert r_free_flags.indices().size() == self.f_obs.indices().size()
    self.r_free_flags = r_free_flags
    self.work = ~r_free_flags.data()
    self.test =  r_free_flags.data()
    if (self.work.count(True) == 0): self.work = ~self.work # XXX BAD
    if (self.test.count(True) == 0): self.test = ~self.test # XXX BAD
    self.work_count_true = self.work.count(True)
    self.test_count_true = self.test.size() - self.work_count_true
    self.test_count_true = self.test.count(True) # XXX BAD
    #assert self.work_count_true != 0 XXX fails some tests
    #assert self.test_count_true != 0

  def determine_n_bins(self,
        free_reflections_per_bin,
        max_n_bins=None,
        min_n_bins=1,
        min_refl_per_bin=100):
    assert free_reflections_per_bin > 0
    n_refl = self.test.size()
    n_free = self.test_count_true
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

  def update(self, f_calc                       = None,
                   f_obs                        = None,
                   f_mask                       = None,
                   f_ordered_solvent            = None,
                   r_free_flags                 = None,
                   b_cart                       = None,
                   k_sol                        = None,
                   b_sol                        = None,
                   sf_and_grads_accuracy_params = None,
                   target_name                  = None,
                   abcd                         = None,
                   alpha_beta_params            = None,
                   xray_structure               = None,
                   mask_params                  = None,
                   overall_scale                = None):
    if(f_calc is not None):
       assert f_calc.indices().all_eq(self.core.f_calc.indices())
       self.update_core(f_calc = f_calc)
    if(mask_params is not None):
       self.mask_params = mask_params
    if(f_obs is not None):
       assert f_obs.data().size() == self.f_obs.data().size()
       self.f_obs = f_obs
       self.f_obs_w = self.f_obs.select(self.work)
       self.f_obs_t = self.f_obs.select(self.test)
    if(f_mask is not None):
      assert f_mask.data().size() == self.f_mask().data().size()
      self.update_core(f_mask = f_mask)
    if(f_ordered_solvent is not None):
       if(self.f_ordered_solvent is not None):
          assert f_ordered_solvent.data().size() == \
                 self.f_ordered_solvent.data().size()
       self.f_ordered_solvent = f_ordered_solvent
    if(r_free_flags is not None):
      self.update_r_free_flags(r_free_flags)
      self.update_core(r_free_flags = r_free_flags)
    if(b_cart is not None):
      try: assert b_cart.size() == 6
      except: assert len(b_cart) == 6
      self.update_core(b_cart = b_cart)
    if(k_sol is not None):
       self.update_core(k_sol = k_sol)
    if(b_sol is not None):
       self.update_core(b_sol = b_sol)
    if(sf_and_grads_accuracy_params is not None):
      self.sfg_params = sf_and_grads_accuracy_params
      self.update_xray_structure(update_f_calc  = True)
    if(target_name is not None):
      self.set_target_name(target_name=target_name)
    if(abcd is not None):
      self.abcd = abcd
    if(alpha_beta_params is not None):
      self.alpha_beta_params = alpha_beta_params
    if(xray_structure is not None):
       self.update_xray_structure(xray_structure = xray_structure,
                                  update_f_mask  = True,
                                  update_f_calc  = True)
    if(overall_scale is not None):
       self.overall_scale = overall_scale
       self.update_core()
    return self

  def f_ordered_solvent_w(self):
    assert self.r_free_flags is not None
    if(self.r_free_flags.data().count(True) > 0):
      return self.f_ordered_solvent.select(~self.r_free_flags.data())
    else:
      return self.f_ordered_solvent

  def f_ordered_solvent_t(self):
    assert self.r_free_flags is not None
    if(self.r_free_flags.data().count(True) > 0):
      return self.f_ordered_solvent.select(self.r_free_flags.data())
    else:
      return self.f_ordered_solvent

  def f_bulk(self):
    return miller.array(miller_set = self.f_obs, data = self.core.core.f_bulk)

  def f_bulk_w(self):
    if(self.r_free_flags.data().count(True) > 0):
      return self.f_bulk().select(self.work)
    else:
      return self.f_bulk()

  def f_bulk_t(self):
    if(self.r_free_flags.data().count(True) > 0):
      return self.f_bulk().select(self.test)
    else:
      return self.f_bulk()

  def fb_cart(self):
    return self.core.core.fb_cart

  def fb_cart_w(self):
    return self.fb_cart().select(self.work)

  def fb_cart_t(self):
    return self.fb_cart().select(self.test)

  def f_obs_scaled_with_k2(self):
    scale_k2 = self.scale_k2()
    f_obs = self.f_obs
    d = f_obs.data() * scale_k2
    s = f_obs.sigmas()
    if (s is not None): s = s * scale_k2
    return f_obs.array(data=d, sigmas=s)

  def f_model(self):
    return self.core.f_model

  def f_model_scaled_with_k1(self):
    return miller.array(miller_set = self.f_obs,
                        data       = self.scale_k1()*self.f_model().data())

  def f_model_scaled_with_k1_t(self):
    return miller.array(miller_set = self.f_obs_t,
                        data       = self.scale_k1_t()*self.f_model_t().data())

  def f_model_scaled_with_k1_w(self):
    return miller.array(miller_set = self.f_obs_w,
                        data       = self.scale_k1_w()*self.f_model_w().data())

  def f_model_w(self):
    return self.core.f_model_w

  def f_model_t(self):
    return self.core.f_model_t

  def f_star_w_star_obj(self):
    #XXX why I use self.f_calc and not f_model ????????????????????????????????
    alpha, beta = self.alpha_beta()
    obj = max_lik.f_star_w_star_mu_nu(
                                 f_obs          = self.f_obs.data(),
                                 f_model        = flex.abs(self.f_calc().data()),
                                 alpha          = alpha.data(),
                                 beta           = beta.data(),
                                 space_group    = self.f_obs.space_group(),
                                 miller_indices = self.f_obs.indices())
    return obj

  def f_star_w_star(self):
    obj = self.f_star_w_star_obj()
    f_star = miller.array(miller_set = self.f_obs,
                          data       = obj.f_star())
    w_star = miller.array(miller_set = self.f_obs,
                          data       = obj.w_star())
    return f_star, w_star

  def f_star_w_star_work(self):
    assert self.r_free_flags is not None
    f_star, w_star = self.f_star_w_star()
    flags = self.r_free_flags.data()
    if(flags.count(True) > 0):
       return f_star.select(~flags), w_star.select(~flags)
    else:
       return f_star, w_star

  def f_star_w_star_test(self):
    assert self.r_free_flags is not None
    f_star, w_star = self.f_star_w_star()
    flags = self.r_free_flags.data()
    if(flags.count(True) > 0):
       return f_star.select(flags), w_star.select(flags)
    else:
       return f_star, w_star

  def b_cart(self):
    return self.core.core.b_cart

  def b_iso(self):
    b_cart = self.b_cart()
    return (b_cart[0]+b_cart[1]+b_cart[2])/3.0

  def f_mask(self):
    return self.core.f_mask

  def f_mask_w(self):
    assert self.r_free_flags is not None
    if(self.r_free_flags.data().count(True) > 0):
      return self.f_mask().select(~self.r_free_flags.data())
    else:
      return self.f_mask()

  def f_mask_t(self):
    assert self.r_free_flags is not None
    if(self.r_free_flags.data().count(True) > 0):
      return self.f_mask().select(self.r_free_flags.data())
    else:
      return self.f_mask()

  def f_calc(self):
    return self.core.f_calc

  def f_calc_w(self):
    assert self.r_free_flags is not None
    if(self.r_free_flags.data().count(True) > 0):
      return self.f_calc().select(~self.r_free_flags.data())
    else:
      return self.f_calc()

  def f_calc_t(self):
    assert self.r_free_flags is not None
    if(self.r_free_flags.data().count(True) > 0):
      return self.f_calc().select(self.r_free_flags.data())
    else:
      return self.f_calc()

  def k_sol(self):
    return self.core.core.k_sol

  def b_sol(self):
    return self.core.core.b_sol

  def k_sol_b_sol(self):
    return self.k_sol(), self.b_sol()

  def alpha_beta(self, f_obs = None, f_model = None):
    global time_alpha_beta
    timer = user_plus_sys_time()
    if(f_obs is None): f_obs = self.f_obs
    if(f_model is None): f_model = self.f_model()
    alpha, beta = None, None
    ab_params = self.alpha_beta_params
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
             flags           = self.r_free_flags.data(),
             interpolation   = self.alpha_beta_params.interpolation) \
               .alpha_beta()
         else:
           p = self.alpha_beta_params.sigmaa_estimator
           alpha, beta = sigmaa_estimator(
             miller_obs  = f_obs,
             miller_calc = f_model,
             r_free_flags=self.r_free_flags,
             kernel_width_free_reflections=p.kernel_width_free_reflections,
             kernel_on_chebyshev_nodes=p.kernel_on_chebyshev_nodes,
             n_sampling_points=p.number_of_sampling_points,
             n_chebyshev_terms=p.number_of_chebyshev_terms,
             use_sampling_sum_weights=p.use_sampling_sum_weights).alpha_beta()
       else:
         check = flex.max(flex.abs(self.f_ordered_solvent_dist.data()))
         if(check < 1.e-9):
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
            alpha, beta = max_like_non_uniform.alpha_beta(
             f_dist                   = self.f_ordered_solvent_dist,
             n_atoms_included         = ab_params.n_atoms_included,
             n_nonwater_atoms_absent  = ab_params.number_of_macromolecule_atoms_absent,
             n_water_atoms_absent     = ab_params.number_of_waters_absent,
             bf_atoms_absent          = ab_params.bf_atoms_absent,
             final_error              = ab_params.final_error,
             absent_atom_type         = ab_params.absent_atom_type)
    else:
       alpha, beta = maxlik.alpha_beta_est_manager(
                                    f_obs           = f_obs,
                                    f_calc          = f_model,
                                    free_reflections_per_bin = 140,
                                    flags           = self.r_free_flags.data(),
                                    interpolation   = False).alpha_beta()
    time_alpha_beta += timer.elapsed()
    return alpha, beta

  def alpha_beta_w(self, only_if_required_by_target=False):
    if (only_if_required_by_target):
      if (self.target_name not in ["ml", "mlhl"]): return None, None
    assert self.r_free_flags is not None
    alpha, beta = self.alpha_beta()
    if(self.r_free_flags.data().count(True) > 0):
      return alpha.select(~self.r_free_flags.data()), \
             beta.select(~self.r_free_flags.data())
    else:
      return alpha, beta

  def alpha_beta_t(self):
    assert self.r_free_flags is not None
    alpha, beta = self.alpha_beta()
    if(self.r_free_flags.data().count(True) > 0):
      return alpha.select(self.r_free_flags.data()), \
             beta.select(self.r_free_flags.data())
    else:
      return alpha, beta

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
    fmat6 = self.resolution_filter(d_max=6.0)
    for fmodel in [
          fmat6,
          self]:
      ss = 1./flex.pow2(fmodel.f_obs.d_spacings().data())
      if (estimation_algorithm == "analytical"):
        try:
          alpha, beta = maxlik.alpha_beta_est_manager(
            f_obs           = fmodel.f_obs,
            f_calc          = fmodel.f_model(),
            free_reflections_per_bin = free_reflections_per_bin,
            flags           = fmodel.r_free_flags.data(),
            interpolation   = True).alpha_beta()
          break
        except KeyboardInterrupt: raise
        except Exception, e: est_exceptions.append(str(e))
      else:
        p = self.alpha_beta_params.sigmaa_estimator
        try:
          alpha, beta = sigmaa_estimator(
            miller_obs=fmodel.f_obs,
            miller_calc=fmodel.f_model(),
            r_free_flags=fmodel.r_free_flags,
            kernel_width_free_reflections=p.kernel_width_free_reflections,
            kernel_on_chebyshev_nodes=p.kernel_on_chebyshev_nodes,
            n_sampling_points=p.number_of_sampling_points,
            n_chebyshev_terms=p.number_of_chebyshev_terms,
            use_sampling_sum_weights=p.use_sampling_sum_weights).alpha_beta()
          break
        except KeyboardInterrupt: raise
        except Exception, e: est_exceptions.append(str(e))
    else:
      raise RuntimeError(
        "Failure estimating alpha, beta coefficients:\n"
        + est_exceptions[0] + "\n"
        + "  " + "-"*77 + "\n"
        + est_exceptions[1])
    save_self_overall_scale = fmodel.overall_scale
    fmodel.overall_scale = fmodel.scale_k3_w()
    fmodel.update_core()
    alpha = fmodel.alpha_beta()[0].data()
    omega = flex.double()
    for ae,ssi in zip(alpha,ss):
      if(ae >  1.0): ae = 1.0
      if(ae <= 0.0): ae = 1.e-6
      coeff = -4./(math.pi**3*ssi)
      omega.append( math.sqrt( math.log(ae) * coeff ) )
    #omega_ma  = miller.array(miller_set= self.f_obs,data= flex.double(omega))
    fmodel.overall_scale = save_self_overall_scale
    fmodel.update_core()
    omega_mean = flex.mean(omega)
    #sel = (omega < omega_mean * 3.0) & (omega > omega_mean / 3.0)
    #if(sel.count(True) > 0):
    #   omega_mean = flex.mean(omega.select(sel))
    return omega_mean
    #return flex.mean(omega), flex.max(omega), flex.min(omega)

  def _r_factor(self, f_obs, f_model, d_min=None, d_max=None, d_spacings=None,
                                                               selection=None):
    global time_r_factors
    if(selection is not None):
       assert [d_min, d_max].count(None) == 2
    if([d_min, d_max].count(None) < 2):
       assert selection is None and d_spacings is not None
    timer = user_plus_sys_time()
    if([d_min, d_max].count(None) == 0):
       keep = flex.bool(d_spacings.size(), True)
       if (d_max): keep &= d_spacings <= d_max
       if (d_min): keep &= d_spacings >= d_min
       f_obs   = f_obs.select(keep)
       f_model = f_model.select(keep)
    if(selection is not None):
       f_obs   = f_obs.select(selection)
       f_model = f_model.select(selection)
    result = bulk_solvent.r_factor(f_obs, f_model)
    time_r_factors += timer.elapsed()
    return result

  def r_work(self, d_min = None, d_max = None, selection = None):
    return self._r_factor(f_obs      = self.f_obs_w.data(),
                          f_model    = self.core.f_model_w.data(),
                          d_min      = d_min,
                          d_max      = d_max,
                          d_spacings = self.d_spacings_w,
                          selection  = selection)

  def r_free(self, d_min = None, d_max = None, selection = None):
    return self._r_factor(f_obs      = self.f_obs_t.data(),
                          f_model    = self.core.f_model_t.data(),
                          d_min      = d_min,
                          d_max      = d_max,
                          d_spacings = self.d_spacings_t,
                          selection  = selection)

  def r_all(self):
    return self._r_factor(f_obs = self.f_obs.data(),
                          f_model = self.core.f_model.data())

  def scale_k1(self, selection = None):
    return _scale_helper(
      num=self.f_obs.data(),
      den=flex.abs(self.core.f_model.data()),
      selection=selection)

  def scale_k1_w(self, selection = None):
    return _scale_helper(
      num=self.f_obs_w.data(),
      den=flex.abs(self.core.f_model_w.data()),
      selection=selection)

  def scale_k1_t(self, selection = None):
    return _scale_helper(
      num=self.f_obs_t.data(),
      den=flex.abs(self.core.f_model_t.data()),
      selection=selection)

  def scale_k2(self, selection = None):
    return _scale_helper(
      num=flex.abs(self.core.f_model.data()),
      den=self.f_obs.data(),
      selection=selection)

  def scale_k2_w(self, selection = None):
    return _scale_helper(
      num=flex.abs(self.core.f_model_w.data()),
      den=self.f_obs_w.data(),
      selection=selection)

  def scale_k2_t(self, selection = None):
    return _scale_helper(
      num=flex.abs(self.core.f_model_t.data()),
      den=self.f_obs_t.data(),
      selection=selection)

  def scale_k3_w(self, selection = None):
    mul_eps_sqrt = flex.sqrt(
        self.f_obs_w.multiplicities().data().as_double()
      / self.f_obs_w.epsilons().data().as_double())
    result = _scale_helper(
      num=self.f_obs_w.data() * mul_eps_sqrt,
      den=flex.abs(self.core.f_model_w.data()) * mul_eps_sqrt,
      selection=selection,
      num_num=True)
    if (result is None): return None
    return result**0.5

  def scale_k3_t(self, selection = None):
    mul_eps_sqrt = flex.sqrt(
        self.f_obs_t.multiplicities().data().as_double()
      / self.f_obs_t.epsilons().data().as_double())
    result = _scale_helper(
      num=self.f_obs_t.data() * mul_eps_sqrt,
      den=flex.abs(self.core.f_model_t.data()) * mul_eps_sqrt,
      selection=selection,
      num_num=True)
    if (result is None): return None
    return result**0.5

  def r_overall_low_high(self, d = 6.0):
    r_work = self.r_work()
    d_max, d_min = self.f_obs_w.d_max_min()
    if(d_max < d): d = d_max
    if(d_min > d): d = d_min
    n_low = self.f_obs_w.resolution_filter(d_min = d, d_max = 999.9).data().size()
    if(n_low > 0):
       r_work_l = self.r_work(d_min = d, d_max = 999.9)
    else:
       r_work_l = None
    n_high = self.f_obs_w.resolution_filter(d_min = 0.0, d_max = d).data().size()
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
    if (self.alpha_beta_params.fix_scale_for_calc_option is None):
      return self.scale_ml()
    return self.alpha_beta_params.fix_scale_for_calc_option

  def scale_ml(self):
    #assert self.alpha_beta_params.method == "calc"
    alpha, beta = self.alpha_beta_w()
    scale_manager = bss.uaniso_ksol_bsol_scaling_minimizer(
               self.f_calc_w(),
               self.f_obs_w,
               self.f_mask_w(),
               k_initial = 0.,
               b_initial = 0.,
               u_initial = [0,0,0,0,0,0],
               scale_initial = self.scale_k3_w(),
               refine_k = False,
               refine_b = False,
               refine_u = False,
               refine_scale = True,
               alpha = alpha.data(),
               beta = beta.data(),
               lbfgs_exception_handling_params = lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True,
                         ignore_line_search_failed_step_at_upper_bound = True,
                         ignore_line_search_failed_maxfev              = True))
    return scale_manager.scale_min

  def figures_of_merit(self):
    alpha, beta = self.alpha_beta()
    global time_foms
    timer = user_plus_sys_time()
    result = max_lik.fom_and_phase_error(
                           f_obs          = self.f_obs.data(),
                           f_model        = flex.abs(self.core.f_model.data()),
                           alpha          = alpha.data(),
                           beta           = beta.data(),
                           space_group    = self.f_obs.space_group(),
                           miller_indices = self.f_obs.indices()).fom()
    time_foms += timer.elapsed()
    return result

  def figures_of_merit_work(self):
    assert self.r_free_flags is not None
    fom = self.figures_of_merit()
    if(self.r_free_flags.data().count(True) > 0):
      return fom.select(~self.r_free_flags.data())
    else:
      return fom

  def phase_errors(self):
    alpha, beta = self.alpha_beta()
    global time_phase_errors
    timer = user_plus_sys_time()
    result = max_lik.fom_and_phase_error(
                           f_obs          = self.f_obs.data(),
                           f_model        = flex.abs(self.core.f_model.data()),
                           alpha          = alpha.data(),
                           beta           = beta.data(),
                           space_group    = self.f_obs.space_group(),
                           miller_indices = self.f_obs.indices()).phase_error()
    time_phase_errors += timer.elapsed()
    return result

  def phase_errors_test(self):
    assert self.r_free_flags is not None
    pher = self.phase_errors()
    if(self.r_free_flags.data().count(True) > 0):
      return pher.select(self.r_free_flags.data())
    else:
      return pher

  def phase_errors_work(self):
    assert self.r_free_flags is not None
    pher = self.phase_errors()
    if(self.r_free_flags.data().count(True) > 0):
      return pher.select(~self.r_free_flags.data())
    else:
      return pher

  def _map_coeff(self, f_obs, f_model, f_obs_scale, f_model_scale):
    d_obs = miller.array(miller_set = f_model,
                         data       = f_obs.data()*f_obs_scale
                        ).phase_transfer(phase_source = f_model)
    return miller.array(miller_set = f_model,
                        data       = d_obs.data()-f_model.data()*f_model_scale)

  def map_coefficients(self, map_type):
    map_name_manager = mmtbx.map_names(map_name_string = map_type)
    if(map_name_manager.anomalous and self.f_obs.anomalous_flag()):
      anom_diff = self.f_obs.anomalous_differences()
      fom = miller.array(miller_set = self.f_obs,
                         data       = self.figures_of_merit())
      fom = fom.as_non_anomalous_array().merge_equivalents().array()
      fom_match_anom_diff = fom.common_set(other = anom_diff)
      fmodel = self.f_model().as_non_anomalous_array().merge_equivalents().array()
      fc_match_anom_diff = fmodel.common_set(other = anom_diff)
      fc_match_anom_diff.indices().all_eq(anom_diff.indices())
      tmp = miller.array(miller_set = anom_diff,
                         data       = anom_diff.data()*fom_match_anom_diff.data()
                         ).phase_transfer(phase_source = fc_match_anom_diff)
      return miller.array(miller_set = tmp,
                          data = tmp.data()/(2j))
    #
    fb_cart  = self.fb_cart()
    scale_k2 = self.scale_k2()
    f_obs_scale   = 1.0 / fb_cart * scale_k2
    f_model_scale = 1.0 / fb_cart
    if(not map_name_manager.ml_map):
       return self._map_coeff(
         f_obs         = self.f_obs,
         f_model       = self.f_model(),
         f_obs_scale   = map_name_manager.k * f_obs_scale,
         f_model_scale = map_name_manager.n * f_model_scale)
    if(map_name_manager.ml_map):
      alpha, beta = self.alpha_beta(
        f_obs   = self.f_obs.array(data = self.f_obs.data() * f_obs_scale),
        f_model = self.f_model().array(data = self.f_model().data()*f_model_scale))
      return self._map_coeff(
        f_obs         = self.f_obs,
        f_model       = self.f_model(),
        f_obs_scale   = map_name_manager.k*self.figures_of_merit()*f_obs_scale,
        f_model_scale = map_name_manager.n*alpha.data() * f_model_scale)
      ####
      #result = miller.array(miller_set = self.f_calc,
      #                      data       = d_obs.data() - d_model)
      #centrics  = result.select_centric()
      #acentrics = result.select_acentric()
      #acentrics_data = acentrics.data() * 2.0
      #centrics_data  = centrics.data()
      #new = acentrics.customized_copy(
      #          indices = acentrics.indices().concatenate(centrics.indices()),
      #          data    = acentrics_data.concatenate(centrics_data) )
      ####
      #return new

  def electron_density_map(self,
                           map_type,
                           resolution_factor = 1/3.,
                           symmetry_flags = None):
    return self.map_coefficients(map_type = map_type).fft_map(
      resolution_factor = resolution_factor,
      symmetry_flags    = symmetry_flags)

  def info(self, free_reflections_per_bin = None, max_number_of_bins = None):
    if(free_reflections_per_bin is None):
      free_reflections_per_bin= self.alpha_beta_params.free_reflections_per_bin
    if(max_number_of_bins is None):
      max_number_of_bins = self.max_number_of_bins
    return mmtbx.f_model.info(
      fmodel                   = self,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins       = max_number_of_bins)

  def fft_vs_direct(self, reflections_per_bin = 250,
                          n_bins              = 0,
                          out                 = None):
    if(out is None): out = sys.stdout
    f_obs_w = self.f_obs_w
    f_obs_w.setup_binner(reflections_per_bin = reflections_per_bin,
                         n_bins              = n_bins)
    fmodel_dc = self.deep_copy()
    fft_p = sf_and_grads_accuracy_params.extract()
    fft_p.algorithm = "fft"
    dir_p = sf_and_grads_accuracy_params.extract()
    dir_p.algorithm = "direct"
    if(self.sfg_params.algorithm == "fft"):
       fmodel_dc.update(sf_and_grads_accuracy_params = dir_p)
    elif(self.sfg_params.algorithm == "direct"):
       fmodel_dc.update(sf_and_grads_accuracy_params = fft_p)
    print >> out, "|"+"-"*77+"|"
    print >> out, "| Bin     Resolution   Compl.  No.       Scale_k1             R-work          |"
    print >> out, "|number     range              Refl.      direct       fft direct fft-direct,%|"
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
           print >> out, format % (i_bin,d_range,compl,n_ref,
                                            scale_k1_2,r_work_1,r_work_2,delta)
        else:
           print >> out, format % (i_bin,d_range,compl,n_ref,
                                            scale_k1_1,r_work_2,r_work_1,delta)
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def explain_members(self, out=None, prefix="", suffix=""):
    if (out is None): out = sys.stdout
    def zero_if_almost_zero(v, eps=1.e-6):
      if (abs(v) < eps): return 0
      return v
    for line in [
          "Fmodel   = scale_k1 * fb_cart * (Fcalc + Fbulk)",
          "Fcalc    = structure factors calculated from atomic model",
          "Fbulk    = k_sol * exp(-b_sol*s**2/4) * Fmask",
          "scale_k1 = SUM(Fobs * Fmodel) / SUM(Fmodel**2)",
          "fb_cart  = exp(-h^t * A^-1 * B_cart * A^-t * h)",
          "A        = orthogonalization matrix",
          "k_sol    = %.6g" % self.k_sol(),
          "b_sol    = %.6g" % zero_if_almost_zero(self.b_sol()),
          "scale_k1 = %.6g" % self.scale_k1(),
          "B_cart = (B11, B22, B33, B12, B13, B23)",
          "       = (%s)" % ", ".join(
            ["%.6g" % zero_if_almost_zero(v) for v in self.b_cart()])]:
      print >> out, prefix + line + suffix

  def export(self, out=None, format="mtz"):
    assert format in ["mtz", "cns"]
    file_name = None
    if (out is None):
      out = sys.stdout
    elif (hasattr(out, "name")):
      file_name = libtbx.path.canonical_path(file_name=out.name)
    warning = [
      "DO NOT USE THIS FILE AS INPUT FOR REFINEMENT!",
      "Resolution and sigma cutoffs may have been applied to FOBS."]
    width = max([len(line) for line in warning])
    warning.insert(0, "*" * width)
    warning.append(warning[0])
    if (format == "cns"):
      for line in warning:
        print >> out, "{ %s%s }" % (line, " "*(width-len(line)))
      print >> out, "{ %s }" % date_and_time()
      if (file_name is not None):
        print >> out, "{ file name: %s }" % os.path.basename(file_name)
        print >> out, "{ directory: %s }" % os.path.dirname(file_name)
      self.explain_members(out=out, prefix="{ ", suffix=" }")
      crystal_symmetry_as_cns_comments(
        crystal_symmetry=self.f_obs, out=out)
      print >> out, "NREFlection=%d" % self.f_obs.indices().size()
      print >> out, "ANOMalous=%s" % {0: "FALSE"}.get(
        int(self.f_obs.anomalous_flag()), "TRUE")
      have_sigmas = self.f_obs.sigmas() is not None
      for n_t in [("FOBS", "REAL"),
                  ("SIGFOBS", "REAL"),
                  ("R_FREE_FLAGS", "INTEGER"),
                  ("FMODEL", "COMPLEX"),
                  ("FCALC", "COMPLEX"),
                  ("FMASK", "COMPLEX"),
                  ("FBULK", "COMPLEX"),
                  ("FB_CART", "REAL"),
                  ("FOM", "REAL"),
                  ("ALPHA", "REAL"),
                  ("BETA", "REAL")]:
        if (not have_sigmas and n_t[0] == "SIGFOBS"): continue
        print >> out, "DECLare NAME=%s DOMAin=RECIprocal TYPE=%s END" % n_t
      f_model            = self.f_model_scaled_with_k1()
      f_model_amplitudes = f_model.amplitudes().data()
      f_model_phases     = f_model.phases(deg=True).data()
      f_calc_amplitudes  = self.f_calc().amplitudes().data()
      f_calc_phases      = self.f_calc().phases(deg=True).data()
      f_mask_amplitudes  = self.f_mask().amplitudes().data()
      f_mask_phases      = self.f_mask().phases(deg=True).data()
      f_bulk_amplitudes  = self.f_bulk().amplitudes().data()
      f_bulk_phases      = self.f_bulk().phases(deg=True).data()
      alpha, beta        = [item.data() for item in self.alpha_beta()]
      arrays = [
        self.f_obs.indices(), self.f_obs.data(), self.f_obs.sigmas(),
        self.r_free_flags.data(),
        f_model_amplitudes, f_model_phases,
        f_calc_amplitudes, f_calc_phases,
        f_mask_amplitudes, f_mask_phases,
        f_bulk_amplitudes, f_bulk_phases,
        self.fb_cart(),
        self.figures_of_merit(),
        alpha, beta]
      if (not have_sigmas):
        del arrays[2]
        i_r_free_flags = 2
      else:
        i_r_free_flags = 3
      for values in zip(*arrays):
        print >> out, "INDE %d %d %d" % values[0]
        print >> out, " FOBS= %.6g" % values[1],
        if (have_sigmas):
          print >> out, " SIGFOBS= %.6g" % values[2],
        print >> out, \
          " R_FREE_FLAGS= %d FMODEL= %.6g %.6g\n" \
          " FCALC= %.6g %.6g FMASK= %.6g %.6g FBULK= %.6g %.6g\n" \
          " FB_CART= %.6g FOM= %.6g ALPHA= %.6g BETA= %.6g"  \
            % values[i_r_free_flags:]
      if (file_name is not None):
        out.close()
    else:
      assert file_name is not None
      mtz_dataset = self.f_obs.as_mtz_dataset(column_root_label="FOBS")
      mtz_dataset.add_miller_array(
        miller_array=self.r_free_flags, column_root_label="R_FREE_FLAGS")
      mtz_dataset.add_miller_array(
        miller_array=self.f_model_scaled_with_k1(), column_root_label="FMODEL")
      mtz_dataset.add_miller_array(
        miller_array=self.f_calc(), column_root_label="FCALC")
      mtz_dataset.add_miller_array(
        miller_array=self.f_mask(), column_root_label="FMASK")
      mtz_dataset.add_miller_array(
        miller_array=self.f_bulk(), column_root_label="FBULK")
      mtz_dataset.add_miller_array(
        miller_array=self.f_obs.array(data=self.fb_cart()),
        column_root_label="FB_CART", column_types="W")
      mtz_dataset.add_miller_array(
        miller_array=self.f_obs.array(data=self.figures_of_merit()),
        column_root_label="FOM", column_types="W")
      alpha, beta = self.alpha_beta()
      mtz_dataset.add_miller_array(
        miller_array=alpha, column_root_label="ALPHA", column_types="W")
      mtz_dataset.add_miller_array(
        miller_array=beta, column_root_label="BETA", column_types="W")
      mtz_history_buffer = flex.std_string(warning)
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

class phaser_sad_target_functor(object):

  def __init__(self, f_obs, r_free_flags, xray_structure, f_calc):
    adopt_init_args(self, locals())
    self.refine_sad_object = phaser.phenix_adaptors.sad_target.data_adaptor(
      f_obs=f_obs,
      r_free_flags=r_free_flags,
      verbose=True).target(xray_structure=xray_structure)
    self.refine_sad_object.set_f_calc(f_calc=f_calc)
    if (not f_obs.space_group().is_centric()):
      self.refine_sad_object.refine_variance_terms()

  def __call__(self, f_calc, compute_gradients):
    rso = self.refine_sad_object
    self.refine_sad_object.set_f_calc(f_calc=f_calc)
    target_work = rso.functional(use_working_set=True)
    gradients_work = rso.gradients()
    target_test = rso.functional(use_working_set=False)
    return xray.targets_common_results(
      target_per_reflection=flex.double(),
      target_work=target_work,
      target_test=target_test,
      gradients_work=gradients_work.data())

class target_functor(object):

  def __init__(self, manager):
    self.manager = manager
    target_name = manager.target_name
    assert target_name is not None
    attr = manager.target_attributes()
    if (target_name == "ml_sad"):
      if (phaser is None):
        raise Sorry(
          "ml_sad target requires phaser extension, which is not available"
          " in this installation.")
      self.core = phaser_sad_target_functor(
        f_obs=manager.f_obs_scaled_with_k2(),
        r_free_flags=manager.r_free_flags,
        xray_structure=manager.xray_structure,
        f_calc=manager.f_model())
    elif (attr.family == "ml"):
      if (attr.requires_experimental_phases()):
        experimental_phases = manager.abcd
      else:
        experimental_phases = None
      self.core = xray.target_functors.max_like(
        f_obs=manager.f_obs,
        r_free_flags=manager.r_free_flags,
        experimental_phases=experimental_phases,
        alpha_beta=manager.alpha_beta(),
        scale_factor=manager.scale_ml_wrapper(),
        integration_step_size=5.0)
    else:
      if (attr.pseudo_ml):
        f_obs, weights = manager.f_star_w_star()
        weights = weights.data()
        if   (target_name == "lsm_k1"):
          scale_factor = 0
        elif (target_name == "lsm_k2"):
          scale_factor = 0
        elif (target_name == "lsm_k1ask3_fixed"):
          scale_factor = manager.scale_k3_w()
        elif (target_name == "lsm_k1_fixed"):
          scale_factor = manager.scale_k1_w()
        elif (target_name == "lsm_kunit"):
          scale_factor = 1.0
        else:
          raise RuntimeError
      else:
        f_obs = manager.f_obs
        if (target_name.startswith("ls_wunit_")):
          weights = flex.double(f_obs.data().size(), 1.0)
          if   (target_name == "ls_wunit_k1"):
            scale_factor = 0
          elif (target_name == "ls_wunit_k1_fixed"):
            scale_factor = manager.scale_k1_w()
          elif (target_name == "ls_wunit_k2"):
            scale_factor = 0
          elif (target_name == "ls_wunit_kunit"):
            scale_factor = 1.0
          elif (target_name == "ls_wunit_k1ask3_fixed"):
            scale_factor = manager.scale_k3_w()
          else:
            raise RuntimeError
        elif (target_name.startswith("ls_wexp_")):
          weights = ls_sigma_weights(f_obs)
          if   (target_name == "ls_wexp_k1"):
            scale_factor = 0
          elif (target_name == "ls_wexp_k2"):
            scale_factor = 0
          elif (target_name == "ls_wexp_kunit"):
            scale_factor = 1.0
          else:
            raise RuntimeError
        elif (target_name.startswith("ls_wff_")):
          weights = ls_ff_weights(f_obs, "N", 25.0)
          if   (target_name == "ls_wff_k1"):
            scale_factor = 0
          elif (target_name == "ls_wff_k1_fixed"):
            scale_factor = manager.scale_k1_w()
          elif (target_name == "ls_wff_k1ask3_fixed"):
            scale_factor = manager.scale_k3_w()
          elif (target_name == "ls_wff_k2"):
            scale_factor = 0
          elif (target_name == "ls_wff_kunit"):
            scale_factor = 1.0
          else:
            raise RuntimeError
        else:
          raise RuntimeError
      self.core = xray.target_functors.least_squares(
        apply_scale_to_f_calc=attr.ls_apply_scale_to_f_calc(),
        compute_scale_using_all_data=False,
        f_obs=f_obs,
        r_free_flags=manager.r_free_flags,
        weights=weights,
        scale_factor=scale_factor)

  def __call__(self, compute_gradients=False):
    return target_result(
      manager=self.manager,
      core_result=self.core(
        f_calc=self.manager.f_model(),
        compute_gradients=compute_gradients))

class target_result_mixin(object):

  def gradients_wrt_atomic_parameters(self,
        selection=None,
        site=False,
        u_iso=False,
        u_aniso=False,
        occupancy=False,
        tan_b_iso_max=None,
        u_iso_refinable_params=None):
    if (tan_b_iso_max is not None and tan_b_iso_max != 0):
      raise RuntimeError("Not implemented:\n"
        "  See CVS revision 1.87, 2007/03/03 01:53:05\n"
        "  method: manager.gradient_wrt_atomic_parameters()")
    global time_gradients_wrt_atomic_parameters
    timer = user_plus_sys_time()
    manager = self.manager
    xray_structure = manager.xray_structure
    if (selection is not None):
      xray_structure = xray_structure.select(selection)
    d_target_d_f_calc = self.d_target_d_f_calc_work()
    result = None
    if (u_aniso):
      result = manager.structure_factor_gradients_w(
        u_iso_refinable_params=None,
        d_target_d_f_calc=d_target_d_f_calc.data(),
        xray_structure=xray_structure,
        n_parameters=0,
        miller_set=d_target_d_f_calc,
        algorithm=manager.sfg_params.algorithm).d_target_d_u_cart()
    elif(u_iso):
      result = manager.structure_factor_gradients_w(
        u_iso_refinable_params=None,
        d_target_d_f_calc=d_target_d_f_calc.data(),
        xray_structure=xray_structure,
        n_parameters=0,
        miller_set=d_target_d_f_calc,
        algorithm=manager.sfg_params.algorithm).d_target_d_u_iso()
    elif(occupancy):
      result = manager.structure_factor_gradients_w(
        u_iso_refinable_params=None,
        d_target_d_f_calc=d_target_d_f_calc.data(),
        xray_structure=xray_structure,
        n_parameters=0,
        miller_set=d_target_d_f_calc,
        algorithm=manager.sfg_params.algorithm).d_target_d_occupancy()
    else:
      result = manager.structure_factor_gradients_w(
        u_iso_refinable_params=u_iso_refinable_params,
        d_target_d_f_calc=d_target_d_f_calc.data(),
        xray_structure=xray_structure,
        n_parameters=xray_structure.n_parameters_XXX(),
        miller_set=d_target_d_f_calc,
        algorithm=manager.sfg_params.algorithm)
    time_gradients_wrt_atomic_parameters += timer.elapsed()
    return result

  def d_target_d_site_cart(self):
    manager = self.manager
    xray.set_scatterer_grad_flags(
      scatterers=manager.xray_structure.scatterers(),
      site=True)
    return flex.vec3_double(
      self.gradients_wrt_atomic_parameters().packed())

class target_result(target_result_mixin):

  def __init__(self, manager, core_result):
    self.manager = manager
    self.core_result = core_result

  def target_per_reflection(self):
    return self.core_result.target_per_reflection()

  def target_work(self):
    return self.core_result.target_work()

  def target_test(self):
    return self.core_result.target_test()

  def d_target_d_f_model_work(self):
    return self.manager.f_obs_w.array(
      data=self.core_result.gradients_work())

  def d_target_d_f_calc_work(self):
    return self.manager.f_obs_w.array(
      data=self.core_result.gradients_work() * self.manager.core.fb_cart_w)


def ls_ff_weights(f_obs, atom, B):
  d_star_sq_data = f_obs.d_star_sq().data()
  table = wk1995(atom).fetch()
  ff = table.at_d_star_sq(d_star_sq_data) * flex.exp(-B/4.0*d_star_sq_data)
  weights = 1.0/flex.pow2(ff)
  return weights

def ls_sigma_weights(f_obs):
  if(f_obs.sigmas() is not None):
     sigmas_squared = flex.pow2(f_obs.sigmas())
  else:
     sigmas_squared = flex.double(f_obs.data().size(), 1.0)
  assert sigmas_squared.all_gt(0)
  weights = 1 / sigmas_squared
  return weights

def kb_range(x_max, x_min, step):
  x_range = []
  x = x_min
  while x <= x_max + 0.0001:
    x_range.append(x)
    x += step
  return x_range

def n_as_s(format, value):
  if (value is None):
    return format_value(format=format, value=value)
  if (isinstance(value, (int, float))):
    return (format % value).strip()
  return [(format % v).strip() for v in value]


### Information summary holding class with formatted printing functionality ###

class resolution_bin(object):
  def __init__(self,
               i_bin         = None,
               d_range       = None,
               completeness  = None,
               alpha_work    = None,
               beta_work     = None,
               r_work        = None,
               r_free        = None,
               target_work   = None,
               target_free   = None,
               n_work        = None,
               n_free        = None,
               mean_f_obs    = None,
               fom_work      = None,
               scale_k1_work = None,
               pher_work     = None,
               pher_free     = None):
    adopt_init_args(self, locals())

class info(object):
  def __init__(self,
               fmodel,
               free_reflections_per_bin = 140,
               max_number_of_bins = 30):
    mp = fmodel.mask_params
    self.target_name = fmodel.target_name
    if(self.target_name == "twin_lsq_f"):
      self.twin_fraction = fmodel.twin_fraction_object.twin_fraction
      self.twin_law = fmodel.twin_law
    else:
      self.twin_fraction = None
      self.twin_law = None
    self.r_work = fmodel.r_work()
    self.r_free = fmodel.r_free()
    self.r_all = fmodel.r_all()
    self.target_work = fmodel.target_w()
    self.target_free = fmodel.target_t()
    self.overall_scale_k1 = fmodel.scale_k1()
    self.number_of_test_reflections = fmodel.f_calc_t().data().size()
    self.number_of_work_reflections = fmodel.f_calc_w().data().size()
    self.number_of_reflections = fmodel.f_obs.data().size()
    self.k_sol = fmodel.k_sol()
    self.b_sol = fmodel.b_sol()
    self.b_cart = fmodel.b_cart()
    self.b_iso = fmodel.b_iso()
    self.mask_solvent_radius = mp.solvent_radius
    self.mask_shrink_radius = mp.shrink_truncation_radius
    self.mask_grid_step_factor = mp.grid_step_factor
    self.ml_phase_error = flex.mean(fmodel.phase_errors())
    self.ml_coordinate_error = fmodel.model_error_ml()
    self.d_max, self.d_min = fmodel.f_obs.resolution_range()
    self.completeness_in_range = fmodel.f_obs.completeness(d_max = self.d_max)
    self.min_f_obs_over_sigma = fmodel.f_obs.min_f_over_sigma()
    self.sf_algorithm = fmodel.sfg_params.algorithm
    alpha_w, beta_w = fmodel.alpha_beta_w()
    self.alpha_work_min, self.alpha_work_max, self.alpha_work_mean = \
      alpha_w.data().min_max_mean().as_tuple()
    self.beta_work_min, self.beta_work_max, self.beta_work_mean = \
      beta_w.data().min_max_mean().as_tuple()
    self.fom_work_min, self.fom_work_max, self.fom_work_mean = \
      fmodel.figures_of_merit_work().min_max_mean().as_tuple()
    self.pher_work_min, self.pher_work_max, self.pher_work_mean = \
      fmodel.phase_errors_work().min_max_mean().as_tuple()
    self.pher_free_min, self.pher_free_max, self.pher_free_mean = \
      fmodel.phase_errors_test().min_max_mean().as_tuple()
    self.bins = self.statistics_in_resolution_bins(
      fmodel = fmodel,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins = max_number_of_bins)

  def statistics_in_resolution_bins(self, fmodel, free_reflections_per_bin,
                                    max_number_of_bins):
    if(self.target_name == "twin_lsq_f"):
      return fmodel.statistics_in_resolution_bins()
    result = []
    target_functor = fmodel.target_functor()
    target_result = target_functor(compute_gradients=False)
    tpr = target_result.target_per_reflection()
    if(tpr.size() != 0):
      tpr_w = tpr.select(fmodel.work)
      tpr_t = tpr.select(fmodel.test)
    fo_t = fmodel.f_obs_t
    fc_t = fmodel.f_model_t()
    fo_w = fmodel.f_obs_w
    fc_w = fmodel.f_model_w()
    alpha_w, beta_w = fmodel.alpha_beta_w()
    alpha_t, beta_t = fmodel.alpha_beta_t()
    pher_w = fmodel.phase_errors_work()
    pher_t = fmodel.phase_errors_test()
    fom = fmodel.figures_of_merit_work()
    fmodel.f_obs.setup_binner(n_bins=fmodel.determine_n_bins(
      free_reflections_per_bin=free_reflections_per_bin,
      max_n_bins=max_number_of_bins))
    fo_t.use_binning_of(fmodel.f_obs)
    fc_t.use_binning_of(fo_t)
    fo_w.use_binning_of(fo_t)
    fc_w.use_binning_of(fo_t)
    alpha_w.use_binning_of(fo_t)
    alpha_t.use_binning_of(fo_t)
    beta_w.use_binning_of(fo_t)
    beta_t.use_binning_of(fo_t)
    for i_bin in fo_t.binner().range_used():
      sel_t = fo_t.binner().selection(i_bin)
      sel_w = fo_w.binner().selection(i_bin)
      sel_all = fmodel.f_obs.binner().selection(i_bin)
      sel_fo_all = fmodel.f_obs.select(sel_all)
      sel_fo_t = fo_t.select(sel_t)
      sel_fc_t = fc_t.select(sel_t)
      sel_fo_w = fo_w.select(sel_w)
      sel_fc_w = fc_w.select(sel_w)
      if (tpr.size() == 0):
        sel_tpr_w = None
        sel_tpr_t = None
      else:
        denom_w = sel_fo_w.data().size()
        denom_t = sel_fo_t.data().size()
        if(denom_w != 0):
           sel_tpr_w = flex.sum(tpr_w.select(sel_w))/denom_w
        else:
           sel_tpr_w = flex.sum(tpr_w.select(sel_w))
        if(denom_t != 0):
           sel_tpr_t = flex.sum(tpr_t.select(sel_t))/denom_t
        else:
           sel_tpr_t = flex.sum(tpr_t.select(sel_t))
      d_max_,d_min_ = sel_fo_all.d_max_min()
      completeness = fmodel.f_obs.resolution_filter(
        d_min= d_min_,d_max= d_max_).completeness(d_max = d_max_)
      d_range = fo_t.binner().bin_legend(
        i_bin=i_bin, show_bin_number=False, show_counts=False)
      s_fo_w_d = sel_fo_w.data()
      s_fc_w_d = sel_fc_w.data()
      s_fc_w_d_a = flex.abs(s_fc_w_d)
      bin = resolution_bin(
        i_bin        = i_bin,
        d_range      = d_range,
        completeness = completeness,
        alpha_work   = flex.mean_default(alpha_w.select(sel_w).data(),None),
        beta_work    = flex.mean_default(beta_w.select(sel_w).data(),None),
        r_work       = bulk_solvent.r_factor(s_fo_w_d, s_fc_w_d),
        r_free       = bulk_solvent.r_factor(sel_fo_t.data(), sel_fc_t.data()),
        target_work  = sel_tpr_w,
        target_free  = sel_tpr_t,
        n_work       = sel_fo_w.data().size(),
        n_free       = sel_fo_t.data().size(),
        scale_k1_work= _scale_helper(num=s_fo_w_d, den=s_fc_w_d_a),
        mean_f_obs   = flex.mean_default(sel_fo_all.data(),None),
        fom_work     = flex.mean_default(fom.select(sel_w),None),
        pher_work    = flex.mean_default(pher_w.select(sel_w),None),
        pher_free    = flex.mean_default(pher_t.select(sel_t),None))
      result.append(bin)
    return result

  def show_remark_3(self, out = None):
    if(out is None): out = sys.stdout
    pr = "REMARK   3  "
    print >> out,pr+"REFINEMENT TARGET : %s"%self.target_name.upper()
    print >> out,pr
    print >> out,pr+"DATA USED IN REFINEMENT."
    print >> out,pr+" RESOLUTION RANGE HIGH (ANGSTROMS) : %s"%format_value("%-8.3f", self.d_min)
    print >> out,pr+" RESOLUTION RANGE LOW  (ANGSTROMS) : %s"%format_value("%-8.3f", self.d_max)
    print >> out,pr+" MIN(FOBS/SIGMA_FOBS)              : %s"%format_value("%-6.2f", self.min_f_obs_over_sigma)
    print >> out,pr+" COMPLETENESS FOR RANGE        (%s) : %-6.2f"%\
      ("%", self.completeness_in_range*100.0)
    print >> out,pr+" NUMBER OF REFLECTIONS             : %-10d"%self.number_of_reflections
    print >> out,pr
    print >> out,pr+"FIT TO DATA USED IN REFINEMENT."
    print >> out,pr+" R VALUE     (WORKING + TEST SET) : %s"%format_value("%-6.4f",self.r_all)
    print >> out,pr+" R VALUE            (WORKING SET) : %s"%format_value("%-6.4f", self.r_work)
    print >> out,pr+" FREE R VALUE                     : %s"%format_value("%-6.4f", self.r_free)
    print >> out,pr+" FREE R VALUE TEST SET SIZE   (%s) : %-6.2f"%("%",
      float(self.number_of_test_reflections)/self.number_of_reflections*100.)
    print >> out,pr+" FREE R VALUE TEST SET COUNT      : %-10d"%self.number_of_test_reflections
    print >> out,pr
    print >> out,pr+"FIT TO DATA USED IN REFINEMENT (IN BINS)."
    print >> out,pr+" BIN  RESOLUTION RANGE  COMPL.    NWORK NFREE   RWORK  RFREE"
    fmt = " %s %s    %s %s %s  %s %s"
    for bin in self.bins:
      print >> out,pr+fmt%(
        format_value("%3d", bin.i_bin),
        format_value("%-17s", bin.d_range),
        format_value("%4.2f", bin.completeness),
        format_value("%8d", bin.n_work),
        format_value("%5d", bin.n_free),
        format_value("%6.4f", bin.r_work),
        format_value("%6.4f", bin.r_free))
    print >> out,pr
    print >> out,pr+"BULK SOLVENT MODELLING."
    print >> out,pr+" METHOD USED        : FLAT BULK SOLVENT MODEL"
    print >> out,pr+" SOLVENT RADIUS     : %s"%format_value("%-8.2f", self.mask_solvent_radius)
    print >> out,pr+" SHRINKAGE RADIUS   : %s"%format_value("%-8.2f", self.mask_shrink_radius)
    print >> out,pr+" GRID STEP FACTOR   : %s"%format_value("%-8.2f", self.mask_grid_step_factor)
    print >> out,pr+" K_SOL              : %s"%format_value("%-8.3f", self.k_sol)
    print >> out,pr+" B_SOL              : %s"%format_value("%-8.3f", self.b_sol)
    print >> out,pr
    if(self.twin_fraction is not None):
      print >> out,pr+"TWINNING INFORMATION."
      print >> out,pr+" FRACTION: %s"%format_value("%-8.3f", self.twin_fraction)
      print >> out,pr+" OPERATOR: %s"%\
        format_value("%-s", sgtbx.change_of_basis_op(self.twin_law).as_hkl())
    print >> out,pr
    print >> out,pr+"ERROR ESTIMATES."
    print >> out,pr+" COORDINATE ERROR (MAXIMUM-LIKELIHOOD BASED)     : %s"%\
      format_value("%-8.2f", self.ml_coordinate_error)
    print >> out,pr+" PHASE ERROR (DEGREES, MAXIMUM-LIKELIHOOD BASED) : %s"%\
      format_value("%-8.2f", self.ml_phase_error)
    print >> out,pr
    print >> out,pr+"OVERALL SCALE FACTORS."
    print >> out,pr+" SCALE = SUM(|F_OBS|*|F_MODEL|)/SUM(|F_MODEL|**2) : %s"%\
      format_value("%-12.4f", self.overall_scale_k1)
    print >> out,pr+" ANISOTROPIC SCALE MATRIX ELEMENTS (IN CARTESIAN BASIS)."
    print >> out,pr+"  B11 : %s"%format_value("%-15.4f", self.b_cart[0])
    print >> out,pr+"  B22 : %s"%format_value("%-15.4f", self.b_cart[1])
    print >> out,pr+"  B33 : %s"%format_value("%-15.4f", self.b_cart[2])
    print >> out,pr+"  B12 : %s"%format_value("%-15.4f", self.b_cart[3])
    print >> out,pr+"  B13 : %s"%format_value("%-15.4f", self.b_cart[4])
    print >> out,pr+"  B23 : %s"%format_value("%-15.4f", self.b_cart[5])
    print >> out,pr
    print >> out,pr+"R FACTOR FORMULA."
    print >> out,pr+" R = SUM(||F_OBS|-SCALE*|F_MODEL||)/SUM(|F_OBS|)"
    print >> out,pr
    print >> out,pr+"TOTAL MODEL STRUCTURE FACTOR (F_MODEL)."
    print >> out,pr+" F_MODEL = FB_CART * (F_CALC_ATOMS + F_BULK)"
    print >> out,pr+"  F_BULK = K_SOL * EXP(-B_SOL * S**2 / 4) * F_MASK"
    print >> out,pr+"  F_CALC_ATOMS = ATOMIC MODEL STRUCTURE FACTORS"
    print >> out,pr+"  FB_CART = EXP(-H(t) * A(-1) * B * A(-1t) * H)"
    print >> out,pr+"   A = orthogonalization matrix, H = MILLER INDEX"
    print >> out,pr+"   (t) = TRANSPOSE, (-1) = INVERSE"
    print >> out,pr
    print >> out,pr+"STRUCTURE FACTORS CALCULATION ALGORITHM : %-s"%\
      self.sf_algorithm.upper()
    out.flush()

  def show_targets(self, out = None, text = ""):
    if(out is None): out = sys.stdout
    part1 = "|-"+text
    part2 = "-|"
    n = 79 - len(part1+part2)
    print >> out, part1 + "-"*n + part2
    part3 = "| target_work(%s) = %s  r_work = %s  r_free = %s" % (
      self.target_name,
      format_value(format="%.6g",  value = self.target_work),
      format_value(format="%6.4f", value = self.r_work),
      format_value(format="%6.4f", value = self.r_free))
    n = 78 - len(str(part3)+"|")
    print >> out, part3, " "*n +"|"
    print >> out, "|" +"-"*77+"|"
    out.flush()

  def _header_resolutions_nreflections(self, header, out):
    if(header is None): header = ""
    line1 = "(resolution: "
    line2 = format_value("%6.2f",self.d_min).strip()
    line3 = format_value("%6.2f",self.d_max).strip()
    line4 = " - "
    line5 = " A; n_refl. = "
    line6 = format_value("%d",self.number_of_reflections).strip()
    tl = header+"-"+line1+line2+line4+line3+line5+line6+")"
    line_len = len("|-"+"|"+tl)
    fill_len = 80-line_len-1
    print >> out, "|-"+tl+"-"*(fill_len)+"|"
    out.flush()

  def _rfactors_and_bulk_solvent_and_scale_params(self, out):
    out.flush()
    r_work = format_value("%6.4f",self.r_work).strip()
    r_free = format_value("%6.4f",self.r_free).strip()
    scale  = format_value("%6.3f",self.overall_scale_k1).strip()
    k_sol  = format_value("%4.2f",self.k_sol).strip()
    b_sol  = format_value("%6.2f",self.b_sol).strip()
    b0,b1,b2,b3,b4,b5 = n_as_s("%7.2f",self.b_cart)
    b_iso  = format_value("%7.2f",self.b_iso).strip()
    line = "| r_work= "+r_work+"   r_free= "+r_free+"   ksol= "+k_sol+\
           "   Bsol= "+b_sol+"   scale= "+scale
    np = 79 - (len(line) + 1)
    if(np < 0): np = 0
    print >> out, line + " "*np + "|"
    print >> out, "| "+"  "*38+"|"
    print >> out, "| overall anisotropic scale matrix (Cartesian basis; B11,B22,B33,B12,B13,B23):|"
    c = ","
    line4 = "| ("+b0+c+b1+c+b2+c+b3+c+b4+c+b5+"); trace/3= "+b_iso
    np = 79 - (len(line4) + 1)
    line4 = line4 + " "*np + "|"
    print >> out, line4
    out.flush()

  def show_rfactors_targets_scales_overall(self, header = None, out=None):
    if(out is None): out = sys.stdout
    out.flush()
    p = " "
    self._header_resolutions_nreflections(header=header, out=out)
    print >> out, "| "+"  "*38+"|"
    self._rfactors_and_bulk_solvent_and_scale_params(out=out)
    err = format_value("%6.2f",self.ml_coordinate_error)
    print >> out, "| "+"  "*38+"|"
    line6="| maximum likelihood estimate for coordinate error: "+err+" A"
    np = 79 - (len(line6) + 1)
    line6 = line6 + " "*np + "|"
    print >> out, line6
    line7="| x-ray target function (%s) for work reflections: %s"% (
      self.target_name, n_as_s("%15.6f",self.target_work))
    np = 79 - (len(line7) + 1)
    line7 = line7 + " "*np + "|"
    print >> out, line7
    if(self.twin_fraction is not None):
      line8="| twin fraction: "+format_value("%-4.2f",self.twin_fraction)+\
        "  twin operator: "+\
        format_value("%-s",sgtbx.change_of_basis_op(self.twin_law).as_hkl())
      np = 79 - (len(line8) + 1)
      line8 = line8 + " "*np + "|"
      print >> out, line8
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_rfactors_targets_in_bins(self, out = None):
    if(out is None): out = sys.stdout
    print >> out, "|"+"-"*77+"|"
    print >> out, "| Bin     Resolution   Compl.  No. Refl.    R-factors          Targets        |"
    print >> out, "|number     range              work test   work   test        work        test|"
    for bin in self.bins:
      print >> out, "|%3d: %-17s %4.2f %6d %4d %s %s %s %s|" % (
        bin.i_bin,
        bin.d_range,
        bin.completeness,
        bin.n_work,
        bin.n_free,
        format_value("%6.4f",  bin.r_work),
        format_value("%6.4f",  bin.r_free),
        format_value("%11.5g", bin.target_work),
        format_value("%11.5g", bin.target_free))
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_fom_pher_alpha_beta_in_bins(self, out = None):
    if(out is None): out = sys.stdout
    print >> out, "|"+"-"*77+"|"
    print >> out, "|R-free likelihood based estimates for figures of merit, absolute phase error,|"
    print >> out, "|and distribution parameters alpha and beta (Acta Cryst. (1995). A51, 880-887)|"
    print >> out, "|"+" "*77+"|"
    print >> out, "| Bin     Resolution      No. Refl.   FOM  Phase Scale    Alpha        Beta   |"
    print >> out, "|  #        range        work  test        error factor                       |"
    for bin in self.bins:
      print >> out, "|%3d: %-17s%6d%6d%s%s%s%s%s|" % (
        bin.i_bin,
        bin.d_range,
        bin.n_work,
        bin.n_free,
        format_value("%6.2f",  bin.fom_work),
        format_value("%7.2f",  bin.pher_work),
        format_value("%7.2f",  bin.scale_k1_work),
        format_value("%9.2f",  bin.alpha_work),
        format_value("%14.2f", bin.beta_work))
    print >>out, "|alpha:            min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.alpha_work_min),
      format_value("%16.2f", self.alpha_work_max),
      format_value("%13.2f", self.alpha_work_mean))
    print >>out, "|beta:             min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.beta_work_min),
      format_value("%16.2f", self.beta_work_max),
      format_value("%13.2f", self.beta_work_mean))
    print >>out, "|figures of merit: min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.fom_work_min),
      format_value("%16.2f", self.fom_work_max),
      format_value("%13.2f", self.fom_work_mean))
    print >>out, "|phase err.(work): min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.pher_work_min),
      format_value("%16.2f", self.pher_work_max),
      format_value("%13.2f", self.pher_work_mean))
    print >>out, "|phase err.(test): min =%s max =%s mean =%s|"%(
      format_value("%12.2f", self.pher_free_min),
      format_value("%16.2f", self.pher_free_max),
      format_value("%13.2f", self.pher_free_mean))
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_all(self, header = "", out = None):
    if(out is None): out = sys.stdout
    self.show_rfactors_targets_scales_overall(header = header, out = out)
    print >> out
    self.show_rfactors_targets_in_bins(out = out)
    print >> out
    self.show_fom_pher_alpha_beta_in_bins(out = out)

def show_histogram(data, n_slots, log):
  hm = flex.histogram(data = data, n_slots = n_slots)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    print >> log, "%10.3f - %-10.3f : %d" % (lc_1, hc_1, n_1)
    lc_1 = hc_1
