from __future__ import division
try: import phaser
except ImportError: phaser = None
else: import phaser.phenix_adaptors.sad_target
from cctbx.array_family import flex
import math, time
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
import sys, random
from cctbx import miller
import cctbx.xray.structure_factors
from cctbx.array_family import flex
from stdlib import math
from cctbx import xray
from cctbx import adptbx
import boost.python
import mmtbx
from libtbx.math_utils import iround
from libtbx.utils import Sorry, user_plus_sys_time
from libtbx.str_utils import show_string

ext = boost.python.import_ext("mmtbx_f_model_ext")

time_bulk_solvent_and_scale         = 0.0
time_mask                           = 0.0
number_mask                         = 0
time_f_calc                         = 0.0
time_alpha_beta                     = 0.0
time_target                         = 0.0
time_gradient_wrt_atomic_parameters = 0.0
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
          time_gradient_wrt_atomic_parameters +\
          time_fmodel_core_data               +\
          time_r_factors                      +\
          time_phase_errors                   +\
          time_foms
  print >> out, "  Micro-tasks:"
  print >> out, "    mask                           = %-7.2f number of calls = %3d " % \
                (time_mask,number_mask)
  print >> out, "    f_calc                         = %-7.2f" % time_f_calc
  print >> out, "    alpha_beta                     = %-7.2f" % time_alpha_beta
  print >> out, "    target                         = %-7.2f" % time_target
  print >> out, "    gradient_wrt_atomic_parameters = %-7.2f" % \
                                            time_gradient_wrt_atomic_parameters
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

class target_attributes(xray.target_functors.target_attributes):

  def validate(self):
    if (self.family == "lsm"):
      self.family = "ls"
      self.pseudo_ml = True
    else:
      self.pseudo_ml = False
    if (xray.target_functors.target_attributes.validate(self)):
      return True
    if (self.family == "ml" and self.specialization == "sad"):
      return True
    return False

class manager(object):

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

  def __init__(self, f_obs                 = None,
                     r_free_flags          = None,
                     b_cart                = [0.,0.,0.,0.,0.,0.],
                     k_sol                 = 0.0,
                     b_sol                 = 0.0,
                     sf_algorithm          = "fft",
                     sf_cos_sin_table      = True,
                     target_name           = None,
                     abcd                  = None,
                     alpha_beta_params     = None,
                     xray_structure        = None,
                     f_mask                = None,
                     f_calc                = None,
                     mask_params           = None,
                     trust_xray_structure  = False,
                     update_xray_structure = True,
                     use_f_model_scaled    = False,
                     overall_scale         = 1.0,
                     f_ordered_solvent     = None,
                     f_ordered_solvent_dist = None,
                     n_ordered_water = 0,
                     b_ordered_water = 0.0):
    assert f_obs is not None
    assert f_obs.is_real_array()
    self.f_obs             = f_obs
    self.r_free_flags      = None
    self.sf_algorithm      = sf_algorithm
    self.sf_cos_sin_table  = sf_cos_sin_table
    self.abcd              = abcd
    self.alpha_beta_params = alpha_beta_params
    self.xray_structure    = xray_structure
    self.use_f_model_scaled= use_f_model_scaled
    self.overall_scale     = overall_scale
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
    self.structure_factor_gradients_w = cctbx.xray.structure_factors.gradients(
                                         miller_set    = self.f_obs_w,
                                         cos_sin_table = self.sf_cos_sin_table)
    self.uc = self.f_obs.unit_cell()
    self.d_spacings = self.f_obs.d_spacings().data()
    self.d_spacings_w = self.d_spacings.select(self.work)
    self.d_spacings_t = self.d_spacings.select(self.test)
    self.ss = 1./flex.pow2(self.d_spacings) / 4.
    if(self.xray_structure is None):
       self.xray_structure_mask_cache = None
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
       self.xray_structure_mask_cache = \
                                     self.xray_structure.deep_copy_scatterers()
       if(not trust_xray_structure):
          assert [f_calc, f_mask].count(None) == 2
       if(update_xray_structure):
          self.update_xray_structure(xray_structure       = self.xray_structure,
                                     update_f_calc        = True,
                                     update_f_mask        = True,
                                     force_update_f_mask  = True,
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
    if(self.sf_algorithm not in ("fft", "direct")):
       raise RuntimeError("Unknown s.f. calculation method: %s"%
                                                             self.sf_algorithm)
    self._setup_target_functors(target_name=target_name)

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

  def deep_copy(self):
    if(self.abcd is not None):
       abcd = self.abcd.deep_copy()
    else:
       abcd = None
    if(self.xray_structure is None):
       xrs = None
    else:
       xrs = self.xray_structure.deep_copy_scatterers()
    return manager(f_obs              = self.f_obs.deep_copy(),
                r_free_flags          = self.r_free_flags.deep_copy(),
                b_cart                = self.b_cart(),
                f_mask                = self.f_mask().deep_copy(),
                k_sol                 = self.k_sol(),
                b_sol                 = self.b_sol(),
                sf_algorithm          = self.sf_algorithm,
                sf_cos_sin_table      = self.sf_cos_sin_table,
                target_name           = self.target_name,
                abcd                  = abcd,
                alpha_beta_params     = self.alpha_beta_params,
                xray_structure        = xrs,
                f_calc                = self.f_calc().deep_copy(),
                mask_params           = self.mask_params,
                trust_xray_structure  = True,
                update_xray_structure = False,
      overall_scale         = self.overall_scale,
      f_ordered_solvent      = self.f_ordered_solvent.deep_copy(),
      f_ordered_solvent_dist = self.f_ordered_solvent_dist.deep_copy(),
      n_ordered_water        = self.n_ordered_water,
      b_ordered_water        = self.b_ordered_water)

  def select(self, selection, update_xray_structure=False):
    if (self.abcd  is not None):
      abcd = self.abcd.select(selection=selection)
    else:
      abcd = None
    return manager(
      f_obs=self.f_obs.select(selection=selection),
      r_free_flags=self.r_free_flags.select(selection=selection),
      b_cart=tuple(self.b_cart()),
      k_sol=self.k_sol(),
      b_sol=self.b_sol(),
      sf_algorithm=self.sf_algorithm,
      sf_cos_sin_table=self.sf_cos_sin_table,
      target_name=self.target_name,
      abcd=abcd,
      alpha_beta_params=self.alpha_beta_params,
      xray_structure=self.xray_structure,
      f_calc=self.f_calc().select(selection=selection),
      f_mask=self.f_mask().select(selection=selection),
      mask_params=self.mask_params,
      trust_xray_structure=True,
      update_xray_structure=update_xray_structure,
      overall_scale=self.overall_scale,
      f_ordered_solvent=self.f_ordered_solvent.select(
        selection=selection),
      f_ordered_solvent_dist=self.f_ordered_solvent_dist.select(
        selection=selection),
      n_ordered_water=self.n_ordered_water,
      b_ordered_water=self.b_ordered_water)

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
    uc = self.xray_structure.unit_cell()
    b_min = adptbx.u_as_b(flex.min(
          self.xray_structure.scatterers().u_cart_eigenvalues(uc).as_double()))
    assert b_min >= 0.0
    b_iso = self.b_iso()
    b_test = b_min+b_iso
    if(b_test < 0.0): b_adj = b_iso + abs(b_test) + 0.001
    else: b_adj = b_iso
    b_cart = self.b_cart()
    b_cart_new = [b_cart[0]-b_adj,b_cart[1]-b_adj,b_cart[2]-b_adj,
                  b_cart[3],      b_cart[4],      b_cart[5]]
    self.update(b_cart = b_cart_new)
    self.update(b_sol = self.k_sol_b_sol()[1] + b_adj)
    self.xray_structure.shift_us(b_shift = b_adj)
    b_min = adptbx.u_as_b(flex.min(
          self.xray_structure.scatterers().u_cart_eigenvalues(uc).as_double()))
    assert b_min >= 0.0
    self.xray_structure.tidy_us(u_min = 1.e-6)
    self.update_xray_structure(xray_structure           = self.xray_structure,
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

  def _get_step(self, update_f_ordered_solvent = False):
    step = self.f_obs.d_min()/self.mask_params.grid_step_factor
    if(step < 0.2): step = 0.2
    step = min(0.8, step)
    if(update_f_ordered_solvent): step = 0.3
    return step

  def _update_f_mask_flag(self, xray_structure, mean_shift):
    if(self.xray_structure_mask_cache is None):
       self.xray_structure_mask_cache = xray_structure.deep_copy_scatterers()
       return True
    else:
       sites_cart_1 = self.xray_structure_mask_cache.sites_cart()
       sites_cart_2 = xray_structure.sites_cart()
       self.xray_structure_mask_cache = xray_structure.deep_copy_scatterers()
       if(sites_cart_1.size() != sites_cart_2.size()): return True
       atom_atom_distances = flex.sqrt((sites_cart_1 - sites_cart_2).dot())
       mean_shift_ = flex.mean(atom_atom_distances)
       update_f_mask = False
       if(mean_shift_ >= mean_shift):
          update_f_mask = True
       return update_f_mask

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
    if(update_f_mask):
       if(force_update_f_mask):
          consider_mask_update = True
       else:
          consider_mask_update = self._update_f_mask_flag(
                  xray_structure = self.xray_structure,
                  mean_shift     = self.mask_params.mean_shift_for_mask_update)
    step = self._get_step(update_f_ordered_solvent = update_f_ordered_solvent)
    f_calc = None
    if(update_f_calc):
       global time_f_calc
       timer = user_plus_sys_time()
       assert self.xray_structure is not None
       f_calc = self.f_obs.structure_factors_from_scatterers(
                               xray_structure = self.xray_structure,
                               algorithm      = self.sf_algorithm,
                               cos_sin_table  = self.sf_cos_sin_table).f_calc()
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
    if(update_f_mask and consider_mask_update):
       global time_mask, number_mask
       number_mask += 1
       timer = user_plus_sys_time()
       if(update_f_ordered_solvent == False): nu_map = None
       bulk_solvent_mask_obj = self.bulk_solvent_mask()
       if (nu_map is not None):
         bulk_solvent_mask_obj.subtract_non_uniform_solvent_region_in_place(
                                                     non_uniform_mask = nu_map)
       if(out is not None and self.mask_params.verbose > 0):
          bulk_solvent_mask_obj.show_summary(out = out)
       f_mask = bulk_solvent_mask_obj.structure_factors(miller_set= self.f_obs)
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

  def bulk_solvent_mask(self):
    step = self._get_step()
    result = masks.bulk_solvent(
          xray_structure           = self.xray_structure,
          grid_step                = step,
          solvent_radius           = self.mask_params.solvent_radius,
          shrink_truncation_radius = self.mask_params.shrink_truncation_radius)
    return result

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
          rw_ = self.r_work()
          if(rw_ < rw):
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
    print >> out, "| r_work= %6.4f     r_free= %6.4f     Rad_solv= %4.2f     Rad_shrink= %4.2f   |"%\
     (self.r_work(), self.r_free(), self.mask_params.solvent_radius,
     self.mask_params.shrink_truncation_radius)
    print >> out, "|"+"-"*77+"|"
    print >> out

  def update_solvent_and_scale(self, params = None, out = None, verbose=None):
    global time_bulk_solvent_and_scale
    timer = user_plus_sys_time()
    if(out is None): out = sys.stdout
    if(params is None):
       params = bss.solvent_and_scale_params()
    else:
       params = bss.solvent_and_scale_params(params = params)
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

  def _setup_target_functors(self, target_name):
    if (target_name is None):
      self._target_name = None
      return
    if (target_name == "ls"): target_name = "ls_wunit_k1"
    self._target_name = target_name
    attr = self.target_attributes()
    if (attr.family == "ml"):
      f_obs = self.f_obs
      weights = None
      scale_factor = 0
    elif (attr.pseudo_ml):
      f_obs, weights = self.f_star_w_star()
      weights = weights.data()
      if   (target_name == "lsm_k1"):
        scale_factor = 0
      elif (target_name == "lsm_k2"):
        scale_factor = 0
      elif (target_name == "lsm_k1ask3_fixed"):
        scale_factor = self.scale_k3_w()
      elif (target_name == "lsm_k1_fixed"):
        scale_factor = self.scale_k1_w()
      elif (target_name == "lsm_kunit"):
        scale_factor = 1.0
      else:
        raise RuntimeError
    else:
      f_obs = self.f_obs
      if (target_name.startswith("ls_wunit_")):
        weights = flex.double(self.f_obs.data().size(), 1.0)
        if   (target_name == "ls_wunit_k1"):
          scale_factor = 0
        elif (target_name == "ls_wunit_k1_fixed"):
          scale_factor = self.scale_k1_w()
        elif (target_name == "ls_wunit_k2"):
          scale_factor = 0
        elif (target_name == "ls_wunit_kunit"):
          scale_factor = 1.0
        elif (target_name == "ls_wunit_k1ask3_fixed"):
          scale_factor = self.scale_k3_w()
        else:
          raise RuntimeError
      elif (target_name.startswith("ls_wexp_")):
        weights = ls_sigma_weights(self.f_obs)
        if   (target_name == "ls_wexp_k1"):
          scale_factor = 0
        elif (target_name == "ls_wexp_k2"):
          scale_factor = 0
        elif (target_name == "ls_wexp_kunit"):
          scale_factor = 1.0
        else:
          raise RuntimeError
      elif (target_name.startswith("ls_wff_")):
        weights = ls_ff_weights(self.f_obs, "N", 25.0)
        if   (target_name == "ls_wff_k1"):
          scale_factor = 0
        elif (target_name == "ls_wff_k1_fixed"):
          scale_factor = self.scale_k1_w()
        elif (target_name == "ls_wff_k1ask3_fixed"):
          scale_factor = self.scale_k3_w()
        elif (target_name == "ls_wff_k2"):
          scale_factor = 0
        elif (target_name == "ls_wff_kunit"):
          scale_factor = 1.0
        else:
          raise RuntimeError
      else:
        raise RuntimeError
    if (attr.specialization == "sad"): # XXX
      self.target_functors = None
      self.target_functor_w = None
      self.target_functor_t = None
    else:
      self.target_functors = xray.target_functors.manager(
        target_attributes=attr,
        f_obs=f_obs,
        r_free_flags=self.r_free_flags,
        experimental_phases=self.abcd,
        weights=weights,
        scale_factor=scale_factor)
      self.target_functor_w = self.target_functors.target_functor_w()
      self.target_functor_t = self.target_functors.target_functor_t()

  def xray_target_functor_result(self, compute_gradients = None,
                                       alpha             = None,
                                       beta              = None,
                                       scale_ml          = None,
                                       flag              = None):
    assert compute_gradients in (True,False)
    assert flag in ("work", "test")
    if(flag == "work"):
       f_model = self.f_model_w()
    else:
       f_model = self.f_model_t()
    if(self.target_name in ("ml","mlhl")):
       if(alpha is None and beta is None):
          if(flag == "work"):
             alpha, beta = self.alpha_beta_w()
          else:
             alpha, beta = self.alpha_beta_t()
       else:
          assert alpha.data().size() == f_model.data().size()
          assert beta.data().size()  == f_model.data().size()
       if (scale_ml is None):
         scale_ml = self.scale_ml_wrapper()
       if(flag == "work"):
          return self.target_functor_w(f_model,
                                       alpha.data(),
                                       beta.data(),
                                       scale_ml,
                                       compute_gradients)
       else:
          return self.target_functor_t(f_model,
                                       alpha.data(),
                                       beta.data(),
                                       scale_ml,
                                       compute_gradients)
    if(self.target_name.startswith("ls")):
       if(flag == "work"):
          return self.target_functor_w(f_model, compute_gradients)
       else:
          return self.target_functor_t(f_model, compute_gradients)
    raise RuntimeError

  def target_w(self, alpha=None, beta=None, scale_ml=None):
    global time_target
    timer = user_plus_sys_time()
    result = self.xray_target_functor_result(
      compute_gradients = False,
      alpha             = alpha,
      beta              = beta,
      scale_ml          = scale_ml,
      flag              = "work").target_work()
    time_target += timer.elapsed()
    return result

  def target_t(self, alpha=None, beta=None, scale_ml=None):
    global time_target
    timer = user_plus_sys_time()
    result = self.xray_target_functor_result(
      compute_gradients = False,
      alpha             = alpha,
      beta              = beta,
      scale_ml          = scale_ml,
      flag              = "test").target_work()
    time_target += timer.elapsed()
    return result

  def gradient_wrt_atomic_parameters(self, selection     = None,
                                           site          = False,
                                           u_iso         = False,
                                           u_aniso       = False,
                                           occupancy     = False,
                                           alpha         = None,
                                           beta          = None,
                                           tan_b_iso_max = None,
                                           u_iso_refinable_params = None):
    global time_gradient_wrt_atomic_parameters
    timer = user_plus_sys_time()
    xrs = self.xray_structure
    if([site, u_iso, u_aniso, occupancy].count(True) > 0):
       tan_u_iso = False
       param = 0
       if(u_iso):
          assert tan_b_iso_max is not None
          if(tan_b_iso_max > 0.0):
             tan_u_iso = True
             param = int(tan_b_iso_max)
       assert [site, u_iso, u_aniso].count(None) == 0
       if(selection is not None):
          xrs = self.xray_structure.deep_copy_scatterers()
       else:
          xrs = self.xray_structure
       ##XXX very inefficient code:
       #xray.set_scatterer_grad_flags(scatterers = xrs.scatterers(),
       #                              site       = site,
       #                              u_iso      = u_iso,
       #                              u_aniso    = u_aniso,
       #                              tan_u_iso  = tan_u_iso,
       #                              param      = param)
    #structure_factor_gradients = cctbx.xray.structure_factors.gradients(
    #                                     miller_set    = self.f_obs_w,
    #                                     cos_sin_table = self.sf_cos_sin_table)
    #XXX clear with target names
    attr = self.target_attributes()
    if(attr.family == "ml" or attr.pseudo_ml):
       if([alpha, beta].count(None) == 2):
          alpha, beta = self.alpha_beta_w()
    else:
       assert [alpha, beta].count(None) == 2
    if(selection is not None):
       xrs = xrs.select(selection)
    #XXX make it general
    xrtfr = self.xray_target_functor_result(compute_gradients = True,
                                            alpha             = alpha,
                                            beta              = beta,
                                            scale_ml          = None,
                                            flag              = "work")
    if(u_iso and u_iso_refinable_params is None):
       # XXX here is not clean too
       if(tan_b_iso_max != 0):
          u_iso_max = adptbx.b_as_u(tan_b_iso_max)
          u_iso_refinable_params = flex.tan(math.pi*
            (self.xray_structure.scatterers().extract_u_iso()/u_iso_max-1./2.))
       if(tan_b_iso_max == 0):
          u_iso_refinable_params = None
    result = None
    if(u_aniso):
       result = self.structure_factor_gradients_w(
         u_iso_refinable_params = None,
         d_target_d_f_calc  = xrtfr.gradients_work() * self.core.fb_cart_w,
         xray_structure     = xrs,
         n_parameters       = 0,
         miller_set         = self.f_obs_w,
         algorithm          = self.sf_algorithm).d_target_d_u_cart()
    else:
       result = self.structure_factor_gradients_w(
         u_iso_refinable_params = u_iso_refinable_params,
         d_target_d_f_calc  = xrtfr.gradients_work() * self.core.fb_cart_w,
         xray_structure     = xrs,
         n_parameters       = xrs.n_parameters_XXX(),
         miller_set         = self.f_obs_w,
         algorithm          = self.sf_algorithm)
    time_gradient_wrt_atomic_parameters += timer.elapsed()
    return result

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

  def update(self, f_calc              = None,
                   f_obs               = None,
                   f_mask              = None,
                   f_ordered_solvent   = None,
                   r_free_flags        = None,
                   b_cart              = None,
                   k_sol               = None,
                   b_sol               = None,
                   sf_algorithm        = None,
                   target_name         = None,
                   abcd                = None,
                   alpha_beta_params   = None,
                   xray_structure      = None,
                   mask_params         = None,
                   overall_scale       = None):
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
    if(sf_algorithm is not None):
      assert sf_algorithm in ("fft", "direct")
      self.sf_algorithm = sf_algorithm
      self.update_xray_structure(update_f_calc  = True)
    if(target_name is not None):
      self._setup_target_functors(target_name=target_name)
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

  def f_model(self):
    return self.core.f_model

  def f_model_scaled_with_k1(self):
    return miller.array(miller_set = self.f_obs,
                        data       = self.scale_k1()*self.f_model().data())

  def f_model_scaled_with_k1_t(self):
    return miller.array(miller_set = self.f_obs_t(),
                        data       = self.scale_k1_t()*self.f_model_t().data())

  def f_model_scaled_with_k1_w(self):
    return miller.array(miller_set = self.f_obs_w(),
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

  def r_work_in_lowest_resolution_bin(self, free_reflections_per_bin=140):
    fo_w = self.f_obs_w
    fc_w = self.f_model_w()
    fo_w.setup_binner(
      n_bins=self.determine_n_bins(
        free_reflections_per_bin=free_reflections_per_bin))
    fo_w.use_binning_of(fo_w)
    fc_w.use_binning_of(fo_w)
    r = []
    for i_bin in fo_w.binner().range_used():
        sel_w = fo_w.binner().selection(i_bin)
        sel_fo_w = fo_w.select(sel_w)
        sel_fc_w = fc_w.select(sel_w)
        r.append(bulk_solvent.r_factor(sel_fo_w.data(), sel_fc_w.data()))
    return r[0]

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
    for fmodel in [
          self.resolution_filter(d_max=6.0),
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

  def scale_k1(self, selection = None):
    fo = self.f_obs.data()
    fc = flex.abs(self.core.f_model.data())
    if(selection is not None):
       fo = fo.select(selection)
       fc = fc.select(selection)
    return flex.sum(fo*fc) / flex.sum(fc*fc)

  def scale_k1_w(self, selection = None):
    fo = self.f_obs_w.data()
    fc = flex.abs(self.core.f_model_w.data())
    if(selection is not None):
       fo = fo.select(selection)
       fc = fc.select(selection)
    return flex.sum(fo*fc) / flex.sum(fc*fc)

  def scale_k1_t(self):
    fo = self.f_obs_t.data()
    fc = flex.abs(self.core.f_model_t.data())
    return flex.sum(fo*fc) / flex.sum(fc*fc)

  def scale_k2(self):
    fo = self.f_obs.data()
    fc = flex.abs(self.core.f_model.data())
    return flex.sum(fo*fc) / flex.sum(fo*fo)

  def scale_k2_w(self):
    fo = self.f_obs_w.data()
    fc = flex.abs(self.core.f_model_w.data())
    return flex.sum(fo*fc) / flex.sum(fo*fo)

  def scale_k2_t(self):
    fo = self.f_obs_t.data()
    fc = flex.abs(self.core.f_model_t.data())
    return flex.sum(fo*fc) / flex.sum(fo*fo)

  def scale_k3_w(self):
    eps = self.f_obs_w.epsilons().data().as_double()
    mul = self.f_obs_w.multiplicities().data().as_double()
    fo = self.f_obs_w.data()
    fc = flex.abs(self.core.f_model_w.data())
    return math.sqrt(flex.sum(fo * fo * mul / eps) / \
                     flex.sum(fc * fc * mul / eps) )

  def scale_k3_t(self):
    eps = self.f_obs_t.epsilons().data().as_double()
    mul = self.f_obs_t.multiplicities().data().as_double()
    fo = self.f_obs_t.data()
    fc = flex.abs(self.core.f_model_t.data())
    return math.sqrt(flex.sum(fo * fo * mul / eps) / \
                     flex.sum(fc * fc * mul / eps) )

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


  def map_coefficients(self,
                       map_type          = None,
                       k                 = None,
                       n                 = None,
                       w1                = None,
                       w2                = None):
    assert map_type in ("k*Fobs-n*Fmodel",
                        "2m*Fobs-D*Fmodel",
                        "m*Fobs-D*Fmodel")
    if(map_type == "k*Fobs-n*Fmodel"):
       if([k,n].count(None) != 0):
          raise Sorry("Map coefficients (k and n) must be provided.")
    fb_cart  = self.fb_cart()
    scale_k2 = self.scale_k2()
    f_obs_scale   = 1.0 / fb_cart * scale_k2
    f_model_scale = 1.0 / fb_cart
    if(map_type == "k*Fobs-n*Fmodel"):
       return self._map_coeff(f_obs         = self.f_obs,
                              f_model       = self.f_model(),
                              f_obs_scale   = k * f_obs_scale,
                              f_model_scale = n * f_model_scale)
    if(map_type == "2m*Fobs-D*Fmodel"):
      alpha, beta = self.alpha_beta(
        f_obs   = self.f_obs.array(data = self.f_obs.data() * f_obs_scale),
        f_model = self.f_model().array(data = self.f_model().data()*f_model_scale))
      return self._map_coeff(
                        f_obs         = self.f_obs,
                        f_model       = self.f_model(),
                        f_obs_scale   = 2.*self.figures_of_merit()*f_obs_scale,
                        f_model_scale = alpha.data() * f_model_scale)
    if(map_type == "m*Fobs-D*Fmodel"):
      alpha, beta = self.alpha_beta(
        f_obs   = self.f_obs.array(data = self.f_obs.data() * f_obs_scale),
        f_model = self.f_model().array(data = self.f_model().data()*f_model_scale))
      return self._map_coeff(
                        f_obs         = self.f_obs,
                        f_model       = self.f_model(),
                        f_obs_scale   = self.figures_of_merit() * f_obs_scale,
                        f_model_scale = alpha.data() * f_model_scale)
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
      #f = open("qq","w")
      #fom = self.figures_of_merit()
      #for i, a, b in zip(self.f_calc.indices(), fom, alpha.data()):
      #    print >> f, "%5d%5d%5d %10.3f %10.3f" % (i[0], i[1], i[2], a, b)
      #return miller.array(miller_set = f_model,
      #                    data       = d_obs.data() - d_model)

  def electron_density_map(self,
                           map_type          = "k*Fobs-n*Fmodel",
                           k                 = 1,
                           n                 = 1,
                           w1                = None,
                           w2                = None,
                           resolution_factor = 1/3.,
                           symmetry_flags = None):
    assert map_type in ("k*Fobs-n*Fmodel",
                        "2m*Fobs-D*Fmodel",
                        "m*Fobs-D*Fmodel",
                        "m*w1*Fobs-n*w2*Fmodel")
    return self.map_coefficients(
                       map_type          = map_type,
                       k                 = k,
                       n                 = n,
                       w1                = w1,
                       w2                = w2).fft_map(
                                         resolution_factor = resolution_factor,
                                         symmetry_flags    = symmetry_flags)

  def show_targets(self, out=None, text=""):
    global time_show
    timer = user_plus_sys_time()
    if(out is None): out = sys.stdout
    part1 = "|-"+text
    part2 = "-|"
    n = 79 - len(part1+part2)
    print >> out, part1 + "-"*n + part2
    part3 = "| target_work(%s"%self.target_name+") = %.6e  r_work = %6.4f  r_free = %6.4f"%\
                                (self.target_w(), self.r_work(), self.r_free())
    n = 78 - len(str(part3)+"|")
    print >> out, part3, " "*n +"|"
    print >> out, "|" +"-"*77+"|"
    out.flush()
    time_show += timer.elapsed()

  def show(self, out=None):
    global time_show
    timer = user_plus_sys_time()
    if(out is None): out = sys.stdout
    print >> out, "f_calc          = ", self.f_calc()
    print >> out, "f_obs           = ", self.f_obs
    print >> out, "f_mask          = ", self.f_mask()
    print >> out, "r_free_flags    = ", self.r_free_flags
    print >> out, "b_cart          = ", self.b_cart()
    print >> out, "k_sol           = ", self.k_sol()
    print >> out, "b_sol           = ", self.b_sol()
    print >> out, "sf_algorithm    = ", self.sf_algorithm
    print >> out, "target_name     = ", self.target_name
    out.flush()
    time_show += timer.elapsed()

  def show_k_sol_b_sol_b_cart_target(self, header=None,target=None,out=None):
    global time_show
    timer = user_plus_sys_time()
    if(out is None): out = sys.stdout
    p = " "
    if(header is None): header = ""
    line_len = len("|-"+"|"+header)
    fill_len = 80-line_len-1
    print >> out, "|-"+header+"-"*(fill_len)+"|"
    k_sol = self.k_sol()
    b_sol = self.b_sol()
    u0,u1,u2,u3,u4,u5 = self.b_cart()
    if(target is None):
       target_w = self.target_w()
    else:
       target_w = target
    alpha, beta = self.alpha_beta_w()
    alpha_d = alpha.data()
    a_mean = flex.mean(alpha_d)
    a_zero = (alpha_d <= 0.0).count(True)
    r_work = self.r_work()
    u_isos = self.xray_structure.extract_u_iso_or_u_equiv()
    b_iso_mean = flex.mean(u_isos * math.pi**2*8)
    print >> out, "| k_sol=%5.2f b_sol=%7.2f target_w =%20.6f r_work=%7.4f" % \
                  (k_sol, b_sol, target_w, r_work) + 5*p+"|"
    print >> out, "| B(11,22,33,12,13,23)=%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f |" % \
                  (u0,u1,u2,u3,u4,u5)
    print >> out, "| trace(B) = (B11 + B22 + B33)/3 = %-10.3f                                 |"%self.b_iso()
    if(flex.mean(flex.abs(self.f_ordered_solvent.data())) > 1.e-6):
       print >> out, "| n_ordered_solv=%6d b_ordered_solv=%7.2f b_mean=%7.2f " \
                  "n_atoms=%7d |" % (self.n_ordered_water,\
                                 self.b_ordered_water,b_iso_mean,u_isos.size())
    print >> out, "| mean alpha:%8.4f  number of alpha <= 0.0:%7d" % \
                  (a_mean, a_zero)+25*p+"|"
    print >> out, "|"+"-"*77+"|"
    out.flush()
    time_show += timer.elapsed()

  def show_essential(self, header = None, out=None):
    global time_show
    timer = user_plus_sys_time()
    if(out is None): out = sys.stdout
    out.flush()
    p = " "
    if(header is None): header = ""
    d_max, d_min = self.f_obs.d_max_min()
    line1 = "---(resolution: "
    line2 = n_as_s("%6.2f",d_min)
    line3 = n_as_s("%6.2f",d_max)
    line4 = " - "
    line5 = " A)"
    tl = header+line1+line2+line4+line3+line5
    line_len = len("|-"+"|"+tl)
    fill_len = 80-line_len-1
    print >> out, "|-"+tl+"-"*(fill_len)+"|"
    print >> out, "| "+"  "*38+"|"
    r_work = n_as_s("%6.4f",self.r_work()    )
    r_free = n_as_s("%6.4f",self.r_free()    )
    scale  = n_as_s("%6.3f",self.scale_k1_w())
    k_sol  = n_as_s("%4.2f",self.k_sol())
    b_sol  = n_as_s("%6.2f",self.b_sol())
    b0,b1,b2,b3,b4,b5 = n_as_s("%7.2f",self.b_cart())
    b_iso  = n_as_s("%7.2f",self.b_iso())
    err    = n_as_s("%6.2f",self.model_error_ml())
    try:    target_work = n_as_s("%.4g",self.target_w())
    except: target_work = str(None)
    line = "| r_work= "+r_work+"   r_free= "+r_free+"   ksol= "+k_sol+\
           "   Bsol= "+b_sol+"   scale= "+scale
    np = 79 - (len(line) + 1)
    if(np < 0): np = 0
    print >> out, line + p*np + "|"
    print >> out, "| "+"  "*38+"|"
    print >> out, "| overall anisotropic scale matrix (Cartesian basis):    "\
                  "                     |"
    c = ","
    line4 = "| (B11,B22,B33,B12,B13,B23)= ("+b0+c+b1+c+b2+c+b3+c+b4+c+b5+")"
    np = 79 - (len(line4) + 1)
    line4 = line4 + " "*np + "|"
    print >> out, line4
    line5 = "| (B11+B22+B33)/3 = "+b_iso
    np = 79 - (len(line5) + 1)
    line5 = line5 + " "*np + "|"
    print >> out, line5
    print >> out, "| "+"  "*38+"|"
    line6="| Target ("+self.target_name+")= "+target_work+\
          " | ML estimate for coordinates error: "+err+" A"
    np = 79 - (len(line6) + 1)
    line6 = line6 + " "*np + "|"
    print >> out, line6
    print >> out, "|"+"-"*77+"|"
    out.flush()
    time_show += timer.elapsed()

  def r_work_scale_k1_completeness_in_bins(self, reflections_per_bin = 500,
                                                 n_bins              = 0,
                                                 prefix              = "",
                                                 out                 = None):
    if(out is None): out = sys.stdout
    f_obs_w = self.f_obs_w
    reflections_per_bin = min(f_obs_w.data().size(), reflections_per_bin)
    f_obs_w.setup_binner(reflections_per_bin = reflections_per_bin,
                         n_bins              = n_bins)
    print >> out, prefix+"Bin     Resolution    Compl.  No.   Scale_k1(work) R-factor(work)"
    print >> out, prefix+"number     range              refl."
    for i_bin in f_obs_w.binner().range_used():
        sel         = f_obs_w.binner().selection(i_bin)
        r_work      = self.r_work(selection = sel)
        scale_k1    = self.scale_k1_w(selection = sel)
        f_obs_sel   = f_obs_w.select(sel)
        d_max,d_min = f_obs_sel.d_max_min()
        compl       = f_obs_sel.completeness(d_max = d_max)
        n_ref       = sel.count(True)
        d_range     = f_obs_w.binner().bin_legend(
                   i_bin = i_bin, show_bin_number = False, show_counts = False)
        format = "%3d: %-17s %5.2f %6d %6.3f     %10.4f"
        print >> out, prefix+format % (i_bin,d_range,compl,n_ref,scale_k1,
                                                                        r_work)
    print >> out, prefix+"where:"
    print >> out, prefix+"   R-factor = SUM(||Fobs|-Scale_k1*|Fmodel||)/SUM(|Fobs|)"
    print >> out, prefix+"   Scale_k1 = SUM(|Fobs|*|Fmodel|)/SUM(|Fmodel|**2)"
    print >> out, prefix+"   Fmodel   = fb_cart * (Fcalc + Fbulk)"
    print >> out, prefix+"   Fbulk    = k_sol * exp(-b_sol*s**2/4) * Fmask"
    print >> out, prefix+"   Fcalc    = structure factors calculated from atomic model"
    print >> out, prefix+"   fb_cart  = exp(-h(t) * A(-1) * B_cart * A(-1t) * h), "
    print >> out, prefix+"   A - orthogonalization matrix"
    out.flush()


  def fft_vs_direct(self, reflections_per_bin = 250,
                          n_bins              = 0,
                          out                 = None):
    if(out is None): out = sys.stdout
    f_obs_w = self.f_obs_w
    f_obs_w.setup_binner(reflections_per_bin = reflections_per_bin,
                         n_bins              = n_bins)
    fmodel_dc = self.deep_copy()
    if(self.sf_algorithm=="fft"):      fmodel_dc.update(sf_algorithm="direct")
    elif(self.sf_algorithm=="direct"): fmodel_dc.update(sf_algorithm="fft")
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
        if(self.sf_algorithm == "fft"):
           print >> out, format % (i_bin,d_range,compl,n_ref,
                                            scale_k1_2,r_work_1,r_work_2,delta)
        else:
           print >> out, format % (i_bin,d_range,compl,n_ref,
                                            scale_k1_1,r_work_2,r_work_1,delta)
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_comprehensive(self, header = "",
                               free_reflections_per_bin = 140,
                               max_number_of_bins  = 30,
                               out=None):
    if(out is None): out = sys.stdout
    self.show_essential(header = header, out = out)
    print >> out
    self.statistics_in_resolution_bins(
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins  = max_number_of_bins,
      out                 = out)
    print >> out
    self.show_fom_phase_error_alpha_beta_in_bins(
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins  = max_number_of_bins,
      out                 = out)

  def statistics_in_resolution_bins(self, free_reflections_per_bin = 140,
                                          max_number_of_bins  = 30,
                                          out=None):
    statistics_in_resolution_bins(
      fmodel          = self,
      target_functors = self.target_functors,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins  = max_number_of_bins,
      out=out)

  def r_factors_in_resolution_bins(self, free_reflections_per_bin = 140,
                                          max_number_of_bins  = 30,
                                          out=None):
    if(out is None): out = sys.stdout
    r_factors_in_resolution_bins(
      fmodel              = self,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins  = max_number_of_bins,
      out=out)

  def show_fom_phase_error_alpha_beta_in_bins(self,
        free_reflections_per_bin = 140,
        max_number_of_bins = 30,
        out=None):
    if(out is None): out = sys.stdout
    show_fom_phase_error_alpha_beta_in_bins(
      fmodel              = self,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins  = max_number_of_bins,
      out=out)

  def dump(self, out = None):
    if(out is None): out = sys.stdout
    f_model            = self.f_model_scaled_with_k1()
    f_model_amplitudes = f_model.amplitudes().data()
    f_model_phases     = f_model.phases(deg=True).data()
    f_calc_amplitudes  = self.f_calc().amplitudes().data()
    f_calc_phases      = self.f_calc().phases(deg=True).data()
    f_mask_amplitudes  = self.f_mask().amplitudes().data()
    f_mask_phases      = self.f_mask().phases(deg=True).data()
    f_bulk_amplitudes  = self.f_bulk().amplitudes().data()
    f_bulk_phases      = self.f_bulk().phases(deg=True).data()
    fom                = self.figures_of_merit()
    alpha, beta        = [item.data() for item in self.alpha_beta()]
    f_obs              = self.f_obs.data()
    flags              = self.r_free_flags.data()
    indices            = self.f_obs.indices()
    resolution         = self.f_obs.d_spacings().data()
    fb_cart            = self.fb_cart()
    scale_k1           = self.scale_k1()
    print >> out, "Fmodel   = scale_k1 * fb_cart * (Fcalc + Fbulk)"
    print >> out, "Fbulk    = k_sol * exp(-b_sol*s**2/4) * Fmask"
    print >> out, "Fcalc    = structure factors calculated from atomic model"
    print >> out, "fb_cart  = exp(-h(t) * A(-1) * B_cart * A(-1t) * h), "+\
                  "A - orthogonalization matrix"
    print >> out, "scale_k1 = SUM(Fobs * Fmodel) / SUM(Fmodel**2)"
    print >> out, "k_sol    = %10.4f"%self.k_sol()
    print >> out, "b_sol    = %10.4f"%self.b_sol()
    print >> out, "scale_k1 = %10.4f"%self.scale_k1()
    b0,b1,b2,b3,b4,b5 = n_as_s("%7.4f",self.b_cart())
    c=","
    print >>out,"(B11,B22,B33,B12,B13,B23) = ("+b0+c+b1+c+b2+c+b3+c+b4+c+b5+")"
    format="inde= %5d%5d%5d Fobs= %11.4f Fcalc= %10.3f %9.3f Fmask= %10.3f %9.3f"\
           " Fmodel= %13.6f %12.6f fom= %5.3f alpha= %8.4f beta= %12.4f"\
           " R-free-flags= %1d resol= %7.3f Fbulk= %10.3f %9.3f fb_cart= %8.4f "
    for ind, f_obs_i, f_calc_amplitudes_i, f_calc_phases_i, f_mask_amplitudes_i,\
        f_mask_phases_i, f_model_amplitudes_i, f_model_phases_i, fom_i, alpha_i,\
        beta_i, flags_i, resolution_i, f_bulk_amplitudes_i, f_bulk_phases_i, fb_cart_i in \
        zip(indices, f_obs, f_calc_amplitudes, \
        f_calc_phases, f_mask_amplitudes, f_mask_phases, f_model_amplitudes,    \
        f_model_phases, fom, alpha, beta, flags, resolution, f_bulk_amplitudes, \
        f_bulk_phases, fb_cart):
        print >> out, format%(ind[0], ind[1], ind[2], f_obs_i,
          f_calc_amplitudes_i, f_calc_phases_i, f_mask_amplitudes_i,
          f_mask_phases_i, f_model_amplitudes_i, f_model_phases_i, fom_i,
          alpha_i, beta_i, flags_i, resolution_i, f_bulk_amplitudes_i,
          f_bulk_phases_i, fb_cart_i)

class phaser_sad_target_functor(object):

  def __init__(self, f_obs, r_free_flags, xray_structure, f_calc):
    adopt_init_args(self, locals())
    self.refine_sad_object = phaser.phenix_adaptors.sad_target.data_adaptor(
      f_obs=f_obs,
      r_free_flags=r_free_flags,
      verbose=False).target(xray_structure=xray_structure)
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
        f_obs=manager.f_obs,
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

class target_result(object):

  def __init__(self, manager, core_result):
    self.manager = manager
    self.core_result = core_result

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

  def d_target_d_site_cart(self):
    manager = self.manager
    xray.set_scatterer_grad_flags(
      scatterers=manager.xray_structure.scatterers(),
      site=True)
    d_t_d_fc = self.d_target_d_f_calc_work()
    return manager.structure_factor_gradients_w(
      u_iso_refinable_params=None,
      d_target_d_f_calc=d_t_d_fc.data(),
      xray_structure=manager.xray_structure,
      n_parameters=0,
      miller_set=d_t_d_fc,
      algorithm=manager.sf_algorithm).d_target_d_site_cart()

def statistics_in_resolution_bins(fmodel,
                                  target_functors,
                                  free_reflections_per_bin,
                                  max_number_of_bins,
                                  out=None):
  global time_show
  timer = user_plus_sys_time()
  if(out is None): out = sys.stdout
  fo_t = fmodel.f_obs_t
  fc_t = fmodel.f_model_t()
  fo_w = fmodel.f_obs_w
  fc_w = fmodel.f_model_w()
  alpha_w, beta_w = fmodel.alpha_beta_w()
  alpha_t, beta_t = fmodel.alpha_beta_t()
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
  print >> out, "|"+"-"*77+"|"
  print >> out, "| Bin     Resolution   Compl.  No. Refl.    R-factors          Targets        |"
  print >> out, "|number     range              work test   work   test        work        test|"
  for i_bin in fo_t.binner().range_used():
    sel_t = fo_t.binner().selection(i_bin)
    sel_w = fo_w.binner().selection(i_bin)
    sel_all = fmodel.f_obs.binner().selection(i_bin)
    sel_fo_all = fmodel.f_obs.select(sel_all)
    sel_fo_t = fo_t.select(sel_t)
    sel_fc_t = fc_t.select(sel_t)
    sel_fo_w = fo_w.select(sel_w)
    sel_fc_w = fc_w.select(sel_w)
    sel_alpha_t = alpha_t.select(sel_t)
    sel_beta_t  = beta_t.select(sel_t)
    sel_alpha_w = alpha_w.select(sel_w)
    sel_beta_w  = beta_w.select(sel_w)
    xray_target_functor_w = target_functors.target_functor_w(selection = sel_w)
    xray_target_functor_t = target_functors.target_functor_t(selection = sel_t)
    d_max_,d_min_ = sel_fo_all.d_max_min()
    ch = fmodel.f_obs.resolution_filter(d_min= d_min_,d_max= d_max_).completeness(d_max = d_max_)
    if(fmodel.target_attributes().family == "ls"):
      target_w = "%11.5g" % xray_target_functor_w(sel_fc_w, False).target_work()
      target_t = "%11.5g" % xray_target_functor_t(sel_fc_t, False).target_work()
    elif(fmodel.target_name in ["ml", "mlhl"]):
      if (sel_fc_w.indices().size() != 0):
        target_w = "%11.5g" % xray_target_functor_w(
          sel_fc_w,
          sel_alpha_w.data(),
          sel_beta_w.data(),
          1.0,
          False).target_work()
      else:
        target_w = "%11s" % "None"
      if (sel_fc_t.indices().size() != 0):
        target_t = "%11.5g" % xray_target_functor_t(
          sel_fc_t,
          sel_alpha_t.data(),
          sel_beta_t.data(),
          1.0,
          False).target_work()
      else:
        target_t = "%11s" % "None"
    nw = sel_fo_w.data().size()
    nt = sel_fo_t.data().size()
    if (nw != 0):
      r_w = "%6.4f" % bulk_solvent.r_factor(sel_fo_w.data(), sel_fc_w.data())
    else:
      r_w = "  None"
    if (nt != 0):
      r_t = "%6.4f" % bulk_solvent.r_factor(sel_fo_t.data(), sel_fc_t.data())
    else:
      r_t = "  None"
    d_range = fo_t.binner().bin_legend(
      i_bin=i_bin, show_bin_number=False, show_counts=False)
    print >> out, "|%3d: %-17s %4.2f %6d %4d %s %s %s %s|" % (
      i_bin, d_range, ch, nw, nt, r_w, r_t, target_w, target_t)
  print >> out, "|"+"-"*77+"|"
  out.flush()
  time_show += timer.elapsed()

def r_factors_in_resolution_bins(fmodel,
                                 free_reflections_per_bin,
                                 max_number_of_bins,
                                 out=None):
  global time_show
  timer = user_plus_sys_time()
  if(out is None): out = sys.stdout
  fo_t = fmodel.f_obs_t
  fc_t = fmodel.f_model_t()
  fo_w = fmodel.f_obs_w
  fc_w = fmodel.f_model_w()
  fo_t.setup_binner(n_bins=fmodel.determine_n_bins(
    free_reflections_per_bin=free_reflections_per_bin,
    max_n_bins=max_number_of_bins))
  fo_w.use_binning_of(fo_t)
  fc_t.use_binning_of(fo_t)
  fc_w.use_binning_of(fo_t)
  print >> out, " Bin     Resolution       No. Refl.      R-factors"
  print >> out, "number     range         work   test     work   test"
  for i_bin in fo_t.binner().range_used():
    sel_t = fo_t.binner().selection(i_bin)
    sel_w = fo_w.binner().selection(i_bin)
    sel_fo_t = fo_t.select(sel_t)
    sel_fc_t = fc_t.select(sel_t)
    sel_fo_w = fo_w.select(sel_w)
    sel_fc_w = fc_w.select(sel_w)
    r_w = bulk_solvent.r_factor(sel_fo_w.data(), sel_fc_w.data())
    r_t = bulk_solvent.r_factor(sel_fo_t.data(), sel_fc_t.data())
    nt = sel_fo_t.data().size()
    nw = sel_fo_w.data().size()
    d_range = fo_t.binner().bin_legend(
      i_bin=i_bin, show_bin_number=False, show_counts=False)
    print >> out, "%3d: %-17s %6d %6d   %6.4f %6.4f" % (
      i_bin, d_range, nw, nt, r_w, r_t)
  out.flush()
  time_show += timer.elapsed()


def show_fom_phase_error_alpha_beta_in_bins(fmodel,
                                            free_reflections_per_bin,
                                            max_number_of_bins,
                                            out=None):
  global time_show
  timer = user_plus_sys_time()
  if(out is None): out = sys.stdout
  mi_fom = fmodel.f_obs.array(data = fmodel.figures_of_merit())
  mi_fom.setup_binner(n_bins=fmodel.determine_n_bins(
    free_reflections_per_bin=free_reflections_per_bin,
    max_n_bins=max_number_of_bins))
  phase_errors_work = fmodel.phase_errors_work()
  phase_errors_test = fmodel.phase_errors_test()
  alpha, beta = fmodel.alpha_beta()
  mi_per_work = fmodel.f_obs_w.array(data = phase_errors_work)
  mi_per_test = fmodel.f_obs_t.array(data = phase_errors_test)
  mi_per_test.use_binning_of(mi_fom)
  mi_per_work.use_binning_of(mi_fom)
  print >> out, "|"+"-"*77+"|"
  print >> out, "|R-free likelihood based estimates for figures of merit," \
                  " absolute phase error,|"
  print >> out, "|and distribution parameters alpha and beta" \
                  " (Acta Cryst. (1995). A51, 880-887)|"
  print >> out, "|"+" "*77+"|"
  print >> out, "| Bin     Resolution      No. Refl.   FOM   Phase error   "\
                " Alpha        Beta  |"
  print >> out, "|  #        range        work  test        work    test  "\
                "                     |"
  for i_bin in mi_fom.binner().range_used():
    sel = mi_fom.binner().selection(i_bin).iselection()
    sel_work = mi_per_work.binner().selection(i_bin).iselection()
    sel_test = mi_per_test.binner().selection(i_bin).iselection()
    assert sel.size() == sel_work.size() + sel_test.size() or \
           2*sel.size() == sel_work.size() + sel_test.size()
    print >> out, "|%3d: %-17s%6d%6d%s%s%s%s%s|" % (
      i_bin,
      mi_fom.binner().bin_legend(
        i_bin=i_bin, show_bin_number=False, show_counts=False),
      sel_work.size(),
      sel_test.size(),
      mi_fom.data().select(sel).format_mean("%6.2f"),
      mi_per_work.data().select(sel_work).format_mean("%7.2f"),
      mi_per_test.data().select(sel_test).format_mean("%7.2f"),
      alpha.data().select(sel).format_mean("%9.2f"),
      beta.data().select(sel).format_mean("%14.2f"))
  alpha_stats = alpha.data().min_max_mean()
  beta_stats = beta.data().min_max_mean()
  print >>out, "|alpha:            min =%12.2f max =%16.2f mean =%13.2f|"%\
    alpha_stats.as_tuple()
  print >>out, "|beta:             min =%12.2f max =%16.2f mean =%13.2f|"%\
    beta_stats.as_tuple()
  print >>out, "|figures of merit: min =%12.2f max =%16.2f mean =%13.2f|"%\
    mi_fom.data().min_max_mean().as_tuple()
  print >>out, "|phase err.(work): min =%12.2f max =%16.2f mean =%13.2f|"%\
    phase_errors_work.min_max_mean().as_tuple()
  print >>out, "|phase err.(test): min =%12.2f max =%16.2f mean =%13.2f|"%\
    phase_errors_test.min_max_mean().as_tuple()
  if(alpha_stats.min <= 0.0):
    print >> out, "| *** f_model warning: there are some alpha <= 0.0 ***" \
      "                        |"
    amz = alpha.data() <= 0.0
    print >> out, "|                      number of alpha <= 0.0: %6d" \
      "                         |" % (amz.count(True))
    bmz = beta.data() <= 0.0
  if(beta_stats.min <= 0.0):
    print >> out, "| *** f_model warning: there are some beta <= 0.0 ***" \
      "                         |"
    bmz = beta.data() <= 0.0
    print >> out, "|   number of beta <= 0.0: %6d |" % (bmz.count(True))
  print >> out, "|"+"-"*77+"|"
  out.flush()
  time_show += timer.elapsed()

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
  vt = type(value).__name__
  if(vt in ["int","float"]):
     return str(format%value).strip()
  else:
     new = []
     for item in value:
       new.append( str(format%item).strip() )
     return new
