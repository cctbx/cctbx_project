import libtbx.load_env
if(not libtbx.env.has_module(name="phaser")):
  phaser = None
else:
  import phaser.phenix_adaptors.sad_target
from cctbx.array_family import flex
from cctbx import xray
import boost.python
from libtbx.utils import Sorry, user_plus_sys_time
from cctbx.eltbx.xray_scattering import wk1995
from cctbx import adptbx
from libtbx import adopt_init_args


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

class target_attributes(object):

  def __init__(self, family, specialization=None):
    adopt_init_args(self, locals())
    assert self.validate()

  def validate(self):
    if (self.family == "lsm"):
      self.family = "ls"
      self.pseudo_ml = True
    else:
      self.pseudo_ml = False
    if (self.family == "ls"):
      return self.specialization is None
    elif (self.family == "ml"):
      return self.specialization in [None, "hl", "sad"]
    return False

  def requires_experimental_phases(self):
    return (self.family == "ml" and self.specialization == "hl")

target_names = {
  "ls_wunit_k1": target_attributes("ls"),
  "ls_wunit_kunit": target_attributes("ls"),
  "ls_wunit_k1_fixed": target_attributes("ls"),
  "ls_wunit_k1ask3_fixed": target_attributes("ls"),
  "ls_wexp_k1": target_attributes("ls"),
  "ls_wexp_kunit": target_attributes("ls"),
  "ls_wff_k1": target_attributes("ls"),
  "ls_wff_kunit": target_attributes("ls"),
  "ls_wff_k1_fixed": target_attributes("ls"),
  "ls_wff_k1ask3_fixed": target_attributes("ls"),
  "lsm_k1": target_attributes("lsm"),
  "lsm_kunit": target_attributes("lsm"),
  "lsm_k1_fixed": target_attributes("lsm"),
  "lsm_k1ask3_fixed": target_attributes("lsm"),
  "ml": target_attributes("ml"),
  "mlhl": target_attributes("ml", "hl"),
  "ml_sad": target_attributes("ml", "sad")}

class phaser_sad_target_functor(object):

  def __init__(self,
        f_obs,
        r_free_flags,
        xray_structure,
        f_calc,
        target_memory):
    self.f_obs = f_obs
    self.r_free_flags = r_free_flags
    self.xray_structure = xray_structure
    self.f_calc = f_calc
    if (target_memory is None): # XXX could be more elegant!
      den = self.f_obs.data()
      num = flex.abs(self.f_calc.data())
      denom = flex.sum(num*den)
      numerator = flex.sum(den*den)
      if (denom == 0):
        raise RuntimeError("Zero denominator in scale calculation.")
      previous_overall_scaleK = numerator/denom
      previous_overall_scaleU = 0.
      previous_variances = None
      adaptor = phaser.phenix_adaptors.sad_target.data_adaptor(
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        verbose=True)
      self.refine_sad_object = adaptor.target(
        xray_structure=xray_structure,
        previous_overall_scaleK=previous_overall_scaleK,
        previous_overall_scaleU=previous_overall_scaleU,
        previous_variances=previous_variances)
      self.refine_sad_object.set_f_calc(f_calc=f_calc)
      target_memory = self.target_memory()

    assert len(target_memory) == 4
    assert target_memory[0] == "ml_sad"
    previous_overall_scaleK = target_memory[1]
    previous_overall_scaleU = target_memory[2]
    previous_variances = target_memory[3]
    adaptor = phaser.phenix_adaptors.sad_target.data_adaptor(
      f_obs=f_obs,
      r_free_flags=r_free_flags,
      verbose=True)
    self.refine_sad_object = adaptor.target(
      xray_structure=xray_structure,
      previous_overall_scaleK=previous_overall_scaleK,
      previous_overall_scaleU=previous_overall_scaleU,
      previous_variances=previous_variances)
    self.refine_sad_object.set_f_calc(f_calc=f_calc)
    self.refine_sad_object.reject_outliers()

  def prepare_for_minimization(self):
    rso = self.refine_sad_object
    rso.refine_variance_terms()
    self.refined_overall_b_iso = adptbx.u_as_b(
          rso.refine_sad_instance.get_refined_scaleU())
    rso.refine_sad_instance.set_scaleU(0.)

  def target_memory(self):
    rsi = self.refine_sad_object.refine_sad_instance
    return ("ml_sad", rsi.get_refined_scaleK(),
                      rsi.get_refined_scaleU(),rsi.get_variance_array())

  def __call__(self, f_calc, compute_gradients):
    self.refine_sad_object.set_f_calc(f_calc=f_calc)
    rso = self.refine_sad_object
    target_work = rso.functional(use_working_set=True)
    da_db, daa_dbb_dab = rso.derivatives(curvs=True)
    target_test = rso.functional(use_working_set=False)
    return xray.targets_common_results(
      target_per_reflection=flex.double(),
      target_work=target_work,
      target_test=target_test,
      gradients_work=da_db.data(),
      curvatures_work=daa_dbb_dab.data())

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
        f_obs=manager.f_obs(),
        r_free_flags=manager.r_free_flags(),
        xray_structure=manager.xray_structure,
        f_calc=manager.f_model(),
        target_memory=manager._target_memory)
      manager._target_memory = self.core.target_memory()
    elif (attr.family == "ml"):
      if (attr.requires_experimental_phases()):
        experimental_phases = manager.hl_coeffs()
      else:
        experimental_phases = None
      self.core = xray.target_functors.max_like(
        f_obs=manager.f_obs(),
        r_free_flags=manager.r_free_flags(),
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
        elif (target_name == "lsm_k1ask3_fixed"):
          scale_factor = manager.scale_k3_w()
        elif (target_name == "lsm_k1_fixed"):
          scale_factor = manager.scale_k1_w()
        elif (target_name == "lsm_kunit"):
          scale_factor = 1.0
        else:
          raise RuntimeError
      else:
        f_obs = manager.f_obs()
        if (target_name.startswith("ls_wunit_")):
          weights = flex.double(f_obs.data().size(), 1.0)
          if   (target_name == "ls_wunit_k1"):
            scale_factor = 0
          elif (target_name == "ls_wunit_k1_fixed"):
            scale_factor = manager.scale_k1_w()
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
          elif (target_name == "ls_wff_kunit"):
            scale_factor = 1.0
          else:
            raise RuntimeError
        else:
          raise RuntimeError
      self.core = xray.target_functors.least_squares(
        compute_scale_using_all_data=False,
        f_obs=f_obs,
        r_free_flags=manager.r_free_flags(),
        weights=weights,
        scale_factor=scale_factor)

  def prepare_for_minimization(self):
    if (self.manager.target_name == "ml_sad"):
      self.core.prepare_for_minimization()
      self.manager.adopt_external_b_iso_adjustments(
        overall_b_iso_shift=self.core.refined_overall_b_iso)

  def target_function_is_invariant_under_allowed_origin_shifts(self):
    return (self.manager.target_name != "mlhl")

  def __call__(self, compute_gradients=False):
    result = target_result(
      manager=self.manager,
      core_result=self.core(
        f_calc=self.manager.f_model(),
        compute_gradients=compute_gradients))
    target_memory = getattr(self.core, "target_memory", None)
    if (target_memory is not None):
      self.manager._target_memory = target_memory()
    return result

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
        n_parameters=xray_structure.n_parameters(),
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
    return self.manager.f_obs_work().array(
      data=self.core_result.gradients_work())

  def d_target_d_f_calc_work(self):
    return self.manager.f_obs_work().array(
      data=self.core_result.gradients_work()
          *self.manager.fb_cart_work())

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
