from scitbx.python_utils.misc import adopt_init_args
from cctbx.xray import ext
from cctbx.array_family import flex

class target_functor_base:

  def __call__(self, f_calc, compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(
           self.f_obs().unit_cell())
    assert f_calc.space_group() == self.f_obs().space_group()
    if (self.weights() is not None):
      return self._target_calculator(self.f_obs().data(),
                                     self.weights(),
                                     f_calc.data(),
                                     compute_derivatives)
    else:
      return self._target_calculator(self.f_obs().data(),
                                     f_calc.data(),
                                     compute_derivatives)

class target_functors_manager:

  def __init__(self, target_name,
                     f_obs,
                     flags,
                     abcd                   = None,
                     weights                = None,
                     use_sigmas_as_weights  = False,
                     scale_factor           = 0):
    adopt_init_args(self, locals())
    assert self.target_name in (
      "ls_wunit_k1","ls_wunit_k2","ls_wunit_kunit","ls_wunit_k1_fixed",
      "ls_wunit_k1ask3_fixed",
      "ls_wexp_k1" ,"ls_wexp_k2" ,"ls_wexp_kunit",
      "ls_wff_k1"  ,"ls_wff_k2"  ,"ls_wff_kunit","ls_wff_k1_fixed",
      "ls_wff_k1ask3_fixed",
      "lsm_k1"     ,"lsm_k2"     ,"lsm_kunit",
      "ml","mlhl")
    assert self.f_obs.data().size() == self.flags.size()
    if(self.flags.count(True) > 0):
      self.f_obs_w = self.f_obs.select(~self.flags)
      self.f_obs_t = self.f_obs.select( self.flags)
    else:
      self.f_obs_w = self.f_obs
      self.f_obs_t = self.f_obs
    if(self.target_name == "mlhl"):
      assert self.abcd is not None
    if(self.abcd is not None):
      if(self.target_name == "mlhl"):
        if(self.flags.count(True) > 0):
          self.abcd_w = self.abcd.select(~self.flags)
          self.abcd_t = self.abcd.select( self.flags)
        else:
          self.abcd_w = self.abcd
          self.abcd_t = self.abcd
      else:
        self.abcd_w, self.abcd_t = None, None
    else:
      self.abcd_w, self.abcd_t = None, None
    if(self.weights is not None):
      assert self.target_name.count("ls") == 1
      if(self.flags.count(True) > 0):
        self.weights_w = self.weights.select(~self.flags)
        self.weights_t = self.weights.select( self.flags)
      else:
        self.weights_w = self.weights
        self.weights_t = self.weights
    else:
      self.weights_w, self.weights_t = None, None

  def target_functor_w(self, selection = None):
    if(selection is None):
      f_obs   = self.f_obs_w
      weights = self.weights_w
      abcd = self.abcd_w
    else:
      assert selection.size() == self.f_obs_w.data().size()
      f_obs   = self.f_obs_w.select(selection)
      if(self.weights_w is not None): weights = self.weights_w.select(selection)
      else:                           weights = self.weights_w
      if(self.abcd_w is not None): abcd = self.abcd_w.select(selection)
      else:                        abcd = self.abcd_w
    if(self.target_name.count("k1") == 1 and self.target_name.count("k1as") == 0):
       assert self.scale_factor == 0
       return ls_k1(f_obs            = f_obs,
                    weights          = weights,
                    scale_factor     = self.scale_factor,
                    fix_scale_factor = False)
    if(self.target_name.count("k2") == 1 and self.target_name.count("k2as") == 0):
       assert self.scale_factor == 0
       return ls_k2(f_obs            = f_obs,
                    weights          = weights,
                    scale_factor     = self.scale_factor,
                    fix_scale_factor = False)
    if(self.target_name.count("kunit") == 1 or self.target_name.count("ask") == 1):
       assert self.scale_factor != 0.0
       return ls_k1(f_obs            = f_obs,
                    weights          = weights,
                    scale_factor     = self.scale_factor,
                    fix_scale_factor = True)
    if(self.target_name == "ml"):
      return maximum_likelihood_criterion(f_obs = f_obs)
    if(self.target_name == "mlhl"):
      return maximum_likelihood_criterion_hl(
                              f_obs                = f_obs,
                              abcd                 = abcd.data())

  def target_functor_t(self, selection = None):
    if(selection is None):
      f_obs   = self.f_obs_t
      weights = self.weights_t
      abcd = self.abcd_t
    else:
      assert selection.size() == self.f_obs_t.data().size()
      f_obs   = self.f_obs_t.select(selection)
      if(self.weights_t is not None): weights = self.weights_t.select(selection)
      else:                           weights = self.weights_t
      if(self.abcd_t is not None): abcd = self.abcd_t.select(selection)
      else:                        abcd = self.abcd_t
    if(self.target_name.count("k1") == 1 and self.target_name.count("k1as") == 0):
       assert self.scale_factor == 0
       return ls_k1(f_obs            = f_obs,
                    weights          = weights,
                    scale_factor     = self.scale_factor,
                    fix_scale_factor = False)
    if(self.target_name.count("k2") == 1 and self.target_name.count("k2as") == 0):
       assert self.scale_factor == 0
       return ls_k2(f_obs            = f_obs,
                    weights          = weights,
                    scale_factor     = self.scale_factor,
                    fix_scale_factor = False)
    if(self.target_name.count("kunit") == 1 or self.target_name.count("ask") == 1):
       assert self.scale_factor != 0.0
       return ls_k1(f_obs            = f_obs,
                    weights          = weights,
                    scale_factor     = self.scale_factor,
                    fix_scale_factor = True)
    if(self.target_name == "ml"):
      return maximum_likelihood_criterion(f_obs = f_obs)
    if(self.target_name == "mlhl"):
      return maximum_likelihood_criterion_hl(
                              f_obs                = f_obs,
                              abcd                 = abcd.data())

class ls_k1:

  def __init__(self, f_obs,
                     weights,
                     scale_factor,
                     fix_scale_factor):
    adopt_init_args(self, locals(), hide=True)
    if(self._fix_scale_factor == True):
       assert self._scale_factor > 0.0

  def f_obs(self):
    return self._f_obs

  def weights(self):
    return self._weights

  def scale_factor(self):
    return self._scale_factor

  def __call__(self, f_calc, compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(
           self.f_obs().unit_cell())
    assert f_calc.space_group() == self.f_obs().space_group()
    assert f_calc.data().size() == self.f_obs().data().size()
    return ext.ls_target_with_scale_k1(
                                  f_obs               = self.f_obs().data(),
                                  weights             = self.weights(),
                                  f_calc              = f_calc.data(),
                                  compute_derivatives = compute_derivatives,
                                  fix_scale           = self._fix_scale_factor,
                                  scale               = self._scale_factor)

class ls_k2:

  def __init__(self, f_obs,
                     weights,
                     scale_factor,
                     fix_scale_factor):
    adopt_init_args(self, locals(), hide=True)
    if(self._fix_scale_factor == True):
       assert self._scale_factor > 0.0

  def f_obs(self):
    return self._f_obs

  def weights(self):
    return self._weights

  def scale_factor(self):
    return self._scale_factor

  def __call__(self, f_calc, compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(
           self.f_obs().unit_cell())
    assert f_calc.space_group() == self.f_obs().space_group()
    assert f_calc.data().size() == self.f_obs().data().size()
    return ext.ls_target_with_scale_k2(
                                  f_obs               = self.f_obs().data(),
                                  weights             = self.weights(),
                                  f_calc              = f_calc.data(),
                                  compute_derivatives = compute_derivatives,
                                  fix_scale           = self._fix_scale_factor,
                                  scale               = self._scale_factor)

class least_squares_residual:

  def __init__(self, f_obs,
                     weights               = None,
                     use_sigmas_as_weights = False,
                     scale_factor          = 0):
    adopt_init_args(self, locals(), hide=True)
    assert self._weights is None or self._use_sigmas_as_weights == False
    if (self._use_sigmas_as_weights):
      sigmas_squared = flex.pow2(self._f_obs.sigmas())
      assert sigmas_squared.all_gt(0)
      self._weights = 1 / sigmas_squared

  def f_obs(self):
    return self._f_obs

  def weights(self):
    return self._weights

  def use_sigmas_as_weights(self):
    return self._use_sigmas_as_weights

  def __call__(self, f_calc, compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(
           self.f_obs().unit_cell())
    assert f_calc.space_group() == self.f_obs().space_group()
    if (self.weights() is not None):
      return ext.targets_least_squares_residual(
        self.f_obs().data(),
        self.weights(),
        f_calc.data(),
        compute_derivatives,
        self._scale_factor)
    else:
      return ext.targets_least_squares_residual(
        self.f_obs().data(),
        f_calc.data(),
        compute_derivatives,
        self._scale_factor)

class intensity_correlation(target_functor_base):

  def __init__(self, f_obs, weights=None,
               use_multiplicities_as_weights=False):
    adopt_init_args(self, locals(), hide=True)
    assert self._weights is None or self._use_multiplicities_as_weights==False
    self._target_calculator = ext.targets_intensity_correlation
    if (self._use_multiplicities_as_weights):
      self._weights = self._f_obs.multiplicities().data()

  def f_obs(self):
    return self._f_obs

  def weights(self):
    return self._weights

  def use_multiplicities_as_weights(self):
    return self._use_multiplicities_as_weights

class maximum_likelihood_criterion:

  def __init__(self, f_obs):
    adopt_init_args(self, locals(), hide=True)
    self.epsilons = f_obs.epsilons().data()
    self.centric_flags = flex.int(flex.to_list(f_obs.centric_flags().data()))

  def f_obs(self):
    return self._f_obs

  def __call__(self, f_calc,
                     alpha,
                     beta,
                     k,
                     compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(self.f_obs().unit_cell())
    assert f_calc.space_group() == self.f_obs().space_group()
    return ext.targets_maximum_likelihood_criterion(
        self.f_obs().data(),
        f_calc.data(),
        alpha,
        beta,
        k,
        self.epsilons,
        self.centric_flags,
        compute_derivatives)

class maximum_likelihood_criterion_hl:

  def __init__(self, f_obs,
                     abcd,
                     step_for_integration = 5.0):
    adopt_init_args(self, locals(), hide=True)

  def f_obs(self):
    return self._f_obs

  def abcd(self):
    return self._abcd

  def step_for_integration(self):
    return self._step_for_integration

  def __call__(self, f_calc,
                     alpha,
                     beta,
                     k,
                     compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(self.f_obs().unit_cell())
    assert f_calc.space_group() == self.f_obs().space_group()
    return ext.targets_maximum_likelihood_criterion_hl(
        self.f_obs().data(),
        f_calc.data(),
        alpha,
        beta,
        f_calc.epsilons().data(),
        flex.int(flex.to_list(f_calc.centric_flags().data())),
        compute_derivatives,
        self.abcd(),
        self.step_for_integration())

def registry():
  return {
    "least_squares_residual": least_squares_residual,
    "intensity_correlation": intensity_correlation,
  }
