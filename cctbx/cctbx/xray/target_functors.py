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
                     step_for_integration   = None,
                     weights                = None,
                     use_sigmas_as_weights  = False,
                     scale_factor           = 0):
    adopt_init_args(self, locals())
    assert self.target_name == "ml" or self.target_name == "ls" or \
           self.target_name == "mlhl"
    assert self.f_obs.data().size() == self.flags.size()
    if(self.flags.count(True) > 0):
      self.f_obs_w = self.f_obs.select(~self.flags)
      self.f_obs_t = self.f_obs.select( self.flags)
    else:
      self.f_obs_w = self.f_obs
      self.f_obs_t = self.f_obs
    if(self.target_name == "ml"):
      self.tf_w = maximum_likelihood_criterion(f_obs = self.f_obs_w)
      self.tf_t = maximum_likelihood_criterion(f_obs = self.f_obs_t)
    if(self.target_name == "mlhl"):
      assert self.abcd is not None and self.step_for_integration is not None
      if(self.flags.count(True) > 0):
        self.abcd_w = self.abcd.select(~self.flags)
        self.abcd_t = self.abcd.select( self.flags)
      else:
        self.abcd_w = self.abcd
        self.abcd_t = self.abcd
      self.tf_w = maximum_likelihood_criterion_hl(
                              f_obs                = self.f_obs_w,
                              abcd                 = self.abcd_w.data(),
                              step_for_integration = self.step_for_integration)
      self.tf_t = maximum_likelihood_criterion_hl(
                              f_obs                = self.f_obs_t,
                              abcd                 = self.abcd_t.data(),
                              step_for_integration = self.step_for_integration)
    if(self.target_name == "ls"):
      if(self.weights is not None):
        if(self.flags.count(True) > 0):
          self.weights_w = self.weights.select(~self.flags)
          self.weights_t = self.weights.select( self.flags)
        else:
          self.weights_w = self.weights
          self.weights_t = self.weights
      else:
        self.weights_w, self.weights_t = None, None
      self.tf_w = least_squares_residual(
                            f_obs                 = self.f_obs_w,
                            weights               = self.weights_w,
                            use_sigmas_as_weights = self.use_sigmas_as_weights,
                            scale_factor          = self.scale_factor)
      self.tf_t = least_squares_residual(
                            f_obs                 = self.f_obs_t,
                            weights               = self.weights_t,
                            use_sigmas_as_weights = self.use_sigmas_as_weights,
                            scale_factor          = self.scale_factor)

  def target_functor_w(self):
    return self.tf_w

  def target_functor_t(self):
    return self.tf_t


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
                     step_for_integration):
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
