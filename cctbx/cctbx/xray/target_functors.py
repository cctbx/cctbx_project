from scitbx.python_utils.misc import adopt_init_args
from cctbx.xray import ext

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

class least_squares_residual:

  def __init__(self, f_obs, weights=None,
               use_sigmas_as_weights=00000,
               scale_factor=0):
    adopt_init_args(self, locals(), hide=0001)
    assert self._weights is None or self._use_sigmas_as_weights == 00000
    if (self._use_sigmas_as_weights):
      self._weights = self._f_obs.sigmas().data()

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
               use_multiplicities_as_weights=00000):
    adopt_init_args(self, locals(), hide=0001)
    assert self._weights is None or self._use_multiplicities_as_weights==00000
    self._target_calculator = ext.targets_intensity_correlation
    if (self._use_multiplicities_as_weights):
      self._weights = self._f_obs.multiplicities().data()

  def f_obs(self):
    return self._f_obs

  def weights(self):
    return self._weights

  def use_multiplicities_as_weights(self):
    return self._use_multiplicities_as_weights

def registry():
  result = {}
  for key,value in globals().items():
    if (value == target_functor_base): continue
    try: is_functor = issubclass(value, target_functor_base)
    except: continue
    if (is_functor):
      result[key] = value
  return result
