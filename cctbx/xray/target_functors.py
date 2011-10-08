from cctbx.xray import ext
from cctbx.xray import weighting_schemes
from cctbx import miller
from cctbx.array_family import flex
from libtbx import adopt_init_args

class least_squares(object):

  def __init__(self,
        compute_scale_using_all_data,
        f_obs,
        weights,
        r_free_flags,
        scale_factor):
    assert scale_factor >= 0
    adopt_init_args(self, locals())

  def __call__(self, f_calc, compute_gradients):
    assert f_calc.unit_cell().is_similar_to(self.f_obs.unit_cell())
    assert f_calc.space_group() == self.f_obs.space_group()
    return ext.targets_least_squares(
      compute_scale_using_all_data=self.compute_scale_using_all_data,
      obs_type="F",
      obs=self.f_obs.data(),
      weights=self.weights,
      r_free_flags=self.r_free_flags.data(),
      f_calc=f_calc.data(),
      derivatives_depth=int(compute_gradients),
      scale_factor=self.scale_factor)

class max_like(object):

  def __init__(self,
        f_obs,
        r_free_flags,
        experimental_phases,
        alpha_beta,
        scale_factor,
        integration_step_size=5.0):
    adopt_init_args(self, locals())
    self.epsilons = f_obs.epsilons().data()
    self.centric_flags = f_obs.centric_flags().data()

  def __call__(self, f_calc, compute_gradients):
    assert f_calc.unit_cell().is_similar_to(self.f_obs.unit_cell())
    assert f_calc.space_group() == self.f_obs.space_group()
    if (self.experimental_phases is None):
      return ext.targets_maximum_likelihood_criterion(
        f_obs=self.f_obs.data(),
        r_free_flags=self.r_free_flags.data(),
        f_calc=f_calc.data(),
        alpha=self.alpha_beta[0].data(),
        beta=self.alpha_beta[1].data(),
        scale_factor=self.scale_factor,
        epsilons=self.epsilons,
        centric_flags=self.centric_flags,
        compute_gradients=compute_gradients)
    return ext.targets_maximum_likelihood_criterion_hl(
      f_obs=self.f_obs.data(),
      r_free_flags=self.r_free_flags.data(),
      experimental_phases=self.experimental_phases.data(),
      f_calc=f_calc.data(),
      alpha=self.alpha_beta[0].data(),
      beta=self.alpha_beta[1].data(),
      epsilons=self.epsilons,
      centric_flags=self.centric_flags,
      integration_step_size=self.integration_step_size,
      compute_gradients=compute_gradients)

class unified_least_squares_residual(object):
  """ A least-square residual functor for refinement against F or F^2. """

  def __init__(self, obs,
               is_amplitude=True,
               weighting=None
               ):
    """
    Construct a least-square residuals

    .. |Sigma|  unicode:: U+003A3 .. GREEK CAPITAL LETTER SIGMA

    |Sigma|:sub:`i`  w[i] ( f_obs.data[i] - k abs(f_calc.data[i]) )^2
    /  |Sigma|:sub:`i` w[i] f_obs.data[i]^2

    or

    |Sigma|:sub:`i` w[i] ( f_obs_square.data[i] - k abs(f_calc.data[i])^2 )^2
    /  |Sigma|:sub:`i` w[i] f_obs_square.data[i]^2

    depending on which of f_obs and f_obs_square is not None.

    Note that
      - the sums are over the indices i of the reflections,
      - f_calc is to be passed to the __call__ method,
      - the weights w and the scale factor k are discussed below.

    :Parameters:

      obs : real miller.array
        the observed reflections, with F and sigma(F) (or F^2 and sigma(F^2))
        respectively in obs.data() and obs.sigmas()

      is_amplitude : bool
        a flag to discriminate the type of data in obs if the latter does not
        spell out whether it contains amplitudes or intensities;
        the default means that amplitudes is the default when data type is
        unknown.

      weighting
        a weighting scheme for the data (c.f. cctbx.xray.weighting for common
        ones) or None, which means no weights
    """
    if not(obs.is_xray_amplitude_array() or obs.is_xray_intensity_array()):
      self._obs = miller.array(obs, data=obs.data(), sigmas=obs.sigmas())
      if is_amplitude:
        self._obs.set_observation_type_xray_amplitude()
      else:
        self._obs.set_observation_type_xray_intensity()
    else:
      self._obs = obs
    assert(self._obs.is_xray_amplitude_array()
           or self._obs.is_xray_intensity_array())
    if self.obs().is_xray_amplitude_array():
      self._ext_ls_residual = ext.targets_least_squares_residual
      default_weighting = weighting_schemes.amplitude_unit_weighting()
    elif self.obs().is_xray_intensity_array():
      self._ext_ls_residual = ext.targets_least_squares_residual_for_intensity
      default_weighting = weighting_schemes.intensity_quasi_unit_weighting()
    if weighting is None: weighting = default_weighting
    self._weighting = weighting
    self._weighting.observed = self._obs
    self._scale_factor = None

  def obs(self):
    """ The obs passed to the constructor """
    return self._obs

  f_obs = obs
  """ For compatibility with the optimiser's interface
  (as minimization.lbfgs) """

  def weighting(self):
    """ The weighting scheme """
    return self._weighting

  def __call__(self, f_calc, compute_derivatives, scale_factor=0):
    """
    Compute the least-squares residual value and perhaps its derivatives
    wrt to the calculated structure factor F_c of the i-th reflection

    :Parameters:

      f_calc : complex `cctbx.miller.array`
        f_calc.data()[i] constains F_c for the i-th reflection in f_obs()
      compute_derivatives : bool
        whether to compute the derivatives of the least square residual or not
      scale_factor : number
        the scale factor k if not null,
        otherwise k will be computed as a by-product of calling this method

    :return:
       An object holding the residual value, derivatives and scale factor
    :rtype:
       Boost.Python binding of
       :doxyclass:`cctbx::xray::targets::least_squares_residual`
    """
    assert f_calc.unit_cell().is_similar_to(self.obs().unit_cell())
    assert f_calc.space_group() == self.obs().space_group()
    if self._scale_factor is None:
      self._scale_factor = 0
    else:
      self._scale_factor = scale_factor
    self.weighting().calculated = f_calc
    if scale_factor == 0: scale_factor = None
    if (isinstance(self.weighting(), weighting_schemes.shelx_weighting) or
        isinstance(self.weighting(), weighting_schemes.simple_shelx_weighting)):
      self.weighting().compute(f_calc, scale_factor)
      self._scale_factor = self.weighting().scale_factor
    else:
      self.weighting().compute()
    if (self.weighting().weights is not None):
      result = self._ext_ls_residual(
        self.obs().data(),
        self.weighting().weights,
        f_calc.data(),
        compute_derivatives,
        self._scale_factor)
    else:
      result = self._ext_ls_residual(
        self.obs().data(),
        f_calc.data(),
        compute_derivatives,
        self._scale_factor)
    self._scale_factor = result.scale_factor()
    return result


class least_squares_residual_for_amplitude(unified_least_squares_residual):
  """ A least-square residual functor for F refinement. """

  statistical_weighting = weighting_schemes.pure_statistical_weighting()

  def __init__(self, f_obs, use_sigmas_as_weights = False):
    if use_sigmas_as_weights:
      weighting = self.statistical_weighting
    else:
      weighting = None
    super(least_squares_residual, self).__init__(f_obs, True, weighting)

  def use_sigmas_as_weights(self):
    return isinstance(self._weighting, statistical_weighting)

least_squares_residual = least_squares_residual_for_amplitude


class intensity_correlation(object):

  def __init__(self, f_obs, weights=None,
               use_multiplicities_as_weights=False):
    adopt_init_args(self, locals(), hide=True)
    assert self._weights is None or not self._use_multiplicities_as_weights
    if (self._use_multiplicities_as_weights):
      self._weights = self._f_obs.multiplicities().data().as_double()

  def f_obs(self):
    return self._f_obs

  def weights(self):
    return self._weights

  def use_multiplicities_as_weights(self):
    return self._use_multiplicities_as_weights

  def __call__(self, f_calc, compute_derivatives):
    assert f_calc.is_similar_symmetry(self.f_obs())
    return ext.targets_correlation(
      obs_type="I",
      obs=flex.pow2(self.f_obs().data()),
      weights=self.weights(),
      r_free_flags=None,
      f_calc=f_calc.data(),
      derivatives_depth=int(compute_derivatives))

def registry():
  return {
    "least_squares_residual": least_squares_residual,
    "unified_least_squares_residual": unified_least_squares_residual,
    "intensity_correlation": intensity_correlation,
  }
