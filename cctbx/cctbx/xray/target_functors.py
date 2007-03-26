from cctbx.xray import ext
from cctbx.array_family import flex
from libtbx import adopt_init_args

class least_squares(object):

  def __init__(self,
        apply_scale_to_f_calc,
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
    return ext.targets_ls_with_scale(
      apply_scale_to_f_calc=self.apply_scale_to_f_calc,
      compute_scale_using_all_data=self.compute_scale_using_all_data,
      f_obs=self.f_obs.data(),
      weights=self.weights,
      r_free_flags=self.r_free_flags.data(),
      f_calc=f_calc.data(),
      compute_gradients=compute_gradients,
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

  def __init__(self, f_obs                 = None,
                     f_obs_square          = None,
                     weights               = None,
                     use_sigmas_as_weights = False,
                     scale_factor          = 0):
    """
    Construct a least-square residuals

    S{sum} w[i] ( f_obs.data[i] - k abs(f_calc.data[i]) )^2
    / S{sum} w[i] f_obs.data[i]^2

    or

    S{sum} w[i] ( f_obs_square.data[i] - k abs(f_calc.data[i])^2 )^2
    / S{sum} w[i] f_obs_square.data[i]^2

    depending on which of f_obs and f_obs_square is not None.

    Note that
      - the sums are over the indices i of the reflections,
      - f_calc is to be passed to the __call__ method,
      - the weights w and the scale factor k are discussed below.

    @type f_obs: real miller.array
    @param f_obs: the observed reflections, with F and sigma(F)
    respectively in f_obs.data() and f_obs.sigmas() or None
    @type f_obs_square: real miller.array
    @param f_obs_square: the observed reflections, with F^2 and sigma(F^2)
    respectively in f_obs_square.data() and f_obs_square.sigmas()
    @type weights: flex.double
    @param weights: the weights w or None in which case w = 1
    @type use_sigmas_as_weights: bool
    @param use_sigmas_as_weights: whether to use w = 1/sigmas()^2
    @type scale_factor: number
    @param scale_factor: the scale factor k if not null, otherwise k will
    be computed as a by-product by the __call__ method
    """
    adopt_init_args(self, locals(), hide=True)
    assert self._f_obs is not None or self._f_obs_square is not None
    assert self._weights is None or self._use_sigmas_as_weights == False
    if (self._use_sigmas_as_weights):
      sigmas_squared = flex.pow2(self._f_obs.sigmas())
      assert sigmas_squared.all_gt(0)
      self._weights = 1 / sigmas_squared

  def f_obs(self):
    """ The f_obs passed to the constructor """
    return self._f_obs

  def f_obs_square(self):
    """ The f_obs_square passed to the constructor """
    return self._f_obs

  def weights(self):
    """ The weights w """
    return self._weights

  def use_sigmas_as_weights(self):
    """ The flag with the same name passed to the constructor """
    return self._use_sigmas_as_weights

  def __call__(self, f_calc, compute_derivatives):
    """
    Compute the least-squares residual value and perhaps its derivatives
    wrt to the calculated structure factor F_c of the i-th reflection
    @type f_calc: complex miller.array
    @param f_calc: f_calc.data()[i] constains F_c for the i-th reflection
    in f_obs()
    @type compute_derivatives: bool
    @param compute_derivatives: whether to compute the derivatives of the
    least square residual or not
    @rtype: Boost.Python binding of
    U{least_squares_residual<DOXYCLASS:cctbx::xray::targets::least_squares_residual>}
    @return: An object holding the residual value, derivatives and scale
    factor
    """
    if(self.f_obs() is not None):
      ext_ls_residual = ext.targets_least_squares_residual
      y_obs = self.f_obs()
    else:
      ext_ls_residual = ext.targets_least_squares_residual_for_F_square
      y_obs = self.f_obs_square()
    assert f_calc.unit_cell().is_similar_to(
           y_obs.unit_cell())
    assert f_calc.space_group() == y_obs.space_group()
    if (self.weights() is not None):
      return ext_ls_residual(
        y_obs().data(),
        self.weights(),
        f_calc.data(),
        compute_derivatives,
        self._scale_factor)
    else:
      return ext_ls_residual(
        y_obs.data(),
        f_calc.data(),
        compute_derivatives,
        self._scale_factor)


class least_squares_residual(unified_least_squares_residual):
  """ A least-square residual functor for F refinement. """

  def __init__(self, f_obs,
                     weights               = None,
                     use_sigmas_as_weights = False,
                     scale_factor          = 0):
    super(least_squares_residual, self).__init__(f_obs,
                                                 None,
                                                 weights,
                                                 use_sigmas_as_weights,
                                                 scale_factor)


class intensity_correlation(object):

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

def registry():
  return {
    "least_squares_residual": least_squares_residual,
    "intensity_correlation": intensity_correlation,
  }
