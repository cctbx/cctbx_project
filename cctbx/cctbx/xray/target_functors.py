from cctbx.xray import ext
from cctbx.array_family import flex
from libtbx import adopt_init_args

class target_attributes_base(object):

  def __init__(self,
        family,
        specialization=None,
        requires_external_scale=False):
    adopt_init_args(self, locals())
    assert self.validate()

  def validate(self):
    raise RuntimeError("Internal error: missing override.")

class target_attributes(target_attributes_base):

  def validate(self):
    if (self.family == "ls"):
      return self.specialization in [None, "k1", "k2"]
    if (self.family == "ml"):
      return self.specialization in [None, "hl"]
    return False

  def requires_experimental_phases(self):
    return (self.family == "ml" and self.specialization == "hl")

  def ls_apply_scale_to_f_calc(self):
    if (self.family == "ls"):
      return (self.specialization == "k1")
    return None

class manager(object):

  target_attributes_validate = target_attributes.validate

  def __init__(self,
        target_attributes,
        f_obs,
        r_free_flags,
        experimental_phases=None,
        weights=None,
        scale_factor=0):
    assert manager.target_attributes_validate(target_attributes)
    assert r_free_flags.data().size() == f_obs.data().size()
    if (experimental_phases is not None):
      assert experimental_phases.data().size() == f_obs.data().size()
    adopt_init_args(self, locals())
    rffd = r_free_flags.data()
    rffd_true = self.r_free_flags.data().count(True)
    if (rffd_true > 0):
      self.f_obs_w = self.f_obs.select(~rffd)
      self.f_obs_t = self.f_obs.select( rffd)
    else:
      self.f_obs_w = self.f_obs
      self.f_obs_t = self.f_obs # XXX very misleading
    self.experimental_phases_w, self.experimental_phases_t = None, None
    if (target_attributes.requires_experimental_phases()):
      assert experimental_phases is not None
      if (rffd_true > 0):
        self.experimental_phases_w = self.experimental_phases.select(~rffd)
        self.experimental_phases_t = self.experimental_phases.select( rffd)
      else:
        self.experimental_phases_w = self.experimental_phases
        # XXX very misleading
        self.experimental_phases_t = self.experimental_phases
    if (self.weights is not None):
      assert self.target_attributes.family == "ls"
      if (rffd_true):
        self.weights_w = self.weights.select(~rffd)
        self.weights_t = self.weights.select( rffd)
      else:
        self.weights_w = self.weights
        self.weights_t = self.weights # XXX very misleading
    else:
      self.weights_w, self.weights_t = None, None
    self.requires_external_scale = target_attributes.requires_external_scale
    self.ls_apply_scale_to_f_calc =target_attributes.ls_apply_scale_to_f_calc()

  def _target_functor(self, f_obs, weights, experimental_phases, selection):
    if (selection is not None):
      assert selection.size() == f_obs.data().size()
      f_obs = f_obs.select(selection)
      if (weights is not None):
        weights = weights.select(selection)
      if (experimental_phases is not None):
        experimental_phases = experimental_phases.select(selection)
    if (self.ls_apply_scale_to_f_calc is not None):
       if (self.requires_external_scale):
         assert self.scale_factor > 0
       return _ls_functor(
         apply_scale_to_f_calc=self.ls_apply_scale_to_f_calc,
         f_obs=f_obs,
         weights=weights,
         scale_factor=self.scale_factor)
    return _ml_functor(
      f_obs=f_obs,
      experimental_phases=experimental_phases)

  def target_functor_w(self, selection=None):
    return self._target_functor(
      f_obs=self.f_obs_w,
      weights=self.weights_w,
      experimental_phases=self.experimental_phases_w,
      selection=selection)

  def target_functor_t(self, selection=None):
    return self._target_functor(
      f_obs=self.f_obs_t,
      weights=self.weights_t,
      experimental_phases=self.experimental_phases_t,
      selection=selection)

class _ls_functor(object):

  def __init__(self, apply_scale_to_f_calc, f_obs, weights, scale_factor):
    assert scale_factor >= 0
    adopt_init_args(self, locals())

  def __call__(self, f_calc, compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(self.f_obs.unit_cell())
    assert f_calc.space_group() == self.f_obs.space_group()
    return ext.ls_with_scale(
      apply_scale_to_f_calc=self.apply_scale_to_f_calc,
      f_obs=self.f_obs.data(),
      weights=self.weights,
      f_calc=f_calc.data(),
      compute_derivatives=compute_derivatives,
      scale_factor=self.scale_factor)

class _ml_functor(object):

  def __init__(self,
        f_obs,
        experimental_phases=None,
        step_for_integration=5.0):
    adopt_init_args(self, locals())
    self.epsilons = f_obs.epsilons().data()
    self.centric_flags = f_obs.centric_flags().data()

  def __call__(self,
        f_calc,
        alpha,
        beta,
        k,
        compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(self.f_obs.unit_cell())
    assert f_calc.space_group() == self.f_obs.space_group()
    if (self.experimental_phases is None):
      return ext.targets_maximum_likelihood_criterion(
        self.f_obs.data(),
        f_calc.data(),
        alpha,
        beta,
        k,
        self.epsilons,
        self.centric_flags,
        compute_derivatives)
    return ext.targets_maximum_likelihood_criterion_hl(
      self.f_obs.data(),
      f_calc.data(),
      alpha,
      beta,
      self.epsilons,
      self.centric_flags,
      compute_derivatives,
      self.experimental_phases.data(),
      self.step_for_integration)

class least_squares_residual(object):
  """ A least-square residual functor. """

  def __init__(self, f_obs,
                     weights               = None,
                     use_sigmas_as_weights = False,
                     scale_factor          = 0):
    """
    Construct a least-square residuals

    S{sum} w[i] ( f_obs.data[i] - k abs(f_calc.data[i]) )^2
    / S{sum} w[i] f_obs.data[i]^2

    where
      - the sums are over the indices i of the reflections,
      - f_calc is to be passed to the __call__ method,
      - the weights w and the scale factor k are discussed below.

    @type f_obs: real miller.array
    @param f_obs: the observed reflections, with F and sigma(F)
    respectively in f_obs.data() and f_obs.sigmas()
    @type weights: flex.double
    @param weights: the weights w or None in which case w = 1
    @type use_sigmas_as_weights: bool
    @param use_sigmas_as_weights: whether to use w = 1/f_obs.sigmas()^2
    @type scale_factor: number
    @param scale_factor: the scale factor k is not null, otherwise k will
    be computed as a by-product by the __call__ method
    """
    adopt_init_args(self, locals(), hide=True)
    assert self._weights is None or self._use_sigmas_as_weights == False
    if (self._use_sigmas_as_weights):
      sigmas_squared = flex.pow2(self._f_obs.sigmas())
      assert sigmas_squared.all_gt(0)
      self._weights = 1 / sigmas_squared

  def f_obs(self):
    """ The f_obs passed to the constructor """
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
