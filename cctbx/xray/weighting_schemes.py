from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import math


class amplitude_unit_weighting(object):
  """ Mere unit weights for F """

  def compute(self):
    assert(self.observed.is_xray_amplitude_array())
    self.weights = None
    self.derivatives_wrt_f_c = None


class intensity_quasi_unit_weighting(object):
  """ Quasi-unit  weights 1/(4 Fo^2) for F^2.

      The weights are replaced by 1/(4 sigma(Fo^2)^2)
      for weak Fo^2, which are defined as Fo^2 < n_sigma * sigma(Fo^2)
  """

  def __init__(self, n_sigma=1.0):
    self.n_sigma = n_sigma

  def compute(self):
    assert(self.observed.is_xray_intensity_array())
    f_sqr = self.observed.data()
    sig_f_sqr = self.observed.sigmas()
    w = flex.double(f_sqr.size())
    if sig_f_sqr is None:
      strongs = flex.bool(f_sqr.size(), True)
    else:
      if self.n_sigma == 1:
        strongs = f_sqr > sig_f_sqr
      else:
        strongs = f_sqr > self.n_sigma * sig_f_sqr
    weaks = ~strongs
    w.set_selected(strongs, 0.25/f_sqr.select(strongs))
    if sig_f_sqr is not None:
      w.set_selected(weaks,   0.25/flex.pow2( sig_f_sqr.select(weaks) ))
    self.weights = w
    self.derivatives_wrt_f_c = None


class pure_statistical_weighting(object):
  """ 1/sigma^2 weights """

  def compute(self):
    assert(self.observed is not None)
    sigmas_squared = flex.pow2(self.observed.sigmas())
    assert sigmas_squared.all_gt(0)
    self.weights = 1 / sigmas_squared
    self.derivatives_wrt_f_c = None


class simple_shelx_weighting(object):
  """ As the WGHT instruction in ShelXL with only a and b terms.
      The implementation is straightforward without any attempt to save
      space or time: mostly a reference for tests.
  """

  def __init__(self, a=0.1, b=0):
    self._params = a,b

  def compute(self, f_calc, scale_factor=None):
    self.calculated = f_calc
    assert(self.observed.is_xray_intensity_array())
    assert(self.calculated.is_complex_array())
    a,b = self._params
    f_c = self.calculated.data()
    if scale_factor is None:
      scale_factor = self.observed.scale_factor(
        self.calculated, cutoff_factor=0.99)
    self.scale_factor = scale_factor
    f_c = f_c * math.sqrt(scale_factor) # don't modify f_c in place
    sigmas_square = flex.pow2(self.observed.sigmas())
    f_obs_square_plus = self.observed.data().deep_copy()
    negatives = self.observed.data() < 0
    f_obs_square_plus.set_selected(negatives, 0)
    p = (f_obs_square_plus + 2*flex.norm(f_c))/3
    w = 1/(sigmas_square + flex.pow2(a*p) + b*p)
    dw_dfc = -(2*a*a*p + b) * flex.pow2(w) * (4./3*f_c)
    self.weights, self.derivatives_wrt_f_c = w, dw_dfc


class shelx_weighting(object):
  """ As the WGHT instruction in ShelXL """

  def __init__(self, a=0.1, b=0, c=0, d=0, e=0, f=1./3, wavelength=None):
    assert(e == 0 or wavelength is not None)
    self._params = a,b,c,d,e,f
    self._wavelength = wavelength
    self._obs_part_dirty = True
    self.computing_derivatives_wrt_f_c = False

  def _propdef():
    def fget(self):
      try: return self._observed
      except Exception: return None
    def fset(self, obs):
      assert(obs.is_xray_intensity_array())
      self._observed = obs
      self._obs_part_dirty = True
    return locals()
  observed = property(**_propdef())

  def _propdef():
    def fget(self):
      try: return self._calculated
      except Exception: return None
    def fset(self, f_calc):
      assert(f_calc.is_complex_array())
      self._calculated = f_calc
    return locals()
  calculated = property(**_propdef())

  def compute(self, f_calc, scale_factor=None):
    self.calculated = f_calc
    a,b,c,d,e,f = self._params
    if self._obs_part_dirty:
      # The part depending only on |F_o|^2
      if c == 0:
        q = None
      else:
        exp_args = self.observed.sin_theta_over_lambda_sq().data()
        if e != 0: k_sqr = exp_args.deep_copy()
        exp_args *= c
        exp_vals = flex.exp(exp_args)
        if c > 0:
          q = exp_vals
        else:
          q = 1 - exp_vals
      self._q = q
      if self.observed.sigmas() is not None:
        self._den_obs = flex.pow2(self.observed.sigmas())
        if d != 0:
          self._den_obs += d
        if e != 0:
          e_times_sin_theta_sq = k_sqr
          e_times_sin_theta_sq *= e * self._wavelength
          self._den_obs += e_times_sin_theta_sq
      else:
        self._den_obs = None
      negatives = self.observed.data() < 0
      self._p_obs = self.observed.data().deep_copy()
      self._p_obs.set_selected(negatives, 0)
      self._p_obs *= f
      #
      self._obs_part_dirty = False
    # The part depending on |F_c|^2 as well
    q = self._q
    f_c = self.calculated.data()
    p = flex.norm(f_c)
    if scale_factor is None:
      scale_factor = self.observed.scale_factor(
        self.calculated, cutoff_factor=0.99)
    self.scale_factor = scale_factor
    p *= scale_factor
    p *= 1 - f
    p += self._p_obs
    den = p.deep_copy()
    den *= a * a
    der = None
    if self.computing_derivatives_wrt_f_c: der = 2*den
    den += b
    if self.computing_derivatives_wrt_f_c: der += b
    den *= p
    if self._den_obs is not None: den += self._den_obs
    if q is None:
      w = 1 / den
    else:
      w =  q / den
    if self.computing_derivatives_wrt_f_c:
      if scale_factor is not None:
        # don't modify f_c in place
        f_c = f_c * math.sqrt(scale_factor)
      der *= -flex.pow2(w)
      der *= 4./3
      der = der * f_c
    self.weights = w
    self.derivatives_wrt_f_c = der
