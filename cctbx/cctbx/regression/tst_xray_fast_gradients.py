from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import xray
from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import adptbx
from cctbx import matrix
from cctbx.array_family import flex
from cctbx.regression.tst_xray_derivatives import linear_regression_test
from cctbx.regression.tst_sampled_model_density import assign_custom_gaussians
from scitbx.python_utils.misc import adopt_init_args
from libtbx.test_utils import approx_equal
from scitbx import fftpack
import random
import sys

def randomize_gradient_flags(gradient_flags, anomalous_flag,
                             thresholds=(2/3., 1/3.)):
  r = random.random()
  if (r >= thresholds[0]):
    gradient_flags = xray.structure_factors.gradient_flags(default=0001)
  elif (r >= thresholds[1]):
    gradient_flags = gradient_flags.copy()
    if (random.random() > 0.5): gradient_flags.site = 0001
    if (random.random() > 0.5): gradient_flags.u_iso = 0001
    if (random.random() > 0.5): gradient_flags.u_aniso = 0001
    if (random.random() > 0.5): gradient_flags.occupancy = 0001
    if (random.random() > 0.5): gradient_flags.fp = 0001
    if (anomalous_flag):
      if (random.random() > 0.5): gradient_flags.fdp = 0001
  return gradient_flags

class resampling(crystal.symmetry):

  def __init__(self, miller_set=None,
                     crystal_symmetry=None,
                     d_min=None,
                     grid_resolution_factor=1/3.,
                     symmetry_flags=maptbx.use_space_group_symmetry,
                     mandatory_grid_factors=None,
                     quality_factor=100000, u_extra=None, b_extra=None,
                     wing_cutoff=1.e-10,
                     exp_table_one_over_step_size=-100,
                     max_prime=5):
    assert miller_set is None or crystal_symmetry is None
    assert [quality_factor, u_extra, b_extra].count(None) == 2
    if (miller_set is None):
      assert crystal_symmetry is not None and d_min is not None
    else:
      crystal_symmetry = miller_set
      if (d_min is None):
        d_min = miller_set.d_min()
      else:
        assert d_min <= miller_set.d_min()
    crystal.symmetry._copy_constructor(self, crystal_symmetry)
    quality_factor = xray.structure_factors.quality_factor_from_any(
      d_min, grid_resolution_factor, quality_factor, u_extra, b_extra)
    del miller_set
    del u_extra
    del b_extra
    adopt_init_args(self, locals(), hide=0001)
    self._crystal_gridding = None
    self._crystal_gridding_tags = None
    self._rfft = None
    self._u_extra = None

  def d_min(self):
    return self._d_min

  def grid_resolution_factor(self):
    return self._grid_resolution_factor

  def symmetry_flags(self):
    return self._symmetry_flags

  def mandatory_grid_factors(self):
    return self._mandatory_grid_factors

  def quality_factor(self):
    return self._quality_factor

  def wing_cutoff(self):
    return self._wing_cutoff

  def exp_table_one_over_step_size(self):
    return self._exp_table_one_over_step_size

  def max_prime(self):
    return self._max_prime

  def crystal_gridding(self, assert_shannon_sampling=0001):
    if (self._crystal_gridding is None):
      self._crystal_gridding = maptbx.crystal_gridding(
        unit_cell=self.unit_cell(),
        d_min=self.d_min(),
        resolution_factor=self.grid_resolution_factor(),
        symmetry_flags=self.symmetry_flags(),
        space_group_info=self.space_group_info(),
        mandatory_factors=self.mandatory_grid_factors(),
        max_prime=self.max_prime(),
        assert_shannon_sampling=assert_shannon_sampling)
    return self._crystal_gridding

  def crystal_gridding_tags(self, assert_shannon_sampling=0001):
    if (self._crystal_gridding_tags is None):
      self._crystal_gridding_tags = self.crystal_gridding(
        assert_shannon_sampling).tags()
    return self._crystal_gridding_tags

  def rfft(self):
    if (self._rfft is None):
      self._rfft = fftpack.real_to_complex_3d(self.crystal_gridding().n_real())
    return self._rfft

  def u_extra(self):
    if (self._u_extra is None):
      self._u_extra = xray.calc_u_extra(
        self.d_min(),
        self.grid_resolution_factor(),
        self.quality_factor())
    return self._u_extra

  def setup_fft(self):
    self.crystal_gridding_tags()
    self.rfft()
    self.u_extra()
    return self

  def ft_dp(self, dp):
    multiplier = (  self.unit_cell().volume()
                  / matrix.row(self.rfft().n_real()).product()
                  * self.space_group().order_z()
                  / dp.multiplicities().data().as_double())
    coeff = dp.deep_copy()
    xray.apply_u_extra(
      self.unit_cell(),
      self.u_extra(),
      coeff.indices(),
      coeff.data(),
      multiplier)
    return miller.fft_map(
      crystal_gridding=self.crystal_gridding(),
      fourier_coefficients=coeff)

  def __call__(self, xray_structure,
                     dp,
                     gradient_flags,
                     n_parameters,
                     electron_density_must_be_positive=0001,
                     tolerance_positive_definite=1.e-5,
                     verbose=0):
    gradient_map = self.ft_dp(dp)
    if (not gradient_map.anomalous_flag()):
      gradient_map_real = gradient_map.real_map()
      gradient_map_complex = flex.complex_double(flex.grid(0,0,0))
    else:
      gradient_map_real = flex.double(flex.grid(0,0,0))
      gradient_map_complex = gradient_map.complex_map()
      assert not gradient_map_complex.is_padded()
      if (0 or verbose):
        gradient_flags.show_summary()
        print "grid:", gradient_map_complex.focus()
        print "ft_dt_map real: %.4g %.4g" % (
          flex.min(flex.real(gradient_map_complex)),
          flex.max(flex.real(gradient_map_complex)))
        print "ft_dt_map imag: %.4g %.4g" % (
          flex.min(flex.imag(gradient_map_complex)),
          flex.max(flex.imag(gradient_map_complex)))
        print
    result = xray.fast_gradients(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      xray_structure.scattering_dict(),
      gradient_map_real,
      gradient_map_complex,
      gradient_flags,
      n_parameters,
      self.u_extra(),
      self.wing_cutoff(),
      self.exp_table_one_over_step_size(),
      electron_density_must_be_positive,
      tolerance_positive_definite)
    if (0 or verbose):
      print "max_shell_radii:", result.max_shell_radii()
      print "exp_table_size:", result.exp_table_size()
      print
    return result

class judge:

  def __init__(self, scatterer, label, reference, other, top):
    label += [" iso", " aniso"][int(scatterer.anisotropic_flag)]
    s = ""
    r = (reference-other)/top
    s += " %.5f " % r + label
    self.is_bad = 00000
    if (abs(r) > 0.03):
      s += " very large mismatch"
      self.is_bad = 0001
    elif (abs(r) > 0.01):
      s += " large mismatch"
    self.s = s.lstrip()

  def __str__(self):
    return self.s

class shifted_site:

  def __init__(self, f_obs, structure, i_scatterer, i_xyz, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    site = list(self.structure_shifted.scatterers()[i_scatterer].site)
    site[i_xyz] += shift
    self.structure_shifted.scatterers()[i_scatterer].site = site
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted, algorithm="direct").f_calc()

def site(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_site(f_obs, structure_ideal, 0, 0, 0.01)
  if (0 or verbose):
    print "site"
    sh.structure_shifted.show_summary().show_scatterers()
    print
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(site=0001),
    f_obs.anomalous_flag())
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    gradient_flags=gradient_flags,
    n_parameters=0,
    verbose=verbose)
  sfd_d_target_d_site_cart = sfd.d_target_d_site_cart()
  top_gradient = None
  for i_scatterer in sh.structure_shifted.scatterers().indices():
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    for i_xyz in (0,1,2):
      direct_summ = sfd_d_target_d_site_cart[i_scatterer][i_xyz]
      if (top_gradient is None): top_gradient = direct_summ
      fast_gradie = map0.d_target_d_site_cart()[i_scatterer][i_xyz]
      match = judge(scatterer, "site", direct_summ, fast_gradie, top_gradient)
      if (0 or verbose):
        print "direct summ[%d][%d]: " % (i_scatterer,i_xyz), direct_summ
        print "fast gradie[%d][%d]: " % (i_scatterer,i_xyz), fast_gradie, match
        print
      assert not match.is_bad
  sys.stdout.flush()

class shifted_u_iso:

  def __init__(self, f_obs, structure, i_scatterer, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    self.structure_shifted.scatterers()[i_scatterer].u_iso += shift
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted).f_calc()

def u_iso(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_u_iso(f_obs, structure_ideal, 0, 0.05)
  if (0 or verbose):
    print "u_iso"
    sh.structure_shifted.show_summary().show_scatterers()
    print
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(u_iso=0001),
    f_obs.anomalous_flag())
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    gradient_flags=gradient_flags,
    n_parameters=0,
    verbose=verbose)
  top_gradient = None
  for i_scatterer in sh.structure_shifted.scatterers().indices():
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    direct_summ = sfd.d_target_d_u_iso()[i_scatterer]
    if (top_gradient is None): top_gradient = direct_summ
    fast_gradie = map0.d_target_d_u_iso()[i_scatterer]
    match = judge(scatterer, "u_iso", direct_summ, fast_gradie, top_gradient)
    if (0 or verbose):
      print "direct summ[%d]: " % i_scatterer, direct_summ
      print "fast gradie[%d]: " % i_scatterer, fast_gradie, match
      print
    assert not match.is_bad
  sys.stdout.flush()

class shifted_u_star:

  def __init__(self, f_obs, structure, i_scatterer, ij, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    scatterer = self.structure_shifted.scatterers()[i_scatterer]
    u_star = list(scatterer.u_star)
    u_star[ij] += shift
    scatterer.u_star = u_star
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted).f_calc()

def u_star(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_u_star(f_obs, structure_ideal, 0, 0, 0.0001)
  if (0 or verbose):
    print "u_star"
    sh.structure_shifted.show_summary().show_scatterers()
    print
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(u_aniso=0001),
    f_obs.anomalous_flag())
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    gradient_flags=gradient_flags,
    n_parameters=0,
    verbose=verbose)
  sfd_d_target_d_u_cart = sfd.d_target_d_u_cart()
  map0_d_target_d_u_cart = map0.d_target_d_u_cart()
  top_gradient = None
  for i_scatterer in sh.structure_shifted.scatterers().indices():
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    sfd_star = sfd.d_target_d_u_star()[i_scatterer]
    sfd_cart = adptbx.grad_u_star_as_u_cart(
      structure_ideal.unit_cell(), sfd_star)
    assert approx_equal(
      sfd_star,
      adptbx.grad_u_cart_as_u_star(structure_ideal.unit_cell(), sfd_cart))
    for ij in xrange(6):
      direct_summ = sfd_d_target_d_u_cart[i_scatterer][ij]
      if (top_gradient is None): top_gradient = direct_summ
      fast_gradie = map0_d_target_d_u_cart[i_scatterer][ij]
      match = judge(scatterer, "u_cart", direct_summ,fast_gradie,top_gradient)
      if (0 or verbose):
        print "direct summ[%d][%d]: " % (i_scatterer, ij), direct_summ
        print "fast gradie[%d][%d]: " % (i_scatterer, ij), fast_gradie, match
        print
      assert not match.is_bad
  sys.stdout.flush()

class shifted_occupancy:

  def __init__(self, f_obs, structure, i_scatterer, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    self.structure_shifted.shift_occupancy(i_scatterer, shift)
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted).f_calc()

def occupancy(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_occupancy(f_obs, structure_ideal, 0, 0.2)
  if (0 or verbose):
    print "occupancy"
    sh.structure_shifted.show_summary().show_scatterers()
    print
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(occupancy=0001),
    f_obs.anomalous_flag())
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    gradient_flags=gradient_flags,
    n_parameters=0,
    verbose=verbose)
  top_gradient = None
  for i_scatterer in sh.structure_shifted.scatterers().indices():
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    direct_summ = sfd.d_target_d_occupancy()[i_scatterer]
    if (top_gradient is None): top_gradient = direct_summ
    fast_gradie = map0.d_target_d_occupancy()[i_scatterer]
    match = judge(scatterer, "occupancy", direct_summ,fast_gradie,top_gradient)
    if (0 or verbose):
      print "direct summ[%d]: " % i_scatterer, direct_summ
      print "fast gradie[%d]: " % i_scatterer, fast_gradie, match
      print
    assert not match.is_bad
  sys.stdout.flush()

class shifted_fp:

  def __init__(self, f_obs, structure, i_scatterer, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    self.structure_shifted.scatterers()[i_scatterer].fp += shift
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted).f_calc()

def fp(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_fp(f_obs, structure_ideal, 0, -0.2)
  if (0 or verbose):
    print "fp"
    sh.structure_shifted.show_summary().show_scatterers()
    print
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(fp=0001),
    f_obs.anomalous_flag())
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    gradient_flags=gradient_flags,
    n_parameters=0,
    verbose=verbose)
  top_gradient = None
  for i_scatterer in sh.structure_shifted.scatterers().indices():
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    direct_summ = sfd.d_target_d_fp()[i_scatterer]
    if (top_gradient is None): top_gradient = direct_summ
    fast_gradie = map0.d_target_d_fp()[i_scatterer]
    match = judge(scatterer, "fp", direct_summ, fast_gradie, top_gradient)
    if (0 or verbose):
      print "direct summ[%d]: " % i_scatterer, direct_summ
      print "fast gradie[%d]: " % i_scatterer, fast_gradie, match
      print
    assert not match.is_bad
  sys.stdout.flush()

class shifted_fdp:

  def __init__(self, f_obs, structure, i_scatterer, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    self.structure_shifted.scatterers()[i_scatterer].fdp += shift
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted).f_calc()

def fdp(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_fdp(f_obs, structure_ideal, 0, 2)
  if (0 or verbose):
    print "fdp"
    sh.structure_shifted.show_summary().show_scatterers()
    print
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(fdp=0001),
    f_obs.anomalous_flag())
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    gradient_flags=gradient_flags,
    n_parameters=0,
    verbose=verbose)
  top_gradient = None
  for i_scatterer in sh.structure_shifted.scatterers().indices():
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    direct_summ = sfd.d_target_d_fdp()[i_scatterer]
    if (top_gradient is None): top_gradient = direct_summ
    fast_gradie = map0.d_target_d_fdp()[i_scatterer]
    match = judge(scatterer, "fdp", direct_summ, fast_gradie, top_gradient)
    if (0 or verbose):
      print "direct summ[%d]: " % i_scatterer, direct_summ
      print "fast gradie[%d]: " % i_scatterer, fast_gradie, match
      print
    assert not match.is_bad
  sys.stdout.flush()

def shift_all(structure_ideal, f_obs, anomalous_flag, anisotropic_flag):
  sh = shifted_site(None, structure_ideal, 0, 0, 0.01)
  if (not anisotropic_flag):
    sh = shifted_u_iso(None, sh.structure_shifted, 0, 0.05)
  else:
    sh = shifted_u_star(None, sh.structure_shifted, 0, 0, 0.0001)
  sh = shifted_occupancy(None, sh.structure_shifted, 0, 0.2)
  if (anomalous_flag):
    sh = shifted_fdp(None, sh.structure_shifted, 0, 2)
  sh = shifted_fp(f_obs, sh.structure_shifted, 0, -0.2)
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  return sh, ls

def exercise_packed(structure_ideal, f_obs,
                    anomalous_flag, anisotropic_flag,
                    verbose=0):
  sh, ls = shift_all(structure_ideal, f_obs, anomalous_flag, anisotropic_flag)
  flag = (random.random() > 0.5)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(site=flag, u=not flag),
    f_obs.anomalous_flag(),
    thresholds=(1/2.,0))
  n_parameters = structure_ideal.n_parameters(gradient_flags)
  assert n_parameters > 0
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=n_parameters)
  assert sfd.packed().size() == n_parameters
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    gradient_flags=gradient_flags,
    n_parameters=n_parameters,
    verbose=verbose)
  assert map0.packed().size() == n_parameters
  correlation = flex.linear_correlation(sfd.packed(), map0.packed())
  assert correlation.is_well_defined()
  assert correlation.coefficient() > 0.999

def exercise_gradient_manager(structure_ideal, f_obs,
                              anomalous_flag, anisotropic_flag,
                              verbose=0):
  sh, ls = shift_all(structure_ideal, f_obs, anomalous_flag, anisotropic_flag)
  grad_manager = xray.structure_factors.gradients(
    miller_set=f_obs,
    quality_factor=100000,
    wing_cutoff=1.e-10)
  gradient_flags=xray.structure_factors.gradient_flags(default=0001)
  if (random.random() > 0.5):
    n_parameters = 0
  else:
    n_parameters = structure_ideal.n_parameters(gradient_flags)
  gd = grad_manager(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=n_parameters,
    algorithm="direct")
  gf = grad_manager(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    gradient_flags=gradient_flags,
    n_parameters=n_parameters,
    algorithm="fft")
  if (n_parameters == 0):
    d = gd.d_target_d_site_frac()
    f = gf.d_target_d_site_frac()
    linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    d = gd.d_target_d_site_cart()
    f = gf.d_target_d_site_cart()
    linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    d = gd.d_target_d_u_iso()
    f = gf.d_target_d_u_iso()
    linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    d = gd.d_target_d_u_cart()
    f = gf.d_target_d_u_cart()
    linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    d = gd.d_target_d_occupancy()
    f = gf.d_target_d_occupancy()
    linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    d = gd.d_target_d_fp()
    f = gf.d_target_d_fp()
    linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    if (anomalous_flag):
      d = gd.d_target_d_fdp()
      f = gf.d_target_d_fdp()
      linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
  else:
    correlation = flex.linear_correlation(gd.packed(), gf.packed())
    assert correlation.is_well_defined()
    assert correlation.coefficient() > 0.999

def run_one(space_group_info, n_elements=3, volume_per_atom=1000, d_min=2,
            anomalous_flag=0, anisotropic_flag=0, verbose=0):
  if (random.random() < 0.5):
    random_f_prime_scale=0.6
  else:
    random_f_prime_scale=0
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=(("C","N","O")*(n_elements/3+1))[:n_elements],
    volume_per_atom=volume_per_atom,
    min_distance=5,
    general_positions_only=1,
    random_f_prime_d_min=d_min-1,
    random_f_prime_scale=random_f_prime_scale,
    random_f_double_prime=anomalous_flag,
    anisotropic_flag=anisotropic_flag,
    random_u_iso=0001,
    random_u_iso_scale=.3,
    random_u_cart_scale=.3,
    random_occupancy=0001)
  if (random.random() < 0.5):
    assign_custom_gaussians(structure_ideal)
  if (0 or verbose):
    structure_ideal.show_summary().show_scatterers()
    if (anisotropic_flag):
      uc = structure_ideal.unit_cell()
      for scatterer in structure_ideal.scatterers():
        print "u_iso:", adptbx.u_star_as_u_iso(uc, scatterer.u_star)
    print
  f_obs = abs(structure_ideal.structure_factors(
    d_min=d_min, anomalous_flag=anomalous_flag, algorithm="direct").f_calc())
  if (1):
    site(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1):
    if (not anisotropic_flag):
      u_iso(structure_ideal, d_min, f_obs, verbose=verbose)
    else:
      u_star(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1):
    occupancy(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1):
    fp(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1 and anomalous_flag):
    fdp(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1):
    exercise_gradient_manager(
      structure_ideal, f_obs, anomalous_flag, anisotropic_flag)
  if (1):
    exercise_packed(
      structure_ideal, f_obs, anomalous_flag, anisotropic_flag)

def run_call_back(flags, space_group_info):
  for anomalous_flag in [0,1]:
    for anisotropic_flag in [0,1]:
      run_one(
        space_group_info=space_group_info,
        anomalous_flag=anomalous_flag,
        anisotropic_flag=anisotropic_flag,
        verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
