from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import xray
import cctbx.xray.structure_factors.global_counters
from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import adptbx
from cctbx.regression.tst_xray_derivatives import linear_regression_test
from cctbx.regression.tst_sampled_model_density import assign_custom_gaussians
from scitbx import fftpack
from scitbx import matrix
import omptbx
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
import libtbx.utils
import libtbx.introspection
import random
import sys, math
from six.moves import range
from six.moves import zip

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def randomize_gradient_flags(gradient_flags, anomalous_flag,
                             thresholds=(2/3., 1/3.)):
  r = random.random()
  if (r >= thresholds[0]):
    gradient_flags = xray.structure_factors.gradient_flags(default=True)
  elif (r >= thresholds[1]):
    gradient_flags = gradient_flags.copy()
    if (random.random() > 0.5): gradient_flags.site = True
    if (random.random() > 0.5): gradient_flags.u_iso = True
    if (random.random() > 0.5): gradient_flags.u_aniso = True
    if (random.random() > 0.5): gradient_flags.occupancy = True
    if (random.random() > 0.5): gradient_flags.fp = True
    if (anomalous_flag):
      if (random.random() > 0.5): gradient_flags.fdp = True
  return gradient_flags

class resampling(crystal.symmetry):

  def __init__(self, miller_set=None,
                     crystal_symmetry=None,
                     d_min=None,
                     grid_resolution_factor=1/3.,
                     symmetry_flags=maptbx.use_space_group_symmetry,
                     mandatory_grid_factors=None,
                     quality_factor=1000000, u_base=None, b_base=None,
                     wing_cutoff=1.e-10,
                     exp_table_one_over_step_size=-100,
                     max_prime=5):
    assert miller_set is None or crystal_symmetry is None
    assert [quality_factor, u_base, b_base].count(None) == 2
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
      d_min, grid_resolution_factor, quality_factor, u_base, b_base)
    del miller_set
    del u_base
    del b_base
    adopt_init_args(self, locals(), hide=True)
    self._crystal_gridding = None
    self._crystal_gridding_tags = None
    self._rfft = None
    self._u_base = None

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

  def crystal_gridding(self, assert_shannon_sampling=True):
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

  def crystal_gridding_tags(self, assert_shannon_sampling=True):
    if (self._crystal_gridding_tags is None):
      self._crystal_gridding_tags = self.crystal_gridding(
        assert_shannon_sampling).tags()
    return self._crystal_gridding_tags

  def rfft(self):
    if (self._rfft is None):
      self._rfft = fftpack.real_to_complex_3d(self.crystal_gridding().n_real())
    return self._rfft

  def u_base(self):
    if (self._u_base is None):
      self._u_base = xray.calc_u_base(
        self.d_min(),
        self.grid_resolution_factor(),
        self.quality_factor())
    return self._u_base

  def setup_fft(self):
    self.crystal_gridding_tags()
    self.rfft()
    self.u_base()
    return self

  def ft_dp(self, dp, u_extra):
    multiplier = (  self.unit_cell().volume()
                  / matrix.row(self.rfft().n_real()).product()
                  * self.space_group().order_z()
                  / dp.multiplicities().data().as_double())
    coeff = dp.deep_copy()
    xray.apply_u_extra(
      self.unit_cell(),
      u_extra,
      coeff.indices(),
      coeff.data())
    coeff_data = coeff.data()
    coeff_data *= flex.polar(multiplier, 0)
    return miller.fft_map(
      crystal_gridding=self.crystal_gridding(),
      fourier_coefficients=coeff)

  def __call__(self, xray_structure,
                     u_iso_refinable_params,
                     dp,
                     n_parameters,
                     verbose=0):
    omptbx.env.num_threads = libtbx.introspection.number_of_processors()
    result = xray.fast_gradients(
      unit_cell=xray_structure.unit_cell(),
      scatterers=xray_structure.scatterers(),
      scattering_type_registry=xray_structure.scattering_type_registry(),
      u_base=self.u_base(),
      wing_cutoff=self.wing_cutoff(),
      exp_table_one_over_step_size=self.exp_table_one_over_step_size(),
      tolerance_positive_definite=1.e-5)
    if (0 or verbose):
      print("u_base:", result.u_base())
      print("u_extra:", result.u_extra())
    gradient_map = self.ft_dp(dp, u_extra=result.u_extra())
    if (not gradient_map.anomalous_flag()):
      gradient_map = gradient_map.real_map()
    else:
      gradient_map = gradient_map.complex_map()
      assert not gradient_map.is_padded()
      if (0 or verbose):
        print("grid:", gradient_map.focus())
        print("ft_dt_map real: %.4g %.4g" % (
          flex.min(flex.real(gradient_map)),
          flex.max(flex.real(gradient_map))))
        print("ft_dt_map imag: %.4g %.4g" % (
          flex.min(flex.imag(gradient_map)),
          flex.max(flex.imag(gradient_map))))
        print()
    result.sampling(
      scatterers=xray_structure.scatterers(),
      u_iso_refinable_params=u_iso_refinable_params,
      scattering_type_registry=xray_structure.scattering_type_registry(),
      site_symmetry_table=xray_structure.site_symmetry_table(),
      ft_d_target_d_f_calc=gradient_map,
      n_parameters=n_parameters,
      sampled_density_must_be_positive=False)
    if (0 or verbose):
      print("max_sampling_box_edges:", result.max_sampling_box_edges())
      print("exp_table_size:", result.exp_table_size())
      print()
    return result

class judge(object):

  def __init__(self, scatterer, label, reference, other, top):
    if(scatterer.flags.use_u_iso()):   label += " iso "
    if(scatterer.flags.use_u_aniso()): label += " aniso "
    s = ""
    r = (reference-other)/max(abs(top), min(abs(reference), abs(other)))
    s += " %.5f " % r + label
    self.is_bad = False
    if (abs(r) > 0.03):
      s += " very large mismatch"
      self.is_bad = True
    elif (abs(r) > 0.01):
      s += " large mismatch"
    self.s = s.lstrip()

  def __str__(self):
    return self.s

class shifted_site(object):

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
    print("site")
    sh.structure_shifted.show_summary().show_scatterers()
    print()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), True, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(site=True),
    f_obs.anomalous_flag())
  xray.set_scatterer_grad_flags(scatterers = sh.structure_shifted.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    n_parameters=0,
    verbose=verbose)
  sfd_d_target_d_site_cart = sfd.d_target_d_site_cart()
  top_gradient = None
  for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
    for i_xyz in (0,1,2):
      direct_summ = sfd_d_target_d_site_cart[i_scatterer][i_xyz]
      if (top_gradient is None): top_gradient = direct_summ
      fast_gradie = map0.d_target_d_site_cart()[i_scatterer][i_xyz]
      match = judge(scatterer, "site", direct_summ, fast_gradie, top_gradient)
      if (0 or verbose):
        print("direct summ[%d][%d]: " % (i_scatterer,i_xyz), direct_summ)
        print("fast gradie[%d][%d]: " % (i_scatterer,i_xyz), fast_gradie, match)
        print()
      assert not match.is_bad
  sys.stdout.flush()

class shifted_u_iso(object):

  def __init__(self, f_obs, structure, i_scatterer, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    if (self.structure_shifted.scatterers()[i_scatterer].flags.use_u_iso()):
      self.structure_shifted.scatterers()[i_scatterer].u_iso += shift
      if (f_obs is not None):
        self.f_calc = f_obs.structure_factors_from_scatterers(
          xray_structure=self.structure_shifted).f_calc()

def u_iso(structure_ideal, d_min, f_obs, tan_u_iso, verbose=0):
  sh = shifted_u_iso(f_obs, structure_ideal, 0, 0.05)
  if (0 or verbose):
    print("u_iso")
    sh.structure_shifted.show_summary().show_scatterers()
    print()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), True, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(u_iso=True),
    f_obs.anomalous_flag())
  if(tan_u_iso):
     u_iso_refinable_params = flex.double()
  else:
     u_iso_refinable_params = None
  for scatterer in sh.structure_shifted.scatterers():
      scatterer.flags.set_grad_site(gradient_flags.site)
      scatterer.flags.set_grad_u_iso(gradient_flags.u_iso)
      scatterer.flags.set_grad_u_aniso(gradient_flags.u_aniso)
      scatterer.flags.set_grad_occupancy(gradient_flags.occupancy)
      scatterer.flags.set_grad_fp(gradient_flags.fp)
      scatterer.flags.set_grad_fdp(gradient_flags.fdp)
      if(tan_u_iso):
         scatterer.flags.set_tan_u_iso(True)
         param = random.randint(90,120)
         scatterer.flags.param= param
         value=math.tan(scatterer.u_iso*math.pi/adptbx.b_as_u(param)-math.pi/2)
         u_iso_refinable_params.append(value)
  if (0):
    print("u_iso")
    print("gradient_flags.site      ", gradient_flags.site)
    print("gradient_flags.u_iso     ", gradient_flags.u_iso)
    print("gradient_flags.u_aniso   ", gradient_flags.u_aniso)
    print("gradient_flags.occupancy ", gradient_flags.occupancy)
    print("gradient_flags.fp        ", gradient_flags.fp)
    print("gradient_flags.fdp       ", gradient_flags.fdp)
    cntr_use_u_iso = 0
    cntr_use_u_aniso = 0
    cntr_grad_u_iso = 0
    cntr_grad_u_aniso = 0
    for scatterer in sh.structure_shifted.scatterers():
      if (scatterer.flags.use_u_iso()):  cntr_use_u_iso += 1
      if (scatterer.flags.use_u_aniso()):  cntr_use_u_aniso += 1
      if (scatterer.flags.grad_u_iso()):  cntr_grad_u_iso += 1
      if (scatterer.flags.grad_u_aniso()):  cntr_grad_u_aniso += 1
    print("use_u_iso                ", cntr_use_u_iso,cntr_grad_u_iso)
    print("use_u_aniso              ", cntr_use_u_aniso,cntr_grad_u_aniso)
  grad_flags_counts = \
            xray.scatterer_grad_flags_counts(sh.structure_shifted.scatterers())
  if(grad_flags_counts.n_parameters() > 0):
     sfd = xray.structure_factors.gradients_direct(
       xray_structure=sh.structure_shifted,
       u_iso_refinable_params=u_iso_refinable_params,
       miller_set=f_obs,
       d_target_d_f_calc=ls.derivatives(),
       n_parameters=0)
     re = resampling(miller_set=f_obs)
     map0 = re(
       xray_structure=sh.structure_shifted,
       u_iso_refinable_params=u_iso_refinable_params,
       dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
       n_parameters=0,
       verbose=verbose)
     if(grad_flags_counts.u_aniso > 0):
        sfd_d_target_d_u_cart = sfd.d_target_d_u_cart()
        map0_d_target_d_u_cart = map0.d_target_d_u_cart()
     top_gradient = None
     gradients_1 = []
     for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
       if(0):
         print("i_scatterer= ", i_scatterer,scatterer.flags.use_u_iso(),\
           scatterer.flags.grad_u_iso(), scatterer.flags.use_u_aniso(),\
           scatterer.flags.grad_u_aniso(), scatterer.u_iso, scatterer.u_star)
       if(scatterer.flags.use_u_iso()): parameter_name = "u_iso"
       if(scatterer.flags.use_u_aniso()): parameter_name = "u_star"
       if(parameter_name == "u_iso" and scatterer.flags.grad_u_iso() and
                                                  scatterer.flags.use_u_iso()):
          direct_summ = sfd.d_target_d_u_iso()[i_scatterer]
          if (top_gradient is None): top_gradient = direct_summ
          fast_gradie = map0.d_target_d_u_iso()[i_scatterer]
          sys.stdout.flush()
          gradients_1.append([direct_summ, fast_gradie])
          match = judge(scatterer, parameter_name, direct_summ, fast_gradie,
                                                                  top_gradient)
          if (0 or verbose):
            print("direct summ[%d]: " % i_scatterer, direct_summ)
            print("fast gradie[%d]: " % i_scatterer, fast_gradie, match)
            print()
          assert not match.is_bad
       if(parameter_name == "u_star" and scatterer.flags.grad_u_aniso() and
                                                scatterer.flags.use_u_aniso()):
        sfd_star = sfd.d_target_d_u_star()[i_scatterer]
        sfd_cart = adptbx.grad_u_star_as_u_cart(
          structure_ideal.unit_cell(), sfd_star)
        assert approx_equal(
          sfd_star,
          adptbx.grad_u_cart_as_u_star(structure_ideal.unit_cell(), sfd_cart))
        for ij in range(6):
          direct_summ = sfd_d_target_d_u_cart[i_scatterer][ij]
          if (top_gradient is None): top_gradient = direct_summ
          fast_gradie = map0_d_target_d_u_cart[i_scatterer][ij]
          gradients_1.append([direct_summ, fast_gradie])
          match =judge(scatterer,"u_star",direct_summ,fast_gradie,top_gradient)
          if (0 or verbose or match.is_bad):
            print("direct summ[%d][%d]: " % (i_scatterer, ij), direct_summ)
            print("fast gradie[%d][%d]: " % (i_scatterer, ij),fast_gradie,match)
            print()
          assert not match.is_bad
     # Making sure that gradients_1 = gradients_2
     for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
         if(not scatterer.flags.use_u_iso()):
            scatterer.u_iso = -12345.0
         if(not scatterer.flags.use_u_aniso()):
            scatterer.u_star =(-999.,-999.,-999.,-999.,-999.,-999.)
     sfd = xray.structure_factors.gradients_direct(
       xray_structure=sh.structure_shifted,
       u_iso_refinable_params=u_iso_refinable_params,
       miller_set=f_obs,
       d_target_d_f_calc=ls.derivatives(),
       n_parameters=0)
     re = resampling(miller_set=f_obs)
     map0 = re(
       xray_structure=sh.structure_shifted,
       u_iso_refinable_params=u_iso_refinable_params,
       dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
       n_parameters=0,
       verbose=verbose)
     grad_flags_counts = \
            xray.scatterer_grad_flags_counts(sh.structure_shifted.scatterers())
     if(grad_flags_counts.u_aniso):
        sfd_d_target_d_u_cart = sfd.d_target_d_u_cart()
        map0_d_target_d_u_cart = map0.d_target_d_u_cart()
     gradients_2 = []
     for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
       if(scatterer.flags.use_u_iso()):   parameter_name = "u_iso"
       if(scatterer.flags.use_u_aniso()): parameter_name = "u_star"
       if(parameter_name == "u_iso" and scatterer.flags.grad_u_iso() and
                                                  scatterer.flags.use_u_iso()):
          direct_summ = sfd.d_target_d_u_iso()[i_scatterer]
          fast_gradie = map0.d_target_d_u_iso()[i_scatterer]
          gradients_2.append([direct_summ, fast_gradie])
       if(parameter_name == "u_star" and scatterer.flags.grad_u_aniso() and
                                                scatterer.flags.use_u_aniso()):
        sfd_star = sfd.d_target_d_u_star()[i_scatterer]
        sfd_cart = adptbx.grad_u_star_as_u_cart(structure_ideal.unit_cell(),
                                                                      sfd_star)
        assert approx_equal(
          sfd_star,
          adptbx.grad_u_cart_as_u_star(structure_ideal.unit_cell(), sfd_cart))
        for ij in range(6):
          direct_summ = sfd_d_target_d_u_cart[i_scatterer][ij]
          fast_gradie = map0_d_target_d_u_cart[i_scatterer][ij]
          gradients_2.append([direct_summ, fast_gradie])
     for g1,g2 in zip(gradients_1, gradients_2):
       assert approx_equal(g1, g2)
     sys.stdout.flush()

class shifted_u_star(object):

  def __init__(self, f_obs, structure, i_scatterer, ij, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    if (self.structure_shifted.scatterers()[i_scatterer].flags.use_u_aniso()):
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
    print("u_star")
    sh.structure_shifted.show_summary().show_scatterers()
    print()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), True, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(u_aniso=True),
    f_obs.anomalous_flag())
  xray.set_scatterer_grad_flags(scatterers = sh.structure_shifted.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  grad_flags_counts = xray.scatterer_grad_flags_counts(sh.structure_shifted.scatterers())
  if(grad_flags_counts.n_parameters() > 0):
     if (0):
       print("u_aniso")
       print("gradient_flags.site      ", gradient_flags.site)
       print("gradient_flags.u_iso     ", gradient_flags.u_iso)
       print("gradient_flags.u_aniso   ", gradient_flags.u_aniso)
       print("gradient_flags.occupancy ", gradient_flags.occupancy)
       print("gradient_flags.fp        ", gradient_flags.fp)
       print("gradient_flags.fdp       ", gradient_flags.fdp)
       cntr_use_u_iso = 0
       cntr_use_u_aniso = 0
       cntr_grad_u_iso = 0
       cntr_grad_u_aniso = 0
       for scatterer in sh.structure_shifted.scatterers():
         if (scatterer.flags.use_u_iso()):  cntr_use_u_iso += 1
         if (scatterer.flags.use_u_aniso()):  cntr_use_u_aniso += 1
         if (scatterer.flags.grad_u_iso()):  cntr_grad_u_iso += 1
         if (scatterer.flags.grad_u_aniso()):  cntr_grad_u_aniso += 1
       print("use_u_iso                ", cntr_use_u_iso,cntr_grad_u_iso)
       print("use_u_aniso              ", cntr_use_u_aniso,cntr_grad_u_aniso)
     sfd = xray.structure_factors.gradients_direct(
       xray_structure=sh.structure_shifted,
       u_iso_refinable_params=None,
       miller_set=f_obs,
       d_target_d_f_calc=ls.derivatives(),
       n_parameters= 0
       )
     re = resampling(miller_set=f_obs)
     map0 = re(
       xray_structure=sh.structure_shifted,
       u_iso_refinable_params=None,
       dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
       n_parameters= 0,
       verbose=verbose)

     grad_flags_counts = xray.scatterer_grad_flags_counts(sh.structure_shifted.scatterers())
     if(grad_flags_counts.u_aniso):
        sfd_d_target_d_u_cart = sfd.d_target_d_u_cart()
        map0_d_target_d_u_cart = map0.d_target_d_u_cart()
     top_gradient = None
     gradients_1 = []
     for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
       if(scatterer.flags.use_u_iso()):   parameter_name = "u_iso"
       if(scatterer.flags.use_u_aniso()): parameter_name = "u_star"
       if(parameter_name == "u_iso" and scatterer.flags.grad_u_iso()):
          direct_summ = sfd.d_target_d_u_iso()[i_scatterer]
          if (top_gradient is None): top_gradient = direct_summ
          fast_gradie = map0.d_target_d_u_iso()[i_scatterer]
          sys.stdout.flush()
          gradients_1.append([direct_summ, fast_gradie])
          match = judge(scatterer, parameter_name, direct_summ, fast_gradie,
                                                                  top_gradient)
          if (0 or verbose):
            print("direct summ[%d]: " % i_scatterer, direct_summ)
            print("fast gradie[%d]: " % i_scatterer, fast_gradie, match)
            print()
          assert not match.is_bad
       if parameter_name == "u_star" and scatterer.flags.grad_u_aniso():
        sfd_star = sfd.d_target_d_u_star()[i_scatterer]
        sfd_cart = adptbx.grad_u_star_as_u_cart(
          structure_ideal.unit_cell(), sfd_star)
        assert approx_equal(
          sfd_star,
          adptbx.grad_u_cart_as_u_star(structure_ideal.unit_cell(), sfd_cart))
        for ij in range(6):
          direct_summ = sfd_d_target_d_u_cart[i_scatterer][ij]
          if (top_gradient is None): top_gradient = direct_summ
          fast_gradie = map0_d_target_d_u_cart[i_scatterer][ij]
          gradients_1.append([direct_summ, fast_gradie])
          match =judge(scatterer,"u_star",direct_summ,fast_gradie,top_gradient)
          if (0 or verbose or match.is_bad):
            print("direct summ[%d][%d]: " % (i_scatterer, ij), direct_summ)
            print("fast gradie[%d][%d]: " % (i_scatterer, ij),fast_gradie,match)
            print()
          assert not match.is_bad
     # Making sure that gradients_1 = gradients_2
     for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
         if(not scatterer.flags.use_u_iso()):
            scatterer.u_iso = -12345.0
         if(not scatterer.flags.use_u_aniso()):
            scatterer.u_star =(-999.,-999.,-999.,-999.,-999.,-999.)
     sfd = xray.structure_factors.gradients_direct(
       xray_structure=sh.structure_shifted,
       u_iso_refinable_params=None,
       miller_set=f_obs,
       d_target_d_f_calc=ls.derivatives(),
       n_parameters= 0
       )
     re = resampling(miller_set=f_obs)
     map0 = re(
       xray_structure=sh.structure_shifted,
       u_iso_refinable_params=None,
       dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
       n_parameters= 0,
       verbose=verbose)

     grad_flags_counts = \
            xray.scatterer_grad_flags_counts(sh.structure_shifted.scatterers())
     if(grad_flags_counts.u_aniso):
        sfd_d_target_d_u_cart = sfd.d_target_d_u_cart()
        map0_d_target_d_u_cart = map0.d_target_d_u_cart()
     gradients_2 = []
     for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
       if(scatterer.flags.use_u_iso()):   parameter_name = "u_iso"
       if(scatterer.flags.use_u_aniso()): parameter_name = "u_star"
       if(parameter_name == "u_iso" and scatterer.flags.grad_u_iso()):
          direct_summ = sfd.d_target_d_u_iso()[i_scatterer]
          fast_gradie = map0.d_target_d_u_iso()[i_scatterer]
          gradients_2.append([direct_summ, fast_gradie])
       if parameter_name == "u_star" and scatterer.flags.grad_u_aniso():
        sfd_star = sfd.d_target_d_u_star()[i_scatterer]
        sfd_cart = adptbx.grad_u_star_as_u_cart(
          structure_ideal.unit_cell(), sfd_star)
        assert approx_equal(
          sfd_star,
          adptbx.grad_u_cart_as_u_star(structure_ideal.unit_cell(), sfd_cart))
        for ij in range(6):
          direct_summ = sfd_d_target_d_u_cart[i_scatterer][ij]
          fast_gradie = map0_d_target_d_u_cart[i_scatterer][ij]
          gradients_2.append([direct_summ, fast_gradie])
     for g1,g2 in zip(gradients_1, gradients_2):
       assert approx_equal(g1, g2)
     sys.stdout.flush()

class shifted_occupancy(object):

  def __init__(self, f_obs, structure, i_scatterer, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    self.structure_shifted.scatterers()[i_scatterer].occupancy += shift
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted).f_calc()

def occupancy(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_occupancy(f_obs, structure_ideal, 0, 0.2)
  if (0 or verbose):
    print("occupancy")
    sh.structure_shifted.show_summary().show_scatterers()
    print()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), True, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(occupancy=True),
    f_obs.anomalous_flag())
  xray.set_scatterer_grad_flags(scatterers = sh.structure_shifted.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    n_parameters=0,
    verbose=verbose)
  top_gradient = None
  for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
    direct_summ = sfd.d_target_d_occupancy()[i_scatterer]
    if (top_gradient is None): top_gradient = direct_summ
    fast_gradie = map0.d_target_d_occupancy()[i_scatterer]
    match = judge(scatterer, "occupancy", direct_summ,fast_gradie,top_gradient)
    if (0 or verbose):
      print("direct summ[%d]: " % i_scatterer, direct_summ)
      print("fast gradie[%d]: " % i_scatterer, fast_gradie, match)
      print()
    assert not match.is_bad
  sys.stdout.flush()

class shifted_fp(object):

  def __init__(self, f_obs, structure, i_scatterer, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    self.structure_shifted.scatterers()[i_scatterer].fp += shift
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted).f_calc()

def fp(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_fp(f_obs, structure_ideal, 0, -0.2)
  if (0 or verbose):
    print("fp")
    sh.structure_shifted.show_summary().show_scatterers()
    print()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), True, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(fp=True),
    f_obs.anomalous_flag())
  xray.set_scatterer_grad_flags(scatterers = sh.structure_shifted.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    n_parameters=0,
    verbose=verbose)
  top_gradient = None
  for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
    direct_summ = sfd.d_target_d_fp()[i_scatterer]
    if (top_gradient is None): top_gradient = direct_summ
    fast_gradie = map0.d_target_d_fp()[i_scatterer]
    match = judge(scatterer, "fp", direct_summ, fast_gradie, top_gradient)
    if (0 or verbose):
      print("direct summ[%d]: " % i_scatterer, direct_summ)
      print("fast gradie[%d]: " % i_scatterer, fast_gradie, match)
      print()
    assert not match.is_bad
  sys.stdout.flush()

class shifted_fdp(object):

  def __init__(self, f_obs, structure, i_scatterer, shift):
    self.structure_shifted = structure.deep_copy_scatterers()
    self.structure_shifted.scatterers()[i_scatterer].fdp += shift
    if (f_obs is not None):
      self.f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure=self.structure_shifted).f_calc()

def fdp(structure_ideal, d_min, f_obs, verbose=0):
  sh = shifted_fdp(f_obs, structure_ideal, 0, 2)
  if (0 or verbose):
    print("fdp")
    sh.structure_shifted.show_summary().show_scatterers()
    print()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), True, 1)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(fdp=True),
    f_obs.anomalous_flag())
  xray.set_scatterer_grad_flags(scatterers = sh.structure_shifted.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  sfd = xray.structure_factors.gradients_direct(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    n_parameters=0)
  re = resampling(miller_set=f_obs)
  map0 = re(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
    n_parameters=0,
    verbose=verbose)
  top_gradient = None
  for i_scatterer,scatterer in enumerate(sh.structure_shifted.scatterers()):
    direct_summ = sfd.d_target_d_fdp()[i_scatterer]
    if (top_gradient is None): top_gradient = direct_summ
    fast_gradie = map0.d_target_d_fdp()[i_scatterer]
    match = judge(scatterer, "fdp", direct_summ, fast_gradie, top_gradient)
    if (0 or verbose):
      print("direct summ[%d]: " % i_scatterer, direct_summ)
      print("fast gradie[%d]: " % i_scatterer, fast_gradie, match)
      print()
    assert not match.is_bad
  sys.stdout.flush()

def shift_all(structure_ideal, f_obs, anomalous_flag):
  sh = shifted_site(None, structure_ideal, 0, 0, 0.01)
  sh = shifted_u_iso(None, sh.structure_shifted, 0, 0.05)
  sh = shifted_u_star(None, sh.structure_shifted, 0, 0, 0.0001)
  sh = shifted_occupancy(None, sh.structure_shifted, 0, 0.2)
  if (anomalous_flag):
    sh = shifted_fdp(None, sh.structure_shifted, 0, 2)
  sh = shifted_fp(f_obs, sh.structure_shifted, 0, -0.2)
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), True, 1)
  return sh, ls

def exercise_packed(structure_ideal, f_obs,
                    anomalous_flag,
                    verbose=0):
  sh, ls = shift_all(structure_ideal, f_obs, anomalous_flag)
  flag = (random.random() > 0.5)
  gradient_flags = randomize_gradient_flags(
    xray.structure_factors.gradient_flags(site=flag, u=not flag),
    f_obs.anomalous_flag(),
    thresholds=(1/2.,0))
  u_iso_refinable_params = flex.double()
  for scatterer in sh.structure_shifted.scatterers():
      scatterer.flags.set_grad_site(gradient_flags.site)
      scatterer.flags.set_grad_u_iso(gradient_flags.u_iso)
      scatterer.flags.set_grad_u_aniso(gradient_flags.u_aniso)
      scatterer.flags.set_grad_occupancy(gradient_flags.occupancy)
      scatterer.flags.set_grad_fp(gradient_flags.fp)
      scatterer.flags.set_grad_fdp(gradient_flags.fdp)
      scatterer.flags.set_tan_u_iso(True)
      param = random.randint(90,120)
      scatterer.flags.param= param
      value = math.tan(scatterer.u_iso*math.pi/adptbx.b_as_u(param)-math.pi/2)
      u_iso_refinable_params.append(value)
  n_parameters = xray.scatterer_grad_flags_counts(
                              sh.structure_shifted.scatterers()).n_parameters()
  assert n_parameters == sh.structure_shifted.n_parameters()
  if (n_parameters > 0):
    sfd = xray.structure_factors.gradients_direct(
      xray_structure=sh.structure_shifted,
      u_iso_refinable_params=u_iso_refinable_params,
      miller_set=f_obs,
      d_target_d_f_calc=ls.derivatives(),
      n_parameters=n_parameters)
    assert sfd.packed().size() == n_parameters
    re = resampling(miller_set=f_obs)
    map0 = re(
      xray_structure=sh.structure_shifted,
      u_iso_refinable_params=u_iso_refinable_params,
      dp=miller.array(miller_set=f_obs, data=ls.derivatives()),
      n_parameters=n_parameters,
      verbose=verbose)
    assert map0.packed().size() == n_parameters
    correlation = flex.linear_correlation(sfd.packed(), map0.packed())
    assert correlation.is_well_defined()
    assert correlation.coefficient() > 0.999

def exercise_gradient_manager(structure_ideal, f_obs,
                              anomalous_flag,
                              verbose=0):
  sh, ls = shift_all(structure_ideal, f_obs, anomalous_flag)
  grad_manager = xray.structure_factors.gradients(
    miller_set=f_obs,
    quality_factor=100000,
    wing_cutoff=1.e-10)
  gradient_flags=xray.structure_factors.gradient_flags(default=True)
  xray.set_scatterer_grad_flags(scatterers = sh.structure_shifted.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  if (0):
    print("exercise_gradient_manager")
    print("gradient_flags.site      ", gradient_flags.site)
    print("gradient_flags.u_iso     ", gradient_flags.u_iso)
    print("gradient_flags.u_aniso   ", gradient_flags.u_aniso)
    print("gradient_flags.occupancy ", gradient_flags.occupancy)
    print("gradient_flags.fp        ", gradient_flags.fp)
    print("gradient_flags.fdp       ", gradient_flags.fdp)
    cntr_use_u_iso = 0
    cntr_use_u_aniso = 0
    cntr_grad_u_iso = 0
    cntr_grad_u_aniso = 0
    for scatterer in sh.structure_shifted.scatterers():
      if (scatterer.flags.use_u_iso()):  cntr_use_u_iso += 1
      if (scatterer.flags.use_u_aniso()):  cntr_use_u_aniso += 1
      if (scatterer.flags.grad_u_iso()):  cntr_grad_u_iso += 1
      if (scatterer.flags.grad_u_aniso()):  cntr_grad_u_aniso += 1
    print("use_u_iso                ", cntr_use_u_iso,cntr_grad_u_iso)
    print("use_u_aniso              ", cntr_use_u_aniso,cntr_grad_u_aniso)
  if (random.random() > 0.5):
    n_parameters = 0
  else:
    n_parameters = xray.scatterer_grad_flags_counts(
                              sh.structure_shifted.scatterers()).n_parameters()
    assert n_parameters == sh.structure_shifted.n_parameters()
  gd = grad_manager(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    n_parameters=n_parameters,
    algorithm="direct")
  gf = grad_manager(
    xray_structure=sh.structure_shifted,
    u_iso_refinable_params=None,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    n_parameters=n_parameters,
    algorithm="fft")
  grad_flags_counts = \
            xray.scatterer_grad_flags_counts(sh.structure_shifted.scatterers())
  if (n_parameters == 0):
    d = gd.d_target_d_site_frac()
    f = gf.d_target_d_site_frac()
    linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    d = gd.d_target_d_site_cart()
    f = gf.d_target_d_site_cart()
    linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    if(grad_flags_counts.u_iso > 0):
       d = gd.d_target_d_u_iso()
       f = gf.d_target_d_u_iso()
       linear_regression_test(d, f, slope_tolerance=1.e-2, verbose=verbose)
    if(grad_flags_counts.u_aniso > 0):
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
    assert correlation.coefficient() > 0.995, correlation.coefficient()

def run_one(space_group_info, n_elements= 9, volume_per_atom=1000, d_min = 2.0,
            anomalous_flag=0,
            verbose=0):
  if (random.random() < 0.5):
    random_f_prime_scale=0.6
  else:
    random_f_prime_scale=0
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=(("O","N","C")*(n_elements))[:n_elements],#(("O","N","C")*(n_elements/3+1))[:n_elements],
    volume_per_atom=volume_per_atom,
    min_distance=5,
    general_positions_only=True,
    random_f_prime_d_min=d_min-1,
    random_f_prime_scale=random_f_prime_scale,
    random_f_double_prime=anomalous_flag,
    use_u_aniso = True,
    use_u_iso = False,
    random_u_iso=True,
    random_u_cart_scale=.3,
    random_occupancy=True)
  random_structure.random_modify_adp_and_adp_flags(
                             scatterers         = structure_ideal.scatterers(),
                             random_u_iso_scale = 0.3,
                             random_u_iso_min   = 0.0)
  if (random.random() < 0.5):
    assign_custom_gaussians(structure_ideal, negative_a=random.random()<0.5)
  if (0 or verbose):
    structure_ideal.show_summary().show_scatterers()
  f_obs = abs(structure_ideal.structure_factors(
    d_min=d_min, anomalous_flag=anomalous_flag, algorithm="direct").f_calc())
  if (1):
    site(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1):
    u_iso(structure_ideal,  d_min, f_obs, tan_u_iso=False, verbose=verbose)
    u_iso(structure_ideal,  d_min, f_obs, tan_u_iso=True, verbose=verbose)
    u_star(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1):
    occupancy(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1):
    fp(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1 and anomalous_flag):
    fdp(structure_ideal, d_min, f_obs, verbose=verbose)
  if (1):
    exercise_gradient_manager(structure_ideal, f_obs, anomalous_flag)
  if (1):
    exercise_packed(structure_ideal, f_obs, anomalous_flag)

def run_call_back(flags, space_group_info):
  for anomalous_flag in [False,True]:
    run_one(
      space_group_info=space_group_info,
      anomalous_flag=anomalous_flag,
      verbose=flags.Verbose)

def run():
  show_times = libtbx.utils.show_times()
  debug_utils.parse_options_loop_space_groups(
    argv=sys.argv[1:],
    call_back=run_call_back,
    show_cpu_times=False)
  xray.structure_factors.global_counters.show()
  show_times()

if (__name__ == "__main__"):
  run()
