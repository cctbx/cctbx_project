from cctbx.development import random_structure
from cctbx import xray
from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.array_family import flex
from cctbx import matrix
from scitbx.python_utils import complex_math
from scitbx.python_utils.misc import adopt_init_args, user_plus_sys_time
from scitbx.test_utils import approx_equal
from scitbx import fftpack
import random
import math

random.seed(0)

class resampling(crystal.symmetry):

  def __init__(self, miller_set=None,
                     crystal_symmetry=None,
                     d_min=None,
                     grid_resolution_factor=1/3.,
                     symmetry_flags=maptbx.use_space_group_symmetry,
                     mandatory_grid_factors=None,
                     quality_factor=100,
                     wing_cutoff=1.e-3,
                     exp_table_one_over_step_size=-100,
                     max_prime=5):
    assert miller_set is None or crystal_symmetry is None
    if (miller_set is None):
      assert crystal_symmetry is not None and d_min is not None
    else:
      crystal_symmetry = miller_set
      if (d_min is None):
        d_min = miller_set.d_min()
      else:
        assert d_min <= miller_set.d_min()
    crystal.symmetry._copy_constructor(self, crystal_symmetry)
    del miller_set
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
    n = self.rfft().n_real()
    norm = self.unit_cell().volume()/(n[0]*n[1]*n[2])
    dpe = dp.deep_copy()
    xray.eliminate_u_extra(
      self.unit_cell(),
      self.u_extra(),
      dpe.indices(),
      dpe.data(),
      norm)
    return miller.fft_map(
      crystal_gridding=self.crystal_gridding(),
      fourier_coefficients=dpe)

  def __call__(self, xray_structure,
                     dp,
                     d_target_d_f_calc=None,
                     derivative_flags=None,
                     force_complex=00000,
                     electron_density_must_be_positive=0001):
    self.setup_fft()
    if (0):
      xray.sampled_model_density(
        xray_structure.unit_cell(),
        xray_structure.scatterers(),
        self.rfft().n_real(),
        self.rfft().m_real(),
        self.u_extra(),
        self.wing_cutoff(),
        self.exp_table_one_over_step_size(),
        force_complex,
        electron_density_must_be_positive)
      print "sampled_model_density OK"
    cmap = self.ft_dp(dp).complex_map()
    print "ft_dt_map real: %.4g %.4g" % (
      flex.min(flex.real(cmap)), flex.max(flex.real(cmap)))
    print "ft_dt_map imag: %.4g %.4g" % (
      flex.min(flex.imag(cmap)), flex.max(flex.imag(cmap)))
    print
    time_sampling = user_plus_sys_time()
    result = xray.fast_gradients(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      self.ft_dp(dp).complex_map(),
      self.rfft().n_real(),
      self.rfft().m_real(),
      self.u_extra(),
      self.wing_cutoff(),
      self.exp_table_one_over_step_size(),
      force_complex,
      electron_density_must_be_positive)
    time_sampling = time_sampling.elapsed()
    return result

def get_gms(xray_structure, f_obs):
  uc = xray_structure.unit_cell()
  gms = []
  for m in xray_structure.scatterers():
    gm = flex.complex_double()
    for i,hkl in f_obs.indices().items():
      if (not m.anisotropic_flag):
        d = adptbx.debye_waller_factor_u_iso(uc, hkl, m.u_iso)
      else:
        d = adptbx.debye_waller_factor_u_star(hkl, m.u_star)
      f = m.caasf.at_d_star_sq(uc.d_star_sq(hkl))
      gm.append(m.weight()*d*(f+m.fp_fdp))
    gms.append(gm)
  return gms

def get_dp0(f_obs, phi, e):
  dp0 = flex.complex_double()
  for i in phi.indices():
    dp0.append(-e[i]*complex_math.polar((1, phi[i])))
  return miller.array(miller_set=f_obs, data=dp0)

class two_p_shifted_site:

  def __init__(self, f_obs, structure, i_scatterer, i_xyz, shift):
    structure_shifted = structure.deep_copy_scatterers()
    site = list(structure_shifted.scatterers()[i_scatterer].site)
    site[i_xyz] += shift
    structure_shifted.scatterers()[i_scatterer].site = site
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_shifted).f_calc()
    f_calc_abs = abs(f_calc)
    e = f_calc_abs.data() - f_obs.data()
    two_p = flex.sum(flex.pow2(e))
    self.structure_shifted = structure_shifted
    self.f_calc = f_calc
    self.e = e
    self.two_p = two_p

def site(structure_ideal, f_obs):
  sum_f_obs_sq = flex.sum(flex.pow2(f_obs.data()))
  sh = two_p_shifted_site(f_obs, structure_ideal, 0, 0, 0.05)
  sh.structure_shifted.show_summary().show_scatterers()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  sfd = xray.structure_factors.from_scatterers_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    derivative_flags=xray.structure_factors.derivative_flags(
      site=0001))
  gms = get_gms(sh.structure_shifted, f_obs)
  phi = flex.arg(sh.f_calc.data())
  dp0 = get_dp0(f_obs, phi, sh.e)
  dps = []
  for i_xyz in (0,1,2):
    dp = flex.complex_double()
    for i,hkl in f_obs.indices().items():
      dp.append(-1j*2*math.pi*hkl[i_xyz]*sh.e[i]
                *complex_math.polar((1, phi[i])))
    dps.append(miller.array(miller_set=f_obs, data=dp))
  print "two_p:", sh.two_p
  print "ls.target():", ls.target() * sum_f_obs_sq
  assert approx_equal(sh.two_p, ls.target() * sum_f_obs_sq)
  print
  re = resampling(miller_set=f_obs)
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  for i_scatterer in (0,1,2):
    for i_xyz in (0,1,2):
      delta = 1.e-6
      pl = two_p_shifted_site(
        f_obs, sh.structure_shifted, i_scatterer,i_xyz,delta)
      mi = two_p_shifted_site(
        f_obs, sh.structure_shifted, i_scatterer,i_xyz,-delta)
      j = (pl.e - mi.e) / (2*delta)
      g = flex.sum(j * sh.e)
      print "  g[%d][%d]: " % (i_scatterer, i_xyz), g
      print "sfd[%d][%d]: " % (i_scatterer, i_xyz), \
            sfd.d_target_d_site()[i_scatterer][i_xyz] * sum_f_obs_sq/2
      rm = matrix.col(sh.structure_shifted.scatterers()[i_scatterer].site)
      gxm = 0
      for i,hkl in f_obs.indices().items():
        p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
        gxm += (dps[i_xyz].data()[i]
                * gms[i_scatterer][i].conjugate()
                * complex_math.polar((1, p)))
      print "gxm[%d][%d]:" % (i_scatterer, i_xyz), gxm
      gl = map0.grad_site()[i_scatterer][i_xyz]
      print " m0[%d][%d]: " % (i_scatterer, i_xyz), gl
      print

class two_p_shifted_u_iso:

  def __init__(self, f_obs, structure, i_scatterer, shift):
    structure_shifted = structure.deep_copy_scatterers()
    structure_shifted.scatterers()[i_scatterer].u_iso += shift
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_shifted).f_calc()
    f_calc_abs = abs(f_calc)
    e = f_calc_abs.data() - f_obs.data()
    two_p = flex.sum(flex.pow2(e))
    self.structure_shifted = structure_shifted
    self.f_calc = f_calc
    self.e = e
    self.two_p = two_p

def u_iso(structure_ideal, f_obs):
  sum_f_obs_sq = flex.sum(flex.pow2(f_obs.data()))
  sh = two_p_shifted_u_iso(f_obs, structure_ideal, 0, 0.05)
  sh.structure_shifted.show_summary().show_scatterers()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  sfd = xray.structure_factors.from_scatterers_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    derivative_flags=xray.structure_factors.derivative_flags(
      u_iso=0001))
  gms = get_gms(sh.structure_shifted, f_obs)
  phi = flex.arg(sh.f_calc.data())
  dp0 = get_dp0(f_obs, phi, sh.e)
  uc = sh.structure_shifted.unit_cell()
  dps = flex.complex_double()
  for i,hkl in f_obs.indices().items():
    s_sq = uc.d_star_sq(hkl)
    dps.append(-s_sq/4*8*math.pi*math.pi*sh.e[i]
               *complex_math.polar((1, phi[i])))
  dps = miller.array(miller_set=f_obs, data=dps)
  print "two_p:", sh.two_p
  print "ls.target():", ls.target() * sum_f_obs_sq
  assert approx_equal(sh.two_p, ls.target() * sum_f_obs_sq)
  print
  re = resampling(miller_set=f_obs)
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  for i_scatterer in (0,1,2):
    delta = 1.e-6
    pl = two_p_shifted_u_iso(
      f_obs, sh.structure_shifted, i_scatterer, delta)
    mi = two_p_shifted_u_iso(
      f_obs, sh.structure_shifted, i_scatterer, -delta)
    j = (pl.e - mi.e) / (2*delta)
    g = flex.sum(j * sh.e)
    print "  g[%d]: " % i_scatterer, g
    print "sfd[%d]: " % i_scatterer, \
          sfd.d_target_d_u_iso()[i_scatterer] * sum_f_obs_sq/2
    rm = matrix.col(sh.structure_shifted.scatterers()[i_scatterer].site)
    gxm = 0
    for i,hkl in f_obs.indices().items():
      p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
      gxm += (dps.data()[i]
              * gms[i_scatterer][i].conjugate()
              * complex_math.polar((1, p)))
    print "gxm[%d]:" % i_scatterer, gxm
    gl = map0.grad_u_iso()[i_scatterer]
    print " m0[%d]: " % i_scatterer, gl
    print

class two_p_shifted_occupancy:

  def __init__(self, f_obs, structure, i_scatterer, shift):
    structure_shifted = structure.deep_copy_scatterers()
    structure_shifted.shift_occupancy(i_scatterer, shift)
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_shifted).f_calc()
    f_calc_abs = abs(f_calc)
    e = f_calc_abs.data() - f_obs.data()
    two_p = flex.sum(flex.pow2(e))
    self.structure_shifted = structure_shifted
    self.f_calc = f_calc
    self.e = e
    self.two_p = two_p

def occupancy(structure_ideal, f_obs):
  sum_f_obs_sq = flex.sum(flex.pow2(f_obs.data()))
  sh = two_p_shifted_occupancy(f_obs, structure_ideal, 0, 0.2)
  sh.structure_shifted.show_summary().show_scatterers()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  sfd = xray.structure_factors.from_scatterers_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    derivative_flags=xray.structure_factors.derivative_flags(
      occupancy=0001))
  gms = get_gms(sh.structure_shifted, f_obs)
  phi = flex.arg(sh.f_calc.data())
  dp0 = get_dp0(f_obs, phi, sh.e)
  dps = flex.complex_double()
  for i,hkl in f_obs.indices().items():
    dps.append(sh.e[i]
               *complex_math.polar((1, phi[i])))
  dps = miller.array(miller_set=f_obs, data=dps)
  print "two_p:", sh.two_p
  print "ls.target():", ls.target() * sum_f_obs_sq
  assert approx_equal(sh.two_p, ls.target() * sum_f_obs_sq)
  print
  re = resampling(miller_set=f_obs)
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  for i_scatterer in (0,1,2):
    delta = 1.e-6
    pl = two_p_shifted_occupancy(
      f_obs, sh.structure_shifted, i_scatterer, delta)
    mi = two_p_shifted_occupancy(
      f_obs, sh.structure_shifted, i_scatterer, -delta)
    j = (pl.e - mi.e) / (2*delta)
    g = flex.sum(j * sh.e)
    print "  g[%d]: " % i_scatterer, g
    print "sfd[%d]: " % i_scatterer, \
          sfd.d_target_d_occupancy()[i_scatterer] * sum_f_obs_sq/2
    m = sh.structure_shifted.scatterers()[i_scatterer]
    rm = matrix.col(m.site)
    gxm = 0
    for i,hkl in f_obs.indices().items():
      p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
      gxm += (dps.data()[i]
              * gms[i_scatterer][i].conjugate()
              * complex_math.polar((1, p)))
    print "gxm[%d]:" % i_scatterer, gxm/m.occupancy
    gl = map0.grad_occupancy()[i_scatterer]
    print " m0[%d]:" % i_scatterer, gl
    print

class two_p_shifted_u_cart:

  def __init__(self, f_obs, structure, i_scatterer, ij, shift):
    structure_shifted = structure.deep_copy_scatterers()
    scatterer = structure_shifted.scatterers()[i_scatterer]
    u_cart = list(
      adptbx.u_star_as_u_cart(structure.unit_cell(), scatterer.u_star))
    u_cart[ij] += shift
    u_star = adptbx.u_cart_as_u_star(structure.unit_cell(), u_cart)
    scatterer.u_star = u_star
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_shifted).f_calc()
    f_calc_abs = abs(f_calc)
    e = f_calc_abs.data() - f_obs.data()
    two_p = flex.sum(flex.pow2(e))
    self.structure_shifted = structure_shifted
    self.f_calc = f_calc
    self.e = e
    self.two_p = two_p

class two_p_shifted_u_star:

  def __init__(self, f_obs, structure, i_scatterer, ij, shift):
    structure_shifted = structure.deep_copy_scatterers()
    scatterer = structure_shifted.scatterers()[i_scatterer]
    u_star = list(scatterer.u_star)
    u_star[ij] += shift
    scatterer.u_star = u_star
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_shifted).f_calc()
    f_calc_abs = abs(f_calc)
    e = f_calc_abs.data() - f_obs.data()
    two_p = flex.sum(flex.pow2(e))
    self.structure_shifted = structure_shifted
    self.f_calc = f_calc
    self.e = e
    self.two_p = two_p

def ij_product(hkl, ij):
  if (ij < 3): return hkl[ij]**2
  if (ij == 3): return 2*hkl[0]*hkl[1]
  if (ij == 4): return 2*hkl[0]*hkl[2]
  if (ij == 5): return 2*hkl[1]*hkl[2]
  raise RuntimeError

def u_star(structure_ideal, f_obs):
  sum_f_obs_sq = flex.sum(flex.pow2(f_obs.data()))
  sh = two_p_shifted_u_star(f_obs, structure_ideal, 0, 0, 0.0001)
  sh.structure_shifted.show_summary().show_scatterers()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  sfd = xray.structure_factors.from_scatterers_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    derivative_flags=xray.structure_factors.derivative_flags(
      u_star=0001))
  gms = get_gms(sh.structure_shifted, f_obs)
  phi = flex.arg(sh.f_calc.data())
  dp0 = get_dp0(f_obs, phi, sh.e)
  dps = []
  for ij in xrange(6):
    dp = flex.complex_double()
    for i,hkl in f_obs.indices().items():
      dp.append(-2*(math.pi**2)*ij_product(hkl,ij)*sh.e[i]
                *complex_math.polar((1, phi[i])))
    dps.append(miller.array(miller_set=f_obs, data=dp))
  print "two_p:", sh.two_p
  print "ls.target():", ls.target() * sum_f_obs_sq
  assert approx_equal(sh.two_p, ls.target() * sum_f_obs_sq)
  print
  re = resampling(miller_set=f_obs)
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  for i_scatterer in (0,1,2):
    sfd_star = [x*sum_f_obs_sq/2 for x in sfd.d_target_d_u_star()[i_scatterer]]
    sfd_cart = adptbx.grad_u_star_as_u_cart(
      structure_ideal.unit_cell(), sfd_star)
    assert approx_equal(
      sfd_star,
      adptbx.grad_u_cart_as_u_star(structure_ideal.unit_cell(), sfd_cart))
    for ij in xrange(6):
      delta = 1.e-6
      pl = two_p_shifted_u_star(
        f_obs, sh.structure_shifted, i_scatterer, ij, delta)
      mi = two_p_shifted_u_star(
        f_obs, sh.structure_shifted, i_scatterer, ij, -delta)
      j = (pl.e - mi.e) / (2*delta)
      g = flex.sum(j * sh.e)
      plc = two_p_shifted_u_cart(
        f_obs, sh.structure_shifted, i_scatterer, ij, delta)
      mic = two_p_shifted_u_cart(
        f_obs, sh.structure_shifted, i_scatterer, ij, -delta)
      jc = (plc.e - mic.e) / (2*delta)
      gc = flex.sum(jc * sh.e)
      print "  g[%d][%d]: " % (i_scatterer, ij), g
      print "sfd[%d][%d]: " % (i_scatterer, ij), \
            sfd.d_target_d_u_star()[i_scatterer][ij] * sum_f_obs_sq/2
      rm = matrix.col(sh.structure_shifted.scatterers()[i_scatterer].site)
      gxm = 0
      for i,hkl in f_obs.indices().items():
        p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
        gxm += (dps[ij].data()[i]
                * gms[i_scatterer][i].conjugate()
                * complex_math.polar((1, p)))
      print "gxm[%d][%d]:" % (i_scatterer, ij), gxm
      gl = map0.grad_u_star()[i_scatterer][ij]
      print " m0[%d][%d]: " % (i_scatterer, ij), gl
      print " gc[%d][%d]: " % (i_scatterer, ij), gc
      print "s2c[%d][%d]: " % (i_scatterer, ij), sfd_cart[ij]
      print

class two_p_shifted_fp:

  def __init__(self, f_obs, structure, i_scatterer, shift):
    structure_shifted = structure.deep_copy_scatterers()
    structure_shifted.scatterers()[i_scatterer].fp_fdp += shift
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_shifted).f_calc()
    f_calc_abs = abs(f_calc)
    e = f_calc_abs.data() - f_obs.data()
    two_p = flex.sum(flex.pow2(e))
    self.structure_shifted = structure_shifted
    self.f_calc = f_calc
    self.e = e
    self.two_p = two_p

def fp(structure_ideal, f_obs):
  sum_f_obs_sq = flex.sum(flex.pow2(f_obs.data()))
  sh = two_p_shifted_fp(f_obs, structure_ideal, 0, -0.2)
  sh.structure_shifted.show_summary().show_scatterers()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  sfd = xray.structure_factors.from_scatterers_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    derivative_flags=xray.structure_factors.derivative_flags(
      fp=0001))
  gms = get_gms(sh.structure_shifted, f_obs)
  phi = flex.arg(sh.f_calc.data())
  dp0 = get_dp0(f_obs, phi, sh.e)
  dps = flex.complex_double()
  for i,hkl in f_obs.indices().items():
    dps.append(sh.e[i]
               *complex_math.polar((1, phi[i])))
  uc = sh.structure_shifted.unit_cell()
  dps = miller.array(miller_set=f_obs, data=dps)
  print "two_p:", sh.two_p
  print "ls.target():", ls.target() * sum_f_obs_sq
  assert approx_equal(sh.two_p, ls.target() * sum_f_obs_sq)
  print
  re = resampling(miller_set=f_obs)
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  for i_scatterer in (0,1,2):
    delta = 1.e-6
    pl = two_p_shifted_fp(
      f_obs, sh.structure_shifted, i_scatterer, delta)
    mi = two_p_shifted_fp(
      f_obs, sh.structure_shifted, i_scatterer, -delta)
    j = (pl.e - mi.e) / (2*delta)
    g = flex.sum(j * sh.e)
    print "  g[%d]: " % i_scatterer, g
    print "sfd[%d]: " % i_scatterer, \
          sfd.d_target_d_fp()[i_scatterer] * sum_f_obs_sq/2
    m = sh.structure_shifted.scatterers()[i_scatterer]
    rm = matrix.col(m.site)
    gxm = 0
    for i,hkl in f_obs.indices().items():
      p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
      f0_fp_fdp = m.caasf.at_d_star_sq(uc.d_star_sq(hkl))+m.fp_fdp
      gxm += (dps.data()[i]
              * (gms[i_scatterer][i]/f0_fp_fdp).conjugate()
              * complex_math.polar((1, p)))
    print "gxm[%d]:" % i_scatterer, gxm
    gl = map0.grad_fp()[i_scatterer]
    print " m0[%d]:" % i_scatterer, gl
    print

class two_p_shifted_fdp:

  def __init__(self, f_obs, structure, i_scatterer, shift):
    structure_shifted = structure.deep_copy_scatterers()
    structure_shifted.scatterers()[i_scatterer].fp_fdp += complex(0,shift)
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_shifted).f_calc()
    f_calc_abs = abs(f_calc)
    e = f_calc_abs.data() - f_obs.data()
    two_p = flex.sum(flex.pow2(e))
    self.structure_shifted = structure_shifted
    self.f_calc = f_calc
    self.e = e
    self.two_p = two_p

def fdp(structure_ideal, f_obs):
  sum_f_obs_sq = flex.sum(flex.pow2(f_obs.data()))
  sh = two_p_shifted_fdp(f_obs, structure_ideal, 0, 2)
  sh.structure_shifted.show_summary().show_scatterers()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  sfd = xray.structure_factors.from_scatterers_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    derivative_flags=xray.structure_factors.derivative_flags(
      fdp=0001))
  gms = get_gms(sh.structure_shifted, f_obs)
  phi = flex.arg(sh.f_calc.data())
  dp0 = get_dp0(f_obs, phi, sh.e)
  dps = flex.complex_double()
  for i,hkl in f_obs.indices().items():
    dps.append(sh.e[i]
               *complex_math.polar((1, phi[i])))
  uc = sh.structure_shifted.unit_cell()
  dps = miller.array(miller_set=f_obs, data=dps)
  print "two_p:", sh.two_p
  print "ls.target():", ls.target() * sum_f_obs_sq
  assert approx_equal(sh.two_p, ls.target() * sum_f_obs_sq)
  print
  re = resampling(miller_set=f_obs)
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  for i_scatterer in (0,1,2):
    delta = 1.e-6
    pl = two_p_shifted_fdp(
      f_obs, sh.structure_shifted, i_scatterer, delta)
    mi = two_p_shifted_fdp(
      f_obs, sh.structure_shifted, i_scatterer, -delta)
    j = (pl.e - mi.e) / (2*delta)
    g = flex.sum(j * sh.e)
    print "  g[%d]: " % i_scatterer, g
    print "sfd[%d]: " % i_scatterer, \
          sfd.d_target_d_fdp()[i_scatterer] * sum_f_obs_sq/2
    m = sh.structure_shifted.scatterers()[i_scatterer]
    rm = matrix.col(m.site)
    gxm = 0
    for i,hkl in f_obs.indices().items():
      p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
      f0_fp_fdp = m.caasf.at_d_star_sq(uc.d_star_sq(hkl))+m.fp_fdp
      gxm += (dps.data()[i]
              * (gms[i_scatterer][i]/(-1j*f0_fp_fdp)).conjugate()
              * complex_math.polar((1, p)))
    print "gxm[%d]:" % i_scatterer, gxm
    gl = map0.grad_fdp()[i_scatterer]
    print " m0[%d]:" % i_scatterer, gl
    print

def run_one(n_elements=3, volume_per_atom=1000, d_min=2,
            fdp_flag=0, anisotropic_flag=0):
  structure_ideal = random_structure.xray_structure(
    sgtbx.space_group_info("P 1"),
    elements=("Se",)*n_elements,
    volume_per_atom=volume_per_atom,
    random_f_prime_d_min=d_min,
    random_f_double_prime=fdp_flag,
    anisotropic_flag=anisotropic_flag,
    random_u_iso=1,
    random_u_cart_scale=.1,
    random_occupancy=1)
  if (0):
    a = structure_ideal.unit_cell().volume()**(1/3.)
    structure_ideal = xray.structure(
      special_position_settings=crystal.special_position_settings(
        crystal_symmetry=crystal.symmetry(
          unit_cell=(1.3*a,a,a/1.3,90,90,90),
          space_group_symbol="P 1")),
      scatterers=structure_ideal.scatterers())
  structure_ideal.show_summary().show_scatterers()
  uc = structure_ideal.unit_cell()
  print "volume:", uc.volume()
  if (anisotropic_flag):
    for scatterer in structure_ideal.scatterers():
      print "u_iso:", adptbx.u_star_as_u_iso(uc, scatterer.u_star)
  print
  f_obs = abs(structure_ideal.structure_factors(
    d_min=d_min, anomalous_flag=0001, direct=0001).f_calc())
  if (1):
    print "site"
    site(structure_ideal, f_obs)
  if (1):
    if (not anisotropic_flag):
      print "u_iso"
      u_iso(structure_ideal, f_obs)
    else:
      print "u_star"
      u_star(structure_ideal, f_obs)
  if (1):
    print "occupancy"
    occupancy(structure_ideal, f_obs)
  if (1):
    print "fp"
    fp(structure_ideal, f_obs)
  if (1):
    print "fdp"
    fdp(structure_ideal, f_obs)

def run():
  for fdp_flag in (0, 1):
    for anisotropic_flag in (0, 1):
      run_one(fdp_flag=fdp_flag, anisotropic_flag=anisotropic_flag)

if (__name__ == "__main__"):
  run()
