from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.development import make_cns_input
from cctbx import xray
from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import adptbx
from cctbx.array_family import flex
from cctbx import matrix
from scitbx.python_utils import complex_math
from scitbx.python_utils.misc import adopt_init_args, user_plus_sys_time
from scitbx.test_utils import approx_equal
from scitbx import fftpack
from scitbx.python_utils import easy_pickle
import random
import math
import sys

random.seed(0)

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
    n = self.rfft().n_real()
    norm = self.unit_cell().volume()/(n[0]*n[1]*n[2])
    dpe = dp.deep_copy()
    xray.eliminate_u_extra(
      self.unit_cell(),
      self.u_extra(),
      dpe.indices(),
      dpe.data(),
      norm)
    print "vol:", self.unit_cell().volume()
    print "n[0]*n[1]*n[2]:", n[0]*n[1]*n[2]
    print "n[0],n[1],n[2]:", n[0],n[1],n[2]
    print "norm:", norm
    print "b_extra:", adptbx.u_as_b(self.u_extra())
    if (0):
      f = open("dtdf_cctbx", "w")
      for i,h in dp.indices().items():
        print >> f, h[0],h[1],h[2], dp.data()[i].real, dp.data()[i].imag
        print "fcalc:", h, dp.data()[i], dpe.data()[i]/2
        #print "fcalc:", h, dpe.data()[i].real/dp.data()[i].real/norm
      f.close()
    dpe = miller.array(dpe, dpe.data() \
                            * flex.polar(dpe.epsilons().data().as_double(),0))
    result = miller.fft_map(
      crystal_gridding=self.crystal_gridding(),
      fourier_coefficients=dpe)
    print "grid:", result.complex_map().focus()
    if (0):
      for i in xrange(result.complex_map().focus()[0]):
        for j in xrange(result.complex_map().focus()[1]):
          for k in xrange(result.complex_map().focus()[2]):
            print "rho:", i,j,k, result.complex_map()[(i,j,k)].real
    if (0):
      cmap = result.complex_map().deep_copy()
      assert not cmap.is_padded()
      cfft = fftpack.complex_to_complex_3d(result.n_real())
      map_bt = cfft.backward(cmap)
      collect_conj = 0
      dp_bt = maptbx.structure_factors.from_map(
      dpe.anomalous_flag(),
      dpe.indices(),
      map_bt,
      collect_conj).data()
      print "dp_bt_begin"
      max_diff = 0
      for i,h in dpe.indices().items():
        x = dpe.data()[i]
        y = dp_bt[i]/map_bt.size()
        print h, x, y
        if (abs(x-y) > max_diff):
          print "DIFF:", abs(x-y)
          max_diff = abs(x-y)
      print "dp_bt_end"
    return result

  def __call__(self, xray_structure,
                     dp,
                     d_target_d_f_calc=None,
                     derivative_flags=None,
                     electron_density_must_be_positive=0001):
    self.setup_fft()
    cmap = self.ft_dp(dp).complex_map()
    assert not cmap.is_padded()
    print "grid:", cmap.focus()
    print "ft_dt_map real: %.4g %.4g" % (
      flex.min(flex.real(cmap)), flex.max(flex.real(cmap)))
    print "ft_dt_map imag: %.4g %.4g" % (
      flex.min(flex.imag(cmap)), flex.max(flex.imag(cmap)))
    print
    if (0):
      rho_cns = easy_pickle.load("rho_cns.pickle")
      assert rho_cns.size() == cmap.size()
      rho_cns.resize(cmap.accessor())
      print "Using rho_cns"
      cmap = flex.polar(rho_cns, 0)
    print "b_extra:", adptbx.u_as_b(self.u_extra())
    time_sampling = user_plus_sys_time()
    result = xray.fast_gradients(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      cmap,
      self.u_extra(),
      self.wing_cutoff(),
      self.exp_table_one_over_step_size(),
      electron_density_must_be_positive)
    time_sampling = time_sampling.elapsed()
    print "max_shell_radii:", result.max_shell_radii()
    print "exp_table_size:", result.exp_table_size()
    print
    return result

def judge(scatterer, label, reference, other,
          agreement_factor=0.99, top=None):
  label += [" iso", " aniso"][int(scatterer.anisotropic_flag)]
  s = ""
  if (top is None):
    if (abs(reference-other) > abs(reference*(1-agreement_factor))):
      s += " BAD " + label
  else:
    r = (reference-other)/top
    s += " %.5f " % r + label
    if (abs(r) > 0.03):
      raise RuntimeError("very large mismatch")
    elif (abs(r) > 0.01):
      s += " large mismatch"
  return s.lstrip()

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
      xray_structure=structure_shifted, direct=0001).f_calc()
    f_calc_abs = abs(f_calc)
    e = f_calc_abs.data() - f_obs.data()
    two_p = flex.sum(flex.pow2(e))
    self.structure_shifted = structure_shifted
    self.f_calc = f_calc
    self.e = e
    self.two_p = two_p

def site(structure_ideal, d_min, f_obs, show_finite):
  sum_f_obs_sq = flex.sum(flex.pow2(f_obs.data()))
  sh = two_p_shifted_site(f_obs, structure_ideal, 0, 0, 0.01)
  sh.structure_shifted.show_summary().show_scatterers()
  if (0 and not sh.structure_shifted.scatterers()[0].anisotropic_flag):
    for file_name,method in (("gradients_direct.inp", "direct"),
                             ("gradients_fft.inp", "fft")):
      make_cns_input.write(file_name,
        make_cns_input.script_xray_gradients(
          d_min, f_obs, sh.structure_shifted, method))
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  ls_derivatives = ls.derivatives().deep_copy()
  if (0):
    centric = f_obs.centric_flags().data()
    n_zero = 0
    for i,c in centric.items():
      if (c):
        ls_derivatives[i] = complex(0)
        n_zero += 1
    print "n_zero:", n_zero
  sfd = xray.structure_factors.from_scatterers_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls_derivatives,
    derivative_flags=xray.structure_factors.derivative_flags(
      site=0001))
  gms = get_gms(sh.structure_shifted, f_obs)
  phi = flex.arg(sh.f_calc.data())
  dp0 = get_dp0(f_obs, phi, sh.e)
  if (f_obs.space_group().order_z() == 1):
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
  assert dp0.indices().all_eq(f_obs.indices())
  dp0 = miller.array(miller_set=f_obs, data=ls_derivatives) # XXX
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  sfd.d_target_d_site_inplace_frac_as_cart(sfd.d_target_d_site())
  sfd.d_target_d_site_inplace_frac_as_cart(map0.d_target_d_site())
  top_gradient = None
  for i_scatterer in (0,1,2):
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    for i_xyz in (0,1,2):
      if (show_finite):
        delta = 1.e-6
        pl = two_p_shifted_site(
          f_obs, sh.structure_shifted, i_scatterer,i_xyz,delta)
        mi = two_p_shifted_site(
          f_obs, sh.structure_shifted, i_scatterer,i_xyz,-delta)
        j = (pl.e - mi.e) / (2*delta)
        g = flex.sum(j * sh.e)
        print "finite diff[%d][%d]: " % (i_scatterer, i_xyz), g
      direct_summ = sfd.d_target_d_site()[i_scatterer][i_xyz]#* sum_f_obs_sq/2
      if (top_gradient is None): top_gradient = direct_summ
      print "direct summ[%d][%d]: " % (i_scatterer, i_xyz), direct_summ,
      if (show_finite):
        print judge(scatterer, "site", g, direct_summ)
      else:
        print
      if (f_obs.space_group().order_z() == 1):
        rm = matrix.col(sh.structure_shifted.scatterers()[i_scatterer].site)
        gxm = 0
        for i,hkl in f_obs.indices().items():
          p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
          gxm += (dps[i_xyz].data()[i]
                  * gms[i_scatterer][i].conjugate()
                  * complex_math.polar((1, p)))
        print "        gxm[%d][%d]:" % (i_scatterer, i_xyz), gxm
      fast_gradie = -map0.d_target_d_site()[i_scatterer][i_xyz] \
                  * f_obs.space_group().n_ltr()
      print "fast gradie[%d][%d]: " % (i_scatterer, i_xyz), fast_gradie,
      print judge(scatterer, "site", direct_summ, fast_gradie,
                  top=top_gradient)
      print
  sys.stdout.flush()
  if (0):
    for i_xyz in (0,1,2):
      for i_scatterer in (0,1,2):
        d = sfd.d_target_d_site()[i_scatterer][i_xyz]
        f = -map0.d_target_d_site()[i_scatterer][i_xyz] \
          * f_obs.space_group().n_ltr()
        print "direct: %12.5g" % d,
        print "  fast: %12.5g" % f,
        print " ratio: %12.5g" % (f/d)
    print
    sys.stdout.flush()

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

def u_iso(structure_ideal, d_min, f_obs, show_finite):
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
  if (f_obs.space_group().order_z() == 1):
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
  top_gradient = None
  for i_scatterer in (0,1,2):
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    if (show_finite):
      delta = 1.e-6
      pl = two_p_shifted_u_iso(
        f_obs, sh.structure_shifted, i_scatterer, delta)
      mi = two_p_shifted_u_iso(
        f_obs, sh.structure_shifted, i_scatterer, -delta)
      j = (pl.e - mi.e) / (2*delta)
      g = flex.sum(j * sh.e)
      print "finite diff[%d]: " % i_scatterer, g
    direct_summ = sfd.d_target_d_u_iso()[i_scatterer] * sum_f_obs_sq/2
    if (top_gradient is None): top_gradient = direct_summ
    print "direct summ[%d]: " % i_scatterer, direct_summ,
    if (show_finite):
      print judge(scatterer, "u_iso", g, direct_summ)
    else:
      print
    if (f_obs.space_group().order_z() == 1):
      rm = matrix.col(sh.structure_shifted.scatterers()[i_scatterer].site)
      gxm = 0
      for i,hkl in f_obs.indices().items():
        p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
        gxm += (dps.data()[i]
                * gms[i_scatterer][i].conjugate()
                * complex_math.polar((1, p)))
      print "        gxm[%d]:" % i_scatterer, gxm
    fast_gradie = map0.d_target_d_u_iso()[i_scatterer] \
                * f_obs.space_group().n_ltr()
    print "fast gradie[%d]: " % i_scatterer, fast_gradie,
    print judge(scatterer, "u_iso", direct_summ, fast_gradie,
                top=top_gradient)
    print
  sys.stdout.flush()

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
    if (0):
      f_calc_fft = f_obs.structure_factors_from_scatterers(
        xray_structure=structure_shifted, fft=0001).f_calc()
      r = xray.targets_least_squares_residual(
        f_calc_abs.data(), f_calc_fft.data(), 0001, 1)
      print "shifted k,r: %.4f, %.4f" % (r.scale_factor(), r.target())

def occupancy(structure_ideal, d_min, f_obs, show_finite):
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
  if (f_obs.space_group().order_z() == 1):
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
  top_gradient = None
  for i_scatterer in (0,1,2):
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    if (show_finite):
      delta = 1.e-6
      pl = two_p_shifted_occupancy(
        f_obs, sh.structure_shifted, i_scatterer, delta)
      mi = two_p_shifted_occupancy(
        f_obs, sh.structure_shifted, i_scatterer, -delta)
      j = (pl.e - mi.e) / (2*delta)
      g = flex.sum(j * sh.e)
      print "finite diff[%d]: " % i_scatterer, g
    direct_summ = sfd.d_target_d_occupancy()[i_scatterer] * sum_f_obs_sq/2
    if (top_gradient is None): top_gradient = direct_summ
    print "direct summ[%d]: " % i_scatterer, direct_summ,
    if (show_finite):
      print judge(scatterer, "occupancy", g, direct_summ)
    else:
      print
    if (f_obs.space_group().order_z() == 1):
      m = sh.structure_shifted.scatterers()[i_scatterer]
      rm = matrix.col(m.site)
      gxm = 0
      for i,hkl in f_obs.indices().items():
        p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
        gxm += (dps.data()[i]
                * gms[i_scatterer][i].conjugate()
                * complex_math.polar((1, p)))
      print "        gxm[%d]:" % i_scatterer, gxm/m.occupancy
    fast_gradie = map0.d_target_d_occupancy()[i_scatterer] \
                * f_obs.space_group().n_ltr()
    print "fast gradie[%d]: " % i_scatterer, fast_gradie,
    print judge(scatterer, "occupancy", direct_summ, fast_gradie,
                top=top_gradient)
    print
  sys.stdout.flush()

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

def u_star(structure_ideal, d_min, f_obs, show_finite):
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
  if (f_obs.space_group().order_z() == 1):
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
  top_gradient = None
  for i_scatterer in (0,1,2):
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    sfd_star = [x*sum_f_obs_sq/2 for x in sfd.d_target_d_u_star()[i_scatterer]]
    sfd_cart = adptbx.grad_u_star_as_u_cart(
      structure_ideal.unit_cell(), sfd_star)
    assert approx_equal(
      sfd_star,
      adptbx.grad_u_cart_as_u_star(structure_ideal.unit_cell(), sfd_cart))
    for ij in xrange(6):
      if (show_finite):
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
        print "finite diff[%d][%d]: " % (i_scatterer, ij), g
      direct_summ = sfd.d_target_d_u_star()[i_scatterer][ij] * sum_f_obs_sq/2
      if (top_gradient is None): top_gradient = direct_summ
      print "direct summ[%d][%d]: " % (i_scatterer, ij), direct_summ,
      if (show_finite):
        print judge(scatterer, "u_star", g, direct_summ)
      else:
        print
      if (f_obs.space_group().order_z() == 1):
        rm = matrix.col(sh.structure_shifted.scatterers()[i_scatterer].site)
        gxm = 0
        for i,hkl in f_obs.indices().items():
          p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
          gxm += (dps[ij].data()[i]
                  * gms[i_scatterer][i].conjugate()
                  * complex_math.polar((1, p)))
        print "        gxm[%d][%d]:" % (i_scatterer, ij), gxm
      fast_gradie = map0.d_target_d_u_star()[i_scatterer][ij] \
                  * f_obs.space_group().n_ltr()
      print "fast gradie[%d][%d]: " % (i_scatterer, ij), fast_gradie,
      print judge(scatterer, "u_star", direct_summ, fast_gradie,
                  top=top_gradient)
      if (0):
        print "         gc[%d][%d]: " % (i_scatterer, ij), gc
        print "        s2c[%d][%d]: " % (i_scatterer, ij), sfd_cart[ij]
      print
  sys.stdout.flush()

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

def fp(structure_ideal, d_min, f_obs, show_finite):
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
  if (f_obs.space_group().order_z() == 1):
    dps = flex.complex_double()
    for i,hkl in f_obs.indices().items():
      dps.append(sh.e[i]
                 *complex_math.polar((1, phi[i])))
    dps = miller.array(miller_set=f_obs, data=dps)
    uc = sh.structure_shifted.unit_cell()
  print "two_p:", sh.two_p
  print "ls.target():", ls.target() * sum_f_obs_sq
  assert approx_equal(sh.two_p, ls.target() * sum_f_obs_sq)
  print
  re = resampling(miller_set=f_obs)
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  top_gradient = None
  for i_scatterer in (0,1,2):
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    if (show_finite):
      delta = 1.e-6
      pl = two_p_shifted_fp(
        f_obs, sh.structure_shifted, i_scatterer, delta)
      mi = two_p_shifted_fp(
        f_obs, sh.structure_shifted, i_scatterer, -delta)
      j = (pl.e - mi.e) / (2*delta)
      g = flex.sum(j * sh.e)
      print "finite diff[%d]: " % i_scatterer, g
    direct_summ = sfd.d_target_d_fp()[i_scatterer] * sum_f_obs_sq/2
    if (top_gradient is None): top_gradient = direct_summ
    print "direct summ[%d]: " % i_scatterer, direct_summ,
    if (show_finite):
      print judge(scatterer, "fp", g, direct_summ)
    else:
      print
    m = sh.structure_shifted.scatterers()[i_scatterer]
    if (f_obs.space_group().order_z() == 1):
      rm = matrix.col(m.site)
      gxm = 0
      for i,hkl in f_obs.indices().items():
        p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
        f0_fp_fdp = m.caasf.at_d_star_sq(uc.d_star_sq(hkl))+m.fp_fdp
        gxm += (dps.data()[i]
                * (gms[i_scatterer][i]/f0_fp_fdp).conjugate()
                * complex_math.polar((1, p)))
      print "        gxm[%d]:" % i_scatterer, gxm
    fast_gradie = map0.d_target_d_fp()[i_scatterer] \
                * f_obs.space_group().n_ltr()
    print "fast gradie[%d]: " % i_scatterer, fast_gradie,
    print judge(scatterer, "fp", direct_summ, fast_gradie,
                top=top_gradient)
    print
  sys.stdout.flush()

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

def fdp(structure_ideal, d_min, f_obs, show_finite):
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
  if (f_obs.space_group().order_z() == 1):
    dps = flex.complex_double()
    for i,hkl in f_obs.indices().items():
      dps.append(sh.e[i]
                 *complex_math.polar((1, phi[i])))
    dps = miller.array(miller_set=f_obs, data=dps)
    uc = sh.structure_shifted.unit_cell()
  print "two_p:", sh.two_p
  print "ls.target():", ls.target() * sum_f_obs_sq
  assert approx_equal(sh.two_p, ls.target() * sum_f_obs_sq)
  print
  re = resampling(miller_set=f_obs)
  map0 = re(xray_structure=sh.structure_shifted, dp=dp0)
  top_gradient = None
  for i_scatterer in (0,1,2):
    scatterer = sh.structure_shifted.scatterers()[i_scatterer]
    if (show_finite):
      delta = 1.e-6
      pl = two_p_shifted_fdp(
        f_obs, sh.structure_shifted, i_scatterer, delta)
      mi = two_p_shifted_fdp(
        f_obs, sh.structure_shifted, i_scatterer, -delta)
      j = (pl.e - mi.e) / (2*delta)
      g = flex.sum(j * sh.e)
      print "finite diff[%d]: " % i_scatterer, g
    direct_summ = sfd.d_target_d_fdp()[i_scatterer] * sum_f_obs_sq/2
    if (top_gradient is None): top_gradient = direct_summ
    print "direct summ[%d]: " % i_scatterer, direct_summ,
    if (show_finite):
      print judge(scatterer, "fdp", g, direct_summ)
    else:
      print
    m = sh.structure_shifted.scatterers()[i_scatterer]
    if (f_obs.space_group().order_z() == 1):
      rm = matrix.col(m.site)
      gxm = 0
      for i,hkl in f_obs.indices().items():
        p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
        f0_fp_fdp = m.caasf.at_d_star_sq(uc.d_star_sq(hkl))+m.fp_fdp
        gxm += (dps.data()[i]
                * (gms[i_scatterer][i]/(-1j*f0_fp_fdp)).conjugate()
                * complex_math.polar((1, p)))
      print "        gxm[%d]:" % i_scatterer, gxm
    fast_gradie = map0.d_target_d_fdp()[i_scatterer] \
                * f_obs.space_group().n_ltr()
    print "fast gradie[%d]: " % i_scatterer, fast_gradie,
    print judge(scatterer, "fdp", direct_summ, fast_gradie,
                top=top_gradient)
    print
  sys.stdout.flush()

def run_one(space_group_info, n_elements=3, volume_per_atom=1000, d_min=2,
            fdp_flag=0, anisotropic_flag=0, show_finite=0):
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=("Se",)*n_elements,
    volume_per_atom=volume_per_atom,
    min_distance=5,
    general_positions_only=1,
    random_f_prime_d_min=d_min-1,
    #random_f_prime_scale=0,# XXX
    random_f_double_prime=fdp_flag,
    anisotropic_flag=anisotropic_flag,
    #random_u_iso=0, # XXX
    #u_iso=adptbx.b_as_u(20), # XXX
    #random_u_iso_scale=0, # XXX
    random_u_cart_scale=.3,
    random_occupancy=0)
  if (0):
    a = structure_ideal.unit_cell().volume()**(1/3.)
    r = 1.0
    structure_ideal = xray.structure(
      special_position_settings=crystal.special_position_settings(
        crystal_symmetry=crystal.symmetry(
          unit_cell=(r*a,a,a/r,90,90,90),
          space_group=structure_ideal.space_group())),
      scatterers=structure_ideal.scatterers())
  structure_ideal.show_summary().show_scatterers()
  uc = structure_ideal.unit_cell()
  print "volume:", uc.volume()
  print "orthogonalization:", uc.orthogonalization_matrix()
  if (anisotropic_flag):
    for scatterer in structure_ideal.scatterers():
      print "u_iso:", adptbx.u_star_as_u_iso(uc, scatterer.u_star)
  print
  f_obs_object = structure_ideal.structure_factors(
    d_min=d_min, anomalous_flag=0001, direct=0001, fft=00000)
  if (0):
    print f_obs_object.manager().rfft().n_real()
    print f_obs_object.manager().u_extra()
  f_obs = abs(f_obs_object.f_calc())
  if (0 and not anisotropic_flag):
    make_cns_input.write("tmp.cns",
      make_cns_input.script_predict_methods_comparison(
        d_min, structure_ideal, f_obs))
  if (0):
    f_obs_fft_object = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_ideal, direct=00000, fft=0001,
      quality_factor=10000, b_extra=None)
    f_obs_fft = f_obs_fft_object.f_calc()
    print f_obs_fft_object.manager().rfft().n_real()
    print f_obs_fft_object.manager().u_extra()
    r = xray.targets_least_squares_residual(
      f_obs.data(), f_obs_fft.data(), 00000, 1)
    print "ideal k,r: %.4f, %.4f" % (r.scale_factor(), r.target())
    from cctbx.utils import phase_error
    from scitbx.python_utils.complex_math import arg
    max_d_ampl = 0
    max_d_phase = 0
    for i,h in f_obs.indices().items():
      x,y = f_obs_object.f_calc().data()[i], f_obs_fft.data()[i]
      print h,x
      s = " "*len(str(h))
      print s,y
      p = phase_error(arg(x,deg=1),arg(y,deg=1),deg=1)
      print s, "d-ampl:", abs(x-y),
      if (abs(x-y) > max_d_ampl):
        max_d_ampl = abs(x-y)
        print "max_d_ampl",
      print
      print s, "d-phase:", p,
      if (abs(x-y) > max_d_phase):
        max_d_phase = abs(x-y)
        print "max_d_phase",
      print
    #return
  if (1):
    print "site"
    site(structure_ideal, d_min, f_obs, show_finite)
  if (1):
    if (not anisotropic_flag):
      print "u_iso"
      u_iso(structure_ideal, d_min, f_obs, show_finite)
    else:
      print "u_star"
      u_star(structure_ideal, d_min, f_obs, show_finite)
  if (1):
    print "occupancy"
    occupancy(structure_ideal, d_min, f_obs, show_finite)
  if (1):
    print "fp"
    fp(structure_ideal, d_min, f_obs, show_finite)
  if (1):
    print "fdp"
    fdp(structure_ideal, d_min, f_obs, show_finite)

def run_call_back(flags, space_group_info):
  for fdp_flag in [0,1]:
    for anisotropic_flag in [0,1]:
      run_one(
        space_group_info=space_group_info,
        fdp_flag=fdp_flag,
        anisotropic_flag=anisotropic_flag)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
