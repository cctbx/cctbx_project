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
                     lifchitz,
                     d_target_d_f_calc=None,
                     derivative_flags=None,
                     force_complex=00000,
                     electron_density_must_be_positive=0001):
    self.setup_fft()
    time_sampling = user_plus_sys_time()
    result = xray.agarwal_1978(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      self.ft_dp(dp).complex_map(),
      lifchitz,
      self.rfft().n_real(),
      self.rfft().m_real(),
      self.u_extra(),
      self.wing_cutoff(),
      self.exp_table_one_over_step_size(),
      force_complex,
      electron_density_must_be_positive)
    time_sampling = time_sampling.elapsed()
    return result

class two_p_shifted:

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

def run(n_elements=3, volume_per_atom=1000, d_min=2):
  structure_ideal = random_structure.xray_structure(
    sgtbx.space_group_info("P 1"),
    elements=("Se",)*n_elements,
    volume_per_atom=volume_per_atom,
    random_f_prime_d_min=0,
    random_f_double_prime=00000,
    anisotropic_flag=00000,
    random_u_iso=00000,
    u_iso=0.4,
    random_occupancy=00000)
  structure_ideal.show_summary().show_scatterers()
  print
  f_obs = abs(structure_ideal.structure_factors(
    d_min=d_min, anomalous_flag=0001, direct=0001).f_calc())
  sum_f_obs_sq = flex.sum(flex.pow2(f_obs.data()))
  sh = two_p_shifted(f_obs, structure_ideal, 0, 0, 0.05)
  sh.structure_shifted.show_summary().show_scatterers()
  ls = xray.targets_least_squares_residual(
    f_obs.data(), sh.f_calc.data(), 0001, 1)
  sfd = xray.structure_factors.from_scatterers_direct(
    xray_structure=sh.structure_shifted,
    miller_set=f_obs,
    d_target_d_f_calc=ls.derivatives(),
    derivative_flags=xray.structure_factors.derivative_flags(
      site=0001))
  phi = flex.arg(sh.f_calc.data())
  uc = sh.structure_shifted.unit_cell()
  gms = []
  for m in sh.structure_shifted.scatterers():
    gm = flex.double()
    for i,hkl in f_obs.indices().items():
      d = adptbx.debye_waller_factor_u_iso(uc, hkl, m.u_iso)
      f = m.caasf.at_d_star_sq(uc.d_star_sq(hkl))
      gm.append(d*f)
    gms.append(gm)
  dp0 = flex.complex_double()
  for i,hkl in f_obs.indices().items():
    dp0.append(-sh.e[i]
               *complex_math.polar((1, phi[i])))
  dp0 = miller.array(miller_set=f_obs, data=dp0)
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
  map0 = re(xray_structure=sh.structure_shifted,
            dp=dp0, lifchitz=0001)
  maps = []
  for i_xyz in (0,1,2):
    map = re(xray_structure=sh.structure_shifted,
             dp=dps[i_xyz], lifchitz=00000)
    maps.append(map)
  for i_scatterer in (0,1,2):
    for i_xyz in (0,1,2):
      delta = 1.e-6
      pl = two_p_shifted(f_obs, sh.structure_shifted, i_scatterer,i_xyz,delta)
      mi = two_p_shifted(f_obs, sh.structure_shifted, i_scatterer,i_xyz,-delta)
      j = (pl.e - mi.e) / (2*delta)
      g = flex.sum(j * sh.e)
      print "  g[%d][%d]: " % (i_scatterer, i_xyz), g
      rm = matrix.col(sh.structure_shifted.scatterers()[i_scatterer].site)
      gxm = 0
      for i,hkl in f_obs.indices().items():
        p = -2*math.pi * (matrix.row(hkl) * rm).elems[0]
        gxm += (dps[i_xyz].data()[i]
                * gms[i_scatterer][i]
                * complex_math.polar((1, p)))
      print "gxm[%d][%d]:" % (i_scatterer, i_xyz), gxm
      print "sfd[%d][%d]: " % (i_scatterer, i_xyz), \
            sfd.d_target_d_site()[i_scatterer][i_xyz] * sum_f_obs_sq/2
      m = maps[i_xyz].grad()[i_scatterer]
      print "map[%d][%d]:" % (i_scatterer, i_xyz), m
      gl = (map0.grad_x, map0.grad_y, map0.grad_z)[i_xyz]()[i_scatterer]
      print " m0[%d][%d]:" % (i_scatterer, i_xyz), gl
      print m.real / g

if (__name__ == "__main__"):
  run()
