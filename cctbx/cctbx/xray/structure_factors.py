from cctbx import miller
from cctbx import maptbx
from cctbx import crystal
from cctbx.array_family import flex
from scitbx import fftpack
from scitbx.python_utils.misc import adopt_init_args, user_plus_sys_time

class derivative_flags:

  def __init__(self, site=00000,
                     u_iso=00000,
                     u_star=00000,
                     occupancy=00000,
                     fp=00000,
                     fdp=00000):
    adopt_init_args(self, locals())

class from_scatterers(crystal.symmetry):

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
    assert miller_set == None or crystal_symmetry == None
    if (miller_set == None):
      assert crystal_symmetry != None and d_min != None
    else:
      crystal_symmetry = miller_set
      if (d_min == None):
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
    self.estimate_time_direct = _estimate_time_direct(
      self.space_group().order_z())
    self.estimate_time_fft = _estimate_time_fft()

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
    if (self._crystal_gridding == None):
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
    if (self._crystal_gridding_tags == None):
      self._crystal_gridding_tags = self.crystal_gridding(
        assert_shannon_sampling).tags()
    return self._crystal_gridding_tags

  def rfft(self):
    if (self._rfft == None):
      self._rfft = fftpack.real_to_complex_3d(self.crystal_gridding().n_real())
    return self._rfft

  def u_extra(self):
    if (self._u_extra == None):
      from cctbx import xray
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

  def __call__(self, xray_structure,
                     miller_set,
                     d_target_d_f_calc=None,
                     derivative_flags=None,
                     direct=00000,
                     fft=00000):
    assert direct == 00000 or fft == 00000
    if (direct == 00000 and fft == 00000):
      n_scatterers = xray_structure.scatterers().size()
      n_miller_indices = miller_set.indices().size()
      if (not self.have_good_timing_estimates()):
        # rough estimate
        if (  n_scatterers * self.space_group().order_z() * n_miller_indices
            < self.crystal_gridding().n_grid_points()):
          direct = 0001
      else:
        if (   self.estimate_time_direct(n_scatterers * n_miller_indices)
            <= self.estimate_time_fft(n_scatterers, n_miller_indices)):
          direct = 0001
    if (direct): f = from_scatterers_direct
    else:        f = from_scatterers_fft
    return f(
      manager=self,
      xray_structure=xray_structure,
      miller_set=miller_set,
      d_target_d_f_calc=d_target_d_f_calc,
      derivative_flags=derivative_flags)

  def have_good_timing_estimates(self):
    return self.estimate_time_direct.have_good_estimate() \
       and self.estimate_time_fft.have_good_estimate()

class _estimate_time_direct:

  def __init__(self, order_z, min_product=100000):
    self.order_z = order_z
    self.min_product = min_product
    self.product = 1
    self.time = 0

  def have_good_estimate(self):
    return self.product > self.min_product / float(self.order_z)

  def register(self, product, time):
    if (product >= self.product):
      self.product = product
      self.time = time

  def __call__(self, product):
    return self.time * product / self.product

def _linear_estimate(x1, x2, y1, y2, x):
  if (y1 == y2): return y1
  if (x1 == x2): return max(y1,y2) * x1 / x
  slope = (y2 - y1) / (x2 - x1)
  if   (x1 == 0): y_intercept = y1
  elif (x2 == 0): y_intercept = y2
  else: y_intercept = y1 - (slope * x1)
  return slope * x + y_intercept

class _estimate_time_fft:

  def __init__(self):
    self.min_n_scatterers = 0
    self.min_time_sampling = 0
    self.max_n_scatterers = 0
    self.max_time_sampling = 0
    self.min_n_miller_indices = 0
    self.min_time_collect = 0
    self.max_n_miller_indices = 0
    self.max_time_collect = 0
    self.time_sampling = 0
    self.time_symmetry_mapping = 0
    self.time_fft = 0
    self.time_collect = 0

  def have_good_estimate(self):
    return self.time_fft > 0

  def register(self, n_scatterers, n_miller_indices,
                     time_sampling, time_symmetry_mapping,
                     time_fft, time_collect):
    if (   self.min_n_scatterers == 0
        or self.min_n_scatterers >= n_scatterers):
      self.min_n_scatterers = n_scatterers
      self.min_time_sampling = time_sampling
    if (self.max_n_scatterers <= n_scatterers):
      self.max_n_scatterers = n_scatterers
      self.max_time_sampling = time_sampling
    if (   self.min_n_miller_indices == 0
        or self.min_n_miller_indices >= n_miller_indices):
      self.min_n_miller_indices = n_miller_indices
      self.min_time_collect = time_collect
    if (self.max_n_miller_indices <= n_miller_indices):
      self.max_n_miller_indices = n_miller_indices
      self.max_time_collect = time_collect
    self.time_sampling = time_sampling
    self.time_symmetry_mapping = time_symmetry_mapping
    self.time_fft = time_fft
    self.time_collect = time_collect

  def __call__(self, n_scatterers, n_miller_indices):
    return _linear_estimate(
             self.min_n_scatterers, self.max_n_scatterers,
             self.min_time_sampling, self.max_time_sampling,
             n_scatterers) \
         + self.time_symmetry_mapping \
         + self.time_fft \
         + _linear_estimate(
             self.min_n_miller_indices, self.max_n_miller_indices,
             self.min_time_collect, self.max_time_collect,
             n_miller_indices)

class _from_scatterers_base:

  def __init__(self, manager, xray_structure, miller_set):
    adopt_init_args(self, locals(), hide=0001)
    assert xray_structure != None and miller_set != None
    assert xray_structure.unit_cell().is_similar_to(miller_set.unit_cell())
    assert xray_structure.space_group() == miller_set.space_group()
    if (manager != None):
      assert xray_structure.unit_cell().is_similar_to(manager.unit_cell())
      assert xray_structure.space_group() == manager.space_group()

  def manager(self):
    return self._manager

  def xray_structure(self):
    return self._xray_structure

  def miller_set(self):
    return self._miller_set

class from_scatterers_direct(_from_scatterers_base):

  def __init__(self, manager=None,
                     xray_structure=None,
                     miller_set=None,
                     d_target_d_f_calc=None,
                     derivative_flags=None):
    _from_scatterers_base.__init__(self, manager, xray_structure, miller_set)
    self._d_target_d_f_calc = d_target_d_f_calc
    if (d_target_d_f_calc == None):
      d_target_d_f_calc = flex.complex_double()
    if (derivative_flags == None):
      derivative_flags = globals()["derivative_flags"]()
    timer = user_plus_sys_time()
    from cctbx import xray
    self._results = xray.structure_factors_direct_with_first_derivatives(
      self._miller_set.unit_cell(),
      self._miller_set.space_group(),
      self._miller_set.indices(),
      self._xray_structure.scatterers(),
      d_target_d_f_calc,
      derivative_flags.site,
      derivative_flags.u_iso,
      derivative_flags.u_star,
      derivative_flags.occupancy,
      derivative_flags.fp,
      derivative_flags.fdp)
    if (manager != None):
      manager.estimate_time_direct.register(
        xray_structure.scatterers().size() * miller_set.indices().size(),
        timer.elapsed())

  def f_calc(self):
    return miller.array(self._miller_set, self._results.f_calc())

  def d_target_d_f_calc(self):
    return self._d_target_d_f_calc

  def d_target_d_site(self):
    d_target_d_site = self._results.d_target_d_site()
    xray_structure = self.xray_structure()
    assert d_target_d_site.size() == xray_structure.scatterers().size()
    return d_target_d_site

  def d_target_d_u_iso(self):
    d_target_d_u_iso = self._results.d_target_d_u_iso()
    xray_structure = self.xray_structure()
    assert d_target_d_u_iso.size() == xray_structure.scatterers().size()
    return d_target_d_u_iso

  def d_target_d_u_star(self):
    d_target_d_u_star = self._results.d_target_d_u_star()
    xray_structure = self.xray_structure()
    assert d_target_d_u_star.size() == xray_structure.scatterers().size()
    return d_target_d_u_star

  def d_target_d_occupancy(self):
    d_target_d_occupancy = self._results.d_target_d_occupancy()
    xray_structure = self.xray_structure()
    assert d_target_d_occupancy.size() == xray_structure.scatterers().size()
    return d_target_d_occupancy

  def d_target_d_fp(self):
    d_target_d_fp = self._results.d_target_d_fp()
    xray_structure = self.xray_structure()
    assert d_target_d_fp.size() == xray_structure.scatterers().size()
    return d_target_d_fp

  def d_target_d_fdp(self):
    d_target_d_fdp = self._results.d_target_d_fdp()
    xray_structure = self.xray_structure()
    assert d_target_d_fdp.size() == xray_structure.scatterers().size()
    return d_target_d_fdp

  def d_target_d_site_inplace_frac_as_cart(self, d_target_d_site):
    from cctbx import xray
    xray.structure_factors_d_target_d_site_in_place_frac_as_cart(
      self.miller_set().unit_cell(), d_target_d_site)

class from_scatterers_fft(_from_scatterers_base):

  def __init__(self, manager,
                     xray_structure,
                     miller_set,
                     d_target_d_f_calc=None,
                     derivative_flags=None,
                     force_complex=00000,
                     electron_density_must_be_positive=0001):
    _from_scatterers_base.__init__(self, manager, xray_structure, miller_set)
    assert manager.symmetry_flags().use_space_group_symmetry()
    assert miller_set.d_min() >= manager.d_min()
    assert d_target_d_f_calc == None, "FFT derivatives not implemented."
    assert derivative_flags == None, "FFT derivatives not implemented."
    manager.setup_fft() # before timing
    time_sampling = user_plus_sys_time()
    from cctbx import xray
    sampled_density = xray.sampled_model_density(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      manager.rfft().n_real(),
      manager.rfft().m_real(),
      manager.u_extra(),
      manager.wing_cutoff(),
      manager.exp_table_one_over_step_size(),
      force_complex,
      electron_density_must_be_positive)
    time_sampling = time_sampling.elapsed()
    time_symmetry_mapping = user_plus_sys_time()
    sampled_density.apply_symmetry(manager.crystal_gridding_tags().tags())
    time_symmetry_mapping = time_symmetry_mapping.elapsed()
    time_fft = user_plus_sys_time()
    if (not sampled_density.anomalous_flag()):
      map = sampled_density.real_map()
      sf_map = manager.rfft().forward(map)
      collect_conj = 1
    else:
      cfft = fftpack.complex_to_complex_3d(manager.rfft().n_real())
      map = sampled_density.complex_map()
      sf_map = cfft.backward(map)
      collect_conj = 0
    time_fft = time_fft.elapsed()
    time_collect = user_plus_sys_time()
    self._f_calc_data = maptbx.structure_factors.from_map(
      sampled_density.anomalous_flag(),
      miller_set.indices(),
      sf_map,
      collect_conj).data()
    sampled_density.eliminate_u_extra_and_normalize(
      miller_set.indices(),
      self._f_calc_data)
    time_collect = time_collect.elapsed()
    manager.estimate_time_fft.register(
      xray_structure.scatterers().size(),
      miller_set.indices().size(),
      time_sampling, time_symmetry_mapping, time_fft, time_collect)

  def f_calc(self):
    return miller.array(self.miller_set(), self._f_calc_data)
