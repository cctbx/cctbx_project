from cctbx.xray.structure_factors.misc import quality_factor_from_any
from cctbx.xray import ext
from cctbx import maptbx
from cctbx import crystal
from scitbx import fftpack
from scitbx.python_utils.misc import adopt_init_args

class manager(crystal.symmetry):

  def __init__(self, miller_set=None,
                     crystal_symmetry=None,
                     d_min=None,
                     grid_resolution_factor=1/3.,
                     symmetry_flags=maptbx.use_space_group_symmetry,
                     mandatory_grid_factors=None,
                     quality_factor=None, u_extra=None, b_extra=None,
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
    quality_factor = quality_factor_from_any(
      d_min, grid_resolution_factor, quality_factor, u_extra, b_extra)
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
      self._u_extra = ext.calc_u_extra(
        self.d_min(),
        self.grid_resolution_factor(),
        self.quality_factor())
    return self._u_extra

  def setup_fft(self):
    self.crystal_gridding_tags()
    self.rfft()
    self.u_extra()
    return self

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

class managed_calculation_base:

  def __init__(self, manager, xray_structure, miller_set):
    adopt_init_args(self, locals(), hide=0001)
    assert xray_structure is not None and miller_set is not None
    assert xray_structure.unit_cell().is_similar_to(miller_set.unit_cell())
    assert xray_structure.space_group() == miller_set.space_group()
    if (manager is not None):
      assert xray_structure.unit_cell().is_similar_to(manager.unit_cell())
      assert xray_structure.space_group() == manager.space_group()

  def manager(self):
    return self._manager

  def xray_structure(self):
    return self._xray_structure

  def miller_set(self):
    return self._miller_set
