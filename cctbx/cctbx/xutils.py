from cctbx_boost.arraytbx import shared
from cctbx_boost import sgtbx
from cctbx_boost import miller
from cctbx_boost import sftbx

# XXX move to uctbx
def are_similar_unit_cells(ucell1, ucell2,
                           length_tolerance=0.01,
                           angle_tolerance=1.):
  p1 = ucell1.getParameters()
  p2 = ucell2.getParameters()
  for i in xrange(3):
    if (min(p1[i], p2[i]) / max(p1[i], p2[i]) > (1 + length_tolerance)):
      return False
  for i in xrange(3, 6):
    if (abs(p1[i] - p2[i]) > angle_tolerance):
      return False
  return True

def space_group_info(space_group_symbol):
  return sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(space_group_symbol)).Info()

class crystal_symmetry:

  def __init__(self, UnitCell, SgInfo, auto_check=0):
    self.UnitCell = UnitCell
    self.SgInfo = SgInfo
    self.SgOps = self.SgInfo.SgOps() # FUTURE: avoid the copy?
    if (auto_check):
      self.check_unit_cell()

  def check_unit_cell(self):
    self.SgOps.CheckUnitCell(self.UnitCell)

  def cell_equivalent_p1(self):
    return crystal_symmetry(self.UnitCell, sgtbx.SpaceGroupInfo())

class miller_set(crystal_symmetry):

  def __init__(self, xsym, H):
    crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.H = H

  def expand_to_p1(self, friedel_flag):
    set_p1_H = shared.miller_Index()
    miller.expand_to_p1(self.SgOps, friedel_flag, self.H, set_p1_H)
    return miller_set(self.cell_equivalent_p1(), set_p1_H)

class reciprocal_space_array(miller_set):

  def __init__(self, miller_set_obj, F, sigmas=0):
    miller_set.__init__(self, miller_set_obj, miller_set_obj.H)
    self.F = F
    self.sigmas = sigmas

  def anomalous_differences(self):
    if (hasattr(self, "SgInfo")):
      jbm = miller.join_bijvoet_mates(self.SgInfo, self.H)
    else:
      jbm = miller.join_bijvoet_mates(self.H)
    h = jbm.miller_indices_in_hemisphere("+")
    f = jbm.minus(self.F)
    s = 0
    if (self.sigmas):
      s = jbm.additive_sigmas(self.sigmas)
    return reciprocal_space_array(miller_set(self, h), f, s)

  def sigma_filter(self, cutoff_factor, negate=0):
    sel = miller.selection(self.H)
    sel.sigma_filter(self.F, self.sigmas, cutoff_factor)
    if (negate): sel.negate()
    h = sel.selected_miller_indices()
    f = sel.selected_data(self.F)
    s = sel.selected_data(self.sigmas)
    return reciprocal_space_array(miller_set(self, h), f, s)

  def rms_filter(self, cutoff_factor, negate=0):
    rms = shared.rms(self.F)
    pass # XXX

  def __add__(self, other):
    if (type(other) != type(self)):
      # add a scalar
      # XXX push to C++
      f = self.F.deep_copy()
      for i in f.indices(): f[i] += other
      s = 0
      if (self.sigmas):
        s = self.sigmas.deep_copy()
        for i in s.indices(): s[i] += other
      return reciprocal_space_array(self, f, s)
    # add arrays
    js = miller.join_sets(self.H, other.H)
    h = js.common_miller_indices()
    f = js.plus(self.F, other.F)
    s = 0
    if (self.sigmas):
      s = js.additive_sigmas(self.sigmas, other.sigmas)
    return reciprocal_space_array(miller_set(self, h), f, s)

class symmetrized_sites(crystal_symmetry):

  def __init__(self, xsym, sites = None,
               MinMateDistance = 0.5,
               Ustar_tolerance = 0.1,
               TestPositiveDefiniteness = 1):
    crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.MinMateDistance = MinMateDistance
    self.Ustar_tolerance = Ustar_tolerance
    self.TestPositiveDefiniteness = TestPositiveDefiniteness
    self.Sites = shared.XrayScatterer()
    self.SpecialPositionOps = shared.RTMx()
    self.SnapParameters = sgtbx.SpecialPositionSnapParameters(
      self.UnitCell, self.SgOps, 1, self.MinMateDistance)
    if (sites):
      self.add_sites(sites)

  def copy_attributes(self):
    return symmetrized_sites(self, sites = None,
      MinMateDistance = self.MinMateDistance,
      Ustar_tolerance = self.Ustar_tolerance,
      TestPositiveDefiniteness = self.TestPositiveDefiniteness)

  def n_special_positions(self):
    n = 0
    for special_op in self.SpecialPositionOps:
      if (not special_op.isUnit()):
        n += 1
    return n

  def get_site_symmetry(self, X):
    return sgtbx.SiteSymmetry(self.SnapParameters, X)

  def add_site(self, site):
    site_symmetry = site.ApplySymmetry(
      self.UnitCell, self.SgOps,
      self.MinMateDistance,
      self.Ustar_tolerance,
      self.TestPositiveDefiniteness)
    special_op = site_symmetry.SpecialOp()
    self.Sites.append(site)
    self.SpecialPositionOps.append(special_op)

  def add_sites(self, sites):
    for site in sites:
      self.add_site(site)

  def __len__(self): return self.size()
  def size(self): return self.Sites.size()

  def __getitem__(self, key):
    return self.Sites[key]

def build_miller_set(xsym, friedel_flag, d_min):
  result = miller_set(
    xsym, miller.BuildIndices(
      xsym.UnitCell, xsym.SgInfo, friedel_flag, d_min))
  result.friedel_flag = friedel_flag
  return result

def calculate_structure_factors(miller_indices, xtal, abs_F = 0):
  F = sftbx.StructureFactorArray(
    miller_indices.UnitCell, miller_indices.SgOps, miller_indices.H,
    xtal.Sites)
  if (abs_F): F = shared.abs(F)
  result = reciprocal_space_array(miller_indices, F)
  result.symmetrized_sites = xtal
  return result

def calculate_structure_factors_fft(miller_indices, xtal, abs_F = 0):
  from cctbx_boost import fftbx
  max_q = xtal.UnitCell.max_Q(miller_indices.H)
  grid_resolution_factor = 1/3.
  max_prime = 5
  mandatory_grid_factors = xtal.SgOps.refine_gridding()
  grid_logical = sftbx.determine_grid(
    xtal.UnitCell,
    max_q, grid_resolution_factor, max_prime, mandatory_grid_factors)
  rfft = fftbx.real_to_complex_3d(grid_logical)
  quality_factor = 100
  u_extra = sftbx.calc_u_extra(max_q, grid_resolution_factor, quality_factor)
  wing_cutoff = 1.e-3
  exp_table_one_over_step_size = -100
  force_complex = 0
  electron_density_must_be_positive = 1
  sampled_density = sftbx.sampled_model_density(
    xtal.UnitCell, xtal.Sites,
    rfft.Nreal(), rfft.Mreal(),
    u_extra, wing_cutoff, exp_table_one_over_step_size,
    force_complex, electron_density_must_be_positive)
  friedel_flag = sampled_density.friedel_flag()
  tags = sftbx.grid_tags(rfft.Nreal())
  sym_flags = sftbx.map_symmetry_flags(1)
  tags.build(xtal.SgInfo, sym_flags)
  sampled_density.apply_symmetry(tags)
  if (friedel_flag):
    map = sampled_density.map_real_as_shared()
    rfft.forward(map)
    collect_conj = 1
    n_complex = rfft.Ncomplex()
  else:
    map = sampled_density.map_complex_as_shared()
    cfft = fftbx.complex_to_complex_3d(rfft.Nreal())
    cfft.backward(map)
    collect_conj = 0
    n_complex = cfft.N()
  fcalc = sftbx.collect_structure_factors(
    friedel_flag, miller_indices.H, map, n_complex, collect_conj)
  sampled_density.eliminate_u_extra_and_normalize(miller_indices.H, fcalc)
  if (abs_F): fcalc = shared.abs(fcalc)
  result = reciprocal_space_array(miller_indices, fcalc)
  result.symmetrized_sites = xtal
  return result

def f_as_ampl_phase(f, deg=1):
  from math import hypot, atan2, pi
  if (abs(f) == 0): return (0., 0.)
  x, y = f.real, f.imag
  if deg:
    return hypot(x, y), atan2(y, x) * 180. / pi
  else:
    return hypot(x, y), atan2(y, x)

def ampl_phase_as_f(ampl_phase, deg=1):
  from math import cos, sin, pi
  a, p = ampl_phase
  if deg:
    p *= pi / 180.0
  return complex(a * cos(p), a * sin(p))
