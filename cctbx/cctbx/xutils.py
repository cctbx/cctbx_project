import math
from cctbx_boost.arraytbx import shared
from cctbx_boost.arraytbx import shared_map
from cctbx_boost import sgtbx
from cctbx_boost import miller
from cctbx_boost import sftbx
from cctbx_boost import fftbx
from cctbx.misc import python_utils

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

def show_binner_summary(binner):
  for i_bin in binner.range_all():
    bin_d_range = binner.bin_d_range(i_bin)
    count = binner.count(i_bin)
    if (i_bin == binner.i_bin_d_too_large()):
      assert bin_d_range[0] == 0
      print "unused:              d > %8.4f: %5d" % (bin_d_range[1], count)
    elif (i_bin == binner.i_bin_d_too_small()):
      assert bin_d_range[1] == 0
      print "unused: %9.4f >  d           : %5d" % (bin_d_range[0], count)
    else:
      print "bin %2d: %9.4f >= d > %8.4f: %5d" % (
        (i_bin,) + bin_d_range + (count,))

def ftor_a_s_sub(lhs, rhs): return lhs.sub(rhs)
def ftor_a_s_div(lhs, rhs): return lhs.div(rhs)

def space_group_info(space_group_symbol):
  return sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(space_group_symbol)).Info()

class crystal_symmetry:

  def __init__(self, UnitCell, SgInfo, auto_check=0):
    self.UnitCell = UnitCell
    self.SgInfo = SgInfo
    self.SgOps = self.SgInfo.SgOps() # FUTURE: avoid the copy?
    if (auto_check):
      self.check_unit_cell()

  def show_summary(self):
    print "Unit cell:", self.UnitCell
    print "Space group symbol:", self.SgInfo.BuildLookupSymbol()

  def check_unit_cell(self):
    self.SgOps.CheckUnitCell(self.UnitCell)

  def cell_equivalent_p1(self):
    return crystal_symmetry(self.UnitCell, sgtbx.SpaceGroupInfo())

class miller_set(crystal_symmetry):

  def __init__(self, xsym, H):
    crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.H = H
    if (hasattr(xsym, "friedel_flag")):
      self.friedel_flag = xsym.friedel_flag

  def friedel_info(self):
    if (not hasattr(self, "friedel_flag")): return "undefined"
    return self.friedel_flag

  def show_summary(self):
    crystal_symmetry.show_summary(self)
    print "Friedel flag:", self.friedel_info()
    print "Number of Miller indices:", self.H.size()

  def set_friedel_flag(self, flag):
    self.friedel_flag = flag
    return self

  def multiplicities(self):
    return reciprocal_space_array(
      self, self.SgOps.multiplicity(self.H, self.friedel_flag))

  def epsilons(self):
    return reciprocal_space_array(
      self, self.SgOps.epsilon(self.H))

  def d_spacings(self):
    return reciprocal_space_array(
      self, self.UnitCell.d(self.H))

  def sin_theta_over_lambda_sq(self):
    return reciprocal_space_array(
      self, self.UnitCell.stol2(self.H))

  def expand_to_p1(self):
    set_p1_H = shared.miller_Index()
    miller.expand_to_p1(self.SgOps, self.friedel_flag, self.H, set_p1_H)
    return miller_set(self.cell_equivalent_p1(), set_p1_H).set_friedel_flag(
      self.friedel_flag)

  def setup_binner(self, d_max=0, d_min=0,
                   auto_binning=0,
                   reflections_per_bin=0,
                   n_bins=0):
    assert auto_binning != 0 or reflections_per_bin != 0 or n_bins != 0
    assert auto_binning != 0 or (reflections_per_bin == 0 or n_bins == 0)
    if (auto_binning):
      if (reflections_per_bin == 0): reflections_per_bin = 200
      if (n_bins == 0): n_bins = 8
      n_per_bin = int(float(self.H.size()) / n_bins + .5)
      if (n_per_bin > reflections_per_bin):
        n_bins = int(self.H.size() / reflections_per_bin + .5)
    elif (reflections_per_bin):
      n_bins = int(self.H.size() / reflections_per_bin + .5)
    assert n_bins > 0
    binning = miller.binning(self.UnitCell, n_bins, self.H, d_max, d_min)
    self.binner = miller.binner(binning, self.H)

class reciprocal_space_array(miller_set):

  def __init__(self, miller_set_obj, F, sigmas=0, info=0):
    miller_set.__init__(self, miller_set_obj, miller_set_obj.H)
    self.F = F
    self.sigmas = sigmas
    self.info = info

  def type_of_data(self):
    if (not hasattr(self, "F")): return "none"
    try:
      return self.F.__class__.__name__ + ", size=%d" % (self.F.size())
    except:
      return "unknown"

  def type_of_sigmas(self):
    if (not hasattr(self, "sigmas")): return "none"
    try:
      return self.sigmas.__class__.__name__ + ", size=%d" % (
        self.sigmas.size())
    except:
      return "unknown"

  def show_summary(self):
    if (self.info != 0): print "Info:", self.info
    miller_set.show_summary(self)
    print "Type of data:", self.type_of_data()
    print "Type of sigmas:", self.type_of_sigmas()

  def map_to_asu(self):
    h = self.H.deep_copy()
    f = self.F.deep_copy()
    miller.map_to_asu(self.SgInfo, self.friedel_flag, h, f)
    return reciprocal_space_array(miller_set(self, h), f, self.sigmas)

  def anomalous_differences(self):
    if (hasattr(self, "SgInfo")):
      asu = self.map_to_asu()
      jbm = miller.join_bijvoet_mates(asu.SgInfo, asu.H)
      f = asu.F
    else:
      jbm = miller.join_bijvoet_mates(self.H)
      f = self.F
    h = jbm.miller_indices_in_hemisphere("+")
    f = jbm.minus(f)
    s = 0
    if (self.sigmas):
      s = jbm.additive_sigmas(self.sigmas)
    return reciprocal_space_array(
      miller_set(self, h), f, s).set_friedel_flag(1)

  def all_selection(self):
    return shared.bool(self.H.size(), 1)

  def apply_selection(self, flags, negate=0):
    if (negate): flags = ~flags
    h = self.H.select(flags)
    f = 0
    if (self.F): f = self.F.select(flags)
    s = 0
    if (self.sigmas): s = self.sigmas.select(flags)
    return reciprocal_space_array(miller_set(self, h), f, s)

  def resolution_range(self):
    min_max_d_star_sq  = self.UnitCell.min_max_Q(self.H)
    result = []
    for d_star_sq in min_max_d_star_sq:
      result.append(1 / math.sqrt(d_star_sq))
    return tuple(result)

  def resolution_filter(self, d_max=0, d_min=0, negate=0):
    d = self.d_spacings().F
    keep = self.all_selection()
    if (d_max): keep &= d <= d_max
    if (d_min): keep &= d >= d_min
    return self.apply_selection(keep, negate)

  def sigma_filter(self, cutoff_factor, negate=0):
    flags = shared.abs(self.F) >= self.sigmas.mul(cutoff_factor)
    return self.apply_selection(flags, negate)

  def generic_binner_action(self, use_binning, use_multiplicities,
                            function,
                            function_weighted):
    if (use_multiplicities):
      mult = self.multiplicities().F.as_double()
    if (not use_binning):
      if (not use_multiplicities):
        result = function(self.F)
      else:
        result = function_weighted(self.F, mult)
    else:
      result = shared.double()
      for i_bin in self.binner.range_used():
        sel = self.binner(i_bin)
        if (sel.count(1) == 0):
          result.append(0)
        else:
          sel_data = self.F.select(sel)
          if (not use_multiplicities):
            result.append(function(sel_data))
          else:
            sel_mult = mult.select(sel)
            result.append(function_weighted(sel_data, sel_mult))
    return result

  def mean(self, use_binning=0, use_multiplicities=0):
    return self.generic_binner_action(use_binning, use_multiplicities,
      shared.mean,
      shared.mean_weighted)

  def mean_sq(self, use_binning=0, use_multiplicities=0):
    return self.generic_binner_action(use_binning, use_multiplicities,
      shared.mean_sq,
      shared.mean_sq_weighted)

  def rms(self, use_binning=0, use_multiplicities=0):
    ms = self.mean_sq(use_binning, use_multiplicities)
    if (not use_binning):
      return math.sqrt(ms)
    else:
      return shared.sqrt(ms)

  def rms_filter(self, cutoff_factor,
                 use_binning=0, use_multiplicities=0, negate=0):
    rms = self.rms(use_binning, use_multiplicities)
    abs_f = shared.abs(self.F)
    if (not use_binning):
      keep = abs_f <= cutoff_factor * rms
    else:
      keep = self.all_selection()
      for i_bin in self.binner.range_used():
        keep &= ~self.binner(i_bin) | (abs_f <= cutoff_factor * rms[i_bin-1])
    return self.apply_selection(keep, negate)

  def statistical_mean(self, use_binning=0):
    if (not use_binning):
      result = miller.statistical_mean(
        self.SgOps, friedel_flag, self.H, self.F)
    else:
      result = shared.double()
      for i_bin in self.binner.range_used():
        sel = self.binner(i_bin)
        if (sel.count(1) == 0):
          result.append(0)
        else:
          result.append(miller.statistical_mean(
            self.SgOps, self.friedel_flag,
            self.H.select(sel),
            self.F.select(sel)))
    return result

  def generic_binner_application(self, binned_values, ftor, result_f):
    result_p = shared.size_t()
    for i_bin in self.binner.range_used():
      result_p.append(self.binner.array_indices(i_bin))
      result_f.append(
        ftor(self.F.select(self.binner(i_bin)),
             binned_values[i_bin-1]))
    return reciprocal_space_array(self, result_f.shuffle(result_p))

  def remove_patterson_origin_peak(self):
    s_mean = self.statistical_mean(use_binning=1)
    result_f = shared.double()
    return self.generic_binner_application(s_mean, ftor_a_s_sub, result_f)

  def normalize_structure_factors(self):
    mean = self.mean(use_binning=1, use_multiplicities=1)
    result_f = shared.double()
    return self.generic_binner_application(mean, ftor_a_s_div, result_f)

  def __add__(self, other):
    if (type(other) != type(self)):
      # add a scalar
      f = 0
      if (self.F): f = self.F.add(other)
      s = 0
      if (self.sigmas): s = self.sigmas.add(other)
      return reciprocal_space_array(self, f, s)
    # add arrays
    js = miller.join_sets(self.H, other.H)
    h = js.paired_miller_indices(0)
    f = js.plus(self.F, other.F)
    s = 0
    if (self.sigmas):
      s = js.additive_sigmas(self.sigmas, other.sigmas)
    return reciprocal_space_array(miller_set(self, h), f, s)

class symmetrized_sites(crystal_symmetry):

  def __init__(self, xsym, sites = None,
               MinMateDistance = 0.5,
               Ustar_tolerance = 0.1,
               TestPositiveDefiniteness = 1,
               keep_special_position_operators = 1):
    crystal_symmetry.__init__(self, xsym.UnitCell, xsym.SgInfo)
    self.MinMateDistance = MinMateDistance
    self.Ustar_tolerance = Ustar_tolerance
    self.TestPositiveDefiniteness = TestPositiveDefiniteness
    self.Sites = shared.XrayScatterer()
    if (keep_special_position_operators):
      self.SpecialPositionOps = shared.RTMx()
    else:
      self.SpecialPositionOps = None
    self.SnapParameters = sgtbx.SpecialPositionSnapParameters(
      self.UnitCell, self.SgOps, 1, self.MinMateDistance)
    if (sites):
      self.add_sites(sites)

  def discard_special_position_info(self):
    self.SpecialPositionOps = None
    self.SnapParameters = None

  def copy_attributes(self):
    return symmetrized_sites(self, sites = None,
      MinMateDistance = self.MinMateDistance,
      Ustar_tolerance = self.Ustar_tolerance,
      TestPositiveDefiniteness = self.TestPositiveDefiniteness)

  def n_special_positions(self):
    assert self.SpecialPositionOps != None
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
    if (self.SpecialPositionOps != None):
      self.SpecialPositionOps.append(special_op)

  def add_sites(self, sites):
    for site in sites:
      self.add_site(site)

  def __len__(self): return self.size()
  def size(self): return self.Sites.size()

  def __getitem__(self, key):
    return self.Sites[key]

  def as_emma_model(self):
    from cctbx import euclidean_model_matching as emma
    positions = []
    for site in self.Sites:
      positions.append(emma.labeled_position(site.Label(), site.Coordinates()))
    return emma.model(self, positions)

def build_miller_set(xsym, friedel_flag, d_min):
  return miller_set(
    xsym, miller.BuildIndices(
      xsym.UnitCell, xsym.SgInfo, friedel_flag, d_min)).set_friedel_flag(
        friedel_flag)

def calculate_structure_factors_direct(miller_set_obj, xtal, abs_F=0):
  F = sftbx.StructureFactorArray(
    miller_set_obj.UnitCell, miller_set_obj.SgOps, miller_set_obj.H,
    xtal.Sites)
  if (abs_F): F = shared.abs(F)
  result = reciprocal_space_array(miller_set_obj, F)
  result.symmetrized_sites = xtal
  return result

def calculate_structure_factors_fft(miller_set_obj, xtal, abs_F=0,
                                    grid_resolution_factor=1/3.,
                                    max_prime=5,
                                    quality_factor=100,
                                    wing_cutoff=1.e-3,
                                    exp_table_one_over_step_size=-100):
  from cctbx_boost import fftbx
  max_q = xtal.UnitCell.max_Q(miller_set_obj.H)
  mandatory_grid_factors = xtal.SgOps.refine_gridding()
  grid_logical = sftbx.determine_grid(
    xtal.UnitCell,
    max_q, grid_resolution_factor, max_prime, mandatory_grid_factors)
  rfft = fftbx.real_to_complex_3d(grid_logical)
  u_extra = sftbx.calc_u_extra(max_q, grid_resolution_factor, quality_factor)
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
    friedel_flag, miller_set_obj.H, map, n_complex, collect_conj)
  sampled_density.eliminate_u_extra_and_normalize(miller_set_obj.H, fcalc)
  if (abs_F): fcalc = shared.abs(fcalc)
  result = reciprocal_space_array(miller_set_obj, fcalc)
  result.symmetrized_sites = xtal
  return result

def calculate_structure_factors(miller_set_obj, xtal, abs_F=0, method=0):
  assert method in (0, "fft", "direct")
  if (method == 0):
    approx_number_of_atoms = xtal.Sites.size() * xtal.SgOps.OrderZ()
    if (approx_number_of_atoms > 30):
      method = "fft"
    else:
      method = "direct"
  if (method == "fft"):
    result = calculate_structure_factors_fft(miller_set_obj, xtal, abs_F)
  else:
    result = calculate_structure_factors_direct(miller_set_obj, xtal, abs_F)
  result.fcalc_method = method
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

class fft_map(crystal_symmetry):

  def __init__(self, coeff_set, grid_resolution_factor=1/3., max_prime=5):
    self.grid_resolution_factor = grid_resolution_factor
    self.max_prime = max_prime
    crystal_symmetry.__init__(
      self, coeff_set.UnitCell, coeff_set.SgInfo)
    self.friedel_flag = coeff_set.friedel_flag
    cf = coeff_set.F
    if (type(cf[0]) != type(0j)):
      if (cf.size() == 0):
        cf = shared.complex_double() # XXX avoid large dummy array
      else:
        ph = shared.double(cf.size(), 0)
        cf = shared.polar(cf, ph)
        del ph
    max_q = self.UnitCell.max_Q(coeff_set.H)
    mandatory_grid_factors = self.SgOps.refine_gridding()
    self.n_real = sftbx.determine_grid(
      self.UnitCell, max_q,
      self.grid_resolution_factor, self.max_prime, mandatory_grid_factors)
    if (self.friedel_flag):
      rfft = fftbx.real_to_complex_3d(self.n_real)
      self.m_real = rfft.Mreal()
      self.n_complex = rfft.Ncomplex()
    else:
      cfft = fftbx.complex_to_complex_3d(self.n_real)
      self.m_real = cfft.N()
      self.n_complex = self.m_real
    conjugate = 0 # XXX correct?
    cmap = sftbx.structure_factor_map(
      self.SgOps, self.friedel_flag,
      coeff_set.H, cf,
      self.n_complex, conjugate)
    if (self.friedel_flag):
      rfft.backward(cmap)
      self.rmap = shared.reinterpret_complex_as_real(cmap)
    else:
      cfft.backward(cmap)
      self.cmap = cmap

  def inplace_unpad(self):
    if (self.friedel_flag):
      shared_map.inplace_unpad(self.rmap, self.n_real, self.m_real)
      self.m_real = self.n_real
      del self.n_complex

  def get_real_map(self):
    if (self.friedel_flag):
      self.inplace_unpad()
      return self.rmap
    else:
      return shared.real(self.cmap)
