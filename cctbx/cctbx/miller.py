import cctbx.sgtbx

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.miller_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from cctbx import crystal
from cctbx import maptbx
from cctbx import uctbx
from cctbx.utils import phase_error
from cctbx.array_family import flex
from scitbx import fftpack
import sys
import math
import types

def _slice_or_none(array, slice_object):
  assert type(slice_object) == types.SliceType
  if (array == None): return None
  return array.__getitem__(slice_object)

class binner(ext.binner):

  def __init__(self, binning, miller_indices):
    ext.binner.__init__(self, binning, miller_indices)

  def show_summary(self, f=sys.stdout):
    for i_bin in self.range_all():
      bin_d_range = self.bin_d_range(i_bin)
      count = self.count(i_bin)
      if (i_bin == self.i_bin_d_too_large()):
        assert bin_d_range[0] == -1
        print >> f, "unused:              d > %8.4f: %5d" % (
          bin_d_range[1], count)
      elif (i_bin == self.i_bin_d_too_small()):
        assert bin_d_range[1] == -1
        print >> f, "unused: %9.4f >  d           : %5d" % (
          bin_d_range[0], count)
      else:
        print >> f, "bin %2d: %9.4f >= d > %8.4f: %5d" % (
          (i_bin,) + bin_d_range + (count,))

def make_lookup_dict(indices): # XXX push to C++
  result = {}
  for i in xrange(len(indices)):
    result[indices[i]] = i
  return result

class set(crystal.symmetry):

  def __init__(self, crystal_symmetry, indices, anomalous_flag=None):
    assert anomalous_flag in (None, 00000, 0001)
    crystal.symmetry._copy_constructor(self, crystal_symmetry)
    self._indices = indices
    self._anomalous_flag = anomalous_flag

  def _copy_constructor(self, other):
    crystal.symmetry._copy_constructor(self, other)
    self._indices = other._indices
    self._anomalous_flag = other._anomalous_flag

  def indices(self):
    return self._indices

  def anomalous_flag(self):
    return self._anomalous_flag

  def deep_copy(self):
    return set(
      crystal_symmetry=crystal.symmetry(
        unit_cell=uctbx.unit_cell(self.unit_cell().parameters()),
        space_group_symbol=str(self.space_group_info())),
      indices=self.indices().deep_copy(),
      anomalous_flag=self.anomalous_flag())

  def __getitem__(self, slice_object):
    assert type(slice_object) == types.SliceType
    assert self.indices() != None
    return set(
      crystal_symmetry=self,
      indices=self.indices().__getitem__(slice_object),
      anomalous_flag=self.anomalous_flag())

  def show_summary(self, f=sys.stdout):
    print >> f, "Number of Miller indices:", len(self.indices())
    print >> f, "Anomalous flag:", self.anomalous_flag()
    crystal.symmetry.show_summary(self, f)
    return self

  def centric_flags(self):
    return array(
      self,
      self.space_group().is_centric(self.indices()))

  def multiplicities(self):
    return array(
      self,
      self.space_group().multiplicity(self.indices(), self.anomalous_flag()))

  def epsilons(self):
    return array(
      self,
      self.space_group().epsilon(self.indices()))

  def d_spacings(self):
    return array(
      self, self.unit_cell().d(self.indices()))

  def sin_theta_over_lambda_sq(self):
    return array(
      self, self.unit_cell().stol_sq(self.indices()))

  def d_min(self):
    return uctbx.d_star_sq_as_d(self.unit_cell().max_d_star_sq(self.indices()))

  def resolution_range(self):
    r = self.unit_cell().min_max_d_star_sq(self.indices())
    return tuple([uctbx.d_star_sq_as_d(x) for x in r])

  def sort(self, by_value="resolution", reverse=00000):
    assert by_value in ("resolution",)
    assert reverse in (00000, 0001)
    p = flex.sort_permutation(
      self.unit_cell().d_star_sq(self.indices()),
      reverse)
    return set(
      crystal_symmetry=self,
      indices=self.indices().shuffle(p),
      anomalous_flag=self.anomalous_flag())

  def change_basis(self, cb_op):
    return set(
      crystal_symmetry=crystal.symmetry.change_basis(self, cb_op),
      indices=cb_op.apply(self.indices()),
      anomalous_flag=self.anomalous_flag())

  def expand_to_p1(self):
    assert self.space_group() != None
    assert self.indices() != None
    assert self.anomalous_flag() != None
    p1 = expand_to_p1(
      self.space_group(), self.anomalous_flag(), self.indices())
    return set(self.cell_equivalent_p1(), p1.indices(), self.anomalous_flag())

  def patterson_symmetry(self):
    assert self.anomalous_flag() == 00000
    return set(
      crystal.symmetry.patterson_symmetry(self),
      self.indices(),
      self.anomalous_flag())

  def crystal_gridding(self, resolution_factor=1/3.,
                             d_min=None,
                             symmetry_flags=None,
                             mandatory_factors=None,
                             max_prime=5,
                             assert_shannon_sampling=0001):
    if (d_min == None): d_min = self.d_min()
    return maptbx.crystal_gridding(
      unit_cell=self.unit_cell(),
      d_min=d_min,
      resolution_factor=resolution_factor,
      symmetry_flags=symmetry_flags,
      space_group_info=self.space_group_info(),
      mandatory_factors=mandatory_factors,
      max_prime=max_prime,
      assert_shannon_sampling=assert_shannon_sampling)

  def structure_factors_from_map(self, map, in_place_fft=00000):
    assert map.focus_size_1d() > 0 and map.nd() == 3 and map.is_0_based()
    assert type(map[0]) in (type(float()), type(complex()))
    assert in_place_fft in (00000, 0001)
    if (type(map[0]) == type(float())):
      fft = fftpack.real_to_complex_3d(map.focus())
      if (not map.is_padded()):
        assert not in_place_fft
        assert map.focus() == fft.n_real()
        map = maptbx.copy(map, flex.grid(fft.m_real()).set_focus(fft.n_real()))
      elif (not in_place_fft):
        map = map.deep_copy()
    else:
      if (not in_place_fft):
        map = map.deep_copy()
      fft = fftpack.complex_to_complex_3d(map.focus())
    map = fft.forward(map)
    conjugate_flag = 0001
    from_map = maptbx.structure_factors.from_map(
      self.anomalous_flag(),
      self.indices(),
      map,
      conjugate_flag)
    return array(miller_set=self, data=from_map.data())

  def setup_binner(self, d_max=0, d_min=0,
                   auto_binning=0,
                   reflections_per_bin=0,
                   n_bins=0):
    assert auto_binning != 0 or reflections_per_bin != 0 or n_bins != 0
    assert auto_binning != 0 or (reflections_per_bin == 0 or n_bins == 0)
    if (auto_binning):
      if (reflections_per_bin == 0): reflections_per_bin = 200
      if (n_bins == 0): n_bins = 8
      n_per_bin = int(float(len(self.indices())) / n_bins + .5)
      if (n_per_bin > reflections_per_bin):
        n_bins = int(len(self.indices()) / reflections_per_bin + .5)
    elif (reflections_per_bin):
      n_bins = int(len(self.indices()) / reflections_per_bin + .5)
    assert n_bins > 0
    assert self.unit_cell() != None
    bng = binning(self.unit_cell(), n_bins, self.indices(), d_max, d_min)
    self._binner = binner(bng, self.indices())
    return self.binner()

  def binner(self):
    return self._binner

  def use_binner_of(self, other):
    self._binner = other._binner

def build_set(crystal_symmetry, anomalous_flag, d_min):
  return set(
    crystal_symmetry,
    index_generator(
      crystal_symmetry.unit_cell(),
      crystal_symmetry.space_group_info().type(),
      anomalous_flag,
      d_min).to_array(),
    anomalous_flag)

def _ftor_a_s_sub(lhs, rhs): return lhs - rhs
def _ftor_a_s_div(lhs, rhs): return lhs / rhs

def _array_info(array):
  if (array == None): return str(None)
  try:
    return array.__class__.__name__ + ", size=%d" % (len(array),)
  except:
    return "Unknown"

class array(set):

  def __init__(self, miller_set, data=None, sigmas=None, info=None):
    set._copy_constructor(self, miller_set)
    self._data = data
    self._sigmas = sigmas
    self._info = info

  def _copy_constructor(self, other):
    set._copy_constructor(self, other)
    self._data = other._data
    self._sigmas = other._sigmas
    self._info = other._info

  def data(self):
    return self._data

  def sigmas(self):
    return self._sigmas

  def info(self):
    return self._info

  def deep_copy(self):
    d = None
    s = None
    if (self.data() != None): d = self.data().deep_copy()
    if (self.sigmas() != None): s = self.sigmas().deep_copy()
    return array(
      miller_set = set.deep_copy(self),
      data=d,
      sigmas=s)

  def __getitem__(self, slice_object):
    return array(
      miller_set=set.__getitem__(self, slice_object),
      data=_slice_or_none(self.data(), slice_object),
      sigmas=_slice_or_none(self.sigmas(), slice_object))

  def show_summary(self, f=sys.stdout):
    print >> f, "Miller array info:", self.info()
    print >> f, "Type of data:", _array_info(self.data())
    print >> f, "Type of sigmas:", _array_info(self.sigmas())
    set.show_summary(self, f)
    return self

  def set_info(self, info):
    self._info = info

  def f_sq_as_f(self, tolerance=1.e-6):
    from cctbx import xray
    if (self.sigmas() != None):
      r = xray.array_f_sq_as_f(self.data(), self.sigmas(), tolerance)
      return array(self, r.f, r.sigma_f)
    return array(self, xray.array_f_sq_as_f(self.data()).f)

  def f_as_f_sq(self):
    from cctbx import xray
    if (self.sigmas() != None):
      r = xray.array_f_as_f_sq(self.data(), self.sigmas())
      return array(self, r.f_sq, r.sigma_f_sq)
    return array(self, xray.array_f_as_f_sq(self.data()).f_sq)

  def map_to_asu(self):
    i = self.indices().deep_copy()
    d = self.data().deep_copy()
    map_to_asu(
      self.space_group_info().type(),
      self.anomalous_flag(),
      i, d)
    return array(set(self, i, self.anomalous_flag()), d, self.sigmas())

  def sort(self, by_value="resolution", reverse=00000):
    assert reverse in (00000, 0001)
    if (by_value == "resolution"):
      p = flex.sort_permutation(
        self.unit_cell().d_star_sq(self.indices()),
        reverse)
    elif (by_value == "data"):
      p = flex.sort_permutation(self.data(), not reverse)
    elif (by_value == "abs"):
      p = flex.sort_permutation(flex.abs(self.data()), not reverse)
    else:
      p = flex.sort_permutation(by_value, not reverse)
    new_set = set(
      crystal_symmetry=self,
      indices=self.indices().shuffle(p),
      anomalous_flag=self.anomalous_flag())
    d = None
    s = None
    if (self.data() != None): d = self.data().shuffle(p)
    if (self.sigmas() != None): s = self.sigmas().shuffle(p)
    return array(new_set, d, s)

  def patterson_symmetry(self):
    data = self.data()
    if (type(data) == type(flex.complex_double())):
      data = flex.abs(self.data())
    return array(
      set.patterson_symmetry(self),
      data,
      self.sigmas())

  def expand_to_p1(self, phase_deg=None):
    assert self.space_group() != None
    assert self.indices() != None
    assert self.anomalous_flag() != None
    assert self.data() != None
    assert self.sigmas() == None, "Not implemented." # XXX
    if (type(self.data()) == type(flex.complex_double())):
      assert phase_deg == None
      p1 = expand_to_p1(
        self.space_group(), self.anomalous_flag(), self.indices(),
        self.data())
      new_data = p1.structure_factors()
    else:
      assert type(self.data()) == type(flex.double())
      assert phase_deg in (None, 00000, 0001)
      if (phase_deg == None):
        p1 = expand_to_p1(
          self.space_group(), self.anomalous_flag(), self.indices(),
          self.data())
        new_data = p1.amplitudes()
      else:
        p1 = expand_to_p1(
          self.space_group(), self.anomalous_flag(), self.indices(),
          self.data(), phase_deg)
        new_data = p1.phases()
    return array(
      set(self.cell_equivalent_p1(), p1.indices(), self.anomalous_flag()),
      data=new_data)

  def change_basis(self, cb_op):
    new_data = None
    new_sigmas = None
    if (self.data() != None):
      assert type(self.data()) == type(flex.double())
      new_data = self.data().deep_copy()
    if (self.sigmas() != None):
      assert type(self.sigmas()) == type(flex.double())
      new_sigmas = self.sigmas().deep_copy()
    return array(
      miller_set=set.change_basis(self, cb_op),
      data=new_data,
      sigmas=new_sigmas)

  def phase_transfer(self, phase_source, epsilon=1.e-10):
    assert self.data() != None
    if (hasattr(phase_source, "data")):
      phase_source = phase_source.data()
    return array(
      miller_set=self,
      data=phase_transfer(
        self.space_group(),
        self.indices(),
        self.data(),
        phase_source,
        epsilon))

  def mean_weighted_phase_error(self, phase_source):
    assert self.data() != None
    if (hasattr(phase_source, "data")):
      assert flex.order(phase_source.indices(), self.indices()) == 0
      phase_source = phase_source.data()
    p1 = flex.arg(self.data())
    p2 = flex.arg(phase_source)
    assert p1.size() == p2.size()
    e = flex.double()
    for i in p1.indices():
      e.append(phase_error(p1[i], p2[i]))
    w = flex.abs(self.data())
    sum_w = flex.sum(w)
    assert sum_w != 0
    sum_we = flex.sum(w * e)
    return sum_we / sum_w * 180/math.pi

  def match_bijvoet_mates(self):
    assert self.anomalous_flag() == 0001
    assert self.indices() != None
    if (self.space_group() != None):
      asu = self.map_to_asu()
      matches = match_bijvoet_mates(
        asu.space_group_info().type(), asu.indices())
    else:
      asu = self
      matches = match_bijvoet_mates(asu.indices())
    return asu, matches

  def anomalous_differences(self):
    assert self.data() != None
    asu, matches = self.match_bijvoet_mates()
    i = matches.miller_indices_in_hemisphere("+")
    d = matches.minus(asu.data())
    s = None
    if (asu.sigmas() != None):
      s = matches.additive_sigmas(asu.sigmas())
    return array(set(asu, i, anomalous_flag=00000), d, s)

  def hemisphere(self, plus_or_minus):
    assert plus_or_minus in ("+", "-")
    assert self.data() != None
    asu, matches = self.match_bijvoet_mates()
    return asu.apply_selection(
      flags=matches.hemisphere_selection(plus_or_minus),
      anomalous_flag=00000)

  def hemispheres(self):
    assert self.data() != None
    asu, matches = self.match_bijvoet_mates()
    return tuple(
      [asu.apply_selection(flags=matches.hemisphere_selection(plus_or_minus),
                           anomalous_flag=00000)
       for plus_or_minus in ("+", "-")])

  def anomalous_signal(self, use_binning=00000):
    "sqrt((<||f_plus|-|f_minus||**2>)/(1/2(<|f_plus|>**2+<|f_minus|>**2)))"
    assert not use_binning, "Not implemented." # XXX
    assert self.data().all_gt(0)
    f_plus, f_minus = self.hemispheres()
    mean_sq_diff = flex.mean(flex.pow2(f_plus.data() - f_minus.data()))
    assert mean_sq_diff >= 0
    mean_sum_sq = flex.mean(flex.pow2(f_plus.data())+flex.pow2(f_minus.data()))
    assert mean_sum_sq > 0
    return math.sqrt(2 * mean_sq_diff / mean_sum_sq)

  def all_selection(self):
    return flex.bool(self.indices().size(), 0001)

  def apply_selection(self, flags, negate=00000, anomalous_flag=None):
    assert self.indices() != None
    if (anomalous_flag == None):
      anomalous_flag = self.anomalous_flag()
    if (negate): flags = ~flags
    i = self.indices().select(flags)
    d = None
    if (self.data() != None): d = self.data().select(flags)
    s = None
    if (self.sigmas() != None): s = self.sigmas().select(flags)
    return array(set(self, i, anomalous_flag), d, s)

  def remove_systematic_absences(self, negate=00000):
    return self.apply_selection(
      flags=self.space_group().is_sys_absent(self.indices()),
      negate=not negate)

  def resolution_filter(self, d_max=0, d_min=0, negate=0):
    d = self.d_spacings().data()
    keep = self.all_selection()
    if (d_max): keep &= d <= d_max
    if (d_min): keep &= d >= d_min
    return self.apply_selection(keep, negate)

  def sigma_filter(self, cutoff_factor, negate=0):
    assert self.data() != None
    assert self.sigmas() != None
    flags = flex.abs(self.data()) >= self.sigmas() * cutoff_factor
    return self.apply_selection(flags, negate)

  def _generic_binner_action(self, use_binning, use_multiplicities,
                             function,
                             function_weighted):
    assert self.indices() != None
    assert self.data() != None
    if (use_multiplicities):
      mult = self.multiplicities().data().as_double()
    if (not use_binning):
      if (not use_multiplicities):
        result = function(self.data())
      else:
        result = function_weighted(self.data(), mult)
    else:
      result = flex.double()
      for i_bin in self.binner().range_used():
        sel = self.binner().selection(i_bin)
        if (sel.count(0001) == 0):
          result.append(0)
        else:
          sel_data = self.data().select(sel)
          if (not use_multiplicities):
            result.append(function(sel_data))
          else:
            sel_mult = mult.select(sel)
            result.append(function_weighted(sel_data, sel_mult))
    return result

  def mean(self, use_binning=0, use_multiplicities=0):
    return self._generic_binner_action(use_binning, use_multiplicities,
      flex.mean,
      flex.mean_weighted)

  def mean_sq(self, use_binning=0, use_multiplicities=0):
    return self._generic_binner_action(use_binning, use_multiplicities,
      flex.mean_sq,
      flex.mean_sq_weighted)

  def rms(self, use_binning=0, use_multiplicities=0):
    ms = self.mean_sq(use_binning, use_multiplicities)
    if (not use_binning):
      return math.sqrt(ms)
    else:
      return flex.sqrt(ms)

  def rms_filter(self, cutoff_factor,
                 use_binning=0, use_multiplicities=0, negate=0):
    rms = self.rms(use_binning, use_multiplicities)
    abs_data = flex.abs(self.data())
    if (not use_binning):
      keep = abs_data <= cutoff_factor * rms
    else:
      keep = self.all_selection()
      for i_bin in self.binner().range_used():
        keep &= ~self.binner().selection(i_bin) \
             | (abs_data <= cutoff_factor * rms[i_bin-1])
    return self.apply_selection(keep, negate)

  def statistical_mean(self, use_binning=0):
    if (not use_binning):
      result = statistical_mean(
        self.space_group(), self.anomalous_flag(), self.indices(), self.data())
    else:
      result = flex.double()
      for i_bin in self.binner().range_used():
        sel = self.binner().selection(i_bin)
        if (sel.count(0001) == 0):
          result.append(0)
        else:
          result.append(statistical_mean(
            self.space_group(), self.anomalous_flag(),
            self.indices().select(sel),
            self.data().select(sel)))
    return result

  def _generic_binner_application(self, binned_values, ftor, result_data):
    result_perm = flex.size_t()
    for i_bin in self.binner().range_used():
      result_perm.append(self.binner().array_indices(i_bin))
      result_data.append(
        ftor(self.data().select(self.binner().selection(i_bin)),
             binned_values[i_bin-1]))
    assert self.indices().size() == result_data.size()
    return array(self, result_data.unshuffle(result_perm))

  def remove_patterson_origin_peak(self):
    s_mean = self.statistical_mean(use_binning=1)
    result_data = flex.double()
    return self._generic_binner_application(s_mean, _ftor_a_s_sub, result_data)

  def normalize_structure_factors(self, quasi=00000, first_moment=00000):
    assert quasi in (00000, 0001)
    assert first_moment in (00000, 0001)
    assert self.binner() != None
    assert self.data().all_ge(0)
    if (first_moment):
      assert quasi == 00000
      mean = self.mean(use_binning=1, use_multiplicities=1)
      result_data = flex.double()
      return self._generic_binner_application(mean, _ftor_a_s_div, result_data)
    f_sq = flex.pow2(self.data())
    epsilons = self.epsilons().data().as_double()
    e = flex.double()
    e_perm = flex.size_t()
    for i_bin in self.binner().range_used():
      sel = self.binner().selection(i_bin)
      sel_f_sq = f_sq.select(sel)
      if (sel_f_sq.size() > 0):
        sel_epsilons = epsilons.select(sel)
        sel_f_sq_over_epsilon = sel_f_sq / sel_epsilons
        mean_f_sq_over_epsilon = flex.mean(sel_f_sq_over_epsilon)
        if (quasi):
          e.append(flex.sqrt(sel_f_sq / mean_f_sq_over_epsilon))
        else:
          e.append(flex.sqrt(sel_f_sq_over_epsilon / mean_f_sq_over_epsilon))
        e_perm.append(self.binner().array_indices(i_bin))
    assert self.indices().size() == e.size()
    return array(self, e.unshuffle(e_perm))

  def __abs__(self):
    return array(self, flex.abs(self.data()), self.sigmas())

  def arg(self, deg=00000):
    return array(self, flex.arg(self.data(), deg))

  def __add__(self, other):
    assert self.indices() != None
    assert self.data() != None
    if (type(other) != type(self)):
      # add a scalar
      return array(self, self.data() + other)
    # add arrays
    assert other.indices() != None
    assert other.data() != None
    match = match_indices(self.indices(), other.indices())
    i = match.paired_miller_indices(0)
    d = match.plus(self.data(), other.data())
    s = None
    if (self.sigmas() != None and other.sigmas() != None):
      s = match.additive_sigmas(self.sigmas(), other.sigmas())
    return array(set(self, i), d, s)

  def show_array(self, f=sys.stdout):
    assert self.data().size() == self.indices().size()
    if (self.sigmas() == None):
      for i,h in self.indices().items():
        print >> f, h, self.data()[i]
    else:
      assert self.indices().size() == self.sigmas().size()
      for i,h in self.indices().items():
        print >> f, h, self.data()[i], self.sigmas()[i]
    return self

  def fft_map(self, resolution_factor=1/3.,
                    d_min=None,
                    symmetry_flags=None,
                    mandatory_factors=None,
                    max_prime=5,
                    assert_shannon_sampling=0001,
                    f_000=None):
    return fft_map(
      crystal_gridding=self.crystal_gridding(
        resolution_factor=resolution_factor,
        d_min=d_min,
        symmetry_flags=symmetry_flags,
        mandatory_factors=mandatory_factors,
        max_prime=max_prime,
        assert_shannon_sampling=assert_shannon_sampling),
      coeff_array=self,
      f_000=f_000)

  def patterson_map(self, resolution_factor=1/3.,
                          d_min=None,
                          symmetry_flags=None,
                          mandatory_factors=None,
                          max_prime=5,
                          assert_shannon_sampling=0001,
                          f_000=None,
                          sharpening=00000,
                          origin_peak_removal=00000):
    self_patt = self.patterson_symmetry()
    return patterson_map(
      crystal_gridding=self_patt.crystal_gridding(
        resolution_factor=resolution_factor,
        d_min=d_min,
        symmetry_flags=symmetry_flags,
        mandatory_factors=mandatory_factors,
        max_prime=max_prime,
        assert_shannon_sampling=assert_shannon_sampling),
      coeff_array=self_patt,
      f_000=f_000,
      sharpening=sharpening,
      origin_peak_removal=origin_peak_removal)

class fft_map(maptbx.crystal_gridding):

  def __init__(self, crystal_gridding, coeff_array, f_000=None):
    maptbx.crystal_gridding._copy_constructor(self, crystal_gridding)
    assert coeff_array.anomalous_flag() in (00000, 0001)
    assert coeff_array.unit_cell().is_similar_to(self.unit_cell())
    assert coeff_array.space_group() == self.space_group()
    self._anomalous_flag = coeff_array.anomalous_flag()
    if (not self.anomalous_flag()):
      rfft = fftpack.real_to_complex_3d(self.n_real())
      n_complex = rfft.n_complex()
    else:
      cfft = fftpack.complex_to_complex_3d(self.n_real())
      n_complex = cfft.n()
    conjugate_flag = 0001
    assert type(coeff_array.data()) == type(flex.complex_double())
    map = maptbx.structure_factors.to_map(
      self.space_group(),
      self.anomalous_flag(),
      coeff_array.indices(),
      coeff_array.data(),
      self.n_real(),
      flex.grid(n_complex),
      conjugate_flag)
    if (f_000 != None):
      assert map.complex_map()[0] == 0j
      map.complex_map()[0] = complex(f_000)
    if (not self.anomalous_flag()):
      self._real_map = rfft.backward(map.complex_map())
    else:
      self._complex_map = cfft.backward(map.complex_map())

  def anomalous_flag(self):
    return self._anomalous_flag

  def real_map(self):
    if (not self.anomalous_flag()):
      return self._real_map
    else:
      return flex.real(self._complex_map)

  def complex_map(self):
    assert self.anomalous_flag()
    return self._complex_map

def patterson_map(crystal_gridding, coeff_array, f_000=None,
                  sharpening=00000,
                  origin_peak_removal=00000):
  assert coeff_array.is_patterson_symmetry()
  if (sharpening):
    coeff_array.setup_binner(auto_binning=1)
    coeff_array = coeff_array.normalize_structure_factors(quasi=0001)
  coeff_array = coeff_array.f_as_f_sq()
  if (origin_peak_removal):
    coeff_array.setup_binner(auto_binning=1)
    coeff_array = coeff_array.remove_patterson_origin_peak()
  coeff_array = array(coeff_array, data=flex.polar(coeff_array.data(), 0))
  if (f_000 != None):
    f_000 = f_000 * f_000
  return fft_map(crystal_gridding, coeff_array, f_000)
