import cctbx.sgtbx

import boost.python
ext = boost.python.import_ext("cctbx_miller_ext")
from cctbx_miller_ext import *

from cctbx import crystal
from cctbx import maptbx
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
from scitbx import fftpack
import scitbx.math
from scitbx.python_utils.misc import store
from libtbx.itertbx import count
from libtbx.utils import Keep
from libtbx import introspection
import sys
import math
import types

def _slice_or_none(array, slice_object):
  assert type(slice_object) == types.SliceType
  if (array is None): return None
  return array.__getitem__(slice_object)

class _binning(boost.python.injector, ext.binning):

  def __getinitargs__(self):
    raise RuntimeError(
      "cctbx.miller.binning instances are not picklable."
      " If this error appears while pickling a cctbx.miller.set"
      " or cctbx.miller.array instance use the clear_binner()"
      " method before pickling.")

class binner(ext.binner):

  def __init__(self, binning, miller_set):
    ext.binner.__init__(self, binning, miller_set.indices())
    self.space_group_info = miller_set.space_group_info()
    self.anomalous_flag = miller_set.anomalous_flag()
    if (miller_set.indices().size() == 0):
      self._completeness_d_min = binning.d_min()
    else:
      self._completeness_d_min = miller_set.d_min()
    self._counts_given = None
    self._counts_complete = None
    self._have_format_strings = False

  def __getinitargs__(self):
    raise RuntimeError(
      "cctbx.miller.binner instances are not picklable."
      " If this error appears while pickling a cctbx.miller.set"
      " or cctbx.miller.array instance use the clear_binner()"
      " method before pickling.")

  def counts_given(self):
    if (self._counts_given is None):
      self._counts_given = []
      for i_bin in self.range_all():
        self._counts_given.append(self.count(i_bin))
    return self._counts_given

  def counts_complete(self,
        include_centric=True,
        include_acentric=True,
        d_min_tolerance=1.e-6):
    if (self._counts_complete is None):
      assert self.anomalous_flag in (False, True)
      complete_set = build_set(
        crystal_symmetry=crystal.symmetry(
          unit_cell=self.unit_cell(),
          space_group_info=self.space_group_info),
        anomalous_flag=self.anomalous_flag,
        d_min=self._completeness_d_min*(1-d_min_tolerance))
      if (not include_centric):
        complete_set = complete_set.select_acentric()
      if (not include_acentric):
        complete_set = complete_set.select_centric()
      binner_complete = binner(binning=self, miller_set=complete_set)
      self._counts_complete = binner_complete.counts_given()
    return self._counts_complete

  def n_bin_d_too_large(self):
    return self.array_indices(self.i_bin_d_too_large()).size()

  def n_bin_d_too_small(self):
    return self.array_indices(self.i_bin_d_too_small()).size()

  def n_bin_d_too_large_or_small(self):
    return self.n_bin_d_too_large() + self.n_bin_d_too_small()

  def _setup_format_strings(self):
    if (not self._have_format_strings):
      self._have_format_strings = True
      n = max(2, len(str(self.n_bins_used())))
      self.fmt_bin = "bin %%%dd:"%n
      self.fmt_unused = " "*(4+n+1-7) + "unused:"
      n = len("%.4f" % self.bin_d_range(1)[0])
      self.fmt_d = "%%%d.4f"%n
      blank_d = " "*n
      self.fmt_d_range_first = " ".join([blank_d,    "-", self.fmt_d])
      self.fmt_d_range_used  = " ".join([self.fmt_d, "-", self.fmt_d])
      self.fmt_d_range_last  = " ".join([self.fmt_d, "-", blank_d])
      self.fmt_counts_given = "%%%dd"%len(str(max(self.counts_given())))
      self.fmt_counts_complete = "%%-%dd"%len(str(max(self.counts_complete())))
      self.fmt_both_counts \
        = "["+self.fmt_counts_given+"/"+self.fmt_counts_complete+"]"

  def bin_legend(self,
        i_bin,
        show_bin_number=True,
        show_d_range=True,
        show_counts=True):
    self._setup_format_strings()
    is_first = (i_bin == self.i_bin_d_too_large())
    is_last = (i_bin == self.i_bin_d_too_small())
    result = []
    if (show_bin_number):
      if (is_first or is_last):
        result.append(self.fmt_unused)
      else:
        result.append(self.fmt_bin % i_bin)
    bin_d_range = self.bin_d_range(i_bin)
    if (show_d_range):
      if (is_first):
        result.append(self.fmt_d_range_first % bin_d_range[1])
      elif (is_last):
        result.append(self.fmt_d_range_last % bin_d_range[0])
      else:
        result.append(self.fmt_d_range_used % bin_d_range)
    if (show_counts):
      result.append(self.fmt_both_counts % (
        self._counts_given[i_bin], self._counts_complete[i_bin]))
    return " ".join(result)

  def show_summary(self,
        show_bin_number=True,
        show_d_range=True,
        show_counts=True,
        f=None,
        prefix=""):
    if (f is None): f = sys.stdout
    for i_bin in self.range_all():
      print >> f, prefix + self.bin_legend(
        i_bin=i_bin,
        show_bin_number=show_bin_number,
        show_d_range=show_d_range,
        show_counts=show_counts)

  def show_data(self,
        data,
        data_fmt=None,
        show_bin_number=True,
        show_d_range=True,
        show_counts=True,
        show_unused=True,
        f=None,
        prefix=""):
    assert len(data) == self.n_bins_all()
    if (f is None): f = sys.stdout
    if (show_unused):
      i_bins = self.range_all()
    else:
      i_bins = self.range_used()
    legend = None
    for i_bin in i_bins:
      legend = self.bin_legend(
        i_bin=i_bin,
        show_bin_number=show_bin_number,
        show_d_range=show_d_range,
        show_counts=show_counts)
      print >> f, prefix + legend,
      if (data[i_bin] is not None):
        if (isinstance(data[i_bin], str) or data_fmt is None):
          print >> f, data[i_bin],
        elif (isinstance(data_fmt, str)):
          print >> f, data_fmt % data[i_bin],
        else:
          s = data_fmt(binner=self, i_bin=i_bin, data=data)
          if (s is not None and len(s) > 0): print >> f, s,
      print >> f
    if (legend is not None): return len(legend)
    return None

class binned_data:

  def __init__(self, binner, data, data_fmt=None):
    self.binner = binner
    self.data = data
    self.data_fmt = data_fmt

  def show(self,
        data_fmt=None,
        show_bin_number=True,
        show_d_range=True,
        show_counts=True,
        show_unused=True,
        f=None,
        prefix=""):
    if (data_fmt is None): data_fmt = self.data_fmt
    return self.binner.show_data(
      data=self.data,
      data_fmt=data_fmt,
      show_bin_number=show_bin_number,
      show_d_range=show_d_range,
      show_counts=show_counts,
      show_unused=show_unused,
      f=f,
      prefix=prefix)

def make_lookup_dict(indices): # XXX push to C++
  result = {}
  for i in xrange(len(indices)):
    result[indices[i]] = i
  return result

class set(crystal.symmetry):

  def __init__(self, crystal_symmetry, indices, anomalous_flag=None):
    assert anomalous_flag in (None, False, True)
    crystal.symmetry._copy_constructor(self, crystal_symmetry)
    self._indices = indices
    self._anomalous_flag = anomalous_flag
    self._binner = None

  def _copy_constructor(self, other):
    crystal.symmetry._copy_constructor(self, other)
    self._indices = other._indices
    self._anomalous_flag = other._anomalous_flag
    self._binner = None

  def indices(self):
    return self._indices

  def anomalous_flag(self):
    return self._anomalous_flag

  def size(self):
    return self.indices().size()

  def copy(self):
    return set(
      crystal_symmetry=self,
      indices=self._indices,
      anomalous_flag=self._anomalous_flag)

  def deep_copy(self):
    return set(
      crystal_symmetry=crystal.symmetry(
        unit_cell=uctbx.unit_cell(self.unit_cell().parameters()),
        space_group_symbol=str(self.space_group_info())),
      indices=self.indices().deep_copy(),
      anomalous_flag=self.anomalous_flag())

  def customized_copy(self,
        crystal_symmetry=Keep,
        indices=Keep,
        anomalous_flag=Keep,
        unit_cell=Keep,
        space_group_info=Keep):
    if (crystal_symmetry is Keep): crystal_symmetry = self
    if (indices is Keep): indices = self.indices()
    if (anomalous_flag is Keep): anomalous_flag = self.anomalous_flag()
    crystal_symmetry = crystal.symmetry.customized_copy(crystal_symmetry,
      unit_cell=unit_cell,
      space_group_info=space_group_info)
    return set(
      crystal_symmetry=crystal_symmetry,
      indices=indices,
      anomalous_flag=anomalous_flag)

  def array(self, data=None, sigmas=None):
    return array(miller_set=self, data=data, sigmas=sigmas)

  def __getitem__(self, slice_object):
    assert type(slice_object) == types.SliceType
    assert self.indices() is not None
    return set(
      crystal_symmetry=self,
      indices=self.indices().__getitem__(slice_object),
      anomalous_flag=self.anomalous_flag())

  def show_summary(self, f=None, prefix=""):
    "Minimal Miller set summary"
    if (f is None): f = sys.stdout
    print >> f, prefix + "Number of Miller indices:", len(self.indices())
    print >> f, prefix + "Anomalous flag:", self.anomalous_flag()
    crystal.symmetry.show_summary(self, f=f, prefix=prefix)
    return self

  def show_comprehensive_summary(self, f=None, prefix=""):
    "Comprehensive Miller set or array summary"
    if (f is None): f = sys.stdout
    self.show_summary(f=f, prefix=prefix)
    no_sys_abs = self.copy()
    if (self.space_group_info() is not None):
      sys_absent_flags = self.sys_absent_flags().data()
      n_sys_abs = sys_absent_flags.count(True)
      print >> f, prefix + "Systematic absences:", n_sys_abs
      if (n_sys_abs != 0):
        no_sys_abs = self.select(selection=~sys_absent_flags)
        print >> f, prefix + "Systematic absences not included in following:"
      n_centric = no_sys_abs.centric_flags().data().count(True)
      print >> f, prefix + "Centric reflections:", n_centric
    if (self.unit_cell() is not None):
      print >> f, prefix + "Resolution range: %.6g %.6g" % (
        no_sys_abs.resolution_range())
      if (self.space_group_info() is not None and self.indices().size() > 0):
        no_sys_abs.setup_binner(n_bins=1)
        completeness_d_max_d_min = no_sys_abs.completeness(use_binning=True)
        binner = completeness_d_max_d_min.binner
        assert binner.counts_given()[0] == 0
        assert binner.counts_given()[2] == 0
        n_obs = binner.counts_given()[1]
        n_complete = binner.counts_complete()[1]
        if (n_complete != 0):
          print >> f, prefix + "Completeness in resolution range: %.6g" % (
            float(n_obs) / n_complete)
        n_complete += binner.counts_complete()[0]
        if (n_complete != 0):
          print >> f, prefix + "Completeness with d_max=infinity: %.6g" % (
            float(n_obs) / n_complete)
    if (self.space_group_info() is not None and no_sys_abs.anomalous_flag()):
      asu, matches = no_sys_abs.match_bijvoet_mates()
      print >> f, prefix + "Bijvoet pairs:", matches.pairs().size()
      print >> f, prefix + "Lone Bijvoet mates:", \
        matches.n_singles() - n_centric
      if (isinstance(self, array) and self.is_real_array()):
        print >> f, prefix + "Anomalous signal: %.4f" % (
          no_sys_abs.anomalous_signal())
    return self

  def sys_absent_flags(self):
    return array(
      self,
      self.space_group().is_sys_absent(self.indices()))

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

  def d_star_sq(self):
    return array(
      self, self.unit_cell().d_star_sq(self.indices()))

  def d_spacings(self):
    return array(
      self, self.unit_cell().d(self.indices()))

  def sin_theta_over_lambda_sq(self):
    return array(
      self, self.unit_cell().stol_sq(self.indices()))

  def two_theta(self, wavelength, deg=False):
    return array(
      self, self.unit_cell().two_theta(self.indices(), wavelength, deg))

  def d_min(self):
    return uctbx.d_star_sq_as_d(self.unit_cell().max_d_star_sq(self.indices()))

  def n_bijvoet_pairs(self):
    asu, matches = self.match_bijvoet_mates()
    return matches.pairs().size()

  def as_non_anomalous_set(self):
    return set(
      crystal_symmetry=self,
      indices=self.indices(),
      anomalous_flag=False)

  def auto_anomalous(self, min_n_bijvoet_pairs=None,
                           min_fraction_bijvoet_pairs=None):
    assert [min_n_bijvoet_pairs, min_fraction_bijvoet_pairs].count(None) > 0
    if (min_fraction_bijvoet_pairs is not None):
      anomalous_flag = (2.*self.n_bijvoet_pairs()/self.indices().size()
                        >= min_fraction_bijvoet_pairs)
    elif (min_n_bijvoet_pairs is not None):
      anomalous_flag = (self.n_bijvoet_pairs() >= min_n_bijvoet_pairs)
    else:
      anomalous_flag = (self.n_bijvoet_pairs() > 0)
    return set(
      crystal_symmetry=self,
      indices=self.indices(),
      anomalous_flag=anomalous_flag)

  def map_to_asu(self):
    i = self.indices().deep_copy()
    anomalous_flag = self.anomalous_flag()
    if (anomalous_flag is None):
      anomalous_flag = True
    map_to_asu(
      self.space_group_info().type(),
      anomalous_flag,
      i)
    return set(self, i, self.anomalous_flag())

  def complete_set(self, d_min_tolerance=1.e-6):
    assert self.anomalous_flag() in (False, True)
    if (self.indices().size() == 0):
      return set(
        crystal_symmetry=self,
        anomalous_flag=self.anomalous_flag(),
        indices=flex.miller_index())
    return build_set(
      crystal_symmetry=self,
      anomalous_flag=self.anomalous_flag(),
      d_min=self.d_min()*(1-d_min_tolerance))

  def completeness(self, use_binning=False, d_min_tolerance=1.e-6):
    if (not use_binning):
      complete_set = self.complete_set(d_min_tolerance=d_min_tolerance)
      return self.indices().size() \
           / float(max(1,complete_set.indices().size()))
    assert self.binner() is not None
    data = []
    for n_given,n_complete in zip(self.binner().counts_given(),
                                  self.binner().counts_complete()):
      if (n_complete == 0): data.append(None)
      else: data.append(float(n_given)/n_complete)
    return binned_data(binner=self.binner(), data=data, data_fmt="%5.3f")

  def all_selection(self):
    return flex.bool(self.indices().size(), True)

  def select(self, selection, negate=False, anomalous_flag=None):
    assert self.indices() is not None
    if (anomalous_flag is None):
      anomalous_flag = self.anomalous_flag()
    if (negate): selection = ~selection
    i = self.indices().select(selection)
    return set(self, i, anomalous_flag)

  def select_acentric(self):
    return self.select(~self.centric_flags().data())

  def select_centric(self):
    return self.select(self.centric_flags().data())

  def remove_systematic_absences(self, negate=False):
    return self.select(
      selection=self.sys_absent_flags().data(),
      negate=not negate)

  def resolution_filter(self, d_max=0, d_min=0, negate=0):
    d = self.d_spacings().data()
    keep = self.all_selection()
    if (d_max): keep &= d <= d_max
    if (d_min): keep &= d >= d_min
    return self.select(keep, negate)

  def apply_scaling(self, target_max=None, factor=None):
    assert [target_max, factor].count(None) == 1
    assert self.data() is not None
    s = None
    if (target_max is not None):
      current_max = flex.max(flex.abs(self.data()))
      if (current_max == 0): return self.deep_copy()
      factor = target_max / current_max
    d = self.data() * factor
    if (self.sigmas() is not None): s = self.sigmas() * factor
    return (array(
      miller_set = set.deep_copy(self),
      data=d,
      sigmas=s)
      .set_info(self.info())
      .set_observation_type(self))

  def match_bijvoet_mates(self):
    assert self.anomalous_flag() in (None, True)
    assert self.indices() is not None
    if (self.space_group_info() is not None):
      asu = self.map_to_asu()
      matches = match_bijvoet_mates(
        asu.space_group_info().type(), asu.indices())
    else:
      asu = self
      matches = match_bijvoet_mates(asu.indices())
    return asu, matches

  def resolution_range(self):
    r = self.unit_cell().min_max_d_star_sq(self.indices())
    return tuple([uctbx.d_star_sq_as_d(x) for x in r])

  def sort_permutation(self, by_value="resolution", reverse=False):
    assert by_value in ["resolution", "packed_indices"]
    assert reverse in [False, True]
    if (by_value == "resolution"):
      return flex.sort_permutation(
        data=self.unit_cell().d_star_sq(self.indices()), reverse=reverse)
    else:
      return flex.sort_permutation(
        data=index_span(self.indices()).pack(self.indices()), reverse=reverse)

  def sort(self, by_value="resolution", reverse=False):
    return self.select(
      self.sort_permutation(by_value=by_value, reverse=reverse))

  def change_basis(self, cb_op):
    if (isinstance(cb_op, str)): cb_op = sgtbx.change_of_basis_op(cb_op)
    return set.customized_copy(self,
      crystal_symmetry=crystal.symmetry.change_basis(self, cb_op),
      indices=cb_op.apply(self.indices()))

  def expand_to_p1(self):
    assert self.space_group_info() is not None
    assert self.indices() is not None
    assert self.anomalous_flag() is not None
    p1 = expand_to_p1(
      self.space_group(), self.anomalous_flag(), self.indices())
    return set(self.cell_equivalent_p1(), p1.indices(), self.anomalous_flag())

  def patterson_symmetry(self):
    assert self.anomalous_flag() == False
    return set.customized_copy(self,
      crystal_symmetry=crystal.symmetry.patterson_symmetry(self))

  def crystal_gridding(self, resolution_factor=1/3.,
                             d_min=None,
                             grid_step=None,
                             symmetry_flags=None,
                             mandatory_factors=None,
                             max_prime=5,
                             assert_shannon_sampling=True):
    if (d_min is None): d_min = self.d_min()
    return maptbx.crystal_gridding(
      unit_cell=self.unit_cell(),
      d_min=d_min,
      resolution_factor=resolution_factor,
      step=grid_step,
      symmetry_flags=symmetry_flags,
      space_group_info=self.space_group_info(),
      mandatory_factors=mandatory_factors,
      max_prime=max_prime,
      assert_shannon_sampling=assert_shannon_sampling)

  def structure_factors_from_map(self, map, in_place_fft=False):
    assert map.focus_size_1d() > 0 and map.nd() == 3 and map.is_0_based()
    assert isinstance(map, flex.double) or isinstance(map, flex.complex_double)
    assert in_place_fft in (False, True)
    if (isinstance(map, flex.double)):
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
    conjugate_flag = True
    from_map = maptbx.structure_factors.from_map(
      self.anomalous_flag(),
      self.indices(),
      map,
      conjugate_flag)
    return array(miller_set=self, data=from_map.data())

  def structure_factors_from_scatterers(self, xray_structure,
                                        algorithm=None,
                                        cos_sin_table=False,
                                        quality_factor=None,
                                        u_base=None,
                                        b_base=None,
                                        wing_cutoff=None,
                                        exp_table_one_over_step_size=None):
    from cctbx import xray
    if (algorithm == "direct"):
      return xray.structure_factors.from_scatterers_direct(
        xray_structure=xray_structure,
        miller_set=self,
        cos_sin_table=cos_sin_table)
    return xray.structure_factors.from_scatterers(
      miller_set=self,
      cos_sin_table=cos_sin_table,
      quality_factor=quality_factor,
      u_base=u_base,
      b_base=b_base,
      wing_cutoff=wing_cutoff,
      exp_table_one_over_step_size=exp_table_one_over_step_size)(
        xray_structure=xray_structure,
        miller_set=self,
        algorithm=algorithm)

  def f_obs_minus_xray_structure_f_calc(self, f_obs_factor, xray_structure,
        structure_factor_algorithm=None,
        cos_sin_table=False,
        quality_factor=None,
        u_base=None,
        b_base=None,
        wing_cutoff=None,
        exp_table_one_over_step_size=None):
    return self.f_obs_minus_f_calc(
      f_obs_factor=f_obs_factor,
      f_calc=self.structure_factors_from_scatterers(
        xray_structure=xray_structure,
        algorithm=structure_factor_algorithm,
        cos_sin_table=cos_sin_table,
        quality_factor=quality_factor,
        u_base=u_base,
        b_base=b_base,
        wing_cutoff=wing_cutoff,
        exp_table_one_over_step_size=exp_table_one_over_step_size).f_calc())

  def setup_binner(self, d_max=0, d_min=0,
                   auto_binning=False,
                   reflections_per_bin=0,
                   n_bins=0):
    assert auto_binning or reflections_per_bin != 0 or n_bins != 0
    assert auto_binning or (reflections_per_bin == 0 or n_bins == 0)
    if (auto_binning):
      if (reflections_per_bin == 0): reflections_per_bin = 200
      if (n_bins == 0): n_bins = 8
      n_per_bin = int(float(len(self.indices())) / n_bins + .5)
      if (n_per_bin > reflections_per_bin):
        n_bins = int(len(self.indices()) / reflections_per_bin + .5)
    elif (reflections_per_bin):
      n_bins = int(len(self.indices()) / reflections_per_bin + .5)
    assert n_bins > 0
    assert self.unit_cell() is not None
    assert self.indices().size() > 0 or d_min > 0
    bng = binning(self.unit_cell(), n_bins, self.indices(), d_max, d_min)
    self._binner = binner(bng, self)
    return self.binner()

  def binner(self):
    return self._binner

  def use_binning_of(self, other):
    self._binner = binner(other.binner(), self)

  def use_binner_of(self, other):
    assert self.indices().all_eq(other.indices())
    self._binner = other._binner

  def clear_binner(self):
    self._binner = None

def build_set(crystal_symmetry, anomalous_flag, d_min):
  return set(
    crystal_symmetry,
    index_generator(
      crystal_symmetry.unit_cell(),
      crystal_symmetry.space_group_info().type(),
      anomalous_flag,
      d_min).to_array(),
    anomalous_flag)

class array_info:

  def __init__(self,
        source=None,
        source_type=None,
        history=None,
        labels=None,
        merged=False):
    introspection.adopt_init_args()

  def customized_copy(self,
        source=Keep,
        source_type=Keep,
        history=Keep,
        labels=Keep,
        merged=Keep):
    if (source is Keep): source = self.source
    if (source_type is Keep): source_type = self.source_type
    if (history is Keep): history = self.history
    if (labels is Keep): labels = self.labels
    if (merged is Keep): merged = self.merged
    return array_info(
      source=source,
      source_type=source_type,
      history=history,
      labels=labels,
      merged=merged)

  def as_string_part_2(self):
    part_2 = []
    if (self.labels is not None):
      part_2.extend(self.labels)
    if (self.merged):
      part_2.append("merged")
    return part_2

  def label_string(self):
    part_2 = self.as_string_part_2()
    if (len(part_2) > 0): return ",".join(part_2)
    return None

  def __str__(self):
    result = []
    if (self.source is not None):
      result.append(str(self.source))
    elif (self.source_type is not None):
      result.append(str(self.source_type))
    part_2 = self.as_string_part_2()
    if (len(part_2) > 0):
      result.append(",".join(part_2))
    if (len(result) == 0):
      return "None"
    return ":".join(result)

def raw_array_summary(array):
  if (array is None): return str(None)
  try:
    return array.__class__.__name__ + ", size=%d" % (len(array),)
  except:
    return "Unknown"

class array(set):

  def __init__(self, miller_set, data=None, sigmas=None):
    set._copy_constructor(self, miller_set)
    self._data = data
    self._sigmas = sigmas
    self._info = None
    self._observation_type = None

  def _copy_constructor(self, other):
    set._copy_constructor(self, other)
    self._data = other._data
    self._sigmas = other._sigmas
    self._info = other._info
    self._observation_type = other._observation_type

  def set_info(self, info):
    self._info = info
    return self

  def set_observation_type(self, observation_type):
    from cctbx.xray import observation_types
    if (isinstance(observation_type, array)):
      observation_type = observation_type.observation_type()
    assert observation_type is None or isinstance(observation_type, observation_types.any)
    self._observation_type = observation_type
    return self

  def set_observation_type_xray_amplitude(self):
    from cctbx.xray import observation_types
    return self.set_observation_type(observation_types.amplitude())

  def set_observation_type_xray_intensity(self):
    from cctbx.xray import observation_types
    return self.set_observation_type(observation_types.intensity())

  def data(self):
    return self._data

  def sigmas(self):
    return self._sigmas

  def info(self):
    return self._info

  def observation_type(self):
    return self._observation_type

  def size(self):
    assert self.indices() is not None
    assert self.data() is not None
    assert self.data().size() == self.indices().size()
    if (self.sigmas() is not None):
      assert self.sigmas().size() == self.indices().size()
    return set.size(self)

  def is_bool_array(self):
    return isinstance(self.data(), flex.bool)

  def is_integer_array(self):
    return isinstance(self.data(), flex.int) \
        or isinstance(self.data(), flex.long)

  def is_real_array(self):
    return isinstance(self.data(), flex.float) \
        or isinstance(self.data(), flex.double)

  def is_complex_array(self):
    return isinstance(self.data(), flex.complex_double)

  def is_hendrickson_lattman_array(self):
    return isinstance(self.data(), flex.hendrickson_lattman)

  def is_xray_amplitude_array(self):
    from cctbx.xray import observation_types
    return isinstance(self.observation_type(), observation_types.amplitude)

  def is_xray_reconstructed_amplitude_array(self):
    from cctbx.xray import observation_types
    return isinstance(
      self.observation_type(), observation_types.reconstructed_amplitude)

  def is_xray_intensity_array(self):
    from cctbx.xray import observation_types
    return isinstance(self.observation_type(), observation_types.intensity)

  def copy(self):
    return (array(
      miller_set=self,
      data=self.data(),
      sigmas=self.sigmas())
      .set_info(self.info())
      .set_observation_type(self))

  def deep_copy(self):
    d = None
    s = None
    if (self.data() is not None): d = self.data().deep_copy()
    if (self.sigmas() is not None): s = self.sigmas().deep_copy()
    return (array(
      miller_set = set.deep_copy(self),
      data=d,
      sigmas=s)
      .set_info(self.info())
      .set_observation_type(self))

  def customized_copy(self,
        miller_set=Keep,
        data=Keep,
        sigmas=Keep,
        crystal_symmetry=Keep,
        indices=Keep,
        anomalous_flag=Keep,
        unit_cell=Keep,
        space_group_info=Keep):
    if (miller_set is Keep): miller_set = self
    if (data is Keep): data = self.data()
    if (sigmas is Keep): sigmas = self.sigmas()
    miller_set = set.customized_copy(miller_set,
      crystal_symmetry=crystal_symmetry,
      indices=indices,
      anomalous_flag=anomalous_flag,
      unit_cell=unit_cell,
      space_group_info=space_group_info)
    return array(miller_set=miller_set, data=data, sigmas=sigmas)

  def discard_sigmas(self):
    return array.customized_copy(self, sigmas=None)

  def conjugate(self):
    assert self.is_complex_array()
    return array.customized_copy(self, data=flex.conj(self.data()))

  def __getitem__(self, slice_object):
    return array(
      miller_set=set.__getitem__(self, slice_object),
      data=_slice_or_none(self.data(), slice_object),
      sigmas=_slice_or_none(self.sigmas(), slice_object))

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, prefix + "Miller %s info:" % (
      self.__class__.__name__), self.info()
    print >> f, prefix + "Observation type:", self.observation_type()
    print >> f, prefix + "Type of data:", raw_array_summary(self.data())
    print >> f, prefix + "Type of sigmas:", raw_array_summary(self.sigmas())
    set.show_summary(self, f=f, prefix=prefix)
    return self

  def f_sq_as_f(self, tolerance=1.e-6):
    from cctbx import xray
    assert self.observation_type() is None or self.is_xray_intensity_array()
    if (self.sigmas() is None):
      result = array(self, xray.array_f_sq_as_f(self.data()).f)
    else:
      r = xray.array_f_sq_as_f(self.data(), self.sigmas(), tolerance)
      result = array(self, r.f, r.sigma_f)
    if (self.is_xray_intensity_array()):
      result.set_observation_type_xray_amplitude()
    return result

  def f_as_f_sq(self):
    from cctbx import xray
    assert self.observation_type() is None or self.is_xray_amplitude_array()
    if (self.sigmas() is None):
      result = array(self, xray.array_f_as_f_sq(self.data()).f_sq)
    else:
      r = xray.array_f_as_f_sq(self.data(), self.sigmas())
      result = array(self, r.f_sq, r.sigma_f_sq)
    if (self.is_xray_amplitude_array()):
      result.set_observation_type_xray_intensity()
    return result

  def as_amplitude_array(self):
    if (self.is_xray_intensity_array()):
      return self.f_sq_as_f()
    return self

  def as_intensity_array(self):
    if (not self.is_xray_intensity_array()):
      return self.f_as_f_sq()
    return self

  def map_to_asu(self):
    i = self.indices().deep_copy()
    d = self.data().deep_copy()
    if (self.is_complex_array() or self.is_hendrickson_lattman_array()):
      map_to_asu(
        self.space_group_info().type(),
        self.anomalous_flag(),
        i, d)
    else:
      map_to_asu(
        self.space_group_info().type(),
        self.anomalous_flag(),
        i)
    return (array(set(self, i, self.anomalous_flag()), d, self.sigmas())
      .set_observation_type(self))

  def eliminate_sys_absent(self):
    sys_absent_flags = self.sys_absent_flags().data()
    if (sys_absent_flags.all_eq(False)): return self
    return self.select(selection=~sys_absent_flags)

  def adopt_set(self, other, assert_is_similar_symmetry=True):
    if (assert_is_similar_symmetry):
      assert self.is_similar_symmetry(other)
    assert self.indices().size() == other.indices().size()
    assert self.anomalous_flag() == other.anomalous_flag()
    p = match_indices(other.indices(), self.indices()).permutation()
    assert flex.order(self.indices().select(p), other.indices()) == 0
    d = self.data()
    s = self.sigmas()
    if (d is not None): d = d.select(p)
    if (s is not None): s = s.select(p)
    return (array(miller_set=other, data=d, sigmas=s)
      .set_observation_type(self))

  def match_indices(self, other, assert_is_similar_symmetry=True):
    if (assert_is_similar_symmetry):
      assert self.is_similar_symmetry(other)
    assert self.anomalous_flag() == other.anomalous_flag()
    return match_indices(self.indices(), other.indices())

  def common_set(self, other, assert_is_similar_symmetry=True):
    pairs = other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry).pairs()
    return self.select(pairs.column(1))

  def common_sets(self, other, assert_is_similar_symmetry=True):
    pairs = other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry).pairs()
    return [self.select(pairs.column(1)),
            other.select(pairs.column(0))]

  def lone_set(self, other, assert_is_similar_symmetry=True):
    return self.select(other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry).singles(1))

  def lone_sets(self, other, assert_is_similar_symmetry=True):
    matches = other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry)
    return [self.select(matches.singles(1)),
            other.select(matches.singles(0))]

  def sort_permutation(self, by_value="resolution", reverse=False):
    assert reverse in (False, True)
    if (by_value in ["resolution", "packed_indices"]):
      result = set.sort_permutation(self,
        by_value=by_value, reverse=reverse)
    elif (by_value == "data"):
      result = flex.sort_permutation(
        data=self.data(), reverse=not reverse)
    elif (by_value == "abs"):
      result = flex.sort_permutation(
        data=flex.abs(self.data()), reverse=not reverse)
    else:
      result = flex.sort_permutation(
        data=by_value, reverse=not reverse)
    return result

  def patterson_symmetry(self):
    data = self.data()
    if (self.is_complex_array()):
      data = flex.abs(self.data())
    return array(
      set.patterson_symmetry(self),
      data,
      self.sigmas())

  def expand_to_p1(self, phase_deg=None):
    assert self.space_group_info() is not None
    assert self.indices() is not None
    assert self.anomalous_flag() is not None
    assert self.data() is not None
    new_sigmas = None
    if (self.is_complex_array()):
      assert phase_deg is None
      p1 = expand_to_p1(
        self.space_group(), self.anomalous_flag(), self.indices(),
        self.data())
      new_data = p1.structure_factors()
    elif (self.is_hendrickson_lattman_array()):
      p1 = expand_to_p1(
        self.space_group(), self.anomalous_flag(), self.indices(),
        self.data())
      new_data = p1.hendrickson_lattman_coefficients()
    else:
      assert isinstance(self.data(), flex.double)
      assert phase_deg in (None, False, True)
      if (phase_deg is None):
        p1 = expand_to_p1(
          self.space_group(), self.anomalous_flag(), self.indices(),
          self.data())
        new_data = p1.amplitudes()
        if (self.sigmas() is not None):
          assert isinstance(self.sigmas(), flex.double)
          p1 = expand_to_p1(
            self.space_group(), self.anomalous_flag(), self.indices(),
            self.sigmas())
          new_sigmas = p1.amplitudes()
      else:
        p1 = expand_to_p1(
          self.space_group(), self.anomalous_flag(), self.indices(),
          self.data(), phase_deg)
        new_data = p1.phases()
    assert self.sigmas() is None or new_sigmas is not None
    return array(
      set(self.cell_equivalent_p1(), p1.indices(), self.anomalous_flag()),
      data=new_data,
      sigmas=new_sigmas).set_observation_type(self)

  def change_basis(self, cb_op, deg=None):
    if (isinstance(cb_op, str)): cb_op = sgtbx.change_of_basis_op(cb_op)
    if (deg is False or deg is True):
      assert self.is_real_array()
      result = change_basis_phases_double(
        cb_op=cb_op,
        indices_in=self.indices(),
        data_in=self.data(),
        deg=deg)
    elif (self.is_complex_array()):
      result = change_basis_complex_double(
        cb_op=cb_op,
        indices_in=self.indices(),
        data_in=self.data())
    elif (   self.is_bool_array()
          or self.is_integer_array()
          or self.is_real_array()):
      result = store(
        indices=cb_op.apply(self.indices()),
        data=self.data().deep_copy())
    elif (self.is_hendrickson_lattman_array()):
      result = change_basis_hendrickson_lattman(
        cb_op=cb_op,
        indices_in=self.indices(),
        data_in=self.data())
    else:
      raise RuntimeError("Unsupported miller.array data type.")
    result_sigmas = None
    if (self.sigmas() is not None):
      assert isinstance(self.sigmas(), flex.double)
      result_sigmas = self.sigmas().deep_copy()
    return array(
      miller_set=set(
        crystal_symmetry=crystal.symmetry.change_basis(self, cb_op),
        indices=result.indices,
        anomalous_flag=self.anomalous_flag()),
      data=result.data,
      sigmas=result_sigmas)

  def f_obs_minus_f_calc(self, f_obs_factor, f_calc):
    assert self.is_real_array()
    assert f_calc.is_complex_array()
    assert self.indices().all_eq(f_calc.indices())
    assert not self.anomalous_flag()
    assert not f_calc.anomalous_flag()
    return array(
      miller_set=self,
      data=f_obs_factor*self.data()-flex.abs(f_calc.data())).phase_transfer(
        phase_source=f_calc)

  def phase_transfer(self, phase_source, epsilon=1.e-10, deg=False,
                           phase_integrator_n_steps=None):
    assert self.data() is not None
    if (hasattr(phase_source, "data")):
      phase_source = phase_source.data()
    assert isinstance(self.data(), flex.complex_double) or isinstance(self.data(), flex.double)
    assert isinstance(phase_source, flex.complex_double) or isinstance(phase_source, flex.double) or isinstance(phase_source, flex.hendrickson_lattman)
    if (isinstance(phase_source, flex.hendrickson_lattman)):
      if (phase_integrator_n_steps is None):
        integrator = phase_integrator()
      else:
        integrator = phase_integrator(n_steps=phase_integrator_n_steps)
      phase_source = integrator(
        space_group=self.space_group(),
        miller_indices=self.indices(),
        hendrickson_lattman_coefficients=phase_source)
    if (isinstance(phase_source, flex.complex_double)):
      return array(
        miller_set=self,
        data=phase_transfer(
          self.space_group(),
          self.indices(),
          self.data(),
          phase_source,
          epsilon))
    return array(
      miller_set=self,
      data=phase_transfer(
        self.space_group(),
        self.indices(),
        self.data(),
        phase_source,
        deg))

  def phase_integrals(self, n_steps=None, integrator=None):
    assert self.is_hendrickson_lattman_array()
    assert n_steps is None or integrator is None
    if (integrator is None):
      if (n_steps is None):
        integrator = phase_integrator()
      else:
        integrator = phase_integrator(n_steps=n_steps)
    return array(
      miller_set=self,
      data=integrator(
        space_group=self.space_group(),
        miller_indices=self.indices(),
        hendrickson_lattman_coefficients=self.data()))

  def mean_weighted_phase_error(self, phase_source):
    assert self.data() is not None
    if (isinstance(phase_source, array)):
      assert flex.order(phase_source.indices(), self.indices()) == 0
      phase_source = phase_source.data()
    p1 = flex.arg(self.data())
    assert isinstance(phase_source, flex.complex_double) or isinstance(phase_source, flex.double)
    if (isinstance(phase_source, flex.complex_double)):
      p2 = flex.arg(phase_source)
    else:
      p2 = phase_source
    e = scitbx.math.phase_error(phi1=p1, phi2=p2)
    w = flex.abs(self.data())
    sum_w = flex.sum(w)
    assert sum_w != 0
    sum_we = flex.sum(w * e)
    return sum_we / sum_w * 180/math.pi

  def anomalous_differences(self):
    assert self.data() is not None
    assert self.observation_type() is None or self.is_xray_amplitude_array()
    asu, matches = self.match_bijvoet_mates()
    i = matches.miller_indices_in_hemisphere("+")
    d = matches.minus(asu.data())
    s = None
    if (asu.sigmas() is not None):
      s = matches.additive_sigmas(asu.sigmas())
    return array(set(asu, i, anomalous_flag=False), d, s)

  def hemisphere(self, plus_or_minus):
    assert plus_or_minus in ("+", "-")
    assert self.data() is not None
    asu, matches = self.match_bijvoet_mates()
    i_column = "+-".index(plus_or_minus)
    return asu.select(
      selection=matches.pairs().column(i_column),
      anomalous_flag=False)

  def hemispheres(self):
    assert self.data() is not None
    asu, matches = self.match_bijvoet_mates()
    return tuple(
      [asu.select(
        selection=matches.pairs().column(i_column),
        anomalous_flag=False)
       for i_column in (0,1)])

  def anomalous_signal(self, use_binning=False):
    "sqrt((<||F(+)|-|F(-)||**2>)/(1/2(<|F(+)|>**2+<|F(-)|>**2)))"
    assert not use_binning or self.binner() is not None
    if (not use_binning):
      obs = self.select(self.data() > 0)
      if (self.is_xray_intensity_array()):
        obs = obs.f_sq_as_f()
      f_plus, f_minus = obs.hemispheres()
      assert f_plus.data().size() == f_minus.data().size()
      if (f_plus.data().size() == 0): return 0
      mean_sq_diff = flex.mean(flex.pow2(f_plus.data() - f_minus.data()))
      assert mean_sq_diff >= 0
      mean_sum_sq = flex.mean(  flex.pow2(f_plus.data())
                              + flex.pow2(f_minus.data()))
      assert mean_sum_sq > 0
      return math.sqrt(2 * mean_sq_diff / mean_sum_sq)
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      results.append(self.select(sel).anomalous_signal())
    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f")

  def measurability(self, use_binning=False, cutoff=3.0):
    ## Peter Zwart 3/4/2005
    """\
Fraction of reflections for which (|delta I|/sigma_dI) > cutoff
            and min(I_plus/sigma_plus,I_min/sigma_min) > cutoff"""
    assert not use_binning or self.binner() is not None
    assert self.sigmas() is not None
    cutoff = float(cutoff)
    if (not use_binning):
      obs = self.select(self.data() > 0 )
      if (self.is_xray_amplitude_array()):
        obs = obs.f_as_f_sq()
      if (obs.data().size() == 0): return 0
      i_plus, i_minus = obs.hemispheres()
      assert i_plus.data().size() == i_minus.data().size()
      ratio = flex.fabs(i_plus.data()-i_minus.data()) / flex.sqrt(
                             (i_plus.sigmas()*i_plus.sigmas())
                           + (i_minus.sigmas()*i_minus.sigmas()) )
      i_plus_sigma = i_plus.data()/i_plus.sigmas()
      i_minus_sigma = i_minus.data()/i_minus.sigmas()
      meas = (  (ratio > cutoff)
              & (i_plus_sigma > cutoff)
              & (i_minus_sigma > cutoff) ).count(True)
      return float(meas)/float(ratio.size())
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      results.append(self.select(sel).measurability(cutoff=cutoff))
    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f")

  def second_moment(self, use_binning=False):
    "<data^2>/(<data>)^2"
    assert not use_binning or self.binner() is not None
    if (not use_binning):
      if (self.indices().size() == 0): return None
      mean_data_sq = flex.mean(self.data())**2
      if (mean_data_sq == 0): return None
      return flex.mean_sq(self.data()) / mean_data_sq
    result = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      result.append(self.select(sel).second_moment())
    return binned_data(binner=self.binner(), data=result, data_fmt="%7.4f")

  def second_moment_of_intensities(self, use_binning=False):
    "<I^2>/(<I>)^2 (2.0 for untwinned, 1.5 for twinned data)"
    if (self.is_xray_intensity_array()):
      a = self
    else:
      a = self.f_as_f_sq()
      if (use_binning):
        a.use_binner_of(self)
    return a.second_moment(use_binning=use_binning)

  def wilson_ratio(self, use_binning=False):
    "(<F>)^2/<F^2> (0.785 for untwinned, 0.885 for twinned data)"
    if (not self.is_xray_intensity_array()):
      a = self
    else:
      a = self.f_sq_as_f()
      if (use_binning):
        a.use_binner_of(self)
    second_moment = a.second_moment(use_binning=use_binning)
    if (second_moment is None): return None
    if (not use_binning): return 1/second_moment
    result = []
    for sm in second_moment.data:
      if (sm is None or sm == 0): result.append(None)
      else: result.append(1/sm)
    return binned_data(binner=a.binner(), data=result, data_fmt="%7.4f")

  def select(self, selection, negate=False, anomalous_flag=None):
    assert self.indices() is not None
    if (anomalous_flag is None):
      anomalous_flag = self.anomalous_flag()
    if (negate): selection = ~selection
    i = self.indices().select(selection)
    d = None
    if (self.data() is not None): d = self.data().select(selection)
    s = None
    if (self.sigmas() is not None): s = self.sigmas().select(selection)
    return array(set(self, i, anomalous_flag), d, s).set_observation_type(self)

  def sigma_filter(self, cutoff_factor, negate=False):
    assert self.data() is not None
    assert self.sigmas() is not None
    flags = flex.abs(self.data()) >= self.sigmas() * cutoff_factor
    return self.select(flags, negate)

  def mean(self,
        use_binning=False,
        use_multiplicities=False,
        squared=False,
        rms=False):
    assert squared is False or rms is False
    if (not use_binning):
      if (self.data().size() == 0): return None
      if (not squared and not rms):
        if (not use_multiplicities):
          return flex.mean(self.data())
        else:
          return flex.mean_weighted(
            self.data(),
            self.multiplicities().data().as_double())
      if (not use_multiplicities):
        result = flex.mean_sq(self.data())
      else:
        result = flex.mean_sq_weighted(
          self.data(),
          self.multiplicities().data().as_double())
      if (rms): return math.sqrt(result)
      return result
    assert self.binner() is not None
    data = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      data.append(self.select(sel).mean(
        use_multiplicities=use_multiplicities,
        squared=squared,
        rms=rms))
    return binned_data(binner=self.binner(), data=data)

  def mean_sq(self, use_binning=False, use_multiplicities=False):
    return self.mean(
      use_binning=use_binning,
      use_multiplicities=use_multiplicities,
      squared=True)

  def rms(self, use_binning=False, use_multiplicities=False):
    return self.mean(
      use_binning=use_binning,
      use_multiplicities=use_multiplicities,
      rms=True)

  def rms_filter(self,
        cutoff_factor,
        use_binning=False,
        use_multiplicities=False,
        negate=False):
    rms = self.rms(
      use_binning=use_binning,
      use_multiplicities=use_multiplicities)
    abs_data = flex.abs(self.data())
    if (not use_binning):
      keep = abs_data <= cutoff_factor * rms
    else:
      keep = self.all_selection()
      for i_bin in self.binner().range_used():
        keep &= ~self.binner().selection(i_bin) \
             | (abs_data <= cutoff_factor * rms.data[i_bin])
    return self.select(keep, negate)

  def statistical_mean(self, use_binning=0):
    if (not use_binning):
      result = statistical_mean(
        self.space_group(), self.anomalous_flag(), self.indices(), self.data())
    else:
      result = flex.double()
      for i_bin in self.binner().range_used():
        sel = self.binner().selection(i_bin)
        if (sel.count(True) == 0):
          result.append(0)
        else:
          result.append(statistical_mean(
            self.space_group(), self.anomalous_flag(),
            self.indices().select(sel),
            self.data().select(sel)))
    return result

  def remove_patterson_origin_peak(self):
    assert self.observation_type() is None or self.is_xray_intensity_array()
    s_mean = self.statistical_mean(use_binning=True)
    result_data = self.data().deep_copy()
    for i_bin in self.binner().range_used():
      sel = self.binner().array_indices(i_bin)
      if (sel.size() > 0):
        result_data.set_selected(
          sel, self.data().select(sel) - s_mean[i_bin-1])
    return array(self, result_data)

  def quasi_normalized_as_normalized(self):
    assert self.observation_type() is None or self.is_xray_amplitude_array()
    return array(
      miller_set=self,
      data=self.data()/flex.sqrt(self.epsilons().data().as_double()))

  def quasi_normalize_structure_factors(self, d_star_power=1):
    assert self.binner() is not None
    assert self.binner().n_bin_d_too_large_or_small() == 0
    assert self.data().all_ge(0)
    assert self.observation_type() is None or self.is_xray_amplitude_array()
    epsilons = self.epsilons().data().as_double()
    mean_f_sq_over_epsilon = flex.double()
    for i_bin in self.binner().range_used():
      sel = self.binner().selection(i_bin)
      sel_f_sq = flex.pow2(self.data().select(sel))
      if (sel_f_sq.size() > 0):
        sel_epsilons = epsilons.select(sel)
        sel_f_sq_over_epsilon = sel_f_sq / sel_epsilons
        mean_f_sq_over_epsilon.append(flex.mean(sel_f_sq_over_epsilon))
      else:
        mean_f_sq_over_epsilon.append(0)
    mean_f_sq_over_epsilon_interp = self.binner().interpolate(
      mean_f_sq_over_epsilon, d_star_power)
    assert mean_f_sq_over_epsilon_interp.all_gt(0)
    f_sq = flex.pow2(self.data())
    q = flex.sqrt(f_sq / mean_f_sq_over_epsilon_interp)
    assert q.all_ge(0)
    return array(self, q)

  def __abs__(self):
    return array(self, flex.abs(self.data()), self.sigmas())

  def arg(self, deg=False):
    return array(self, flex.arg(self.data(), deg))

  def amplitudes(self):
    assert isinstance(self.data(), flex.complex_double)
    assert self.sigmas() is None
    return abs(self)

  def phases(self, deg=False):
    assert isinstance(self.data(), flex.complex_double)
    assert self.sigmas() is None
    return self.arg(deg)

  def merge_equivalents(self):
    return merge_equivalents(self)

  def as_non_anomalous_array(self):
    return array(
      miller_set=self.as_non_anomalous_set(),
      data=self.data(),
      sigmas=self.sigmas())

  def average_bijvoet_mates(self):
    assert self.observation_type() is None or self.is_xray_amplitude_array()
    return self.as_non_anomalous_array().merge_equivalents().array()

  def __add__(self, other):
    assert self.indices() is not None
    assert self.data() is not None
    if (type(other) != type(self)):
      # add a scalar
      return array(self, self.data() + other)
    # add arrays
    assert other.indices() is not None
    assert other.data() is not None
    match = match_indices(self.indices(), other.indices())
    i = match.paired_miller_indices(0)
    d = match.plus(self.data(), other.data())
    s = None
    if (self.sigmas() is not None and other.sigmas() is not None):
      s = match.additive_sigmas(self.sigmas(), other.sigmas())
    return array(set(self, i), d, s)

  def as_anomalous(self):
    if (self.anomalous_flag()): return self
    sel = ~self.centric_flags().data()
    indices = self.indices().deep_copy()
    indices.extend(-indices.select(sel))
    data = None
    sigmas = None
    if (self.data() is not None):
      data = self.data().deep_copy()
      if (self.is_complex_array()):
        data.extend(flex.conj(data.select(sel)))
      else:
        data.extend(data.select(sel))
    if (self.sigmas() is not None):
      sigmas = self.sigmas().deep_copy()
      sigmas.extend(sigmas.select(sel))
    return array(
      miller_set=set(
        crystal_symmetry=self,
        indices=indices,
        anomalous_flag=True),
      data=data,
      sigmas=sigmas)

  def correlation(self,
        other,
        use_binning=False,
        assert_is_similar_symmetry=True):
    if (assert_is_similar_symmetry):
      assert self.is_similar_symmetry(other)
    assert self.is_real_array()
    assert other.is_real_array()
    assert not use_binning or self.binner() is not None
    lhs = self
    if (lhs.anomalous_flag() and not other.anomalous_flag()):
      other = other.as_anomalous()
    elif (not lhs.anomalous_flag() and other.anomalous_flag()):
      lhs = lhs.as_anomalous()
    lhs, other = lhs.common_sets(
      other=other, assert_is_similar_symmetry=assert_is_similar_symmetry)
    if (not use_binning):
      return flex.linear_correlation(lhs.data(), other.data())
    lhs.use_binning_of(self)
    data = []
    for i_bin in self.binner().range_all():
      sel = lhs.binner().selection(i_bin)
      correlation = flex.linear_correlation(
        lhs.data().select(sel),
        other.data().select(sel))
      if (not correlation.is_well_defined()): data.append(None)
      else: data.append(correlation.coefficient())
    return binned_data(binner=lhs.binner(), data=data, data_fmt="%6.3f")

  def show_array(self, f=None, prefix=""):
    "Listing of Miller indices and data"
    if (f is None): f = sys.stdout
    assert self.data().size() == self.indices().size()
    if (self.sigmas() is None):
      for h,d in zip(self.indices(), self.data()):
        print >> f, prefix + str(h), d
    else:
      assert self.sigmas().size() == self.indices().size()
      for h,d,s in zip(self.indices(), self.data(), self.sigmas()):
        print >> f, prefix + str(h), d, s
    return self

  def fft_map(self, resolution_factor=1/3.,
                    d_min=None,
                    grid_step=None,
                    symmetry_flags=None,
                    mandatory_factors=None,
                    max_prime=5,
                    assert_shannon_sampling=True,
                    f_000=None):
    return fft_map(
      crystal_gridding=self.crystal_gridding(
        d_min=d_min,
        resolution_factor=resolution_factor,
        grid_step=grid_step,
        symmetry_flags=symmetry_flags,
        mandatory_factors=mandatory_factors,
        max_prime=max_prime,
        assert_shannon_sampling=assert_shannon_sampling),
      fourier_coefficients=self,
      f_000=f_000)

  def patterson_map(self, resolution_factor=1/3.,
                          d_min=None,
                          symmetry_flags=None,
                          mandatory_factors=None,
                          max_prime=5,
                          assert_shannon_sampling=True,
                          f_000=None,
                          sharpening=False,
                          origin_peak_removal=False):
    f_patt = self.patterson_symmetry()
    return patterson_map(
      crystal_gridding=f_patt.crystal_gridding(
        resolution_factor=resolution_factor,
        d_min=d_min,
        symmetry_flags=symmetry_flags,
        mandatory_factors=mandatory_factors,
        max_prime=max_prime,
        assert_shannon_sampling=assert_shannon_sampling),
      f_patt=f_patt,
      f_000=f_000,
      sharpening=sharpening,
      origin_peak_removal=origin_peak_removal)

  def as_mtz_dataset(self,
        column_root_label,
        column_types=None,
        column_label_decorator=None,
        title=None,
        crystal_name="crystal",
        project_name="project",
        dataset_name="dataset",
        wavelength=1.0):
    import iotbx.mtz
    return iotbx.mtz.miller_array_as_mtz_dataset(self,
      column_root_label=column_root_label,
      column_types=column_types,
      column_label_decorator=column_label_decorator,
      title=title,
      crystal_name=crystal_name,
      project_name=project_name,
      dataset_name=dataset_name,
      wavelength=wavelength)

  def as_phases_phs(self, out):
    import iotbx.phases
    iotbx.phases.miller_array_as_phases_phs(self=self, out=out)

class merge_equivalents:

  def __init__(self, miller_array):
    assert isinstance(miller_array.data(), flex.double)
    if (miller_array.sigmas() is not None):
      assert isinstance(miller_array.sigmas(), flex.double)
      sel = (miller_array.sigmas() <= 0) & (miller_array.data() == 0)
      if (sel.count(True) > 0):
        miller_array = miller_array.select(~sel)
    asu_set = set.map_to_asu(miller_array)
    perm = asu_set.sort_permutation(by_value="packed_indices")
    if (miller_array.sigmas() is not None):
      sigmas_squared = flex.pow2(miller_array.sigmas().select(perm))
      assert flex.min(sigmas_squared) > 0
      merge_ext = ext.merge_equivalents(
        asu_set.indices().select(perm),
        miller_array.data().select(perm),
        1./sigmas_squared)
      sigmas = merge_ext.sigmas()
    else:
      merge_ext = ext.merge_equivalents(
        asu_set.indices().select(perm),
        miller_array.data().select(perm))
      sigmas = None
    self._array = array(
      miller_set=set(
        crystal_symmetry=miller_array,
        indices=merge_ext.indices(),
        anomalous_flag=miller_array.anomalous_flag()),
      data=merge_ext.data(),
      sigmas=sigmas).set_observation_type(miller_array)
    self._redundancies = merge_ext.redundancies()

  def array(self):
    return self._array

  def redundancies(self):
    return self._redundancies

class fft_map(maptbx.crystal_gridding):

  def __init__(self, crystal_gridding, fourier_coefficients, f_000=None):
    maptbx.crystal_gridding._copy_constructor(self, crystal_gridding)
    assert fourier_coefficients.anomalous_flag() in (False, True)
    assert fourier_coefficients.unit_cell().is_similar_to(self.unit_cell())
    assert fourier_coefficients.space_group() == self.space_group()
    assert isinstance(fourier_coefficients.data(), flex.complex_double)
    self._anomalous_flag = fourier_coefficients.anomalous_flag()
    if (not self.anomalous_flag()):
      rfft = fftpack.real_to_complex_3d(self.n_real())
      n_complex = rfft.n_complex()
    else:
      cfft = fftpack.complex_to_complex_3d(self.n_real())
      n_complex = cfft.n()
    conjugate_flag = True
    map = maptbx.structure_factors.to_map(
      self.space_group(),
      self.anomalous_flag(),
      fourier_coefficients.indices(),
      fourier_coefficients.data(),
      self.n_real(),
      flex.grid(n_complex),
      conjugate_flag)
    if (f_000 is not None):
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

  def real_map_unpadded(self):
    result = self.real_map()
    if (not result.is_padded()): return result
    return maptbx.copy(result, flex.grid(result.focus()))

  def complex_map(self):
    assert self.anomalous_flag()
    return self._complex_map

  def apply_sigma_scaling(self):
    if (not self.anomalous_flag()):
      statistics = maptbx.statistics(self._real_map)
      if (statistics.sigma() != 0):
        self._real_map /= statistics.sigma()
    else:
      statistics = maptbx.statistics(self.real_map())
      if (statistics.sigma() != 0):
        self._complex_map /= complex(statistics.sigma())
    return self

def patterson_map(crystal_gridding, f_patt, f_000=None,
                  sharpening=False,
                  origin_peak_removal=False):
  assert f_patt.is_patterson_symmetry()
  if (sharpening):
    f_patt.setup_binner(auto_binning=True)
    f_patt = f_patt.quasi_normalize_structure_factors()
  i_patt = f_patt.f_as_f_sq()
  if (origin_peak_removal):
    i_patt.setup_binner(auto_binning=True)
    i_patt = i_patt.remove_patterson_origin_peak()
  i_patt = array(i_patt, data=flex.polar(i_patt.data(), 0))
  if (f_000 is not None):
    f_000 = f_000 * f_000
  return fft_map(crystal_gridding, i_patt, f_000)
