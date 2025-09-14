# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import cctbx.sgtbx

import boost_adaptbx.boost.python as bp
from six.moves import range
from six.moves import zip
ext = bp.import_ext("cctbx_miller_ext")
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
from cctbx_miller_ext import *

from cctbx import crystal
from cctbx import maptbx
from cctbx import sgtbx
from cctbx.sgtbx import lattice_symmetry
from cctbx import uctbx
from cctbx import r_free_utils
from cctbx.array_family import flex
from scitbx import fftpack
from scitbx.math import distributions
import scitbx.math
from libtbx.math_utils import iround
from libtbx import complex_math
from scitbx.python_utils.misc import store
from libtbx import adopt_init_args
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, Keep, plural_s
from libtbx import group_args, Auto
import libtbx.table_utils
from itertools import count
import warnings
import random
import math
import time
import sys
from collections import namedtuple
from scitbx import matrix

fp_eps_double = scitbx.math.floating_point_epsilon_double_get()

generate_r_free_params_str = r_free_utils.generate_r_free_params_str

def _slice_or_none(array, slice_object):
  assert type(slice_object) == slice
  if (array is None): return None
  return array.__getitem__(slice_object)

class binner(ext.binner):

  def __init__(self, binning, miller_set):
    ext.binner.__init__(self, binning, miller_set.indices())
    self.__getstate_manages_dict__ = True
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
    return (
      binning(
        self.unit_cell(),
        self.limits()),
      set(
        crystal_symmetry=crystal.symmetry(
          unit_cell=self.unit_cell(),
          space_group_info=self.space_group_info),
        indices=self.miller_indices(),
        anomalous_flag=self.anomalous_flag))

  def require_all_bins_have_data(self,min_counts=1,
     error_string=""):
    for i_bin in self.range_used():
      if self.count(i_bin)<min_counts:
        for j_bin in self.range_used():
          print("BIN: %s  COUNTS: %s " %(j_bin,self.count(j_bin)))
        raise Sorry(
         "Bin %s required to have %s or more counts has %s\n%s" %(
            i_bin,min_counts,self.count(i_bin),error_string))


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
      self.fmt_bin_range = "%%%d.4f"%n
      blank_d = " "*n
      self.fmt_bin_range_first = " ".join([blank_d,    "-", self.fmt_bin_range])
      self.fmt_bin_range_used  = " ".join([self.fmt_bin_range, "-", self.fmt_bin_range])
      self.fmt_bin_range_last  = " ".join([self.fmt_bin_range, "-", blank_d])
      self.fmt_counts_given = "%%%dd"%len(str(max(self.counts_given())))
      self.fmt_counts_complete = "%%-%dd"%len(str(max(self.counts_complete())))
      self.fmt_both_counts \
        = "["+self.fmt_counts_given+"/"+self.fmt_counts_complete+"]"

  def bin_legend(self,
        i_bin,
        show_bin_number=True,
        show_bin_range=True,
        show_d_range=None,
        show_counts=True,
        bin_range_as="d",
        wavelength=None,
        join_string=True):
    if show_d_range is not None:
      # XXX backward compatibility 2009-11-16
      show_bin_range = show_d_range
      bin_range_as = "d"
    assert bin_range_as in ("d", "d_star_sq", "stol", "stol_sq", "two_theta")
    if bin_range_as == "two_theta":
      assert wavelength is not None
    self._setup_format_strings()
    is_first = (i_bin == self.i_bin_d_too_large())
    is_last = (i_bin == self.i_bin_d_too_small())
    result = []
    if (show_bin_number):
      if (is_first or is_last):
        result.append(self.fmt_unused)
      else:
        result.append(self.fmt_bin % i_bin)
    bin_range = self.bin_d_range(i_bin)
    if (show_bin_range):
      if bin_range_as == "d_star_sq":
        bin_range = [uctbx.d_as_d_star_sq(bin_d) for bin_d in bin_range]
      elif bin_range_as == "stol":
        bin_range = [uctbx.d_star_sq_as_stol(uctbx.d_as_d_star_sq(bin_d))
                     for bin_d in bin_range]
      elif bin_range_as == "stol_sq":
        bin_range = [uctbx.d_star_sq_as_stol_sq(uctbx.d_as_d_star_sq(bin_d))
                     for bin_d in bin_range]
      elif bin_range_as == "two_theta":
        bin_range = [uctbx.d_star_sq_as_two_theta(
          uctbx.d_as_d_star_sq(bin_d), wavelength=wavelength, deg=True)
                     for bin_d in bin_range]
      if (is_first):
        result.append(self.fmt_bin_range_first % bin_range[1])
      elif (is_last):
        result.append(self.fmt_bin_range_last % bin_range[0])
      else:
        result.append(self.fmt_bin_range_used % tuple(bin_range))
    if (show_counts):
      result.append(self.fmt_both_counts % (
        self._counts_given[i_bin], self._counts_complete[i_bin]))
    if (join_string):
      return " ".join(result)
    else :
      return result

  def show_summary(self,
        show_bin_number=True,
        show_bin_range=True,
        show_d_range=None,
        show_counts=True,
        bin_range_as="d",
        wavelength=None,
        f=None,
        prefix=""):
    if (f is None): f = sys.stdout
    for i_bin in self.range_all():
      print(prefix + self.bin_legend(
        i_bin=i_bin,
        show_bin_number=show_bin_number,
        show_bin_range=show_bin_range,
        show_d_range=show_d_range,
        show_counts=show_counts,
        bin_range_as=bin_range_as,
        wavelength=wavelength), file=f)

  def show_data(self,
        data,
        data_fmt=None,
        show_bin_number=True,
        show_bin_range=True,
        show_d_range=None,
        show_counts=True,
        show_unused=True,
        bin_range_as="d",
        wavelength=None,
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
        show_bin_range=show_bin_range,
        show_d_range=show_d_range,
        show_counts=show_counts,
        bin_range_as=bin_range_as,
        wavelength=wavelength)
      parts = [prefix + legend]
      if (data[i_bin] is not None):
        if (isinstance(data[i_bin], str) or data_fmt is None):
          parts.append(str(data[i_bin]))
        elif (isinstance(data_fmt, str)):
          parts.append(data_fmt % data[i_bin])
        else:
          s = data_fmt(binner=self, i_bin=i_bin, data=data)
          if (s is not None and len(s) > 0): parts.append(str(s))
      print(" ".join(parts), file=f)
    if (legend is not None): return len(legend)
    return None

  def as_simple_table(self,
      data,
      data_label,
      data_fmt=None,
      replace_none_with="---",
      show_bin_number=False,
      show_unused=False):
    """
    Export table rows for display elsewhere.
    """
    assert len(data) == self.n_bins_all()
    if (show_unused):
      i_bins = self.range_all()
    else:
      i_bins = self.range_used()
    table_rows = []
    for i_bin in i_bins:
      legend = self.bin_legend(
        i_bin=i_bin,
        show_bin_number=show_bin_number,
        show_bin_range=True,
        show_d_range=None,
        show_counts=True,
        bin_range_as="d",
        wavelength=None,
        join_string=False)#wavelength)
      row = legend
      if (data[i_bin] is not None):
        if (isinstance(data[i_bin], str) or data_fmt is None):
          row.append(data[i_bin])
        elif (isinstance(data_fmt, str)):
          row.append(data_fmt % data[i_bin])
        else:
          s = data_fmt(binner=self, i_bin=i_bin, data=data)
          row.append(s)
      else :
        row.append(replace_none_with)
      table_rows.append(row)
    labels = ["Resolution range", "N(obs)/N(possible)", data_label ]
    return libtbx.table_utils.simple_table(
      column_headers=labels,
      table_rows=table_rows)


AnomalousProbabilityPlotResult = namedtuple("AnomalousProbabilityPlotResult", [
  "slope", "intercept", "n_pairs", "expected_delta"])

def bins_mixed(d_spacings_sorted, bin_width_max=0.005):
  """
  d_spacings_sorted: d spacings data sorted from lowest to hightest resolution (
  from largest to smallerst d).
  Return bins as list of pairs, each pair is start:end indices to access data.
  """
  size=d_spacings_sorted.size()
  cntr=0
  i_seqs = flex.int()
  n_left = size
  reached_end = False
  bins = []
  for i, d in enumerate(d_spacings_sorted):
    i_seqs.append(i)
    cntr+=1
    n_left-=1
    if cntr==1:           dl    = d
    if   dl > 10:         n_min = 100
    elif dl<=10 and dl>5: n_min = 250
    else:                 n_min = 500
    if cntr<n_min and not (n_left<n_min): continue
    if cntr>n_min or abs(dl-d)>0.005:
      if n_left<n_min:
        i_seqs.append(size-1)
        cntr+=n_left
        reached_end=True
      bins.append([flex.min(i_seqs), 1+flex.max(i_seqs)])
      i_seqs = flex.int()
      cntr=0
    if reached_end: break
  return bins

class binned_data(object):

  def __init__(self, binner, data, data_fmt=None):
    self.binner = binner
    self.data = data
    self.data_fmt = data_fmt

  def show(self,
        data_fmt=None,
        show_bin_number=True,
        show_bin_range=True,
        show_d_range=None,
        show_counts=True,
        show_unused=True,
        bin_range_as="d",
        wavelength=None,
        f=None,
        prefix=""):
    if (data_fmt is None): data_fmt = self.data_fmt
    return self.binner.show_data(
      data=self.data,
      data_fmt=data_fmt,
      show_bin_number=show_bin_number,
      show_bin_range=show_bin_range,
      show_d_range=show_d_range,
      show_counts=show_counts,
      show_unused=show_unused,
      bin_range_as=bin_range_as,
      wavelength=wavelength,
      f=f,
      prefix=prefix)

  def as_simple_table(self,
        data_label,
        data_fmt=None):
    """
    Export table rows for display elsewhere.  (Used in Xtriage)
    """
    if (data_fmt is None): data_fmt = self.data_fmt
    return self.binner.as_simple_table(
      data=self.data,
      data_label=data_label,
      data_fmt=data_fmt)

def make_lookup_dict(indices): # XXX push to C++
  result = {}
  for i in range(len(indices)):
    result[indices[i]] = i
  return result

def get_offset_list(layers=1, include_origin=False):
  # Ready to run
  if layers==-1:
    return [ [0,0,0] ]

  elif layers==0:
    offset_list=[
     [-1,0,0],[1,0,0],
     [0,-1,0],[0,1,0],
     [0,0,-1],[0,0,1],
      ]
    if include_origin: offset_list.append([0,0,0])
    return offset_list

  offset_list=[]
  for ix in range(-layers,layers+1):
    for iy in range(-layers,layers+1):
      for iz in range(-layers,layers+1):
        if include_origin or ix!=0 or iy!=0 or iz!=0:
          offset_list.append((ix,iy,iz))
  return offset_list

def offset_indices(array,offset=None):
  from scitbx.matrix import col
  offset=col(offset)
  indices=flex.miller_index(
      list(flex.vec3_int((
         array.indices().as_vec3_double()+offset).iround().as_int())))
  return array.customized_copy(indices=indices,data=array.data())

class set(crystal.symmetry):
  """
  Basic class for handling sets of Miller indices (h,k,l), including sorting
  and matching functions, symmetry handling, generation of R-free flags, and
  extraction of associated statistics.  Does not actually contain data, but
  this can be added using the array(...) method.
  """

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
    """
    Return a reference to the internal array of indices.

    :returns: a cctbx.array_family.flex.miller_index array
    """
    return self._indices

  def anomalous_flag(self):
    """Indicate whether the set or array is anomalous or not."""
    return self._anomalous_flag

  def size(self):
    """Return the number of reflections in the set or array."""
    return self.indices().size()

  def copy(self):
    """
    Create a copy of the set, keeping references to the same crystal symmetry
    and indices.
    """
    return set(
      crystal_symmetry=self,
      indices=self._indices,
      anomalous_flag=self._anomalous_flag)

  def deep_copy(self):
    """
    Create a copy of the set, also copying the crystal symmetry and indices.

    :returns: a set object with all-new attributes.
    """
    unit_cell = self.unit_cell()
    if (unit_cell is not None):
      unit_cell = uctbx.unit_cell(parameters=unit_cell.parameters())
    if (self.space_group_info() is None):
      space_group_symbol = None
    else:
      space_group_symbol = str(self.space_group_info())
    return set(
      crystal_symmetry=crystal.symmetry(
        unit_cell=unit_cell,
        space_group_symbol=space_group_symbol),
      indices=self.indices().deep_copy(),
      anomalous_flag=self.anomalous_flag())

  def customized_copy(self,
        crystal_symmetry=Keep,
        indices=Keep,
        anomalous_flag=Keep,
        unit_cell=Keep,
        space_group_info=Keep):
    """
    Create a copy of the set, optionally changing the symmetry, indices,
    and/or anomalous flag (default = keep all unmodified).
    """
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
    """
    Create an array object, given data and/or sigma arrays of identical
    dimensions to the indices array.

    :param data: a flex array (any format) or None
    :param sigmas: a flex array (any format, but almost always double) or None
    """
    if (data is not None):
      assert data.size() == self._indices.size()
    if (sigmas is not None):
      assert sigmas.size() == self._indices.size()
    return array(miller_set=self, data=data, sigmas=sigmas)

  def __getitem__(self, slice_object):
    assert type(slice_object) == slice
    assert self.indices() is not None
    return set(
      crystal_symmetry=self,
      indices=self.indices().__getitem__(slice_object),
      anomalous_flag=self.anomalous_flag())

  def show_summary(self, f=None, prefix=""):
    """Minimal Miller set summary"""
    if (f is None): f = sys.stdout
    print(prefix + "Number of Miller indices:", len(self.indices()), file=f)
    print(prefix + "Anomalous flag:", self.anomalous_flag(), file=f)
    crystal.symmetry.show_summary(self, f=f, prefix=prefix)
    return self

  def miller_indices_as_pdb_file(self, file_name=None, expand_to_p1=False):
    """
    Write out Miller indices as pseudo-waters for visualization.  Note that
    this treats the indices as literal coordinates (times a scale factor),
    and the distances between points will not correspond to the distances
    in reciprocal space.

    See cctbx/miller/display.py and crys3d/hklview for an alternative (but
    less lightweight) approach.
    """
    assert file_name is not None
    uc = self.unit_cell()
    if(expand_to_p1): indices = self.expand_to_p1().indices()
    else: indices = self.indices()
    h = flex.int()
    k = flex.int()
    l = flex.int()
    for i_mi, mi in enumerate(indices):
      h.append(mi[0])
      k.append(mi[1])
      l.append(mi[2])
    scale = 100
    sh,sk,sl = flex.max(flex.abs(h))*scale,flex.max(flex.abs(k))*scale,\
      flex.max(flex.abs(l))*scale
    rp = self.unit_cell().reciprocal_parameters()
    c1fmt = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P1        "
    fmt = "HETATM%5d  O   HOH %5d    %8.3f%8.3f%8.3f  1.00  1.00           O  "
    with open(file_name, "w") as of:
      for i_mi, mi in enumerate(indices):
        rsv = uc.reciprocal_space_vector(mi)
        print(fmt%(i_mi, i_mi, rsv[0]*scale, rsv[1]*scale, rsv[2]*scale), file=of)

  def show_comprehensive_summary(self, f=None, prefix=""):
    """Display comprehensive Miller set or array summary"""
    if (f is None): f = sys.stdout
    self.show_summary(f=f, prefix=prefix)
    no_sys_abs = self.copy()
    if (self.space_group_info() is not None):
      is_unique_set_under_symmetry = no_sys_abs.is_unique_set_under_symmetry()
      sys_absent_flags = self.sys_absent_flags().data()
      n_sys_abs = sys_absent_flags.count(True)
      print(prefix + "Systematic absences:", n_sys_abs, file=f)
      if (n_sys_abs != 0):
        no_sys_abs = self.select(selection=~sys_absent_flags)
        print(prefix + "Systematic absences not included in following:", file=f)
      n_centric = no_sys_abs.centric_flags().data().count(True)
      print(prefix + "Centric reflections:", n_centric, file=f)
    if (self.unit_cell() is not None):
      d_max_min = no_sys_abs.resolution_range()
      print(prefix + "Resolution range: %.6g %.6g" % d_max_min, file=f)
      if (self.space_group_info() is not None
          and self.indices().size() > 0
          and is_unique_set_under_symmetry):
        no_sys_abs.setup_binner(n_bins=1)
        completeness_d_max_d_min = no_sys_abs.completeness(use_binning=True)
        binner = completeness_d_max_d_min.binner
        assert binner.counts_given()[0] == 0
        assert binner.counts_given()[2] == 0
        n_obs = binner.counts_given()[1]
        n_complete = binner.counts_complete()[1]
        if (n_complete != 0 and d_max_min[0] != d_max_min[1]):
          print(prefix + "Completeness in resolution range: %.6g" % (
            n_obs / n_complete), file=f)
        n_complete += binner.counts_complete()[0]
        if (n_complete != 0):
          print(prefix + "Completeness with d_max=infinity: %.6g" % (
            n_obs / n_complete), file=f)
        if (self.anomalous_flag()) and (self.is_xray_intensity_array() or
            self.is_xray_amplitude_array()):
          print(prefix + \
            "Anomalous completeness in resolution range: %.6g" % \
            self.anomalous_completeness(), file=f)
    if (self.space_group_info() is not None
        and no_sys_abs.anomalous_flag()
        and is_unique_set_under_symmetry):
      asu, matches = no_sys_abs.match_bijvoet_mates()
      print(prefix + "Bijvoet pairs:", matches.pairs().size(), file=f)
      print(prefix + "Lone Bijvoet mates:", \
        matches.n_singles() - n_centric, file=f)
      if (isinstance(self, array) and self.is_real_array()):
        print(prefix + "Mean anomalous difference: %.4f" % (
          no_sys_abs.anomalous_signal()), file=f)
    return self

  def show_completeness(self, reflections_per_bin = 500, out = None):
    """
    Display the completeness in resolution bins.
    """
    if(out is None): out = sys.stdout
    self.setup_binner(reflections_per_bin = reflections_per_bin)
    for i_bin in self.binner().range_used():
      sel         = self.binner().selection(i_bin)
      self_sel    = self.select(sel)
      d_max,d_min = self_sel.d_max_min()
      compl       = self_sel.completeness(d_max = d_max)
      n_ref       = sel.count(True)
      d_range     = self.binner().bin_legend(
                 i_bin = i_bin, show_bin_number = False, show_counts = False)
      fmt = "%3d: %-17s %4.2f %6d"
      print(fmt % (i_bin,d_range,compl,n_ref), file=out)
      out.flush()

  def reflection_intensity_symmetry(self):
    assert self.anomalous_flag() is False or self.anomalous_flag() is True
    return self.customized_copy(
      crystal_symmetry=crystal.symmetry.reflection_intensity_symmetry(self,
        anomalous_flag=self.anomalous_flag()))

  def sys_absent_flags(self, integral_only=False):
    """
    Generate a boolean Miller array flagging those reflections which are
    systematically absent under the current symmetry.
    """
    effective_group = self.space_group()
    if (integral_only):
      effective_group = effective_group \
        .build_derived_reflection_intensity_group(
          anomalous_flag=self.anomalous_flag())
    return self.array(data=effective_group.is_sys_absent(self.indices()))

  def centric_flags(self):
    """
    Generate a boolean Miller array flagging centric reflections.
    """
    return array(
      self,
      self.space_group().is_centric(self.indices()))

  def multiplicities(self):
    """
    Generate a Miller array (with integer data) indicating the multiplicity of
    each unique reflection.  (If the set or array is already symmetry-unique,
    the multiplicity will be 1 for every reflection.)

    :returns: array object with flex.int data
    """
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

  def d_star_cubed(self):
    return array(
      self, flex.pow(self.unit_cell().d_star_sq(self.indices()), 3/2))

  def d_spacings(self):
    """
    Generate a double Miller array containing the resolution d of each
    index.
    Note: if f000 term is present, the resolution of that term will be -1
     and not infinity.
    """
    return array(
      self, self.unit_cell().d(self.indices()))

  def sin_theta_over_lambda_sq(self):
    return array(
      self, self.unit_cell().stol_sq(self.indices()))

  def two_theta(self, wavelength, deg=False):
    """
    Generate a double Miller array containing the scattering angle of each
    index.
    """
    return array(
      self, self.unit_cell().two_theta(self.indices(), wavelength, deg))

  def d_min(self):
    """
    High-resolution limit.
    :returns: Python float
    """
    return uctbx.d_star_sq_as_d(self.unit_cell().max_d_star_sq(self.indices()))

  def min_max_d_star_sq(self):
    return self.unit_cell().min_max_d_star_sq(self.indices())

  def d_max_min(self, d_max_is_highest_defined_if_infinite = False):
    """
    Low- and high-resolution limits.
    :returns: Python tuple of floats
    Modified 2020-10-02 to allow return of maximum defined instead of -1
        if F000 present
    """
    if d_max_is_highest_defined_if_infinite:
      (d_max,d_min) = tuple([uctbx.d_star_sq_as_d(d_star_sq)
        for d_star_sq in self.min_max_d_star_sq()])
      if d_max < 0:  # (0,0,0) is present
        indices_copy = list(self.indices())
        index = indices_copy.index((0,0,0))
        new_indices = flex.miller_index(
          indices_copy[:index] + indices_copy[index+1:])
        d_max_d_star_sq,d_min_d_star_sq= self.unit_cell(
             ).min_max_d_star_sq(new_indices)
        (d_max, d_min )= (
          uctbx.d_star_sq_as_d(d_max_d_star_sq),
          uctbx.d_star_sq_as_d(d_min_d_star_sq))
      return (d_max, d_min)
    else: # usual
      return tuple([uctbx.d_star_sq_as_d(d_star_sq)
        for d_star_sq in self.min_max_d_star_sq()])

  def index_span(self):
    return index_span(self.indices())

  def min_max_indices(self):
    """Return the range of h,k,l indices"""
    span = self.index_span()
    return (span.min(), span.max())

  def at_first_index(self, ary, miller_index):
    """
    Returns the element `ary` coresponding to the `miller_index` if
    `miller_index exists, otherwise returns None.

    :param miller_index: Miller index as a 3-tuple
    :type miller_index: tuple
    :param ary: any array (e.g. self.data(), self.sigmas())
    :type ary: sequence (list, array, etc)
    :returns: type of contents of `ary`, or None
    """
    first = self._indices.first_index(miller_index)
    if first is None:
      result = None
    else:
      result = ary[first]
    return result

  def first_index(self, miller_index):
    """
    Returns the first index of the item  matching
    `miller_index`. If the `miller_index` is not found in `self`,
    then returns ``None``.

    :param miller_index: Miller index as a 3-tuple
    :type miller_index: tuple
    :returns: int, None -- index of first occurrence of
              `miller_index` or None
    """
    return self._indices.first_index(miller_index)

  def d_min_along_a_b_c_star(self):
    """
    Returns the effective resolution limits along the reciprocal space axes.
    """
    min_mi, max_mi = self.min_max_indices()
    max_h = max(abs(min_mi[0]), max_mi[0])
    max_k = max(abs(min_mi[1]), max_mi[1])
    max_l = max(abs(min_mi[2]), max_mi[2])
    assert (not None in [max_h, max_k, max_l])
    ast,bst,cst = self.unit_cell().reciprocal_parameters()[:3]
    return (1./(max_h*ast), 1./(max_k*bst), 1./(max_l*cst))

  def minimum_wavelength_based_on_d_min(self, tolerance=1e-2):
    return 2 * self.d_min() * (1-tolerance)

  def resolution_range(self):
    """Synonym for d_max_min()."""
    return self.d_max_min()

  def debye_waller_factors(self,
        u_iso=None,
        b_iso=None,
        u_cart=None,
        b_cart=None,
        u_cif=None,
        u_star=None,
        exp_arg_limit=50,
        truncate_exp_arg=False):
    """
    Given an isotropic or anisotropic displacement or B-factor, alculate
    resolution-dependent scale factors corresponding to the indices.
    (Note that to simply apply one of the input parameters to an existing
    Miller array, you can call array.apply_debye_waller_factors)

    :param u_iso: Isotropic displacement (in Angstroms)
    :param b_iso: Isotropic B-factor (8*pi^2*u_iso^2)
    :param u_cart: Anisotropic displacement tensor
    :param b_cart: Anisotropic B-factor
    :param u_star: Anisotropic displacement tensor in fractional space
    :param u_cif: Anisotropic displacement tensor, dimensionless basis
    :returns: cctbx.miller.array object
    """
    return self.array(data=self.unit_cell().debye_waller_factors(
      miller_indices=self.indices(),
      u_iso=u_iso, b_iso=b_iso,
      u_cart=u_cart, b_cart=b_cart,
      u_cif=u_cif, u_star=u_star,
      exp_arg_limit=exp_arg_limit,
      truncate_exp_arg=truncate_exp_arg))

  def n_bijvoet_pairs(self):
    """Return the number of Bijvoet pairs."""
    asu, matches = self.match_bijvoet_mates(
      assert_is_unique_set_under_symmetry=False)
    return matches.pairs().size()

  def as_non_anomalous_set(self):
    """
    Return a copy of the set using the same indices but with the anomalous flag
    set to false.
    """
    return set(
      crystal_symmetry=self,
      indices=self.indices(),
      anomalous_flag=False)

  def as_anomalous_set(self):
    """
    Return a copy of the set using the same indices but with the anomalous flag
    set to true.
    """
    return set(
      crystal_symmetry=self,
      indices=self.indices(),
      anomalous_flag=True)

  def auto_anomalous(self, min_n_bijvoet_pairs=None,
                           min_fraction_bijvoet_pairs=None):
    """
    Set the anomalous flag automatically depending on whether the data
    contain Bijvoet pairs (optionally given minimum cutoffs).

    :returns: a copy of the set with (maybe) a new anomalous flag
    """
    assert [min_n_bijvoet_pairs, min_fraction_bijvoet_pairs].count(None) > 0
    if (self.indices().size() == 0):
      anomalous_flag = False
    elif (min_fraction_bijvoet_pairs is not None):
      anomalous_flag = (2*self.n_bijvoet_pairs()/self.indices().size()
                        >= min_fraction_bijvoet_pairs)
    elif (min_n_bijvoet_pairs is not None):
      anomalous_flag = (self.n_bijvoet_pairs() >= min_n_bijvoet_pairs)
    else:
      anomalous_flag = (self.n_bijvoet_pairs() > 0)
    return set(
      crystal_symmetry=self,
      indices=self.indices(),
      anomalous_flag=anomalous_flag)

  def is_unique_set_under_symmetry(self):
    """
    Determine whether the indices in the set are symmetry-unique.
    """
    return ext.is_unique_set_under_symmetry(
      space_group_type=self.space_group_info().type(),
      anomalous_flag=self.anomalous_flag(),
      miller_indices=self.indices())

  def unique_under_symmetry_selection(self):
    return ext.unique_under_symmetry_selection(
      space_group_type=self.space_group_info().type(),
      anomalous_flag=self.anomalous_flag(),
      miller_indices=self.indices())

  def unique_under_symmetry(self):
    sel = self.unique_under_symmetry_selection()
    if (sel.size() == self.indices().size()): return self
    return self.select(sel)

  def is_in_asu(self):
    """
    Indicate whether the array is entirely contained within the reciprocal
    space asymmetric unit (ASU).  Warning: this calls map_to_asu internally,
    so it is better to simply call map_to_asu without checking in many cases.
    """
    #XXX could be made more efficient
    return self.map_to_asu().indices().all_eq(self.indices())

  def map_to_asu(self):
    """
    Convert all indices to lie within the canonical asymmetric unit for the
    current space group (while preserving anomalous flag).  Required for many
    downstream steps.
    """
    i = self.indices().deep_copy()
    anomalous_flag = self.anomalous_flag()
    if (anomalous_flag is None):
      anomalous_flag = True
    map_to_asu(
      self.space_group_info().type(),
      anomalous_flag,
      i)
    return set(self, i, self.anomalous_flag())

  def complete_set(self, d_min_tolerance=1.e-6, d_min=None, d_max=None,
                   max_index=None):
    """
    Generate the complete set of Miller indices expected for the current
    symmetry, excepting systematic absences.

    :param d_min_tolerance: tolerance factor for d_min (avoid precision errors)
    :param d_min: High-resolution limit (default = d_min of current set)
    :param d_max: Low-resolution limit (default = d_max of current set)
    """
    assert self.anomalous_flag() in (False, True)
    if (self.indices().size() == 0):
      return set(
        crystal_symmetry=self,
        anomalous_flag=self.anomalous_flag(),
        indices=flex.miller_index())
    if(max_index is not None):
      return build_set(
        crystal_symmetry=self,
        anomalous_flag=self.anomalous_flag(),
        max_index=max_index)
    else:
      d_min_was_none = (d_min is None)
      if (d_min_was_none): d_min = self.d_min()
      d_min_exact = d_min
      if (d_min_tolerance is not None): d_min *= (1-d_min_tolerance)
      result = build_set(
        crystal_symmetry=self,
        anomalous_flag=self.anomalous_flag(),
        d_min=d_min,
        d_max=d_max)
      if (d_min_was_none):
        result = result.select(
          result.d_spacings().data() >= d_min_exact*(1-fp_eps_double))
      return result

  def completeness(self, use_binning=False, d_min_tolerance=1.e-6,
                   return_fail=None, d_max = None, multiplier=1,
                   as_non_anomalous_array=None):
    """
    Calculate the (fractional) completeness of the array relative to the
    theoretical complete set, either overall or in resolution bins.  By default
    the current low-resolution limit will be used.

    :param d_min_tolerance: tolerance factor for d_min (avoid precision errors)
    :param d_max: Low-resolution limit (default = d_max of current set)
    :param multiplier: Factor to multiply the result by (usually 1 or 100)
    :param as_non_anomalous_array: Report values for non-anomalous array
    """
    assert (multiplier > 0)
    if as_non_anomalous_array and self.anomalous_flag():
      # report values for non-anomalous version of array
      merged = self.as_non_anomalous_array().merge_equivalents().array()
      merged.array().set_observation_type(self)
      if use_binning:
        merged.use_binning_of(self)
      return merged.completeness(use_binning=use_binning,
        d_min_tolerance=d_min_tolerance,
        return_fail=return_fail,
        d_max = d_max,
        multiplier = multiplier,
        as_non_anomalous_array=False)

    if (not use_binning):
      complete_set = self.complete_set(d_min_tolerance=d_min_tolerance,
                                       d_max = d_max)
      return min(self.indices().size() / max(1, complete_set.indices().size()),
                 1.0) * multiplier
    assert self.binner() is not None
    data = []
    for n_given,n_complete in zip(self.binner().counts_given(),
                                  self.binner().counts_complete()):
      if (n_complete == 0): data.append(return_fail)
      else: data.append(multiplier*n_given/n_complete)
    return binned_data(binner=self.binner(), data=data, data_fmt="%5.3f")

  def all_selection(self):
    return flex.bool(self.indices().size(), True)

  def select(self, selection, negate=False, anomalous_flag=None):
    """
    Select a subset of reflections.

    :param selection: flex.bool or flex.size_t selection
    :param negate: select the inverse of the selection array
    :param anomalous_flag: anomalous flag for the new set
    :returns: a new set with a subset of indices
    """
    assert self.indices() is not None
    if (anomalous_flag is None):
      anomalous_flag = self.anomalous_flag()
    if (negate): selection = ~selection
    i = self.indices().select(selection)
    return set(self, i, anomalous_flag)

  def select_acentric(self):
    """
    Select only acentric reflections.

    :returns: A Miller set or array (depending on object type).
    """
    return self.select(~self.centric_flags().data())

  def select_centric(self):
    """
    Select only centric reflections.

    :returns: A Miller set or array (depending on object type).
    """
    return self.select(self.centric_flags().data())

  def remove_systematic_absences(self, negate=False):
    return self.select(
      selection=self.sys_absent_flags().data(),
      negate=not negate)

  def resolution_filter_selection(self, d_max=None, d_min=None):
    """
    Obtain the selection (flex.bool array) corresponding to the specified
    resolution range.
    """
    result = self.all_selection()
    d_star_sq = self.d_star_sq().data()
    assert d_star_sq.all_ge(0)
    d_star = flex.sqrt(d_star_sq)
    if (d_max is not None and d_max > 0): result &= (d_star >= 1/d_max)
    if (d_min is not None and d_min > 0): result &= (d_star <= 1/d_min)
    return result

  def resolution_filter(self, d_max=0, d_min=0, negate=0):
    """
    Select a subset within the indicated resolution range.
    Returns a new miller array (does not change existing array)

    :param d_max: Low-resolution cutoff
    :param d_min: High-resolution cutoff
    :param negate: Select outside this range instead
    :returns: set or array depending on object type
    """
    return self.select(
      selection=self.resolution_filter_selection(d_max=d_max, d_min=d_min),
      negate=negate)

  def match_indices(self, other, assert_is_similar_symmetry=True):
    if (assert_is_similar_symmetry):
      assert self.is_similar_symmetry(other)
    assert self.anomalous_flag() == other.anomalous_flag()
    return match_indices(self.indices(), other.indices())

  def common_set(self, other, assert_is_similar_symmetry=True):
    """
    Match the indices in the current set and another set, and return a set
    (or array) containing only those reflections present in both.  Assumes that
    both sets are already in the asymmetric unit (ASU).
    """
    pairs = other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry).pairs()
    return self.select(pairs.column(1))

  def common_sets(self,
        other,
        assert_is_similar_symmetry=True,
        assert_no_singles=False):
    """
    Like common_set(other), but returns a tuple containing matching copies of
    both sets (or arrays).
    """
    matches = other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry)
    if (assert_no_singles):
      assert not matches.have_singles()
    pairs = matches.pairs()
    return [self.select(pairs.column(1)),
            other.select(pairs.column(0))]

  def lone_set(self, other, assert_is_similar_symmetry=True):
    """
    Match the indices in the current set and another set, and return a set
    (or array) containing reflections which are present only in the current
    set.  Assumes that both sets are already in the asymmetric unit.
    """
    return self.select(other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry).singles(1))

  def lone_sets(self, other, assert_is_similar_symmetry=True):
    """
    Like lone_set(other), but returns a tuple containing the reflections
    unique to each set (or array).
    """
    matches = other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry)
    return [self.select(matches.singles(1)),
            other.select(matches.singles(0))]

  def match_bijvoet_mates(self, assert_is_unique_set_under_symmetry=True):
    """
    Group Bijvoet mates (or Friedel mates) together, returning an object that
    allows enumeration over matching pairs and/or singletons.
    """
    assert self.anomalous_flag() in (None, True)
    assert self.indices() is not None
    if (self.space_group_info() is not None):
      asu = self.map_to_asu()
      matches = match_bijvoet_mates(
        asu.space_group_info().type(), asu.indices(),
        assert_is_unique_set_under_symmetry=assert_is_unique_set_under_symmetry)
    else:
      asu = self
      matches = match_bijvoet_mates(
        asu.indices(),
        assert_is_unique_set_under_symmetry=assert_is_unique_set_under_symmetry)
    return asu, matches

  def sort_permutation(self, by_value="resolution", reverse=False):
    """
    Generate the selection array (flex.size_t object) to reorder the array
    by resolution or Miller index.

    :param by_value: sort type, must be "resolution" or "packed_indices"
    :param reverse: invert order
    :returns: flex.size_t object
    """
    assert by_value in ["resolution", "packed_indices","asu_indices"]
    assert reverse in [False, True]
    if (by_value == "resolution"):
      return flex.sort_permutation(
        data=self.unit_cell().d_star_sq(self.indices()),
        reverse=reverse, stable=True)
    elif (by_value == "asu_indices"):
      # sort on asu_indices, keeping anom pairs together, group by Friedel mates
      asu_indices=self.indices().deep_copy()
      map_to_asu(self.space_group_info().type(),False,asu_indices)
      asu_indices_anom=self.indices().deep_copy()
      map_to_asu(self.space_group_info().type(),True,asu_indices_anom)

      centric_flags = self.centric_flags().data()
      friedel_mate_flags = ~(asu_indices == asu_indices_anom)
      offset=flex.double(self.indices().size(), 0)
      offset.set_selected(friedel_mate_flags.iselection(),0.01)
      data=index_span(asu_indices).pack(asu_indices).as_double()+offset
      return flex.sort_permutation(data=data, reverse=reverse, stable=True)

    else:
      return flex.sort_permutation(
        data=index_span(self.indices()).pack(self.indices()),
        reverse=reverse, stable=True)

  def sort(self, by_value="resolution", reverse=False):
    """
    Reorder reflections by resolution or Miller index.

    :param by_value: 'resolution' or 'packed_indices'
    """
    return self.select(
      self.sort_permutation(by_value=by_value, reverse=reverse))

  def generate_r_free_flags(self,
        fraction=0.1,
        max_free=2000,
        lattice_symmetry_max_delta=5.0,
        use_lattice_symmetry=True,
        use_dataman_shells=False,
        n_shells=20,
        format="cns"):
    """
    Create an array of R-free flags for the current set, keeping anomalous
    pairs together.  Requires that the set already be unique under symmetry,
    and generally assumes that the set is in the ASU.

    :param fraction: fraction of reflections to flag for the test set
    :param max_free: limit on size of test set, overrides fraction
    :param lattice_symmetry_max_delta: limit on lattice symmetry calculation
    :param use_lattice_symmetry: given the current symmetry, determine the \
    highest possible lattice symmetry and generate flags for that symmetry, \
    then expand to the current (lower) symmetry if necessary.  This is almost \
    always a good idea.
    :param use_dataman_shells: generate flags in thin resolution shells to \
    avoid bias due to non-crystallographic symmetry.
    :param n_shells: number of resolution shells if use_dataman_shells=True
    :param format: convention of the resulting flags.  'cns' will return a \
    boolean array (True = free), 'ccp4' will return an integer array from \
    0 to X (0 = free, X dependent on fraction), 'shelx' will return an \
    integer array with values 1 (work) or -1 (free).

    :returns: a boolean or integer Miller array, depending on format.
    """
    assert (format in ["cns", "ccp4", "shelx"])
    if use_lattice_symmetry:
      assert lattice_symmetry_max_delta>=0

    if not use_lattice_symmetry:
      return self.generate_r_free_flags_basic(fraction=fraction,
        max_free=max_free,
        use_dataman_shells=use_dataman_shells,
        n_shells=n_shells,
        format=format)
    else:
      return self.generate_r_free_flags_on_lattice_symmetry(fraction=fraction,
        max_free=max_free,
        max_delta=lattice_symmetry_max_delta,
        use_dataman_shells=use_dataman_shells,
        n_shells=n_shells,
        format=format)

  def crystal_symmetry(self):
    """Get crystal symmetry of the miller set

    :returns: a new crystal.symmetry object
    :rtype: cctbx.crystal.symmetry
    """
    return crystal.symmetry(
      unit_cell = self.unit_cell(),
      space_group_info = self.space_group_info())

  def generate_r_free_flags_on_lattice_symmetry(self,
        fraction=0.10,
        max_free=2000,
        max_delta=5.0,
        return_integer_array=False,
        n_partitions=None,
        use_dataman_shells=False,
        n_shells=20,
        format="cns"):
    """
    Generate R-free flags by converting to the highest possible lattice
    symmetry (regardless of intensity symmetry), creating flags, and expanding
    back to the original symmetry.  This is a safeguard against reflections
    that are correlated due to twinning being split between the work and test
    sets.

    This method should usually not be called directly, but rather through
    set.generate_r_free_flags(...).
    """
    # the max_number of reflections is wrst the non anomalous set
    n_original = self.indices().size()
    if n_original<=0:
      raise Sorry("An array of size zero is given for Free R flag generation")
    n_non_ano = n_original
    if self.anomalous_flag():
      matches = self.match_bijvoet_mates()[1]
      n_non_ano = matches.pairs().size() + matches.n_singles()
    assert self.is_unique_set_under_symmetry()
    if fraction is not None:
      if not (fraction > 0 and fraction < 0.5):
        raise Sorry("R-free flags fraction must be greater than 0 and less "+
          "than 0.5.")
      assert n_partitions is None
      assert return_integer_array is False
    if max_free is not None:
      assert fraction is not None
      fraction = min( n_non_ano*fraction,
                      max_free )/float(n_non_ano)
      assert n_partitions is None
      assert return_integer_array is False
    if return_integer_array:
      assert not use_dataman_shells
      assert fraction is None
      assert max_free is None
      assert n_partitions > 1
    #first make a set of temporary flags
    cb_op_to_niggli = self.change_of_basis_op_to_niggli_cell()
    tmp_ma = self.change_basis( cb_op_to_niggli )
    # please get the lattice symmetry of th niggli cell
    lattice_group = lattice_symmetry.group(
      tmp_ma.unit_cell(),
      max_delta=max_delta)
    lattice_xs = crystal.symmetry(unit_cell=tmp_ma.unit_cell(),
                                  space_group=lattice_group,
                                  assert_is_compatible_unit_cell=False)
    # make some flags, and insert lattice symmetry
    tmp_flags = tmp_ma.array().customized_copy(
      crystal_symmetry = lattice_xs,
      data = flex.double( tmp_ma.indices().size(), 0 ) ).map_to_asu()
    # Carry out the merging please
    tmp_flags = tmp_flags.merge_equivalents().array()
    tmp_flags = tmp_flags.average_bijvoet_mates()
    n_non_ano_lat_sym = tmp_flags.indices().size()
    # now we can do the free r assignement in the lattice symmetry
    n = tmp_flags.indices().size()
    result = None
    if not return_integer_array:
      if use_dataman_shells :
        if (format == "ccp4"):
          raise Sorry("CCP4 convention not available when generating R-free "+
            "flags in resolution shells.")
        result = r_free_utils.assign_r_free_flags_by_shells(
          n_refl=n,
          fraction_free=fraction,
          n_bins=n_shells)
      else :
        result = r_free_utils.assign_random_r_free_flags(n, fraction,
          format=format)
    else:
      result = flex.size_t()
      n_times = int(n/n_partitions)+1
      for ii in range(n_times):
        tmp = flex.random_double( n_partitions )
        tmp = flex.sort_permutation( tmp )
        result.extend( tmp )
      result = flex.int( list(result[0:n]) )

    # please sort the reflections by resolution
    indices = tmp_flags.indices()
    result = result.select(
      indices=flex.sort_permutation(
        data=tmp_flags.unit_cell().d_star_sq(indices), reverse=True),
        reverse=True)

    tmp_flags = tmp_flags.customized_copy(
      data=result, sigmas=None)
    # now expand to p1 please
    tmp_flags = tmp_flags.expand_to_p1()
    # now make it into the proper symmetry please
    tmp_flags = tmp_flags.customized_copy(
      crystal_symmetry = crystal.symmetry( unit_cell=tmp_ma.unit_cell(),
                                           space_group=tmp_ma.space_group() )
      )
    tmp_flags = tmp_flags.merge_equivalents().array()
    if self.anomalous_flag():
      tmp_flags = tmp_flags.generate_bijvoet_mates()
    tmp_flags = tmp_flags.change_basis( cb_op_to_niggli.inverse() ).map_to_asu()
    tmp_flags = tmp_flags.common_set( self.map_to_asu() )
    tmp_flags = tmp_flags.customized_copy(
      indices = self.indices(),
      data = tmp_flags.data() )
    tmp_flags = tmp_flags.common_set( self )
    assert tmp_flags.indices().all_eq( self.indices() )
    return tmp_flags

  def generate_r_free_flags_basic(self,
        fraction=0.1,
        max_free=2000,
        use_dataman_shells=False,
        n_shells=20,
        format="cns"):
    """
    Generate R-free flags, without taking lattice symmetry into account (not
    recommended).  Should not normally be called directly - use
    generate_r_free_flags(...) instead.
    """
    if not (fraction > 0 and fraction < 0.5):
      raise Sorry("R-free flags fraction must be greater than 0 and less "+
          "than 0.5.")
    if (use_dataman_shells and n_shells < 5):
      raise Sorry("You must use at least 5 resolution shells when assigning "+
        "R-free flags this way.")
    if (max_free is not None) and (max_free <= 0):
      raise Sorry("The maximum number of free reflections must either be "+
        "None, or a positive number.")
    if (self.anomalous_flag()):
      matches = self.match_bijvoet_mates()[1]
      sel_pp = matches.pairs_hemisphere_selection("+")
      sel_pm = matches.pairs_hemisphere_selection("-")
      sel_sp = matches.singles("+")
      sel_sm = matches.singles("-")
      n = matches.pairs().size() + matches.n_singles()
      del matches
    else:
      assert self.is_unique_set_under_symmetry()
      n = self.indices().size()
    if (max_free is not None):
      fraction = min(fraction, max_free/max(1,n))
    if use_dataman_shells :
      if (format == "ccp4"):
        raise Sorry("CCP4 convention not available when generating R-free "+
          "flags in resolution shells.")
      result = r_free_utils.assign_r_free_flags_by_shells(
        n_refl=n,
        fraction_free=fraction,
        n_bins=n_shells)
    else :
      result = r_free_utils.assign_random_r_free_flags(
        n_refl=n,
        fraction_free=fraction,
        format=format)
    if (not self.anomalous_flag()):
      indices = self.indices()
    else:
      indices = self.indices().select(sel_pp)
      indices.extend(self.indices().select(sel_sp))
      indices.extend(self.indices().select(sel_sm))
      assert indices.size() == n
    result = result.select(
      indices=flex.sort_permutation(
        data=self.unit_cell().d_star_sq(indices), reverse=True),
        reverse=True)
    if (not self.anomalous_flag()):
      return self.array(data=result)
    del indices
    if (format == "ccp4"):
      test_flag_value = 0
      # XXX are we sure this is the right thing to do?
      result_full = flex.int(self.indices().size(), 1)
    elif (format == "shelx"):
      test_flag_value = -1
      result_full = flex.int(self.indices().size(), 1)
    else :
      test_flag_value = True
      result_full = flex.bool(self.indices().size(), False)
    i_pp = sel_pp.size()
    i_pp_sp = i_pp + sel_sp.size()
    r_pp = result[:i_pp]
    result_full.set_selected(sel_pp, r_pp)
    assert result_full.count(test_flag_value) == r_pp.count(test_flag_value)
    result_full.set_selected(sel_pm, r_pp)
    assert result_full.count(test_flag_value) == 2*r_pp.count(test_flag_value)
    del r_pp
    del sel_pm
    del sel_pp
    result_full.set_selected(sel_sp, result[i_pp:i_pp_sp])
    del sel_sp
    result_full.set_selected(sel_sm, result[i_pp_sp:])
    del sel_sm
    return self.array(data=result_full)

  def random_phases_compatible_with_phase_restrictions(self, deg=False):
    random_phases = flex.random_double(size=self.size())-0.5
    if (deg): random_phases *= 360
    else:     random_phases *= 2*math.pi
    return self.array(data=self.space_group().nearest_valid_phases(
      miller_indices=self.indices(),
      phases=random_phases,
      deg=deg))

  def change_basis(self, cb_op):
    """Get a transformation of the miller set with a new basis specified by cb_op

    :param cb_op: object describing the desired transformation of the basis
    :type cb_op: string or sgtbx.change_of_basis_operator

    :returns: a new miller set with the new basis
    :rtype: cctbx.miller.set
    """
    if (isinstance(cb_op, str)): cb_op = sgtbx.change_of_basis_op(cb_op)
    return set.customized_copy(self,
      crystal_symmetry=crystal.symmetry.change_basis(self, cb_op),
      indices=cb_op.apply(self.indices()))

  def expand_to_p1_iselection(self, build_iselection=True):
    assert self.space_group_info() is not None
    assert self.indices() is not None
    assert self.anomalous_flag() is not None
    return expand_to_p1_iselection(
      space_group=self.space_group(),
      anomalous_flag=self.anomalous_flag(),
      indices=self.indices(),
      build_iselection=build_iselection)

  def expand_to_p1(self, return_iselection=False):
    """Get a transformation of the miller set to spacegroup P1

    :returns: a new set of parameters (symmetry, miller indices, anomalous_flag) in spacegroup P1
    :rtype: set(cctbx.crystal.symmetry, cctbx.miller.indices, boolean)
    """
    proxy = self.expand_to_p1_iselection(build_iselection=return_iselection)
    result = set(
      crystal_symmetry=self.cell_equivalent_p1(),
      indices=proxy.indices,
      anomalous_flag=self.anomalous_flag())
    if (return_iselection):
      return result, proxy.iselection
    return result

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
    """
    Calculate real-space grid for FFT given array crystal symmetry, d_min, and
    desired resolution-dependent spacing.  The actual grid dimensions will be
    adjusted to suit the needs of the FFT algorithm.
    """
    if (d_min is None and grid_step is None): d_min = self.d_min()
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

  def structure_factors_from_asu_map(self, asu_map_data, n_real):
    asu_m = asu_map_ext.asymmetric_map(
      self.space_group().type(), asu_map_data, n_real)
    indices = self.indices()
    return self.customized_copy(
      indices = indices,
      data    = asu_m.structure_factors(indices)*self.unit_cell().volume())

  def structure_factors_from_map(self, map, in_place_fft=False,
      use_scale=False, anomalous_flag=None, use_sg=False):
    """
    Run FFT on a real-space map to calculate structure factors corresponding to
    the current set of Miller indices.

    :param map: flex.double map with 3D flex.grid accessor
    :param in_place_fft: perform the FFT in place instead of creating a copy of
      the map first
    :param use_scale: perform volume-scaling using current unit cell
    :param anomalous_flag: determines anomalous_flag for output array
    :param use_sg: use space-group symmetry
    :returns: array with same Miller indices and complex_double data
    """
    assert map.focus_size_1d() > 0 and map.nd() == 3 and map.is_0_based()
    if(isinstance(map, flex.int)): map = map.as_double()
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
    if(use_scale):
      scale = self.unit_cell().volume() \
        / matrix.col(fft.n_real()).product()
      map *= scale
    if(anomalous_flag is None):
      anomalous_flag = self.anomalous_flag()
    if(use_sg):
      from_map = maptbx.structure_factors.from_map(
        space_group=self.space_group(),
        anomalous_flag=anomalous_flag,
        miller_indices=self.indices(),
        complex_map=map,
        conjugate_flag=True)
    else:
      from_map = maptbx.structure_factors.from_map(
        anomalous_flag=anomalous_flag,
        miller_indices=self.indices(),
        complex_map=map,
        conjugate_flag=True)
    data = from_map.data()
    return array(miller_set=self, data=data)

  def structure_factors_from_scatterers(self, xray_structure,
                                        algorithm=None,
                                        cos_sin_table=False,
                                        grid_resolution_factor=1/3.,
                                        quality_factor=None,
                                        u_base=None,
                                        b_base=None,
                                        wing_cutoff=None,
                                        exp_table_one_over_step_size=None,
                                        extra_params=None):
    """
    Calculate structure factors for an :py:class:`cctbx.xray.structure` object
    corresponding to the current set of Miller indices.  Can use either FFT
    or direct summation.

    :param xray_structure: :py:class:`cctbx.xray.structure` object
    :param algorithm: switch method to calculate structure factors - can be
      'direct' or 'fft'
    :returns: array with same Miller indices and complex_double data
    """
    from cctbx import xray
    if (algorithm == "direct"):
      return xray.structure_factors.from_scatterers_direct(
        xray_structure=xray_structure,
        miller_set=self,
        cos_sin_table=cos_sin_table)
    return xray.structure_factors.from_scatterers(
      miller_set=self,
      cos_sin_table=cos_sin_table,
      grid_resolution_factor=grid_resolution_factor,
      quality_factor=quality_factor,
      u_base=u_base,
      b_base=b_base,
      wing_cutoff=wing_cutoff,
      exp_table_one_over_step_size=exp_table_one_over_step_size)(
        xray_structure=xray_structure,
        miller_set=self,
        algorithm=algorithm,
        extra_params=extra_params)

  def amplitude_normalisations(self, asu_contents, wilson_plot):
    """ A miller.array whose data N(h) are the normalisations to convert
    between E's and F's:
    E(h) = F(h) / N(h)
    The argument wilson_plot shall feature attributes
    - wilson_intensity_scale_factor
    - wilson_b
    """
    from cctbx import eltbx
    multiplicities = flex.double()
    gaussians = shared_gaussian_form_factors()
    for chemical_type, multiplicy in asu_contents.items():
      gaussians.append(eltbx.xray_scattering.wk1995(
        chemical_type).fetch())
      multiplicities.append(multiplicy)
    data = ext.amplitude_normalisation(
      form_factors=gaussians,
      multiplicities=multiplicities,
      wilson_intensity_scale_factor=wilson_plot.wilson_intensity_scale_factor,
      wilson_b=wilson_plot.wilson_b,
      indices=self.indices(),
      unit_cell=self.unit_cell(),
      space_group=self.space_group(),
    ).normalisations
    return array(self, data)

  def f_obs_minus_xray_structure_f_calc(self, f_obs_factor, xray_structure,
        structure_factor_algorithm=None,
        cos_sin_table=False,
        quality_factor=None,
        u_base=None,
        b_base=None,
        wing_cutoff=None,
        exp_table_one_over_step_size=None,
        extra_params=None):
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
        exp_table_one_over_step_size=exp_table_one_over_step_size,
        extra_params=extra_params).f_calc())

  def setup_binner(self, d_max=0, d_min=0,
                   auto_binning=False,
                   reflections_per_bin=0,
                   n_bins=0):
    """
    Create internal resolution binner object; required for many other methods
    to work.
    """
    assert auto_binning or reflections_per_bin != 0 or n_bins != 0
    assert auto_binning or (reflections_per_bin == 0 or n_bins == 0)
    if (auto_binning):
      if (reflections_per_bin == 0): reflections_per_bin = 200
      if (n_bins == 0): n_bins = 8
      n_per_bin = int(len(self.indices()) / n_bins + .5)
      if (n_per_bin > reflections_per_bin):
        n_bins = int(len(self.indices()) / reflections_per_bin + .5)
    elif (reflections_per_bin):
      n_bins = int(len(self.indices()) / reflections_per_bin + .5)
    assert n_bins > 0
    assert self.unit_cell() is not None
    assert self.indices().size() > 0 or d_min > 0
    self.use_binning(
      binning=binning(self.unit_cell(), n_bins, self.indices(), d_max, d_min))
    return self.binner()

  def safe_setup_binner(self, n_bins = None, d_max = None, d_min = None,
     min_in_bin = 10):
    # purpose: make sure there are data in all bins
    # Returns: actual n_bins created

    while n_bins >= 1:
      self.setup_binner(n_bins = n_bins, d_max = d_max, d_min = d_min)
      ok = True
      for i_bin in self.binner().range_used():
        sel = self.binner().selection(i_bin)
        if sel.count(True)<min_in_bin:
          n_bins = min(n_bins-1, max(1, int(0.5+n_bins*0.8)))
          ok = False
          break
      if ok: break  # found it
    if n_bins < 1:
      raise Sorry("Unable to set up bins...perhaps no data?")
    return n_bins  # binner is also set up in self.

  def log_binning(self, n_reflections_in_lowest_resolution_bin=100, eps=1.e-4,
                  max_number_of_bins = 30, min_reflections_in_bin=50):
    """
    Create resolution bins on a logarithmic scale.  See Urzhumtsev et al.
    (2009) Acta Crystallogr D Biol Crystallogr. 65:1283-91.
    """
    #XXX Move entire implementation into C++: can be ~100 times faster.
    assert max_number_of_bins > 1
    if(n_reflections_in_lowest_resolution_bin   >= self.indices().size() or
       n_reflections_in_lowest_resolution_bin*2 >= self.indices().size()):
      return [flex.bool(self.indices().size(), True)]
    d_spacings_sorted = self.sort().d_spacings().data()
    ss = 1./flex.pow2(d_spacings_sorted)
    lnss = flex.log(ss)
    s0 = lnss[0]
    s1 = lnss[n_reflections_in_lowest_resolution_bin]
    s2 = lnss[n_reflections_in_lowest_resolution_bin*2]
    d_spacings = self.d_spacings().data()
    step = s2-s1
    def get_limits():
      limits = [math.sqrt(1/math.exp(s0))]
      lnss_min = s1
      while lnss_min <= flex.max(lnss):
        d = math.sqrt(1/math.exp(lnss_min))
        limits.append(d)
        lnss_min += step
      lnss_min = min(flex.max(lnss), lnss_min)
      d = math.sqrt(1/math.exp(lnss_min))
      limits.append(d)
      return limits
    limits = get_limits()
    if(len(limits) > max_number_of_bins):
      # (max_number_of_bins-1) is because adding both ends; see get_limits().
      step = (lnss[len(lnss)-1]-s1)/(max_number_of_bins-1)
      limits = get_limits()
    pairs = []
    for i, d in enumerate(limits):
      if(i<len(limits)-2):
        d1,d2 = d, limits[i+1]
        pairs.append([d1,d2])
      elif(i == len(limits)-1):
        d1,d2=limits[i-1], limits[i]
        pairs.append([d1,d2])
    selections = []
    for i, p in enumerate(pairs):
      d_max, d_min = p
      if(i == 0):
        sel  = d_spacings <= d_max + eps
        sel &= d_spacings >  d_min
      elif(i == len(pairs)-1):
        sel  = d_spacings <= d_max
        sel &= d_spacings >  d_min - eps
      else:
        sel  = d_spacings <= d_max
        sel &= d_spacings >  d_min
      selections.append(sel)
    sl = selections[len(selections)-1]
    sp = selections[len(selections)-2]
    if(sl.count(True) < sp.count(True)/4):
      new_selections = selections[:len(selections)-2]
      new_selections.append(sp | sl)
      selections = new_selections
    new_sel = []
    for i,si in enumerate(selections):
      if(si.count(True)<min_reflections_in_bin and len(new_sel)-1>=0):
        new_sel[len(new_sel)-1] = new_sel[len(new_sel)-1] | si
      else:
        new_sel.append(si)
    selections = new_sel
    return selections

  def setup_binner_d_star_sq_step(self,
                                  auto_binning=True,
                                  d_max=None,
                                  d_min=None,
                                  d_star_sq_step=None):
    assert auto_binning or ( d_min is not None )
    assert auto_binning or ( d_max is not None )
    if d_star_sq_step:
      assert d_star_sq_step > 0 or (d_star_sq_step is None)
    if auto_binning:
      d_spacings = self.d_spacings().data()
      d_max=flex.max(d_spacings)
      d_min=flex.min(d_spacings)
      del d_spacings
      if d_star_sq_step is None:
        d_star_sq_step = 0.004
    assert (d_star_sq_step>0.0)
    d_min, d_max = sorted((d_min, d_max))
    return self.use_binning(binning=binning(self.unit_cell(),
      self.indices(),
      d_max,
      d_min,
      d_star_sq_step))

  def setup_binner_d_star_sq_bin_size(self,
        reflections_per_bin=1000,
        min_bins=6,
        max_bins=50,
        d_max=None,
        d_min=None,
        d_tolerance=1.e-5):
    """
    Set up bins of equal width in d_star_sq by target for mean number per bin.
    """
    (d_max_data, d_min_data) = self.d_max_min(
        d_max_is_highest_defined_if_infinite=True)
    if d_max is not None:
      d_max_work = min(d_max,d_max_data)
    else:
      d_max_work = d_max_data
    if d_min is not None:
      d_min_work = max(d_min,d_min_data)
    else:
      d_min_work = d_min_data
    if (d_max_work < d_max_data or d_min_work > d_min_data):
      working_data = self.d_spacings().resolution_filter(
          d_max=d_max_work,d_min=d_min_work)
      n_ref_work = working_data.size()
      this_d_tol = 0. # No boundary tolerance if explicit resolution limits
    else:
      n_ref_work = self.size()
      this_d_tol = d_tolerance # Avoid losing reflections with rounding errors
    n_bins = iround(max(min(n_ref_work/reflections_per_bin, max_bins), min_bins))
    limits = flex.double()
    d_star_sq_min = 1. / (d_max_work*d_max_work)
    d_star_sq_max = 1. / (d_min_work*d_min_work)
    d_star_sq_step = (d_star_sq_max - d_star_sq_min) / n_bins
    # Avoid losing reflections to rounding error on bin boundaries
    limits.append(d_star_sq_min - this_d_tol*d_star_sq_step)
    for i_bin in range(1,n_bins):
      this_limit = d_star_sq_min + i_bin*d_star_sq_step
      limits.append(this_limit)
    limits.append(d_star_sq_max + this_d_tol*d_star_sq_step)
    return self.use_binning(binning=binning(self.unit_cell(), limits))

  def setup_binner_counting_sorted(self,
        d_max=0,
        d_min=0,
        reflections_per_bin=None,
        n_bins=None,
        d_tolerance=1.e-10):
    assert d_max >= 0
    assert d_min >= 0
    assert isinstance(reflections_per_bin, int) or isinstance(n_bins, int)
    if reflections_per_bin is not None:
      assert reflections_per_bin > 0
    if n_bins is not None:
      assert n_bins > 0
    assert d_tolerance > 0
    assert d_tolerance < 0.5
    d_star_sq = self.d_star_sq().data()
    d_star_sq = d_star_sq.select(flex.sort_permutation(d_star_sq))

    # Select unique values of d_star_sq (within tolerance)
    delta = d_star_sq[1:] - d_star_sq[:-1]
    isel = (delta > d_tolerance).iselection()
    isel += 1
    isel.insert(0, 0)
    d_star_sq = d_star_sq.select(isel)

    if (d_max > 0):
      d_star_sq = d_star_sq.select(d_star_sq >= 1./d_max**2)
    if (d_min > 0):
      d_star_sq = d_star_sq.select(d_star_sq < 1./d_min**2)
    assert d_star_sq.size() > 0
    if reflections_per_bin:
      n_bins_calc = max(1, iround(d_star_sq.size() / float(reflections_per_bin)))
      n_bins = min(n_bins, n_bins_calc) if n_bins else n_bins_calc
    assert n_bins <= d_star_sq.size(), (
      "n_bins (%i) must be <= number of unique d-spacings (%i)" %(
        n_bins, d_star_sq.size())
      )
    reflections_per_bin = d_star_sq.size() / float(n_bins)
    limits = flex.double()
    limits.reserve(n_bins+1)
    if (d_max > 0):
      limits.append(1./d_max**2 * (1-d_tolerance))
    else:
      limits.append(max(0, d_star_sq[0] * (1-d_tolerance)))
    m = d_star_sq.size()-1
    for i_bin in range(1, n_bins):
      i = iround(i_bin * reflections_per_bin)
      limits.append(d_star_sq[i] * (1-d_tolerance))
      if (i == m): break
    if (d_min > 0):
      limits.append(1./d_min**2 * (1+d_tolerance))
    else:
      limits.append(d_star_sq[-1] * (1+d_tolerance))
    return self.use_binning(binning=binning(self.unit_cell(), limits))

  def binner(self):
    """
    Return a reference to the current resolution binner (or None if undefined).
    """
    return self._binner

  def use_binning(self, binning):
    """
    Use the resolution binning of the specified binner object (does not need
    to be from an identically sized set).
    """
    self._binner = binner(binning, self)
    return self._binner

  def use_binning_of(self, other):
    """
    Use the resolution binning of the specified set (does not need to be an
    identical set of indices).
    """
    return self.use_binning(binning=other.binner())

  def use_binner_of(self, other):
    """
    Use the exact binner of another set, which must have identical indices.
    """
    assert self.indices().all_eq(other.indices())
    self._binner = other._binner
    return self._binner

  def clear_binner(self):
    self._binner = None

  def concatenate(self, other, assert_is_similar_symmetry=True):
    """
    Combine two Miller sets.  Both must have the same anomalous flag, and
    similar symmetry is also assumed.
    """
    if (assert_is_similar_symmetry):
      assert self.is_similar_symmetry(other)
    assert self.anomalous_flag() == other.anomalous_flag()
    return set(
      crystal_symmetry = self,
      indices          = self._indices.concatenate(other.indices()),
      anomalous_flag   = self.anomalous_flag())

  def slice(self,
      axis=None,
      axis_index=None,
      slice_index=None,
      slice_start=None,
      slice_end=None):
    if (axis is not None):
      assert (axis_index is None)
      axis_index = ["h","k","l"].index(axis)
    if (slice_index is not None):
      assert (slice_start is None) and (slice_end is None)
      selection = simple_slice(
        indices=self.indices(),
        slice_axis=axis_index,
        slice_index=slice_index)
    else :
      assert (not None in [slice_start, slice_end])
      selection = multi_slice(
        indices=self.indices(),
        slice_axis=axis_index,
        slice_start=slice_start,
        slice_end=slice_end)
    return self.select(selection)

  def delete_index(self, hkl):
    """
    Remove all reflections with the specified Miller index.
    """
    assert (len(hkl) == 3)
    sele = (self.indices() != hkl)
    return self.select(sele)

  def delete_indices(self, other):
    """
    Delete multiple reflections, as specified by the Miller indices of
    another set.
    """
    # XXX inefficient - should probably use match_indices but it seems
    # to not like unmerged data
    sele = flex.bool(self.indices().size(), True)
    for hkl in other.indices():
      sele &= (self.indices() != hkl)
    return self.select(sele)

  def generate_bivoet_mates(self):
    """
    If the array is not already anomalous, expand the miller indices to generate
    anomalous pairs.
    """
    if (self.anomalous_flag()): return self
    sel = ~self.centric_flags().data()
    indices = self.indices().deep_copy()
    indices.extend(-indices.select(sel))
    return set(
        crystal_symmetry=self,
        indices=indices,
        anomalous_flag=True,
    )


def build_set(crystal_symmetry, anomalous_flag, d_min=None, d_max=None,
              max_index=None):
  """
  Given crystal symmetry, anomalous flag, and resolution limits, create a
  complete set object.

  :param crystal_symmetry: cctbx.crystal.symmetry object
  :param anomalous_flag: Boolean, indicates whether to generate anomalous indices
  :param d_min: High-resolution limit (optional if max_index is specified)
  :param d_max: Low-resolution limit (optional)
  :param max_index: highest-resolution Miller index
  :returns: a set object
  """
  if(max_index is not None):
    assert [d_min, d_max].count(None) == 2
    result = set(
      crystal_symmetry,
      index_generator(
        crystal_symmetry.space_group_info().type(),
        anomalous_flag,
        max_index).to_array(),
      anomalous_flag)
  else:
    result = set(
      crystal_symmetry,
      index_generator(
        crystal_symmetry.unit_cell(),
        crystal_symmetry.space_group_info().type(),
        anomalous_flag,
        d_min).to_array(),
      anomalous_flag)
    if(d_max is not None): result = result.resolution_filter(d_max = d_max)
  return result

def union_of_sets(miller_sets):
  assert len(miller_sets) != 0
  uoi = union_of_indices_registry()
  for ms in miller_sets:
    uoi.update(ms.indices())
  return miller_sets[0].customized_copy(indices=uoi.as_array())

class array_info(object):
  """
  Container for metadata associated with a Miller array, including labels
  read from a data file.
  """
  def __init__(self,
        source=None,
        source_type=None,
        history=None,
        labels=None,
        merged=False,
        systematic_absences_eliminated=False,
        crystal_symmetry_from_file=None,
        type_hints_from_file=None,
        wavelength=None):
    adopt_init_args(self, locals())

  def __setstate__(self, state):
    self.type_hints_from_file = None # backward compatibility
    self.__dict__.update(state)
    if (not hasattr(self, "wavelength")) : # backward compatibility
      self.wavelength = None

  def customized_copy(self,
        source=Keep,
        source_type=Keep,
        history=Keep,
        labels=Keep,
        merged=Keep,
        systematic_absences_eliminated=Keep,
        crystal_symmetry_from_file=Keep,
        type_hints_from_file=Keep,
        wavelength=Keep):
    """
    Create a modified copy of the array_info object, keeping all attributes
    that are not explicitly modified.
    """
    if (source is Keep): source = self.source
    if (source_type is Keep): source_type = self.source_type
    if (history is Keep): history = self.history
    if (labels is Keep): labels = self.labels
    if (merged is Keep): merged = self.merged
    if (systematic_absences_eliminated is Keep):
      systematic_absences_eliminated = self.systematic_absences_eliminated
    if (crystal_symmetry_from_file is Keep):
      crystal_symmetry_from_file = self.crystal_symmetry_from_file
    if (type_hints_from_file is Keep):
      type_hints_from_file = self.type_hints_from_file
    if (wavelength is Keep):
      wavelength = self.wavelength
    return array_info(
      source=source,
      source_type=source_type,
      history=history,
      labels=labels,
      merged=merged,
      systematic_absences_eliminated=systematic_absences_eliminated,
      crystal_symmetry_from_file=crystal_symmetry_from_file,
      type_hints_from_file=type_hints_from_file,
      wavelength=wavelength)

  def as_string_part_2(self):
    part_2 = []
    if (self.labels is not None):
      part_2.extend(self.labels)
    if (self.merged):
      part_2.append("merged")
    if (self.systematic_absences_eliminated):
      part_2.append("systematic_absences_eliminated")
    return part_2

  def label_string(self):
    """
    A combined representation of the data labels extracted from the input file.
    This is generally how downstream programs will identify and select Miller
    arrays.
    """
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
  except Exception:
    return "Unknown"

class array(set):
  """
  Extension of the set class with addition of data and (optional) sigmas
  flex arrays, plus an optional array_info object and an optional flag for
  the observation type (amplitude, intensity, or reconstructed amplitude).
  """
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
    # FIXME this should really be restricted to array_info objects or None,
    # but some higher-level apps (e.h. HySS) need to change for this to work.
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
    """
    Flag the array as X-ray amplitudes (F).
    """
    from cctbx.xray import observation_types
    return self.set_observation_type(observation_types.amplitude())

  def set_observation_type_xray_intensity(self):
    """
    Flag the array as X-ray intensities (I).
    """
    from cctbx.xray import observation_types
    return self.set_observation_type(observation_types.intensity())

  def data(self):
    return self._data

  def sigmas(self):
    return self._sigmas

  def set_sigmas(self, sigmas):
    if sigmas is not None:
      assert sigmas.size() == self.indices().size()
    self._sigmas = sigmas

  def __iter__(self):
    if self.sigmas() is not None:
      for item in zip(self.indices(), self.data(), self.sigmas()):
        yield item
    else:
      for item in zip(self.indices(), self.data()):
        yield item

  def info(self):
    """
    Return the associated info object, or None if undefined.
    """
    return self._info

  def show_comprehensive_summary(self, f=None, prefix=""):
    set.show_comprehensive_summary(self, f=f, prefix=prefix)
    if (self.info() is not None) and isinstance(self.info(), array_info):
      wavelength = self.info().wavelength
      if (wavelength is not None):
        print(prefix + "Wavelength: %.4f" % wavelength, file=f)
    return self

  def observation_type(self):
    """
    Return the (experimental) data type, if defined.  See the module
    cctbx.xray.observation_types for details.

    :returns: an object from cctbx.xray.observation_types
    """
    return self._observation_type

  def size(self):
    assert self.indices() is not None
    assert self.data() is not None
    assert self.data().size() == self.indices().size()
    if (self.sigmas() is not None):
      assert self.sigmas().size() == self.indices().size()
    return set.size(self)

  def is_string_array(self):
    return isinstance(self.data(), flex.std_string)

  def is_bool_array(self):
    return isinstance(self.data(), flex.bool)

  def is_integer_array(self):
    return isinstance(self.data(), flex.int) \
        or isinstance(self.data(), flex.long) \
        or isinstance(self.data(), flex.size_t)

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

  def is_xray_data_array(self):
    return (self.is_xray_amplitude_array() or self.is_xray_intensity_array())

  def copy(self):
    """
    Create a new array object using references to internal objects.
    """
    return (array(
      miller_set=self,
      data=self.data(),
      sigmas=self.sigmas())
      .set_info(self.info())
      .set_observation_type(self))

  def deep_copy(self):
    """
    Clone the array, making copies of all internal array objects.
    """
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
        space_group_info=Keep,
        observation_type=Keep,
        info=None):
    if (miller_set is Keep): miller_set = self
    if (data is Keep): data = self.data()
    if (sigmas is Keep): sigmas = self.sigmas()
    if (info is Keep) : info = self.info()
    if observation_type is Keep: observation_type = self.observation_type()
    miller_set = set.customized_copy(miller_set,
      crystal_symmetry=crystal_symmetry,
      indices=indices,
      anomalous_flag=anomalous_flag,
      unit_cell=unit_cell,
      space_group_info=space_group_info)
    return array(miller_set=miller_set, data=data, sigmas=sigmas)\
           .set_observation_type(observation_type).set_info(info)

  def concatenate(self, other, assert_is_similar_symmetry=True):
    if([self.sigmas(), other.sigmas()].count(None) == 0):
      return self.set().concatenate(
        other=other.set(),
        assert_is_similar_symmetry=assert_is_similar_symmetry).array(
          data=self.data().concatenate(other.data()),
          sigmas=self.sigmas().concatenate(other.sigmas()))
    else:
      return self.set().concatenate(
        other=other.set(),
        assert_is_similar_symmetry=assert_is_similar_symmetry).array(
          data=self.data().concatenate(other.data()))

  def set(self,
        crystal_symmetry=Keep,
        indices=Keep,
        anomalous_flag=Keep,
        unit_cell=Keep,
        space_group_info=Keep):
    """
    Return the basic cctbx.miller.set object for the array.
    """
    return set.customized_copy(self,
      crystal_symmetry=crystal_symmetry,
      indices=indices,
      anomalous_flag=anomalous_flag,
      unit_cell=unit_cell,
      space_group_info=space_group_info)

  def discard_sigmas(self):
    """
    Create a copy of the array without sigmas.
    """
    return array.customized_copy(self, sigmas=None)

  def conjugate(self):
    assert self.is_complex_array()
    return array.customized_copy(self, data=flex.conj(self.data()))

  def regularize(self):
    """
    A series of conversions that are required for many downstream tests, such as
    refinement, map calculation, etc.
    """
    result = self.deep_copy()
    info = result.info()
    result = result.eliminate_sys_absent()
    if info is not None:
      info = info.customized_copy(systematic_absences_eliminated = True)
    if(not result.is_unique_set_under_symmetry()):
      merged = result.merge_equivalents()
      result = merged.array()
      if info is not None:
        info = info.customized_copy(merged=True)
    result = result.map_to_asu()
    if(not result.sigmas_are_sensible()):
      result = result.customized_copy(
        indices=result.indices(),
        data=result.data(),
        sigmas=None).set_observation_type(result)
    sel = result.indices()==(0,0,0)
    if(not sel.all_eq(False)):
      result = result.select(~sel)
    sigmas = result.sigmas()
    if(sigmas is not None):
      selection  = result.sigmas() != 0
      selection &= result.data() != 0
      result = result.select(selection)
    if(result.is_xray_amplitude_array()):
      selection_positive = result.data() >= 0
      result = result.select(selection_positive)
    return result.set_info(info)

  def g_function(self, R, s=None, volume_scale=False):
    if s is None:
      s = 1./self.d_spacings().data() # Note f000 term will get s = -1
    return scitbx.math.g_function(s, R, volume_scale)

  def as_double(self):
    """
    Create a copy of the array with the data converted to a flex.double type.
    Will fail for incompatible arrays.
    """
    return self.array(data=self.data().as_double())

  def __getitem__(self, slice_object):
    return array(
      miller_set=set.__getitem__(self, slice_object),
      data=_slice_or_none(self.data(), slice_object),
      sigmas=_slice_or_none(self.sigmas(), slice_object))

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print(prefix + "Miller %s info:" % (
      self.__class__.__name__), self.info(), file=f)
    print(prefix + "Observation type:", self.observation_type(), file=f)
    print(prefix + "Type of data:", raw_array_summary(self.data()), file=f)
    print(prefix + "Type of sigmas:", raw_array_summary(self.sigmas()), file=f)
    set.show_summary(self, f=f, prefix=prefix)
    return self

  def make_up_hl_coeffs(self, k_blur, b_blur):
    assert isinstance(self.data(), flex.complex_double)
    phases = self.phases().data()
    sin_phases = flex.sin(phases)
    cos_phases = flex.cos(phases)
    ss = 1./flex.pow2(self.d_spacings().data()) / 4.
    t = 2*k_blur * flex.exp(-b_blur*ss)
    return self.customized_copy(
      data = flex.hendrickson_lattman(a = t * cos_phases, b = t * sin_phases))

  def disagreeable_reflections(self, f_calc_sq, n_reflections=20):
    assert f_calc_sq.is_xray_intensity_array()
    assert self.is_xray_intensity_array()
    assert self.sigmas() is not None
    assert self.size() == f_calc_sq.size()
    n_reflections = min(n_reflections, self.size())
    fo2 = self
    fc2 = f_calc_sq
    fc = f_calc_sq.as_amplitude_array()
    delta_f_sq = flex.abs(fo2.data() - fc2.data())/fo2.sigmas()
    fc_over_fc_max = fc.data()/flex.max(fc.data())
    perm = flex.sort_permutation(delta_f_sq, reverse=True)[:n_reflections]
    fo2 = fo2.select(perm)
    fc2 = fc2.select(perm)
    delta_f_sq_over_sigma = delta_f_sq.select(perm)
    fc_over_fc_max = fc_over_fc_max.select(perm)
    indices = fo2.indices()
    d_spacings = fo2.d_spacings().data()
    return group_args(indices=indices,
                      fo_sq=fo2,
                      fc_sq=fc2,
                      delta_f_sq_over_sigma=delta_f_sq_over_sigma,
                      fc_over_fc_max=fc_over_fc_max,
                      d_spacings=d_spacings)

  def show_disagreeable_reflections(self, f_calc_sq, n_reflections=20, out=None):
    if out is None: out = sys.stdout
    result = self.disagreeable_reflections(f_calc_sq, n_reflections=n_reflections)
    print("  h   k   l       Fo^2      Fc^2   |Fo^2-Fc^2|/sig(F^2)   Fc/max(Fc)  d spacing(A)", file=out)
    for i in range(result.fo_sq.size()):
      print("%3i %3i %3i" %result.indices[i], end=' ', file=out)
      print(" %9.2f %9.2f        %9.2f         %9.2f     %9.2f" %(
        result.fo_sq.data()[i], result.fc_sq.data()[i],
        result.delta_f_sq_over_sigma[i],
        result.fc_over_fc_max[i], result.d_spacings[i]), file=out)
    return result

  def crystal_symmetry_is_compatible_with_symmetry_from_file(self,
        unit_cell_relative_length_tolerance=0.02,
        unit_cell_absolute_angle_tolerance=3.,
        working_point_group=None):
    return crystal_symmetry_is_compatible_with_symmetry_from_file(
      miller_array=self,
      unit_cell_relative_length_tolerance=unit_cell_relative_length_tolerance,
      unit_cell_absolute_angle_tolerance=unit_cell_absolute_angle_tolerance,
      working_point_group=working_point_group)

  def sigmas_are_sensible(self, critical_ratio=0.75, epsilon=1e-6):
    result=None
    if self.sigmas() is not None and self.sigmas().size() != 0:
      result=True
      suspected = ( self.sigmas() <= epsilon ).count(True)
      all = self.sigmas().size()
      ratio = float(suspected)/float(all)
      if ratio>critical_ratio:
        result = False
    return result

  def enforce_positive_sigmas(self):
    if (self.sigmas() is None):
      return self
    else :
      return self.select(self.sigmas() > 0)

  def enforce_positive_amplitudes(self,i_sig_level=-4.0):
    """
    Takes in an intensity array (including negatives) and spits out amplitudes.
    The basic assumption is that
    P(Itrue) \\propto exp(-(Itrue-Iobs)**2/(2*s))
    where Itrue>=0 (positivity constraint on error free amplitudes).
    For amplitudes, this results in
    P(Ftrue) \\propto 2 Ftrue exp( -(Ftrue**2-Iobs)**2/(2s) )
    A Gaussian approximation is fitted to the Mode of this distribution.
    An analytical solution exists and is implemented below.
    This method does not require any Wilson statistics assumptions.
    """

    assert self.is_xray_intensity_array()
    assert self.sigmas() is not None
    assert self.sigmas_are_sensible()

    self = self.select( self.sigmas() > 0 )
    i_sigi = self.data()/self.sigmas()
    self = self.select( i_sigi > i_sig_level )
    det = flex.sqrt( self.data()*self.data() + 2.0*self.sigmas()*self.sigmas())
    f_saddle = flex.sqrt( (self.data()+det) / 2.0)
    s_saddle = (1.0/(f_saddle*f_saddle)) + (self.data() + 3.0*det) \
             / (self.sigmas()*self.sigmas() )
    s_saddle = flex.sqrt( 1.0/s_saddle )

    result = self.customized_copy(
      data=f_saddle,
      sigmas=s_saddle).set_observation_type( self.as_amplitude_array() )
    return result

  def f_sq_as_f(self, algorithm="xtal_3_7", tolerance=1.e-6):
    """
    Given an intensity/F^2 array (or undefined observation type), return the
    equivalent amplitudes.  Note that negative intensities will be discarded;
    for French-Wilson treatment, call the separate array.french_wilson()
    method.
    """
    from cctbx import xray
    assert self.observation_type() is None or  self.is_xray_intensity_array()
    converters = {
      "xtal_3_7": xray.array_f_sq_as_f_xtal_3_7,
      "crystals": xray.array_f_sq_as_f_crystals}
    converter = converters.get(algorithm)
    if (converter is None):
      raise RuntimeError(
        "Unknown f_sq_as_f algorithm=%s\n" % show_string(algorithm)
        + "  Possible choices are: "
        + ", ".join([show_string(s) for s in converters.keys()]))
    if (self.sigmas() is None):
      result = array(self, converter(self.data()).f)
    else:
      r = converter(self.data(), self.sigmas(), tolerance)
      result = array(self, r.f, r.sigma_f)
    if (self.is_xray_intensity_array()):
      result.set_observation_type_xray_amplitude()
    return result

  def f_as_f_sq(self, algorithm="simple"):
    """
    Convert amplitudes (and associated sigmas, if present) to intensities.
    """
    assert algorithm in ["simple", "shelxl"]
    from cctbx import xray
    assert self.observation_type() is None or self.is_xray_amplitude_array()
    if (self.sigmas() is None):
      result = array(self, xray.array_f_as_f_sq(self.data()).f_sq)
    elif (algorithm == "simple"):
      r = xray.array_f_as_f_sq(self.data(), self.sigmas())
      result = array(self, r.f_sq, r.sigma_f_sq)
    else:
      # shelx97-2.980324 shelxl.f line 6247:
      #   S=2.*AMAX1(.01,S)*ABS(AMAX1(.01,ABS(T),S))
      #   T=T**2
      f_sq = flex.double()
      s_sq = flex.double()
      for f,s in zip(self.data(), self.sigmas()):
        f_sq.append(f**2)
        s_sq.append(2 * max(.01, s) * abs(max(.01, abs(f), s)))
      result = array(miller_set=self, data=f_sq, sigmas=s_sq)
    if (self.is_xray_amplitude_array()):
      result.set_observation_type_xray_intensity()
    return result

  def as_amplitude_array(self, algorithm="xtal_3_7"):
    """
    Convert the array to simple amplitudes if not already in that format.
    Only valid for complex (i.e. F,PHI), intensity, or amplitude arrays.
    """
    if (self.is_complex_array()):
      return array(
        miller_set=self, data=flex.abs(self.data()), sigmas=self.sigmas()) \
          .set_observation_type_xray_amplitude()
    assert self.is_real_array()
    if (self.is_xray_intensity_array()):
      return self.f_sq_as_f(algorithm=algorithm)
    return self

  def as_intensity_array(self, algorithm="simple"):
    """
    Convert the array to intensities if not already in that format.  Only valid
    for complex (F,PHI), amplitude, or intensity arrays.
    """
    if (self.is_complex_array()):
      return self.as_amplitude_array().f_as_f_sq(algorithm=algorithm)
    assert self.is_real_array()
    if (not self.is_xray_intensity_array()):
      return self.f_as_f_sq(algorithm=algorithm)
    return self

  def translational_shift(self, shift_frac, deg=None):
    """
    Adjust a complex array (map coefficients) or phase array
     corresponding to a shift of all coordinates by
     new_xyz_frac = old_xyz_frac + shift_frac.

    If a phase array, must specify whether it is in degrees.
    Only makes sense in P1

    F'=exp( i 2 pi h.(x+shift_frac)) -> F' = F exp ( i 2 pi h.shift_frac)
    phase_shift = 2 * pi * h . shift_frac
    """
    assert self.is_complex_array() or self.is_real_array()
    assert self.space_group().type().number() in [0,1]

    from scitbx.matrix import col

    if self.is_complex_array(): # split into phases and amplitudes, shift phases

      amplitudes = self.amplitudes()
      phases_rad = self.phases(deg=False)
      new_phases_rad = phases_rad.translational_shift(shift_frac,deg=False)

      result = amplitudes.phase_transfer(
        phase_source = new_phases_rad,deg=False)

    else: # phase array

      assert deg is not None

      if deg: # in degrees
        c=360.
      else:  # in radians
        c=2*3.14159

      phase_shift = c * self.indices().as_vec3_double().dot(
         col(shift_frac))
      result = array(self, self.data()+phase_shift)

    return result

  def map_to_asu(self, deg=None):
    """
    Convert all indices to lie within the canonical asymmetric unit for the
    current space group (while preserving anomalous flag).  Required for many
    downstream steps.
    """
    i = self.indices().deep_copy()
    d = self.data().deep_copy()
    if (self.is_complex_array() or self.is_hendrickson_lattman_array()):
      map_to_asu(
        self.space_group_info().type(),
        self.anomalous_flag(),
        i, d)
    elif (deg is None):
      map_to_asu(
        self.space_group_info().type(),
        self.anomalous_flag(),
        i)
    else:
      map_to_asu(
        self.space_group_info().type(),
        self.anomalous_flag(),
        i, d, deg)
    return (array(set(self, i, self.anomalous_flag()), d, self.sigmas())
      .set_observation_type(self))

  def adopt_set(self, other, assert_is_similar_symmetry=True):
    if (assert_is_similar_symmetry):
      assert self.is_similar_symmetry(other)
    assert self.indices().size() == other.indices().size()
    assert self.anomalous_flag() == other.anomalous_flag()
    p = match_indices(other.indices(), self.indices()).permutation()
    assert self.indices().select(p).all_eq(other.indices())
    d = self.data()
    s = self.sigmas()
    if (d is not None): d = d.select(p)
    if (s is not None): s = s.select(p)
    return (array(miller_set=other, data=d, sigmas=s)
      .set_observation_type(self))

  def matching_set(self,
        other,
        data_substitute,
        sigmas_substitute=None,
        assert_is_similar_symmetry=True):
    assert self.data().size() == self.indices().size()
    if (self.sigmas() is not None):
      assert self.sigmas().size() == self.indices().size()
      assert sigmas_substitute is not None
    pairs = other.match_indices(
      other=self,
      assert_is_similar_symmetry=assert_is_similar_symmetry).pairs()
    if isinstance(data_substitute, self.data().__class__):
      assert data_substitute.size() == other.size()
      data = data_substitute.deep_copy()
    else:
      data = self.data().__class__(
        other.indices().size(), data_substitute)
    data.set_selected(
      pairs.column(0), self.data().select(pairs.column(1)))
    if (self.sigmas() is None):
      sigmas = None
    else:
      if isinstance(sigmas_substitute, self.sigmas().__class__):
        assert sigmas_substitute.size() == other.size()
        sigmas = sigmas_substitute.deep_copy()
      else:
        sigmas = self.sigmas().__class__(
          other.indices().size(), sigmas_substitute)
      sigmas.set_selected(
        pairs.column(0), self.sigmas().select(pairs.column(1)))
    return other.array(data=data, sigmas=sigmas)

  def scale(self, other, resolution_dependent=False):
    s, o = self.common_sets(other)
    if(resolution_dependent): ss = 1./flex.pow2(s.d_spacings().data()) / 4.
    s, o = abs(s).data(), abs(o).data()
    scale_factor = 1
    if(resolution_dependent):
      r = scitbx.math.gaussian_fit_1d_analytical(x=flex.sqrt(ss), y=s, z=o)
      ss_other = 1./flex.pow2(other.d_spacings().data()) / 4.
      scale_factor = r.a*flex.exp(-ss_other*r.b)
    else:
      den = flex.sum(o*o)
      if(den != 0):
        scale_factor = flex.sum(s*o)/den
    return other.customized_copy(data = other.data()*scale_factor)

  def complete_with(self, other, scale=False, replace_phases=False):
    s, o = self, other
    if(scale): o = s.scale(other = o)
    if(replace_phases):
      assert isinstance(o.data(), flex.complex_double)
      s_c, c_o = s.common_sets(o)
      s_c = s_c.phase_transfer(phase_source = c_o)
      s_l, c_l = s.lone_sets(o)
      s = s_c.concatenate(s_l)
    ol = o.lone_set(s)
    d_new = s.data().concatenate(ol.data())
    i_new = s.indices().concatenate(ol.indices())
    sigmas_new = None
    if(s.sigmas() is not None):
      sigmas_new = s.sigmas().concatenate(ol.sigmas())
    return self.customized_copy(data = d_new, indices = i_new, sigmas = sigmas_new)

  def combine(self, other, scale = True, scale_for_lones = 1,
        scale_for_matches=1):
    assert self.anomalous_flag() == other.anomalous_flag()
    assert self.sigmas() is None # not implemented
    f1_c, f2_c = self.common_sets(other = other)
    f1_l, f2_l = self.lone_sets(other = other)
    scale_k1 = scale_for_matches
    if(scale):
      den = flex.sum(flex.abs(f2_c.data())*flex.abs(f2_c.data()))
      if(den != 0):
        scale_k1 = flex.sum(flex.abs(f1_c.data())*flex.abs(f2_c.data())) / den
    result_data = f1_c.data() + f2_c.data()*scale_k1
    result_data.extend(f1_l.data()*scale_for_lones)
    result_data.extend(f2_l.data()*scale_k1*scale_for_lones)
    result_indices = f1_c.indices()
    result_indices.extend(f1_l.indices())
    result_indices.extend(f2_l.indices())
    ms = set(
      crystal_symmetry=self.crystal_symmetry(),
      indices=result_indices,
      anomalous_flag=self.anomalous_flag())
    return ms.array(data = result_data)

  def complete_array(self,
        d_min_tolerance=1.e-6,
        d_min=None,
        d_max=None,
        new_data_value=-1,
        new_sigmas_value=-1):
    cs = self.complete_set(
      d_min_tolerance=d_min_tolerance, d_min=d_min, d_max=d_max)
    matches = match_indices(self.indices(), cs.indices())
    # don't assert no singles here, for the cases when
    # d_min > self.d_min() or d_max < self.d_max()
    i = self.indices()
    d = self.data()
    if (d is not None): d = d.deep_copy()
    s = self.sigmas()
    if (s is not None): s = s.deep_copy()
    ms = matches.singles(1)
    n = ms.size()
    if (n == 0):
      i = i.deep_copy()
    else:
      i = i.concatenate(cs.indices().select(ms))
      if (d is not None): d.resize(d.size()+n, new_data_value)
      if (s is not None): s.resize(s.size()+n, new_sigmas_value)
    return self.customized_copy(indices=i, data=d, sigmas=s)

  def sort_permutation(self, by_value="resolution", reverse=False):
    """
    Generate the selection array (flex.size_t object) to reorder the array
    by resolution, Miller index, data values, or absolute data values.

    :param by_value: sort type, must be "resolution", "packed_indices", "data",
                     or "abs"
    :param reverse: invert order
    :returns: flex.size_t object
    """
    assert reverse in (False, True)
    if (by_value in ["resolution", "packed_indices","asu_indices"]):
      result = set.sort_permutation(self,
        by_value=by_value, reverse=reverse)
    elif (by_value == "data"):
      result = flex.sort_permutation(
        data=self.data(), reverse=not reverse, stable=True)
    elif (by_value == "abs"):
      result = flex.sort_permutation(
        data=flex.abs(self.data()), reverse=not reverse, stable=True)
    elif (isinstance(by_value, str)):
      raise ValueError("Unknown: by_value=%s" % by_value)
    else:
      result = flex.sort_permutation(
        data=by_value, reverse=not reverse, stable=True)
    return result

  def patterson_symmetry(self):
    data = self.data()
    if (self.is_complex_array()):
      data = flex.abs(self.data())
    return array(
      set.patterson_symmetry(self),
      data,
      self.sigmas())

  def expand_to_p1(self, phase_deg=None, return_iselection=False):
    """
    Generate the equivalent P1 dataset.
    """
    assert self.space_group_info() is not None
    assert self.indices() is not None
    assert self.anomalous_flag() is not None
    assert self.data() is not None
    p1_sigmas = None
    iselection = None
    expand_type = {
      "complex_double": expand_to_p1_complex,
      "hendrickson_lattman": expand_to_p1_hendrickson_lattman,
    }.get(self.data().__class__.__name__, None)
    if (expand_type is not None):
      assert phase_deg is None
      assert self.sigmas() is None
      p1 = expand_type(
        space_group=self.space_group(),
        anomalous_flag=self.anomalous_flag(),
        indices=self.indices(),
        data=self.data())
      p1_data = p1.data
    elif (phase_deg is not None):
      assert phase_deg in [False, True]
      assert self.sigmas() is None
      assert isinstance(self.data(), flex.double)
      p1 = expand_to_p1_phases(
        space_group=self.space_group(),
        anomalous_flag=self.anomalous_flag(),
        indices=self.indices(),
        data=self.data(),
        deg=phase_deg)
      p1_data = p1.data
    else:
      assert phase_deg is None
      p1 = self.expand_to_p1_iselection()
      p1_data = self.data().select(p1.iselection)
      if (self.sigmas() is not None):
        p1_sigmas = self.sigmas().select(p1.iselection)
      iselection = p1.iselection
    result = set(
      crystal_symmetry=self.cell_equivalent_p1(),
      indices=p1.indices,
      anomalous_flag=self.anomalous_flag()).array(
        data=p1_data,
        sigmas=p1_sigmas).set_observation_type(self)
    if (return_iselection):
      assert iselection is not None # not implemented
      return result, iselection
    return result

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
      sigmas=result_sigmas).set_observation_type(self.observation_type())

  def symmetry_agreement_factor(self, op, assert_is_similar_symmetry=True):
    """ The factor phi_{sym} quantifying whether complex structure factors
    are invariant under the given symmetry operator, as used in Superflip.
    Ref: J. Appl. Cryst. (2008). 41, 975-984
    """
    assert self.is_complex_array()
    f, op_f = self.common_sets(
      self.change_basis(op),
      assert_is_similar_symmetry=assert_is_similar_symmetry)
    assert f.size() == op_f.size()
    f, op_f = f.data(), op_f.data()
    cc_sf = f * flex.conj(op_f) # structure factors of cross-correlation
    weights = flex.abs(cc_sf)
    delta = flex.pow(flex.arg(cc_sf), 2) # flex.arg uses atan2,
                                         # thus already in [-pi, pi]
    # This is eq. (7) from the reference
    return 3/math.pi**2 * flex.mean_weighted(delta, weights)

  def f_obs_minus_f_calc(self, f_obs_factor, f_calc):
    assert f_calc.is_complex_array()
    assert self.indices().all_eq(f_calc.indices())
    assert self.anomalous_flag() is f_calc.anomalous_flag()
    if self.is_complex_array():
      return array(
        miller_set=self,
        data=f_obs_factor*self.data()-f_calc.data())
    else:
      return array(
        miller_set=self,
        data=f_obs_factor*self.data()-flex.abs(f_calc.data())).phase_transfer(
          phase_source=f_calc)

  def phase_transfer(self, phase_source, epsilon=1.e-10, deg=False,
                           phase_integrator_n_steps=None):
    """Combines phases of phase_source with self's data if real (keeping
    the sign of self's data) or with self's amplitudes if complex.

    Centric reflections are forced to be compatible with the phase restrictions.

    phase_source can be a miller.array or a plain flex array.

    epsilon is only used when phase_source is a complex array. If both the
    real and the imaginary part of phase_source[i] < epsilon the phase is
    assumed to be 0.

    deg is only used if phase_source is an array of doubles.
    deg=True indicates that the phases are given in degrees,
    deg=False indicates phases are given in radians.

    phase_integrator_n_steps is only used if phase_source is an
    array of Hendrickson-Lattman coefficients. The centroid
    phases are determined on the fly using the given step size.
    """
    assert self.data() is not None
    # XXX good to enable: assert self.indices().all_eq(phase_source.indices())
    if (hasattr(phase_source, "data")):
      phase_source = phase_source.data()
    assert (   isinstance(self.data(), flex.complex_double)
            or isinstance(self.data(), flex.double))
    assert (   isinstance(phase_source, flex.complex_double)
            or isinstance(phase_source, flex.double)
            or isinstance(phase_source, flex.hendrickson_lattman))
    if (isinstance(phase_source, flex.hendrickson_lattman)):
      if (phase_integrator_n_steps is None):
        integrator = phase_integrator()
      else:
        integrator = phase_integrator(n_steps=phase_integrator_n_steps)
      phase_source = integrator(
        space_group=self.space_group(),
        miller_indices=self.indices(),
        hendrickson_lattman_coefficients=phase_source)
    if isinstance(self.data(), flex.double):
      data = self.data()
    else:
      data = flex.abs(self.data())
    assert self.space_group() is not None
    if (isinstance(phase_source, flex.complex_double)):
      return array(
        miller_set=self,
        data=phase_transfer(
          self.space_group(),
          self.indices(),
          data,
          phase_source,
          epsilon),
        sigmas=self.sigmas())
    return array(
      miller_set=self,
      data=phase_transfer(
        self.space_group(),
        self.indices(),
        data,
        phase_source,
        deg),
      sigmas=self.sigmas())

  def randomize_phases(self):
    random_phases = (2*math.pi)*flex.random_double(self.data().size())
    return self.phase_transfer(random_phases)

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

  def mean_phase_error(self, phase_source):
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
    e = flex.mean(scitbx.math.phase_error(phi1=p1, phi2=p2))
    return e * 180/math.pi

  def randomize_amplitude_and_phase(self, amplitude_error,
        phase_error_deg, selection=None, random_seed=None):
    """
    Add random error to reflections.
    """
    assert self.is_complex_array()
    if(selection is None): selection = flex.bool(self.indices().size(), True)
    import random
    if(random_seed is None): random_seed = random.randint(0, 1000000)
    new_data = ext.randomize_amplitude_and_phase(
      data=self.data(),
      selection=selection,
      amplitude_error=amplitude_error,
      phase_error_deg=phase_error_deg,
      random_seed=random_seed)
    return self.customized_copy(data = new_data)

  def anomalous_differences(self, enforce_positive_sigmas=False):
    """
    Returns an array object with DANO (i.e. F(+) - F(-)) as data, and
    optionally SIGDANO as sigmas.
    """
    assert self.data() is not None
    tmp_array = self
    if (enforce_positive_sigmas):
      tmp_array = tmp_array.enforce_positive_sigmas()
    asu, matches = tmp_array.match_bijvoet_mates()
    i = matches.miller_indices_in_hemisphere("+")
    d = matches.minus(asu.data())
    s = None
    if (asu.sigmas() is not None):
      s = matches.additive_sigmas(asu.sigmas())
    return array(set(asu, i, anomalous_flag=False), d, s)

  def hemisphere_acentrics(self, plus_or_minus):
    assert plus_or_minus in ("+", "-")
    assert self.data() is not None
    asu, matches = self.match_bijvoet_mates()
    i_column = "+-".index(plus_or_minus)
    return asu.select(
      selection=matches.pairs().column(i_column),
      anomalous_flag=False)

  def hemispheres_acentrics(self):
    assert self.data() is not None
    asu, matches = self.match_bijvoet_mates()
    return tuple(
      [asu.select(
        selection=matches.pairs().column(i_column),
        anomalous_flag=False)
       for i_column in (0,1)])

  def anomalous_completeness(self, use_binning=False, d_min_tolerance=1.e-6,
      d_max=None, d_min=None, relative_to_complete_set=True):
    """
    Return the percent of acentric reflections with both h,k,l and -h,-k,-l
    observed (only meaningful for amplitude and intensity arrays).  By default
    this is calculated relative to the complete set.
    """
    assert self.anomalous_flag() and self.is_real_array()
    if (not use_binning):
      merged = self.average_bijvoet_mates()
      if (relative_to_complete_set):
        merged = merged.complete_set(d_max=d_max, d_min=d_min,
          d_min_tolerance=d_min_tolerance)
      centric_flags = merged.centric_flags().data()
      if (centric_flags.count(False) == 0):
        return 0
      merged_acentric = merged.select(~centric_flags)
      n_acentric = merged_acentric.size()
      anom_diffs = self.anomalous_differences()
      return min(anom_diffs.size() / n_acentric, 1.0)
    assert self.binner() is not None
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      array_sel = self.select(sel)
      d_max_bin, d_min_bin = self.binner().bin_d_range(i_bin)
      results.append(array_sel.anomalous_completeness(
        d_max=d_max_bin, d_min=d_min_bin,
        relative_to_complete_set=relative_to_complete_set))
    return binned_data(binner=self.binner(), data=results, data_fmt="%5.3f")

  def convert_to_non_anomalous_if_ratio_pairs_lone_less_than(self, threshold):
    """
    Convert anomalous array into nonanomalous if the number of Bijvoet pairs is
    too small compared to the number of lone Bijvoet mates.
    """
    if(not self.anomalous_flag()): return self
    no_sys_abs = self.copy()
    if (self.space_group_info() is not None):
      is_unique_set_under_symmetry = no_sys_abs.is_unique_set_under_symmetry()
      sys_absent_flags = self.sys_absent_flags().data()
      n_sys_abs = sys_absent_flags.count(True)
      if (n_sys_abs != 0):
        no_sys_abs = self.select(selection=~sys_absent_flags)
      n_centric = no_sys_abs.centric_flags().data().count(True)
    if (self.space_group_info() is not None
        and no_sys_abs.anomalous_flag()
        and is_unique_set_under_symmetry):
      asu, matches = no_sys_abs.match_bijvoet_mates()
      n_pairs = matches.pairs().size()
      n_lone_mates = matches.n_singles() - n_centric
      if(n_lone_mates != 0 and n_pairs*1./n_lone_mates < threshold):
        merged = self.as_non_anomalous_array().merge_equivalents()
        self = merged.array().set_observation_type(self)
    return self

  def anomalous_signal(self, use_binning=False):
    """Get the anomalous signal according to this formula:

    .. math::
       \\sqrt{\\dfrac{<||F(+)|-|F(-)||^2>}{\\frac{1}{2} (<|F(+)|>^2 + <|F(-)|>^2)}}

    :param use_binning: If 'True' the anomalous signal will be calculated for \
    each bin of the data set individually
    :type use_binning: boolean
    :returns: the anomalous signal
    :rtype: float or cctbx.miller.binned_data
    """
    assert not use_binning or self.binner() is not None
    if (not use_binning):
      obs = self.select(self.data() > 0)
      if (self.is_xray_intensity_array()):
        obs = obs.f_sq_as_f()
      f_plus, f_minus = obs.hemispheres_acentrics()
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

  def anomalous_probability_plot(self, expected_delta=None):
    assert self.is_unique_set_under_symmetry()
    assert self.anomalous_flag()

    dI = self.anomalous_differences()
    if not dI.size():
      return AnomalousProbabilityPlotResult(None, None, None, expected_delta)

    y = dI.data() / dI.sigmas()
    perm = flex.sort_permutation(y)
    y = y.select(perm)
    distribution = distributions.normal_distribution()

    x = distribution.quantiles(y.size())
    if expected_delta is not None:
      sel = flex.abs(x) < expected_delta
      x = x.select(sel)
      y = y.select(sel)

    fit = flex.linear_regression(x, y)
    if fit.is_well_defined():
      return AnomalousProbabilityPlotResult(fit.slope(), fit.y_intercept(), x.size(), expected_delta)
    return AnomalousProbabilityPlotResult(None, None, None, expected_delta)

  def phase_entropy(self, exponentiate=False, return_binned_data=False,
                          return_mean=False):
    """Get phase entropy as measured in terms of an base-360 entropy
    (base-2 for centrics).

    An entropy of 0, indicates that the phase uncertainity is as low as possible
    An entropy of 1 however, indicates that the uncertainty is maximal:
    all phases are equally likely!

    :param return_binned_data: if 'True' you receive a binned object rather \
    then a raw array
    :type return_binned_data: boolean
    :param exponentiate: whether or not to exponentiate the entropy. This will \
    return a phase uncertainty in degrees (or the 'alphabet size')
    :type exponentiate: boolean
    """
    assert ([return_binned_data,return_mean]).count(True)!=2

    if self.is_hendrickson_lattman_array():
      integrator = phase_entropy(n_steps=360)
      result = integrator.relative_entropy(self.space_group(), self.indices(), self.data() )
      if exponentiate:
        centric_flags = self.centric_flags().data()
        cen           = flex.exp( math.log(2.0)*result )
        cen           = cen.set_selected( ~centric_flags, 0 )
        acen          = flex.exp( math.log(2.0)*result )
        acen          = acen.set_selected(  centric_flags, 0 )
        result        = cen + acen
      if not return_binned_data:
        if return_mean:
          return flex.mean( result )
        else:
          result  = self.array( data=result )
          return result
      else:
        assert self.binner() is not None
        binned_results = []
        for i_bin in self.binner().range_all():
          sel = self.binner().selection(i_bin)
          sel_data = result.select(sel)
          mean = 0
          if sel_data.size() >0:
            mean = flex.mean( sel_data )
          binned_results.append( mean )
        return binned_data(binner=self.binner(), data=binned_results, data_fmt="%7.4f")
    else:
      return None


  def measurability(self, use_binning=False, cutoff=3.0, return_fail=None):
    ## Peter Zwart 2005-Mar-04
    """Fraction of reflections for which
    (:math:`\\dfrac{|\\Delta I|}{\\sigma_{dI}}` > cutoff and
    :math:`min(\\dfrac{I_{+}}{\\sigma_{+}},\\dfrac{I_{-}}{\\sigma_{-}})` > cutoff
    """
    assert self.anomalous_flag()
    assert not use_binning or self.binner() is not None
    assert self.sigmas() is not None
    cutoff = float(cutoff)
    if (not use_binning):
      obs = self.select(self.data() > 0 )
      if (self.is_xray_amplitude_array()):
        obs = obs.f_as_f_sq()
      if (obs.data().size() == 0): return return_fail
      i_plus, i_minus = obs.hemispheres_acentrics()
      assert i_plus.data().size() == i_minus.data().size()
      top = flex.fabs(i_plus.data()-i_minus.data())
      bottom = flex.sqrt( (i_plus.sigmas()*i_plus.sigmas()) + (i_minus.sigmas()*i_minus.sigmas()) )
      zeros = flex.bool( bottom <= 0 ).iselection()
      bottom = bottom.set_selected( zeros, 1 )
      ratio = top/bottom
      bottom = i_plus.sigmas().set_selected( flex.bool(i_plus.sigmas()<=0).iselection(), 1 )
      i_plus_sigma = i_plus.data()/bottom
      bottom = i_minus.sigmas().set_selected( flex.bool(i_minus.sigmas()<=0).iselection(), 1 )
      i_minus_sigma = i_minus.data()/bottom
      meas = (  (ratio > cutoff)
              & (i_plus_sigma > cutoff)
              & (i_minus_sigma > cutoff) ).count(True)
      if ratio.size()>0:
        return meas/ratio.size()
      else:
        return return_fail
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      results.append(self.select(sel).measurability(cutoff=cutoff,
                                                    return_fail=return_fail))
    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f")

  def bijvoet_ratios(self, obs_type="intensity", measurable_only=True,
      cutoff=3.0):
    assert self.anomalous_flag()
    obs = self.select(self.data() > 0 )
    if (obs_type == "amplitude") and (self.is_xray_intensity_array()):
      obs = obs.f_sq_as_f()
    elif (obs_type == "intensity") and (self.is_xray_amplitude_array()):
      obs = obs.f_as_f_sq()
    i_plus, i_minus = obs.hemispheres_acentrics()
    assert i_plus.data().size() == i_minus.data().size()
    i_mean = (i_plus.data() + i_minus.data()) / 2
    d_ano = flex.fabs(i_plus.data() - i_minus.data())
    if measurable_only :
      assert self.sigmas() is not None
      top = flex.fabs(i_plus.data()-i_minus.data())
      bottom = flex.sqrt( (i_plus.sigmas()*i_plus.sigmas()) +
                          (i_minus.sigmas()*i_minus.sigmas()) )
      zeros = flex.bool( bottom <= 0 ).iselection()
      bottom = bottom.set_selected( zeros, 1 )
      ratio = top/bottom
      bottom = i_plus.sigmas().set_selected(
        flex.bool(i_plus.sigmas()<=0).iselection(), 1 )
      i_plus_sigma = i_plus.data()/bottom
      bottom = i_minus.sigmas().set_selected(
        flex.bool(i_minus.sigmas()<=0).iselection(), 1 )
      i_minus_sigma = i_minus.data()/bottom
      meas = (  (ratio > cutoff)
              & (i_plus_sigma > cutoff)
              & (i_minus_sigma > cutoff) )
      i_mean = i_mean.select(meas)
      d_ano = d_ano.select(meas)
    assert i_mean.size() == d_ano.size()
    non_zero_sele = i_mean > 0
    d_ano = d_ano.select(non_zero_sele)
    i_mean = i_mean.select(non_zero_sele)
    return d_ano / i_mean

  def second_moment(self, use_binning=False):
    """<data^2>/(<data>)^2"""
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

  def wilson_plot(self, use_binning=False):
    """<F^2>"""
    assert not use_binning or self.binner() is not None
    obs = None
    if self.is_xray_amplitude_array():
      obs = self.f_as_f_sq()
      obs.use_binning_of(self)
    else :
      obs = self
    if (not use_binning):
      if (obs.indices().size() == 0): return None
      mean_data = flex.mean(obs.data())
      if (mean_data == 0): return None
      return mean_data
    result = []
    for i_bin in obs.binner().range_all():
      sel = obs.binner().selection(i_bin)
      result.append(obs.select(sel).wilson_plot())
    return binned_data(binner=obs.binner(), data=result, data_fmt="%7.4f")

  def i_over_sig_i(self, use_binning=False,return_fail=None):
    """<I/sigma_I>"""
    assert not use_binning or self.binner is not None

    if (not use_binning):
      if (self.indices().size() == 0): return return_fail
      obs = None
      if self.is_real_array():
        if self.is_xray_amplitude_array():
          obs = self.f_as_f_sq()
        if self.is_xray_intensity_array():
          obs = self
      if obs is not None:
        obs = obs.select(obs.sigmas()>0)
        if (obs.indices().size() == 0): return return_fail
        i_sig_i = flex.mean( obs.data()/obs.sigmas() )
        return i_sig_i
      else:
        return 0
    result = []
    for i_bin in self.binner().range_all():
      sel =  self.binner().selection(i_bin)
      result.append(self.select(sel).i_over_sig_i(return_fail=return_fail) )

    return binned_data(binner=self.binner(),
                       data=result,
                       data_fmt="%7.4f")

  def mean_of_intensity_divided_by_epsilon(self,
                                           use_binning=False,
                                           return_fail=None):
    """ <I/epsilon> """
    assert not use_binning or self.binner is not None
    if (not use_binning):
      if (self.indices().size() == 0): return return_fail
      obs = None
      if self.is_real_array():
        if self.is_xray_amplitude_array():
          obs = self.f_as_f_sq()
        if self.is_xray_intensity_array():
          obs = self
      if obs is not None:
        weighted_mean = flex.mean( obs.data() /
                                   obs.epsilons().data().as_double() )
        return weighted_mean
      else:
        return return_fail
    result = []
    for i_bin in self.binner().range_all():
      sel =  self.binner().selection(i_bin)
      result.append(self.select(sel).mean_of_intensity_divided_by_epsilon(
        return_fail=return_fail) )
    return binned_data(binner=self.binner(),
                       data=result,
                       data_fmt="%7.4f")

  def mean_of_squared_sigma_divided_by_epsilon(self,
                                              use_binning=False,
                                              return_fail=None):
    """ <sigma^2/epsilon> """
    assert not use_binning or self.binner is not None
    if (not use_binning):
      if (self.sigmas().size() == 0): return return_fail
      weighted_mean = flex.mean( self.sigmas() /
                                 self.epsilons().data().as_double() )
      return weighted_mean

    result = []
    for i_bin in self.binner().range_all():
      sel =  self.binner().selection(i_bin)
      result.append(self.select(sel).mean_of_squared_sigma_divided_by_epsilon(
        return_fail=return_fail) )
    return binned_data(binner=self.binner(),
                       data=result,
                       data_fmt="%7.4f")


  def second_moment_of_intensities(self, use_binning=False):
    """<I^2>/(<I>)^2 (2.0 for untwinned, 1.5 for twinned data)"""
    if (self.is_xray_intensity_array()):
      a = self
    else:
      a = self.f_as_f_sq()
      if (use_binning):
        a.use_binner_of(self)
    return a.second_moment(use_binning=use_binning)

  def wilson_ratio(self, use_binning=False):
    """(<F>)^2/<F^2> (0.785 for untwinned, 0.885 for twinned data)"""
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

  def show_r_free_flags_info(self,
        n_bins=10,
        binner_range="used",
        out=None,
        prefix=""):
    assert self.is_bool_array()
    assert binner_range in ["used", "all"]
    print(prefix + "Number of work/free reflections by resolution:", file=out)
    if (n_bins is not None):
      self.setup_binner(n_bins=n_bins)
    else:
      assert self.binner() is not None
    n_works = []
    n_frees = []
    fmt = None
    for i_bin in getattr(self.binner(), "range_"+binner_range)():
      sel = self.binner().selection(i_bin)
      flags = self.data().select(sel)
      n_free = flags.count(True)
      n_work = flags.size() - n_free
      n_works.append(n_work)
      n_frees.append(n_free)
      legend = self.binner().bin_legend(i_bin)
      if (fmt is None):
        width = max(4, len(str(self.indices().size())))
        fmt = "%%%ds" % width
        print(prefix + " ", " "*len(legend), \
          fmt%"work", fmt%"free", " %free", file=out)
        fmt = "%%%dd" % width
      print(prefix + " ", legend, fmt%n_work, fmt%n_free, \
        "%5.1f%%" % (100.*n_free/max(1,n_work+n_free)), file=out)
    n_free = self.data().count(True)
    n_work = self.data().size() - n_free
    print(prefix + " ", (
      "%%%ds"%len(legend))%"overall", fmt%n_work, fmt%n_free, \
      "%5.1f%%" % (100.*n_free/max(1,n_work+n_free)), file=out)
    return n_works, n_frees

  def r_free_flags_accumulation(self):
    assert self.is_bool_array()
    rc = flex.size_t(1, 0)
    ff = flex.double(1, 0)
    d = self.data()
    n = d.size()
    c = 0
    for i,f in enumerate(d):
      if (f):
        c += 1
        rc.append(i+1)
        ff.append(c/n)
    return group_args(reflection_counts=rc, free_fractions=ff)

  def r1_factor(self, other, scale_factor=None, assume_index_matching=False,
                use_binning=False, emulate_sftools=False):
    """Get the R1 factor according to this formula

    .. math::
       R1 = \\dfrac{\\sum{||F| - k|F'||}}{\\sum{|F|}}

    where F is self.data() and F' is other.data() and
    k is the factor to put F' on the same scale as F

    :param other: array object with the same observation type
    :param scale_factor: optional scale factor to be applied to 'other'; if
      Auto, will be determined automatically
    :param assume_index_matching: skips calling self.common_sets(other)
    :param use_binning: divide by resolution shells
    :param emulate_sftools: copies behavior of SFTOOLS in CCP4: instead of
      the denominator being sum(self.data()), it will be 0.5*sum(self.data()+
      other.data())
    :returns: a Python float (if use_binning=False), or a binned_data object
    """
    assert not use_binning or self.binner() is not None
    assert (self.observation_type() is None
            or self.is_complex_array() or self.is_xray_amplitude_array())
    assert (other.observation_type() is None
            or other.is_complex_array() or other.is_xray_amplitude_array())
    assert other.indices().size() == self.indices().size()
    if not use_binning:
      if self.data().size() == 0: return None
      if (assume_index_matching):
        o, c = self, other
      else:
        o, c = self.common_sets(other=other, assert_no_singles=True)
      o  = flex.abs(o.data())
      c = flex.abs(c.data())
      if (scale_factor is Auto):
        den = flex.sum(c * c)
        if (den != 0):
          c *= (flex.sum(o * c) / den)
      elif (scale_factor is not None):
        c *= scale_factor
      if (emulate_sftools):
        return 2 * flex.sum(flex.abs(o - c)) / flex.sum(o+c)
      else :
        return flex.sum(flex.abs(o - c)) / flex.sum(o)
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      results.append(self.select(sel).r1_factor(
        other.select(sel), scale_factor, assume_index_matching))
    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f")

  def select(self, selection, negate=False, anomalous_flag=None):
    """
    Select a sub-array.

    :param selection: flex.bool or flex.size_t selection
    :param negate: select the inverse of the selection array
    :param anomalous_flag: anomalous flag for the new set
    :returns: a new array with a subset of indices and data/sigmas
    """
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

  def select_indices(self, indices, map_indices_to_asu=False, negate=False):
    if map_indices_to_asu:
      indices = indices.deep_copy() # map_to_asu changes indices in place
      map_to_asu(self.space_group().type(), True, indices)
      matched_indices = match_indices(self.map_to_asu().indices(), indices)
    else:
      matched_indices = match_indices(self.indices(), indices)
    try:
      return self.select(matched_indices.pair_selection(0), negate=negate)
    except RuntimeError as e:
      if ('CCTBX_ASSERT(miller_indices_[1].size() == size_processed(1))'
          in str(e)):
        raise RuntimeError(
          "cctbx.miller.array.select_indices(): "
          "This method can only be used reliably on a merged array")
      else:
        raise RuntimeError(e)

  def sigma_filter(self, cutoff_factor=None, negate=False):
    """
    Return a copy of the array filtered to remove reflections whose value is
    less than cutoff_factor*sigma (or the reverse, if negate=True).
    """
    if(cutoff_factor is None): return self
    if(self.data() is None): return self
    if(self.sigmas() is None): return self
    flags = flex.abs(self.data()) >= self.sigmas() * cutoff_factor
    return self.select(flags, negate)

  def min_f_over_sigma(self, return_none_if_zero_sigmas=False):
    result = None
    sigmas = self.sigmas()
    if(sigmas is not None and sigmas.size() != 0):
      if(flex.min(sigmas) == 0.0):
        result = 0.0
      else:
        sigmas_all_not_equal_zero = sigmas.all_ne(0)
        if (not return_none_if_zero_sigmas):
          assert sigmas_all_not_equal_zero
        if (sigmas_all_not_equal_zero):
          result = flex.min(self.data() / sigmas)
    return result

  def apply_scaling(self, target_max=None, factor=None):
    """
    Apply a scale factor to the data (and optionally sigmas).

    :param target_max: target maximum value for the scaled data - the scale
                       factor will be determined automatically
    :param factor: explicit scale factor
    :returns: custumozed copy with scaled data and sigmas
    """
    assert [target_max, factor].count(None) == 1
    assert self.data() is not None
    s = None
    if (target_max is not None):
      current_max = flex.max(flex.abs(self.data()))
      if (current_max == 0): return self.deep_copy()
      factor = target_max / current_max
    d = self.data() * factor
    if (self.sigmas() is not None): s = self.sigmas() * factor
    return self.customized_copy(data=d, sigmas=s) \
      .set_info(self.info()) \
      .set_observation_type(self)

  def multiscale(self, other, use_exp_scale=False, reflections_per_bin=None):
    if use_exp_scale:
      import libtbx.load_env
      if (not libtbx.env.has_module("mmtbx")):
        raise ImportError("mmtbx is required for this functionality.")
      from mmtbx import bulk_solvent
      ss = 1./flex.pow2(self.d_spacings().data()) / 4.
      b_range = flex.double(range(-1000,1000,1))
      if isinstance(self.data(), flex.complex_double):
        scaler_func = bulk_solvent.complex_f_kb_scaled
      else:
        scaler_func = bulk_solvent.f_kb_scaled
    if(reflections_per_bin is None):
      reflections_per_bin = other.indices().size()
    assert self.indices().all_eq(other.indices())
    assert self.is_similar_symmetry(other)
    assert type(self.data()) == type(other.data())

    self.setup_binner(reflections_per_bin = reflections_per_bin)
    other.use_binning_of(self)
    other_data = other.data()*0
    for i_bin in self.binner().range_used():
      sel = self.binner().selection(i_bin)
      f1  = self.select(sel)
      f2  = other.select(sel)
      if not use_exp_scale: # use scalar scale
        scale_ = 1.0
        den = flex.sum(flex.abs(f2.data())*flex.abs(f2.data()))
        if(den != 0):
          scale_ = flex.sum(flex.abs(f1.data())*flex.abs(f2.data())) / den
        other_data.set_selected(sel, f2.data()*scale_)
      else:
        ss_ = ss.select(sel)
        o = scaler_func(
          f1      = f1.data(),
          f2      = f2.data(),
          b_range = b_range,
          ss      = ss_)
        other_data.set_selected(sel, f2.data()*o.k()*flex.exp(-o.b()*ss_))
    return other.array(data = other_data)

  def apply_debye_waller_factors(self,
        u_iso=None,
        b_iso=None,
        u_cart=None,
        b_cart=None,
        u_cif=None,
        u_star=None,
        apply_to_sigmas=True,
        exp_arg_limit=50,
        truncate_exp_arg=False):
    """
    Given an isotropic or anisotropic displacement or B-factor, apply
    resolution-dependent scale factors to the data (and optionally sigmas).

    :param u_iso: Isotropic displacement (in Angstroms)
    :param b_iso: Isotropic B-factor (8*pi^2*u_iso^2)
    :param u_cart: Anisotropic displacement tensor
    :param b_cart: Anisotropic B-factor
    :param u_star: Anisotropic displacement tensor in fractional space
    :param u_cif: Anisotropic displacement tensor, dimensionless basis
    :param apply_to_sigmas: Also scale sigma values (if present)
    :returns: cctbx.miller.array object with scaled data
    """
    dws = self.debye_waller_factors(
      u_iso=u_iso, b_iso=b_iso,
      u_cart=u_cart, b_cart=b_cart,
      u_cif=u_cif, u_star=u_star,
      exp_arg_limit=exp_arg_limit, truncate_exp_arg=truncate_exp_arg).data()
    d = self.data() * dws
    s = self.sigmas()
    if (s is not None and apply_to_sigmas):
      s = s * dws
    return self.customized_copy(data=d, sigmas=s)

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

  def sum(self, use_binning=False, use_multiplicities=False, squared=False):
    mean = self.mean(
      use_binning=use_binning,
      use_multiplicities=use_multiplicities,
      squared=squared)
    if not use_binning:
      return self.size() * mean
    else:
      counts = self.binner().counts()
      for i in self.binner().range_used():
        if (mean.data[i] is not None):
          mean.data[i] *= counts[i]
      return mean

  def sum_sq(self, use_binning=False, use_multiplicities=False):
    return self.sum(
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

  def count_and_fraction_in_bins(self,
        data_value_to_count,
        count_not_equal=False):
    assert self.binner() is not None
    assert self.data().size() == self.indices().size()
    max_n = 0
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      data_sel = self.data().select(sel)
      n = data_sel.count(data_value_to_count)
      if (count_not_equal): n = data_sel.size() - n
      max_n = max(max_n, n)
      results.append((n, n/(max(1,data_sel.size()))))
    data_fmt = "%%%dd %%6.4f" % len("%d"%max_n)
    return binned_data(binner=self.binner(), data=results, data_fmt=data_fmt)

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

  def normalize(self, reflections_per_bin=150, eps_fac=None):
    """
    Compute E-values: E = (F/eps**0.5) / rms of (F/eps**0.5)
    This is 'Karle' approach, that is not using overall B from Wilson plot.
    """
    f_obs = self
    if(self.is_xray_intensity_array()):
      f_obs = self.f_sq_as_f()
    else:
      assert self.is_xray_amplitude_array()
    eps = f_obs.epsilons().data().as_double()
    if(eps_fac is not None):
      eps = eps * eps_fac
    f_obs.setup_binner(reflections_per_bin = reflections_per_bin)
    E = flex.double(f_obs.data().size(), 0)
    for i_bin in f_obs.binner().range_used():
      bin_sel = f_obs.binner().selection(i_bin)
      fo = f_obs.data().select(bin_sel)
      if(fo.size()==0): continue
      e = eps.select(bin_sel)
      fo_eps = fo/flex.sqrt(e)
      E_bin = fo_eps/(flex.sum(fo_eps**2)/fo_eps.size())**0.5
      E = E.set_selected(bin_sel, E_bin)
    return self.array(data = E)

  def second_moments_centric_acentric(self, reflections_per_bin=150,
        eps_fac=None):
    """
    Compute <E**4>/<E**2>**2 for centric and acentric reflections.
    """
    e = self.normalize(reflections_per_bin=reflections_per_bin,
      eps_fac=eps_fac).data()
    centrics_selection  = self.centric_flags().data()
    e_sq = e**2
    result = []
    for prefix, sel in [("centric",  centrics_selection),
                        ("acentric",~centrics_selection)]:
      if(sel.count(True)==0): continue
      e_ = e.select(sel)
      e_sq_ = e_sq.select(sel)
      result.append((prefix, flex.mean(e_sq_*e_sq_)/flex.mean(e_sq_)**2))
    return result

  def amplitude_quasi_normalisations(self, d_star_power=1):
    """ A miller.array whose data N(h) are the normalisations to convert
    between locally normalised E's and F's:
    E(h) = F(h) / N(h)

    self features the F's, which are then binned with the current binner
    and N(h) is the average of F's in the bin h belongs to.
    """
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
    return array(self, flex.sqrt(mean_f_sq_over_epsilon_interp))

  def intensity_quasi_normalisations(self, d_star_power=1):
    """ A miller.array whose data N(h) are the normalisations to convert
      between locally normalised E^2's and I's:
      E^2(h) = I(h) / N(h)

      Intensities are binned with the current binner
      and N(h) is the average of I's in the bin h belongs to.
      """

    # see also cctbx.miller.array.amplitude_quasi_normalisations()

    assert self.binner() is not None
    assert self.binner().n_bin_d_too_large_or_small() == 0
    assert self.data().all_ge(0)
    assert self.observation_type() is None or self.is_xray_intensity_array()

    epsilons = self.epsilons().data().as_double()
    mean_f_sq_over_epsilon = flex.double()
    for i_bin in self.binner().range_used():
      sel = self.binner().selection(i_bin)
      sel_f_sq = self.data().select(sel)
      if (sel_f_sq.size() > 0):
        sel_epsilons = epsilons.select(sel)
        sel_f_sq_over_epsilon = sel_f_sq / sel_epsilons
        mean_f_sq_over_epsilon.append(flex.mean(sel_f_sq_over_epsilon))
      else:
        mean_f_sq_over_epsilon.append(0)
    mean_f_sq_over_epsilon_interp = self.binner().interpolate(
      mean_f_sq_over_epsilon, d_star_power)
    return array(self, mean_f_sq_over_epsilon_interp)

  def quasi_normalize_structure_factors(self, d_star_power=1):
    normalisations = self.amplitude_quasi_normalisations(d_star_power)
    q = self.data() / normalisations.data()
    return array(self, q)

  def phased_translation_function_coeff(self, phase_source, f_calc, fom=None):
    assert self.indices().size() == f_calc.indices().size()
    assert self.indices().size() == phase_source.indices().size()
    assert self.is_real_array()
    assert f_calc.is_complex_array()
    if(fom is not None):
      assert fom.is_real_array()
      assert self.indices().size() == fom.indices().size()
    f1 = self.array(data = self.data()).phase_transfer(
      phase_source = phase_source).data()
    f2 = flex.conj(f_calc.data())
    coeff = f1 * f2
    if(fom is not None):
      coeff = coeff * fom.data()
    if(fom is None):
      f3 = flex.sum(flex.abs(self.data()) * flex.abs(self.data()))
    else:
      f3 = flex.sum(fom.data() * fom.data() * flex.abs(self.data()) *
        flex.abs(self.data()))
    f4 = flex.sum(flex.abs(f_calc.data()) * flex.abs(f_calc.data()))
    den = math.sqrt(f3 * f4)
    assert den != 0
    return self.array(data = coeff/den)

  def __repr__(self):
    """
    Emit a string for debugging of the labels, type of data
    and sigmas array present within this miller_array.
    """
    mstr = self.crystal_symmetry().__repr__()
    if self._info:
      mstr = mstr + "\n" + str(self._info.label_string())
    mstr = mstr + "\n" + self._data.__repr__()
    if self._sigmas:
      mstr = mstr + "\n" + self._sigmas.__repr__()
    mstr = mstr + "\nsize: %d"  %self._data.size()
    return mstr + "\n"

  def __abs__(self):
    """
    Return a copy of the array with data replaced by absolute values, i.e.
    complex arrays will be converted to amplitudes.  Enables abs(array).
    """
    return array(self, flex.abs(self.data()), self.sigmas())

  def norm(self):
    assert isinstance(self.data(), flex.complex_double)
    return array(self, flex.norm(self.data()))

  def arg(self, deg=False):
    return array(self, flex.arg(self.data(), deg))

  def amplitudes(self):
    """
    For a complex array, return array of absolute values.
    """
    assert isinstance(self.data(), flex.complex_double)
    assert self.sigmas() is None
    return abs(self)

  def intensities(self):
    assert isinstance(self.data(), flex.complex_double)
    assert self.sigmas() is None
    return self.norm()

  def phases(self, deg=False):
    """
    For a complex array, return the array of its phases (in radians by default).
    """
    assert isinstance(self.data(), flex.complex_double)
    assert self.sigmas() is None
    return self.arg(deg)

  def merge_equivalents(self, algorithm="gaussian",
                        incompatible_flags_replacement=None,
                        use_internal_variance=True):
    """
    Given a non-unique array, merge the symmetry-related reflections (keeping
    anomalous flag).

    :returns: a merge_equivalents object, from which the merged array may \
    be extracted by calling the array() method.
    """
    return merge_equivalents(
      self, algorithm,
      incompatible_flags_replacement=incompatible_flags_replacement,
      use_internal_variance=use_internal_variance)

  def as_non_anomalous_array(self):
    """
    Return a copy of the array with identical contents (keeping original
    flex arrays) but with the anomalous flag set to false.
    """
    return array(
      miller_set=self.as_non_anomalous_set(),
      data=self.data(),
      sigmas=self.sigmas()).set_observation_type(self)

  def as_anomalous_array(self):
    """
    Return a copy of the array with identical contents (keeping original
    flex arrays) but with the anomalous flag set to true.
    """
    return array(
      miller_set=self.as_anomalous_set(),
      data=self.data(),
      sigmas=self.sigmas()).set_observation_type(self)

  def average_bijvoet_mates(self):
    """
    Given an anomalous array, merge the anomalous pairs and return the
    non-anomalous average.
    """
    if (self.is_complex_array() or self.is_hendrickson_lattman_array()):
      # centrics need special attention
      # very inefficient but simple implementation
      return self \
        .expand_to_p1() \
        .as_non_anomalous_array() \
        .merge_equivalents().array() \
        .customized_copy(crystal_symmetry=self) \
        .merge_equivalents().array()
    else:
      return self.as_non_anomalous_array().merge_equivalents().array()

  def eliminate_sys_absent(self, integral_only=False, log=None, prefix=""):
    """
    Remove all reflections which should be systematically absent in the current
    space group.
    """
    sys_absent_flags = self.sys_absent_flags(
      integral_only=integral_only).data()
    n = sys_absent_flags.count(True)
    if (n == 0): return self
    if (log is not None):
      if (integral_only): q = "integral "
      else: q = ""
      if (   isinstance(self.data(), flex.double)
          or isinstance(self.data(), flex.complex_double)):
        data_abs = flex.abs(self.data())
        c = ":"
      else:
        data_abs = None
        c = "."
      print(prefix + "Removing %d %ssystematic absence%s%s" % (
        n, q, plural_s(n)[1], c), file=log)
    result = self.select(selection=~sys_absent_flags)
    if (log is not None):
      if (data_abs is not None):
        print(prefix + "  Average absolute value of:", file=log)
        mean_absences = flex.mean(data_abs.select(sys_absent_flags))
        print(prefix + "    Absences: %.6g" % mean_absences, file=log)
        if (n != data_abs.size()):
          mean_others = flex.mean(data_abs.select(~sys_absent_flags))
          print(prefix + "      Others: %.6g" % mean_others, file=log)
          if (mean_others != 0 and mean_others > mean_absences * 1.e-20):
            print(prefix + "       Ratio: %.6g" % (
              mean_absences / mean_others), file=log)
      print(file=log)
    return result

  def select_sys_absent(self, integral_only=False):
    return self.select(selection=self.sys_absent_flags(
      integral_only=integral_only).data())

  def __add__(self, other):
    """Overload the '+' operator."""
    assert self.indices() is not None
    assert self.data() is not None
    if (type(other) != type(self)):
      # add a scalar or compatible flex array
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

  def __imul__(self, other):
    assert self.indices() is not None
    assert self.data() is not None
    data = self.data()
    data *= other
    sigmas = self.sigmas()
    if sigmas is not None: sigmas *= other
    return self

  def __mul__(self, other):
    if isinstance(other, array):
      assert self.indices().all_eq(other.indices())
      return self.customized_copy(data=self.data() * other.data())
    result = self.deep_copy()
    result *= other
    return result

  def __itruediv__(self, other):
    assert self.indices() is not None
    assert self.data() is not None
    data = self.data()
    data /= other
    sigmas = self.sigmas()
    if sigmas is not None: sigmas /= other
    return self

  def __truediv__(self, other):
    """ This requires from __future__ import division """
    result = self.deep_copy()
    result /= other
    return result

  def value_at_index(self, hkl):
    """Extract the value of the array for the specified reflection h,k,l"""
    hkl_array = self.select_indices([hkl])
    assert (len(hkl_array.data()) == 1)
    return hkl_array.data()[0]

  def data_at_first_index(self, miller_index):
    """
    Returns the value of data of the first index matching
    `miller_index`. If the `miller_index` is not found in `self`,
    then returns ``None``.

    :param miller_index: Miller index as a 3-tuple
    :type miller_index: tuple
    :returns: int, float, complex, None -- data value or None
    """
    return self.at_first_index(self._data, miller_index)

  def sigma_at_first_index(self, miller_index):
    """
    Returns the value of sigmas of the first index matching
    `miller_index`. If the `miller_index` is not found in `self`,
    then returns ``None``.

    :param miller_index: Miller index as a 3-tuple
    :type miller_index: tuple
    :returns: int, float, complex, None -- sigmas value or None
    """
    assert self._sigmas is not None
    return self.at_first_index(self._sigmas, miller_index)

  def generate_bijvoet_mates(self):
    """
    If the array is not already anomalous, expand to generate anomalous pairs
    (without changing data).
    """
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
      elif (self.is_hendrickson_lattman_array()):
        data.extend(data.select(sel).conj())
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
    """
    Calculate correlation coefficient between two arrays (either globally or
    binned).

    :param other: another array of real numbers
    :param use_binning: calculate CC in resolution bins (default = calculate
                         a single global value)
    :param assert_is_similar_symmetry: check that arrays have compatible
                                       crystal symmetry
    :returns: a Python float (if use_binning=False), or a binned_data object
    """
    if (assert_is_similar_symmetry):
      assert self.is_similar_symmetry(other)
    assert self.is_real_array()
    assert other.is_real_array()
    assert not use_binning or self.binner() is not None
    lhs = self
    if (lhs.anomalous_flag() and not other.anomalous_flag()):
      other = other.generate_bijvoet_mates()
    elif (not lhs.anomalous_flag() and other.anomalous_flag()):
      lhs = lhs.generate_bijvoet_mates()
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

  def show_array(self, f=None, prefix="", deg=None):
    """Listing of Miller indices and data"""
    if (f is None): f = sys.stdout
    assert self.data().size() == self.indices().size()
    if (self.is_complex_array() and deg is not None):
      for h,d in zip(self.indices(), self.data()):
        print(prefix + str(h), d, complex_math.abs_arg(d, deg=deg), file=f)
    elif (self.sigmas() is None):
      for h,d in zip(self.indices(), self.data()):
        print(prefix + str(h), d, file=f)
    else:
      assert self.sigmas().size() == self.indices().size()
      for h,d,s in zip(self.indices(), self.data(), self.sigmas()):
        print(prefix + str(h), d, s, file=f)
    return self

  def fsc(self, other, smooth=False):
    """
    Compute Fourier Shell Correlation (FSC)
    """
    f1, f2 = self, other
    assert f1.indices().all_eq(f2.indices())
    # Get data and order
    ds = f1.d_spacings().data()
    d1 = f1.data()
    d2 = f2.data()
    s = flex.sort_permutation(ds)
    ds = ds.select(s)
    d1 = d1.select(s)
    d2 = d2.select(s)
    # Get bins
    bins = bins_mixed(d_spacings_sorted=ds)
    # Compute FSC
    fsc = flex.double()
    d   = flex.double()
    for bin in bins:
      x = d1[bin[0]:bin[1]]
      y = d2[bin[0]:bin[1]]
      z = ds[bin[0]:bin[1]]
      fsc.append(maptbx.cc_complex_complex(f_1 = x, f_2 = y))
      d.append(flex.mean(z))
    d_inv = 1/d
    # Smooth FSC curve
    # Smoothing can change the input array size causing an assertion crash!
    sel = d>15.
    d_const = None

    if d.size() > 10:

      if sel.count(True)>0:
        d_const   = d    .select(sel)
        fsc_const = fsc  .select(sel)
        d         = d    .select(~sel)
        fsc       = fsc  .select(~sel)
      half_window=10
      from scitbx import smoothing
      d, fsc = smoothing.savitzky_golay_filter(
        x=d,  y=fsc, half_window=half_window, degree=2)
      if d_const is not None:
        d  .extend(d_const)
        fsc.extend(fsc_const)
      d_inv = 1./d

    s = flex.sort_permutation(d_inv)
    return group_args(
      d=d.select(s), d_inv=d_inv.select(s), fsc=fsc.select(s))

  def d_min_from_fsc(self, other=None, fsc_curve=None, fsc_cutoff=0.143):
    """
    Compute Fourier Shell Correlation (FSC) and derive resolution based on
    specified cutoff.
    """
    if(fsc_curve is None):
      assert other is not None
      fsc_curve = self.fsc(other=other)
      if not fsc_curve: return group_args(fsc=None, d_min=None)
    else:
      assert other is None
    i_mid = None
    for i in range(fsc_curve.fsc.size()):
      if(fsc_curve.fsc[i]<fsc_cutoff):
        i_mid = i
        break
    #print "i_mid, cc:", fsc_curve.d[i_mid], fsc_curve.d_inv[i_mid]
    d_min = None
    d_mid = None
    if(i_mid is not None):
      d_mid = fsc_curve.d[i_mid]
      if(i_mid is not None):
        i_min = i_mid-5
        if i_min < 0:
          i_min = 0
        i_max = i_mid+6
        if i_max >= len(fsc_curve.fsc):
          i_max = len(fsc_curve.fsc) - 1
        on_slope=True
        if(fsc_cutoff>0.): # does not have to be on slope around fsc_cutoff=0
          on_slope = [
            fsc_curve.fsc[i_min]>fsc_cutoff,
            fsc_curve.fsc[i_max]<fsc_cutoff].count(True)==2
        if(on_slope or fsc_cutoff):
          x = fsc_curve.d_inv[i_min:i_max]
          y = fsc_curve.fsc[i_min:i_max]
          from scitbx.math import curve_fitting
          c,b,a = curve_fitting.univariate_polynomial_fit(x_obs=x, y_obs=y,
            degree=2, number_of_cycles=5).params
          c = c-fsc_cutoff
          det = b**2-4*a*c
          well_defined = [det >= 0., a != 0.]
          if(well_defined.count(True)==2):
            x1 = (-b+math.sqrt(det))/(2*a)
            x2 = (-b-math.sqrt(det))/(2*a)
            #print "x1,x2", 1./x1,1/x2
            if(x1*x2<0.):
              d_min = 1./max(x1,x2)
            elif(x1>0 and x2>0):
              d1,d2 = 1./x1, 1./x2
              diff1 = abs(d1-d_mid)
              diff2 = abs(d2-d_mid)
              if(diff1<diff2): d_min = d1
              else:            d_min = d2
              if(abs(d_mid-d_min)>0.25): d_min = None
    if(d_min is None): d_min = d_mid
    return group_args(fsc=fsc_curve, d_min=d_min)

  def map_correlation(self, other):
    # Can be replace with a faster version:
    #   return maptbx.cc_complex_complex(f_1 = x, f_2 = y)
    d1 = flex.abs(self.data())
    d2 = flex.abs(other.data())
    p1 = self.phases().data()
    p2 = other.phases().data()
    factor = math.sqrt( flex.sum_sq(d1) * flex.sum_sq(d2) )
    if (factor > 0):
      return flex.sum( d1 * d2 * flex.cos(p2 - p1) ) / factor
    return None

  def as_map_manager(self,
                     resolution_factor=1/4.,
                     crystal_gridding=None,
                     grid_step=None,
                     d_min=None,
                     d_max=None,
                     apply_sigma_scaling=True,
                     apply_volume_scaling=False,
                     wrapping=True):
    assert isinstance(self.data(), flex.complex_double)
    assert [apply_sigma_scaling, apply_volume_scaling].count(True) in [0,1]
    mc = self
    if([d_max, d_min].count(None)>0):
      mc = self.resolution_filter(d_min=d_min, d_max=d_max)
    fft_map_ = mc.fft_map(
      resolution_factor = resolution_factor,
      crystal_gridding  = crystal_gridding,
      symmetry_flags = maptbx.use_space_group_symmetry,
      grid_step         = grid_step)
    if(apply_sigma_scaling):  fft_map_.apply_sigma_scaling()
    if(apply_volume_scaling): fft_map_.apply_volume_scaling()
    return fft_map_.as_map_manager(wrapping=wrapping)

  def fft_map(self, resolution_factor=1/3,
                    d_min=None,
                    grid_step=None,
                    crystal_gridding=None,
                    symmetry_flags=None,
                    mandatory_factors=None,
                    max_prime=5,
                    assert_shannon_sampling=True,
                    f_000=None):
    """
    Calculate the FFT for the array, assuming the data are complex doubles.

    :param resolution_factor: when multiplied times the resolution limit, gives
                              the approximate grid spacing of the map.
    :param d_min: High-resolution cutoff
    :param crystal_gridding: optional gridding to use (overrides automatic
                             gridding)
    :param symmetry_flags: specifies how the grid should be constructed to
                           handle symmetry
    :param f_000: Optional F(0,0,0) value (scalar added to entire map)
    :returns: an fft_map object
    """
    if(crystal_gridding is not None):
      return fft_map(
        crystal_gridding=crystal_gridding,
        fourier_coefficients=self,
        f_000=f_000)
    else:
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

  def direct_summation_at_point(self, site_frac, sigma=None):
    """
    Calculates the exact map value at the specified fractional coordinate
    using direct Fourier summation.  Relatively slow but avoids interpolation
    errors.
    """
    assert (self.is_complex_array())
    map_coeffs = self
    if (self.space_group_info().type().number() != 1):
      map_coeffs = map_coeffs.expand_to_p1()
    if (not map_coeffs.anomalous_flag()):
      map_coeffs = map_coeffs.generate_bijvoet_mates()
    sum = maptbx.direct_summation_at_point(
      miller_indices=map_coeffs.indices(),
      data=map_coeffs.data(),
      site_frac=site_frac)
    if (sigma is not None):
      return sum / sigma
    else :
      return sum / self.unit_cell().volume()

  def complete_with_bin_average(self, reflections_per_bin=100):
    assert isinstance(self.data(), flex.double)
    cs = self.complete_set()
    ls = cs.lone_set(self)
    self.setup_binner(reflections_per_bin = reflections_per_bin)
    result = []
    for i_bin in self.binner().range_used():
      sel = self.binner().selection(i_bin)
      d_range = self.binner().bin_legend(
        i_bin=i_bin, show_bin_number=False, show_counts=False)
      ssel = self.select(selection=sel)
      d_max, d_min = ssel.d_max_min()
      data_mean = flex.mean(ssel.data())
      result.append([d_max, d_min, data_mean])
    data_lone = flex.double()
    indices = flex.miller_index()
    for d, mi in zip(ls.d_spacings().data(), ls.indices()):
      for r in result:
        if(d>=r[1] and d<=r[0]):
          data_lone.append(r[2])
          indices.append(mi)
          break
    lms = set(self, indices, anomalous_flag=False)
    la = array(lms, data_lone)
    return self.complete_with(other=la)

  def hoppe_gassmann_modification(self, mean_scale, n_iterations,
        resolution_factor=0.25, d_min=None):
    assert self.is_complex_array()
    fft_map = self.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_volume_scaling()
    map_data = fft_map.real_map_unpadded()
    maptbx.hoppe_gassman_modification(data=map_data, mean_scale=mean_scale,
      n_iterations=n_iterations)
    miller_array=self
    if(d_min is not None): miller_array = self.complete_set(d_min=d_min)
    return miller_array.structure_factors_from_map(
      map            = map_data,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)

  def double_step_filtration(
        self,
        complete_set=None,
        vol_cutoff_plus_percent=5.0,
        vol_cutoff_minus_percent=5.0,
        resolution_factor=0.25,
        scale_to=None):
    assert self.is_complex_array()
    fft_map = self.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    n_real = fft_map.n_real()
    del fft_map
    value_min = flex.min(map_data.as_1d())
    value_max = flex.max(map_data.as_1d())
    if (value_min == value_max):
      cutoffp = cutoffm = value_max
    else :
      # XXX this avoids the huge memory overhead of sorting the entire map to
      # determine cutoffs, but still not optimal.  to be revisited...
      hist = maptbx.histogram(map_data, n_bins=min(1000, map_data.size()))
      cutoffp, cutoffm = hist.get_percentile_cutoffs(
        map=map_data,
        vol_cutoff_plus_percent=vol_cutoff_plus_percent,
        vol_cutoff_minus_percent=vol_cutoff_minus_percent)
    map_data = maptbx.denmod_simple(
      map_data = map_data,
      n_real   = n_real,
      cutoffp  = cutoffp,
      cutoffm  = cutoffm)
    if(complete_set is None): complete_set = self.complete_set()
    sf = complete_set.structure_factors_from_map(
      map            = map_data,
      use_scale      = True,
      anomalous_flag = self.anomalous_flag(),
      use_sg         = True)
    if(scale_to is None): scale_to = self
    def scale(f1,f2,reflections_per_bin=500):
      assert f2.data().size() >= f1.data().size()
      f1_,f2_ = f1.common_sets(f2)
      reflections_per_bin = min(500, f1_.data().size())
      f1_.setup_binner(reflections_per_bin = reflections_per_bin)
      f2_.use_binning_of(f1_)
      ss = 1./flex.pow2(f1_.d_spacings().data()) / 4.
      scale = flex.double()
      ss_bin = flex.double()
      for i_bin in f1_.binner().range_used():
        sel = f1_.binner().selection(i_bin)
        sel_f1 = abs(f1_).select(sel).data()
        sel_f2 = abs(f2_).select(sel).data()
        n = flex.sum(sel_f1*sel_f2)
        d = flex.sum(sel_f2*sel_f2)
        if(d == 0 or n == 0): return None
        scale.append(n/d)
        ss_bin.append(flex.mean(ss.select(sel)))
      if(scale.size()>1):
        r = scitbx.math.gaussian_fit_1d_analytical(x=flex.sqrt(ss_bin), y=scale)
        ss = 1./flex.pow2(f2.d_spacings().data()) / 4.
        k_isotropic = r.a*flex.exp(-ss*r.b)
      else: k_isotropic = scale[0]
      return f2.customized_copy(data = f2.data()*k_isotropic)
    result = scale(f1=scale_to, f2=sf)
    if(result is None): result = sf
    return result

  def local_standard_deviation_map(self, radius,
                                         mean_solvent_density=0,
                                         resolution_factor=1/3,
                                         d_min=None,
                                         grid_step=None,
                                         symmetry_flags=None,
                                         mandatory_factors=None,
                                         max_prime=5,
                                         assert_shannon_sampling=True,
                                         f_000=None):
    # J. P. Abrahams and A. G. W. Leslie, Acta Cryst. (1996). D52, 30-42
    # This should really have been called "local_variance_map" because the
    # square root is not taken after local averaging of density-squared
    fft = self.fft_map(
      resolution_factor=resolution_factor,
      d_min=d_min,
      grid_step=grid_step,
      symmetry_flags=symmetry_flags,
      mandatory_factors=mandatory_factors,
      max_prime=max_prime,
      assert_shannon_sampling=assert_shannon_sampling,
      f_000=f_000)
    fft.apply_volume_scaling()
    from cctbx import miller
    temp = miller.structure_factor_box_from_map(
       map           = flex.pow2(fft.real_map_unpadded()-mean_solvent_density),
       crystal_symmetry = self.crystal_symmetry(),
       d_min = self.d_min(),
       include_000      = True)

    sphere_reciprocal = self.g_function(R=radius, s=1./temp.d_spacings().data())
    fourier_coeff = temp.array(data=temp.data()*sphere_reciprocal)
    fft = fft_map(
      crystal_gridding=self.crystal_gridding(
        d_min=d_min,
        resolution_factor=resolution_factor,
        grid_step=grid_step,
        symmetry_flags=symmetry_flags,
        mandatory_factors=mandatory_factors,
        max_prime=max_prime,
        assert_shannon_sampling=assert_shannon_sampling),
      fourier_coefficients=fourier_coeff).apply_volume_scaling()
    return fft

  def local_overlap_map(self, other, radius,
                                         resolution_factor=1/3,
                                         d_min=None,
                                         grid_step=None,
                                         symmetry_flags=None,
                                         mandatory_factors=None,
                                         max_prime=5,
                                         assert_shannon_sampling=True,
                                         f_000=None,
                                         include_000=False):
    # Based on local_standard_deviation_map above
    assert self.crystal_symmetry().unit_cell().is_similar_to(
        other.crystal_symmetry().unit_cell())
    complete_set = self.complete_set()
    if(include_000 or f_000 is not None):
      indices = complete_set.indices()
      indices.append((0,0,0))
      complete_set = complete_set.customized_copy(indices = indices)
    sphere_reciprocal = self.g_function(R=radius, s=1./complete_set.d_spacings().data())
    if d_min is None:
      d_min=self.d_min()
    fft = self.fft_map(
      resolution_factor=resolution_factor,
      d_min=d_min,
      grid_step=grid_step,
      symmetry_flags=symmetry_flags,
      mandatory_factors=mandatory_factors,
      max_prime=max_prime,
      assert_shannon_sampling=assert_shannon_sampling,
      f_000=f_000)
    fft.apply_sigma_scaling()
    other_fft = other.fft_map(
      resolution_factor=resolution_factor,
      d_min=d_min,
      grid_step=grid_step,
      symmetry_flags=symmetry_flags,
      mandatory_factors=mandatory_factors,
      max_prime=max_prime,
      assert_shannon_sampling=assert_shannon_sampling,
      f_000=f_000)
    other_fft.apply_sigma_scaling()

    map_data=fft.real_map_unpadded()
    map_data_mean=map_data.as_1d().min_max_mean().mean
    map_data=map_data-map_data_mean
    other_map_data=other_fft.real_map_unpadded()
    other_map_data_mean=other_map_data.as_1d().min_max_mean().mean
    other_map_data=other_map_data-other_map_data_mean

    # get local overlap. They are normalized and have means of zero, so
    # overlap is similar to correlation coefficient:
    # CC is similar to:  <x1 * x2>   (not exactly as <x1**2> is not constant)

    overlap_map_data = map_data * other_map_data

    temp = complete_set.structure_factors_from_map(overlap_map_data)

    fourier_coeff = complete_set.array(data=temp.data()*sphere_reciprocal)
    fft = fft_map(
      crystal_gridding=self.crystal_gridding(
        d_min=d_min,
        resolution_factor=resolution_factor,
        grid_step=grid_step,
        symmetry_flags=symmetry_flags,
        mandatory_factors=mandatory_factors,
        max_prime=max_prime,
        assert_shannon_sampling=assert_shannon_sampling),
      fourier_coefficients=fourier_coeff).apply_volume_scaling()
    return fft


  def patterson_map(self, resolution_factor=1/3,
                          d_min=None,
                          symmetry_flags=None,
                          mandatory_factors=None,
                          max_prime=5,
                          assert_shannon_sampling=True,
                          f_000=None,
                          sharpening=False,
                          origin_peak_removal=False):
    """
    Calculate an unphased Patterson map.
    """
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

  def export_as_shelx_hklf(self,
        file_object=None,
        normalise_if_format_overflow=False,
        full_dynamic_range=False,
        scale_range=None):
    """
    Write reflections to a SHELX-format .hkl file.
    """
    from iotbx.shelx import hklf
    hklf.miller_array_export_as_shelx_hklf(
      miller_array=self,
      file_object=file_object,
      normalise_if_format_overflow=normalise_if_format_overflow,
      full_dynamic_range=full_dynamic_range,
      scale_range=scale_range)

  def export_as_cns_hkl(self,
        file_object,
        file_name=None,
        info=[],
        array_names=None,
        r_free_flags=None):
    """
    Write reflections to a CNS-format file.
    """
    from iotbx.cns.miller_array import export_as_cns_hkl as implementation
    implementation(self,
      file_object=file_object,
      file_name=file_name,
      info=info,
      array_names=array_names,
      r_free_flags=r_free_flags)

  def export_as_scalepack_unmerged(self,
      file_object=None,
      file_name=None,
      batch_numbers=None,
      spindle_flags=None,
      scale_intensities_for_scalepack_merge=False):
    """
    Write data in unmerged scalepack format.

    :param file_object: filehandle-like object
    :param file_name: output file to write
    :param batch_numbers: integer array indicating the batch (image) numbers
                          corresponding to the indices (optional)
    :param spindle_flags: integer array indicating the position of the
                          reflections on the detector (optional)
    """
    from iotbx.scalepack.no_merge_original_index import writer
    writer(i_obs=self,
      file_object=file_object,
      file_name=file_name,
      batch_numbers=batch_numbers,
      spindle_flags=spindle_flags,
      scale_intensities_for_scalepack_merge=
        scale_intensities_for_scalepack_merge)

  def write_mtz(self,
        file_name=None,
        column_root_label=None,
        column_types=None,
        label_decorator=None,
        title=None,
        crystal_name="crystal",
        project_name="project",
        dataset_name="dataset",
        wavelength=None):

    """
     Simple version of write_mtz for a single array, typically map coefficients
    """

    if column_root_label is None:
      if self.info() and self.info().labels:
        column_root_label=self.info().labels[0]
      else:
        column_root_label="F"

    mtz_dataset = self.as_mtz_dataset(
        column_root_label=column_root_label,
        column_types=column_types,
        label_decorator=label_decorator,
        title=title,
        crystal_name=crystal_name,
        project_name=project_name,
        dataset_name=dataset_name,
        wavelength=wavelength)

    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = file_name)

  def as_mtz_dataset(self,
        column_root_label,
        column_types=None,
        label_decorator=None,
        title=None,
        crystal_name="crystal",
        project_name="project",
        dataset_name="dataset",
        wavelength=None):
    """
    Generate an iotbx.mtz.dataset object for the array, which can be extended
    with additional arrays and eventually written to an MTZ file.
    """
    import iotbx.mtz
    if (wavelength is None):
      info = self.info()
      if (info is not None):
        wavelength = info.wavelength
    if (wavelength is None):
      wavelength = 1.0
    return iotbx.mtz.miller_array_as_mtz_dataset(self,
      column_root_label=column_root_label,
      column_types=column_types,
      label_decorator=label_decorator,
      title=title,
      crystal_name=crystal_name,
      project_name=project_name,
      dataset_name=dataset_name,
      wavelength=wavelength)

  def as_cif_block(self, array_type):
    # array_type is 'meas' or 'calc'
    import iotbx.cif
    return iotbx.cif.miller_arrays_as_cif_block(self, array_type).cif_block

  def as_cif_simple(self, array_type, out=None, data_name="global"):
    # array_type is 'meas' or 'calc'
    if out is None: out = sys.stdout
    import iotbx.cif
    cif = iotbx.cif.model.cif()
    cif[data_name] = self.as_cif_block(array_type=array_type)
    print(cif, file=out)

  def as_phases_phs(self,
        out,
        scale_amplitudes=True,
        phases=None,
        phases_deg=None,
        figures_of_merit=None):
    """
    Write phases to .phs file.
    """
    import iotbx.phases
    iotbx.phases.miller_array_as_phases_phs(
      self=self,
      out=out,
      scale_amplitudes=scale_amplitudes,
      phases=phases,
      phases_deg=phases_deg,
      figures_of_merit=figures_of_merit)

  def amplitude_normalisations(self, asu_contents, wilson_plot=None):
    """ Overriden version of set.amplitude_normalisation which computes
    the Wilson parameters from the array data if wilson_plot is None. """
    if wilson_plot is None:
      from cctbx import statistics
      f_obs = self.as_amplitude_array()
      f_obs.setup_binner(n_bins=20)
      wilson_plot = statistics.wilson_plot(f_obs, asu_contents)
    return set.amplitude_normalisations(self, asu_contents, wilson_plot)

  def normalised_amplitudes(self,
                            asu_contents,
                            wilson_plot=None):
    return normalised_amplitudes(self, asu_contents, wilson_plot)

  def scale_factor(self, f_calc, weights=None, cutoff_factor=None,
                   use_binning=False):
    """
    The analytical expression for the least squares scale factor.

    K = sum(w * yo * yc) / sum(w * yc^2)

    If the optional cutoff_factor argument is provided, only the reflections
    whose magnitudes are greater than cutoff_factor * max(yo) will be included
    in the calculation.
    """
    if self.data() is None: return None
    assert not use_binning or self.binner() is not None
    if use_binning: assert cutoff_factor is None
    if weights is not None:
      assert weights.size() == self.data().size()
    assert f_calc.is_complex_array()
    assert f_calc.size() == self.data().size()
    if not use_binning:
      if self.data().size() == 0: return None
      obs = self.data()
      if self.is_xray_intensity_array():
        calc = f_calc.norm().data()
      else:
        calc = flex.abs(f_calc.data())
      if cutoff_factor is not None:
        assert cutoff_factor < 1
        sel = obs >= flex.max(self.data()) * cutoff_factor
        obs = obs.select(sel)
        calc = calc.select(sel)
        if weights is not None:
          weights = weights.select(sel)
      if weights is None:
        return flex.sum(obs*calc) / flex.sum(flex.pow2(calc))
      else:
        return flex.sum(weights * obs * calc) \
             / flex.sum(weights * flex.pow2(calc))
    results = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      weights_sel = None
      if weights is not None:
        weights_sel = weights.select(sel)
      results.append(
        self.select(sel).scale_factor(f_calc.select(sel), weights_sel))
    return binned_data(binner=self.binner(), data=results, data_fmt="%7.4f")

  def from_cif(cls, file_object=None, file_path=None, data_block_name=None):
    """
    Class method for building an array from a CIF file (or filehandle).
    Depends on iotbx.cif.
    """
    import iotbx.cif
    from iotbx.cif import builders
    arrays = iotbx.cif.cctbx_data_structures_from_cif(
      file_object=file_object, file_path=file_path,
      data_block_name=data_block_name,
      data_structure_builder=builders.miller_array_builder).miller_arrays
    if data_block_name is not None:
      return arrays[data_block_name]
    else:
      return arrays
  from_cif = classmethod(from_cif)

  def shelxl_extinction_correction(self, x, wavelength):
    """
    Extinction parameter x, where Fc is multiplied by:
      k[1 + 0.001 x Fc^2 wavelength^3 / sin(2theta)]^(-1/4)

    See SHELX-97 manual, page 7-7 for more information.

    Note: The scale factor, k, is not applied nor calculated by
          this function. The scale factor should be calculated
          and applied ***AFTER*** the application of the extinction
          corrections.
    """
    assert self.is_complex_array()
    fc2 = self.as_intensity_array().data()
    sin_2_theta = self.unit_cell().sin_two_theta(self.indices(), wavelength)
    correction = 0.001 * x * fc2 * math.pow(wavelength, 3) / sin_2_theta
    correction += 1
    correction = flex.pow(correction, -1./4)
    return correction

  def apply_shelxl_extinction_correction(self, x, wavelength):
    correction = self.shelxl_extinction_correction(x, wavelength)
    return self.customized_copy(data=self.data() * correction)

  def f_obs_f_calc_fan_outlier_selection(self,
        f_calc,
        offset_low=0.05,
        offset_high=0.10,
        also_return_x_and_y=False):
    """\
      Preconditions (not checked explicitly):
        self is amplitude array,
        f_calc is complex array or amplitude array.
    """
    assert f_calc.indices().all_eq(self.indices())
    if (f_calc.is_complex_array()):
      x = flex.abs(f_calc.data())
    else:
      x = f_calc.data().deep_copy()
    y = self.data().deep_copy()
    if (flex.min(y) < 0):
      return None
    sum_xx = flex.sum_sq(x)
    sum_xy = flex.sum(x * y)
    if (sum_xx == 0):
      return None
    x *= (sum_xy / sum_xx)
    s = max(flex.max(x), flex.max(y))
    if (s == 0):
      return None
    x *= (1/s)
    y *= (1/s)
    m_low = (1-offset_high) / (1-offset_low)
    b_low = -m_low * offset_low
    m_high = 1/m_low
    b_high = offset_low
    result = (
        (y < m_low  * x + b_low)
      | (y > m_high * x + b_high))
    if (also_return_x_and_y):
      return result, x, y
    return result

  def as_xray_observations(self, scale_indices=None, twin_fractions=None,
                           twin_components=None):
    assert self.observation_type() is None or (
           self.is_xray_amplitude_array() or self.is_xray_intensity_array())
    assert self.is_real_array()
    assert self.sigmas() is not None
    from cctbx.xray import observations
    tw_cmps = ()
    if twin_components is not None:
      tw_cmps = twin_components
    if scale_indices is not None: # HKLF 5
      assert scale_indices.size() == self.indices().size()
      assert not (twin_fractions is None or len(twin_fractions) == 0)
      assert not twin_components
      result = observations.observations(
        self.indices(), self.data(), self.sigmas(),
        scale_indices, twin_fractions)
      result.fo_sq = array(
        miller_set=set(
          crystal_symmetry=self,
          indices=result.indices,
          anomalous_flag=self.anomalous_flag()),
        data=result.data,
        sigmas=result.sigmas).set_observation_type(self)
    else: #HKLF 4
      result = observations.observations(
        self.indices(), self.data(), self.sigmas(), tw_cmps)
      result.fo_sq = self
    # synchronise the life-time of the reference objects
    result.ref_twin_fractions = twin_fractions
    result.ref_twin_components = twin_components
    return result

  def french_wilson(self, **kwds):
    """
    Perform French-Wilson treatment of X-ray intensities to estimate the "true"
    intensities, replacing very weak and/or negative values, and takes the
    square root to obtain amplitudes.

    :returns: an array of all-positive X-ray amplitudes
    """
    assert self.is_xray_intensity_array()
    from cctbx import french_wilson
    kwds = dict(kwds)
    kwds['miller_array'] = self
    return french_wilson.french_wilson_scale(**kwds)

  def remove_cone(self, fraction_percent, vertex=(0,0,0), axis_point_1=(0,0,0),
        axis_point_2=(0,0,1), negate=False):
    """
    Remove reflections corresponding to a cone shape in reciprocal space with
    the apex at the origin.  Used to simulate incomplete data due to poor
    alignment of the crystal with the goniometer axis.
    """
    # single cone equation:
    # cos(half_opening_angle)*|R - VERTEX|*|AXIS| = (R-VERTEX,AXIS)
    # where R is any point on cone surface
    # double-cone requires AXIS*(-1)
    import scitbx.matrix
    fm = self.unit_cell().fractionalization_matrix()
    fm = scitbx.matrix.sqr(fm)
    fm = fm.transpose()
    axis = flex.double(
      fm * list(flex.double(axis_point_2)-flex.double(axis_point_1)))
    axis_length = math.sqrt(axis.dot(axis))
    vertex = flex.double(vertex)
    opening_angles = flex.double()
    for point in self.indices():
      point_minus_vertex = (flex.double(point)-vertex)
      point_minus_vertex = flex.double(fm * list(point_minus_vertex))
      point_minus_vertex_length = math.sqrt(
        point_minus_vertex.dot(point_minus_vertex))
      numerator = point_minus_vertex.dot(axis)
      denominator = point_minus_vertex_length*axis_length
      opening_angle_deg_for_point = 0
      if(point_minus_vertex_length>0):
        assert denominator != 0
        ratio = numerator / denominator
        if(abs(1.-abs(ratio))<1.e-3): ratio = 1.0
        opening_angle_deg_for_point = 180/math.pi * math.acos(ratio)
      opening_angles.append(opening_angle_deg_for_point)
    # smaller than 1 step will make selection accuracy higher
    for oa in range(1,180):
      sel = opening_angles<=oa
      if(100.*sel.count(True)/sel.size() > fraction_percent): break
    if(negate): return self.select(selection = sel)
    return self.select(selection = ~sel)

  def ellipsoidal_resolutions_and_indices_by_sigma(self, sigma_cutoff=3):
    if(self.sigmas() is None): return None
    sorted = self.sort(reverse=True)
    h_cut,k_cut,l_cut = None,None,None
    for mi, f, s, d in zip(sorted.indices(), sorted.data(), sorted.sigmas(),
                        sorted.d_spacings().data()):
      rsv = self.unit_cell().reciprocal_space_vector(mi)
      if(mi[1]==0 and mi[2]==0 and f/s>sigma_cutoff and h_cut is None):
        h_cut = (d, mi[0])
      if(mi[0]==0 and mi[2]==0 and f/s>sigma_cutoff and k_cut is None):
        k_cut = (d, mi[1])
      if(mi[0]==0 and mi[1]==0 and f/s>sigma_cutoff and l_cut is None):
        l_cut = (d, mi[2])
      if([h_cut,k_cut,l_cut].count(None)==0): break
    if([h_cut,k_cut,l_cut].count(None)>0 or
       [h_cut[1],k_cut[1],l_cut[1]].count(0)>0):
      min_mi, max_mi = self.min_max_indices()
      def helper(a,b):
        if abs(a)>b: return a
        if abs(a)<b: return b
        if abs(a)==b: return b
      hm = helper(min_mi[0], max_mi[0])
      km = helper(min_mi[1], max_mi[1])
      lm = helper(min_mi[2], max_mi[2])
      dh,dk,dl = self.d_min_along_a_b_c_star()
      if(h_cut is None or h_cut[1]==0): h_cut = (dh,hm)
      if(k_cut is None or k_cut[1]==0): k_cut = (dk,km)
      if(l_cut is None or l_cut[1]==0): l_cut = (dl,lm)
    return (h_cut,k_cut,l_cut)

  def show_mean_data_over_sigma_along_a_b_c_star(self):
    if(self.sigmas() is None): return None
    sorted = self.sort(reverse=True)
    h,k,l = [],[],[]
    for mi, f, s, d in zip(sorted.indices(), sorted.data(), sorted.sigmas(),
                        sorted.d_spacings().data()):
      data_over_sigma = f/s
      if(mi[1]==0 and mi[2]==0): h.append((d, mi[0], data_over_sigma))
      if(mi[0]==0 and mi[2]==0): k.append((d, mi[1], data_over_sigma))
      if(mi[0]==0 and mi[1]==0): l.append((d, mi[2], data_over_sigma))
    print("a*                      b*                    c*")
    print("d         h   F/sigma   d        k  F/sigma   d        l  F/sigma")
    for i in range( max(len(h),max(len(k),len(l))) ):
      blanc = " "*21
      try: ast = "%7.3f %4d %8.3f"%h[i]
      except Exception: ast = blanc
      try: bst = "%7.3f %4d %8.3f"%k[i]
      except Exception:
        if(len(ast.strip())==0): bst = blanc
        else: bst = blanc
      try: cst = "%7.3f %4d %8.3f"%l[i]
      except Exception: cst = ""
      print(ast,bst,cst)

  def ellipsoidal_truncation_by_sigma(self, sigma_cutoff=3):
    h_cut,k_cut,l_cut = self.ellipsoidal_resolutions_and_indices_by_sigma(
      sigma_cutoff = sigma_cutoff)
    rp = self.unit_cell().reciprocal_parameters()
    ehm = abs(h_cut[1]*rp[0])
    ekm = abs(k_cut[1]*rp[1])
    elm = abs(l_cut[1]*rp[2])
    if([h_cut,k_cut,l_cut].count(None)>0): return self.deep_copy()
    selection = flex.bool(self.indices().size(), False)
    #selection = flex.bool(self.indices().size(), True)
    data = self.data()
    sigmas = self.sigmas()
    for i_mi, mi in enumerate(self.indices()):
      rsv = self.unit_cell().reciprocal_space_vector(mi)
      r = math.sqrt((rsv[0]/ehm)**2 + (rsv[1]/ekm)**2 + (rsv[2]/elm)**2)
      if(r<=1): selection[i_mi] = True
      #if(r>1 and data[i_mi]/sigmas[i_mi]<sigma_cutoff): selection[i_mi] = False
    return self.select(selection=selection)

  def permute_d_range(self, d_max, d_min):
    """
    Randomly shuffle reflections within a given resolution range.  Used for
    control refinements to validate the information content of a dataset.
    """
    d_spacings = self.d_spacings().data()
    all_isel = flex.bool(d_spacings.size(), True).iselection()
    range_selection = (d_spacings <= d_max) & (d_spacings >= d_min)
    range_isel = range_selection.iselection()
    perm = flex.random_permutation(range_isel.size())
    range_isel_permuted = range_isel.select(perm)
    all_isel.set_selected(range_selection, range_isel_permuted)
    data = self.data().select(all_isel)
    sigmas = None
    if (self.sigmas() is not None):
      sigmas = self.sigmas().select(all_isel)
    return self.customized_copy(data=data, sigmas=sigmas)

  def is_unmerged_intensity_array(self):
    """
    Determine whether the array contains unmerged experimental observations
    or not.  In some files only the centric reflections will appear to be
    unmerged, so we specifically check the acentrics (if present).
    """
    if (not self.is_xray_intensity_array()) : return False
    centric_flags = self.centric_flags().data()
    centrics = self.select(centric_flags)
    acentric_flags = ~centric_flags
    if (acentric_flags.count(True) == 0):
      return (not centrics.is_unique_set_under_symmetry())
    acentrics = self.select(acentric_flags)
    return (not acentrics.is_unique_set_under_symmetry())

  # this is tested as part of phenix.merging_statistics (note that the exact
  # values are not reproducible)
  def cc_one_half(self, use_binning=False, n_trials=1, anomalous_flag=False,
                   return_n_refl=False):
    """
    Calculate the correlation between two randomly assigned pools of unmerged
    data ("CC 1/2").  If desired the mean over multiple trials can be taken.
    See Karplus PA & Diederichs K (2012) Science 336:1030-3 for motivation.
    This method assumes that the reflections still have the original indices
    and maps them to the ASU first; the underlying compute_cc_one_half
    function skips this method.
    """
    assert self.is_xray_intensity_array()
    if (not use_binning):
      tmp_array = self.customized_copy(
        anomalous_flag=anomalous_flag).map_to_asu()
      tmp_array = tmp_array.sort("packed_indices")
      return compute_cc_one_half(tmp_array, n_trials=n_trials,
                                 return_n_refl=return_n_refl)
    assert self.binner() is not None
    data = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      bin_array = self.select(sel)
      if (bin_array.size() == 0):
        data.append(None)
      else :
        data.append(bin_array.cc_one_half(
          n_trials=n_trials, return_n_refl=return_n_refl))
    return binned_data(binner=self.binner(), data=data, data_fmt="%6.3f")

  def cc_one_half_sigma_tau(self, use_binning=False, return_n_refl=False):
    """
    Calculation of CC1/2 by the 'sigma-tau' method, avoiding the random
    assignment into half-datasets of the above method.
    See Assmann et al., J. Appl. Cryst. (2016). 49, 1021–1028.

    var_y is the variange of the average intensities across the unique
    reflections of a resolution shell.
    var_e is the average variance of the observations contributing to the
    merged intensities
    """
    assert self.is_xray_intensity_array()
    if (not use_binning):
      # set the sigmas to 1, and calculate the mean intensities and internal variances
      intensities_copy = self.customized_copy(
        sigmas=flex.double(self.size(), 1))
      merging_internal = intensities_copy.merge_equivalents(
        use_internal_variance=True)
      merged = merging_internal.array().select(
        merging_internal.redundancies().data() > 1
      )
      if merged.size() <= 1:
        cc_one_half = 0
      else:
        internal_variances = flex.pow2(merged.sigmas())
        mav = flex.mean_and_variance(merged.data())
        var_y = mav.unweighted_sample_variance()
        var_e = 2 * flex.mean(internal_variances)
        cc_one_half = (var_y - 0.5 * var_e)/(var_y + 0.5 * var_e)
      if return_n_refl:
        return cc_one_half, merged.size()
      return cc_one_half
    assert self.binner() is not None
    data = []
    for i_bin in self.binner().range_all():
      sel = self.binner().selection(i_bin)
      bin_array = self.select(sel)
      if (bin_array.size() == 0):
        data.append(None)
      else :
        data.append(bin_array.cc_one_half_sigma_tau(return_n_refl=return_n_refl))
    return binned_data(binner=self.binner(), data=data, data_fmt="%6.3f")

  def half_dataset_anomalous_correlation(self, use_binning=False, return_n_pairs=False, return_split_datasets=False):
    """
    Calculate the correlation of the anomalous differences of two randomly
    assigned half-datasets (starting from unmerged data).
    """
    assert self.is_xray_intensity_array()
    tmp_array = self.customized_copy(anomalous_flag=True).map_to_asu()
    tmp_array = tmp_array.sort("packed_indices")
    if (not use_binning):
      split_datasets = split_unmerged(
        unmerged_indices=tmp_array.indices(),
        unmerged_data=tmp_array.data(),
        unmerged_sigmas=tmp_array.sigmas(),
        weighted=False)
      base_set = set(
        crystal_symmetry=self,
        indices=split_datasets.indices,
        anomalous_flag=True)
      assert base_set.is_unique_set_under_symmetry()
      array1 = array(
        miller_set=base_set,
        data=split_datasets.data_1)
      dano1 = array1.anomalous_differences()
      array2 = array(
        miller_set=base_set,
        data=split_datasets.data_2)
      dano2 = array2.anomalous_differences()
      assert dano1.indices().all_eq(dano2.indices())
      if return_split_datasets:
        return array1,array2
      if return_n_pairs:
        return dano1.correlation(other=dano2, use_binning=False).coefficient(), dano1.size()
      return dano1.correlation(other=dano2, use_binning=False).coefficient()
    assert self.binner() is not None
    tmp_array.use_binning_of(other=self)
    data = []
    for i_bin in tmp_array.binner().range_all():
      sel = tmp_array.binner().selection(i_bin)
      bin_array = tmp_array.select(sel)
      if (bin_array.size() == 0):
        data.append(None)
      else :
        data.append(bin_array.half_dataset_anomalous_correlation())
    return binned_data(binner=tmp_array.binner(), data=data, data_fmt="%6.3f")

  def cc_anom(self, *args, **kwds):
    """
    Alias for array.half_dataset_anomalous_correlation()
    """
    return self.half_dataset_anomalous_correlation(*args, **kwds)

  def r_anom(self):
    """Calculate R_anom, which measures the agreement between Friedel mates.
    Unlike CC_anom and various other R-factors (such as R_pim, which it is
    usually compared to), this requires merged data.

    .. math::
      R_{anom} = \\dfrac{\\sum_{hkl}{|I_{hkl} - I_{-h,-k,-l}|}}{\\sum_{hkl}{\\left \\langle I_{hkl} \\right \\rangle}}
    """
    assert self.is_xray_intensity_array()
    tmp_array = self.customized_copy(
      anomalous_flag=True).map_to_asu()
    if (not tmp_array.is_unique_set_under_symmetry()):
      tmp_array = tmp_array.merge_equivalents().array()
    d_ano = tmp_array.anomalous_differences()
    # XXX this appears to be consistent with the published descriptions, but
    # I need to confirm that this is what other programs do too
    i_mean = tmp_array.average_bijvoet_mates().common_set(other=d_ano)
    numerator = flex.sum(flex.abs(d_ano.data()))
    denominator = flex.sum(flex.abs(i_mean.data()))
    assert (denominator > 0)
    return numerator / denominator

  def twin_data(self, twin_law, alpha):
    """
    Apply a twin law to the data, returning an array of the same original type.

    :params twin_law: a valid twin law expressed as h,k,l operations
    :params alpha: predicted twin fraction (0 to 0.5)
    :returns: a new array with synthetically twinned data
    """
    assert (alpha is not None) and (0.0 <= alpha <= 0.5)
    assert (twin_law is not None)
    tmp_array = self.map_to_asu()
    if (tmp_array.is_xray_amplitude_array()):
      tmp_array = tmp_array.f_as_f_sq()
    assert tmp_array.is_xray_intensity_array()
    twin_law = sgtbx.rt_mx(twin_law, r_den=24, t_den=288)
    if twin_law.r().determinant() != 1:
      raise Sorry(
        "The determinant of the provided twin law is not equal to unity")
    cb_op = sgtbx.change_of_basis_op(twin_law)
    new_array = tmp_array.change_basis( cb_op ).map_to_asu()
    xa, xb = tmp_array.common_sets(other=new_array,
      assert_is_similar_symmetry=False)
    new_data = (1.0-alpha)*xa.data() + alpha*xb.data()
    new_sigmas = new_data / 100.0 # FIXME this seems wrong...
    twinned = xa.customized_copy(data=new_data,
      sigmas=new_sigmas).set_observation_type_xray_intensity()
    if (self.is_xray_amplitude_array()):
      return twinned.f_sq_as_f()
    return twinned

  def detwin_data(self, twin_law, alpha):
    """
    Detwin data using a known twin fraction, returning an array with the same
    original data type.

    :params twin_law: a valid twin law expressed as h,k,l operations
    :params alpha: predicted twin fraction (0 to 0.5)
    :returns: a new array with detwinned data
    """
    assert (alpha is not None) and (0.0 <= alpha <= 0.5)
    assert (twin_law is not None)
    tmp_array = self.map_to_asu()
    if (tmp_array.is_xray_amplitude_array()):
      tmp_array = tmp_array.f_as_f_sq()
    assert tmp_array.is_xray_intensity_array()
    twin_law = sgtbx.rt_mx(twin_law, r_den=24, t_den=288)
    if twin_law.r().determinant() != 1:
      raise Sorry(
        "The determinant of the provided twin law is not equal to unity")
    cb_op = sgtbx.change_of_basis_op(twin_law)
    new_array = tmp_array.change_basis( cb_op ).map_to_asu()
    # calling common_sets() automatically reorders the indices to be the same
    # between the two arrays, avoiding any further selections
    tmp_array, new_array = tmp_array.common_sets(other=new_array,
      assert_is_similar_symmetry=False)
    assert (tmp_array.size() == new_array.size())
    # adapted from mmtbx/scaling/twinning.h, but the use of common_sets()
    # and other boosted methods simplifies this greatly
    i_old1 = tmp_array.data()
    i_old2 = new_array.data()
    s_old1 = tmp_array.sigmas()
    s_old2 = new_array.sigmas()
    new_data = ((1.0 - alpha)*i_old1 - alpha*i_old2) / (1 - 2.0*alpha)
    new_sigmas = None
    tmp_mult = math.sqrt( 1-2*alpha +2*alpha*alpha)/(1-2.0*alpha)
    if (s_old1 is not None):
      new_sigmas = tmp_mult * flex.sqrt((flex.pow(s_old1,2) + s_old1*s_old2)/2)
    detwinned = tmp_array.customized_copy(
      data=new_data,
      sigmas=new_sigmas).common_set(other=self)
    if (self.is_xray_amplitude_array()):
      return detwinned.f_sq_as_f()
    return detwinned

  # TODO separate test - currently tested implicitly as part of change_symmetry
  def apply_change_of_basis(self,
        change_of_basis,
        eliminate_invalid_indices=True,
        out=None):
    """
    Encapsulates a variety of reindexing operations, including handling for a
    variety of corner cases.

    :param change_of_basis: Python str for change-of-basis operator
    :param eliminate_invalid_indices: remove reflections with non-integral
      indices
    :returns: new Miller array
    """
    if (out is None):
      out = sys.stdout
    miller_array = self
    print("Change of basis:", file=out)
    if   (change_of_basis == "to_reference_setting"):
      cb_op = miller_array.change_of_basis_op_to_reference_setting()
    elif (change_of_basis == "to_primitive_setting"):
      cb_op = miller_array.change_of_basis_op_to_primitive_setting()
    elif (change_of_basis == "to_niggli_cell"):
      cb_op = miller_array.change_of_basis_op_to_niggli_cell()
    elif (change_of_basis == "to_inverse_hand"):
      cb_op = miller_array.change_of_basis_op_to_inverse_hand()
    else:
      try :
        cb_op = sgtbx.change_of_basis_op(change_of_basis)
      except ValueError as e :
        raise Sorry(("The change-of-basis operator '%s' is invalid "+
          "(original error: %s)") % (change_of_basis, str(e)))
    if (cb_op.c_inv().t().is_zero()):
      print("  Change of basis operator in both h,k,l and x,y,z notation:", file=out)
      print("   ", cb_op.as_hkl(), file=out)
    else:
      print("  Change of basis operator in x,y,z notation:", file=out)
    print("    %s [Inverse: %s]" % (cb_op.as_xyz(),
      cb_op.inverse().as_xyz()), file=out)
    d = cb_op.c().r().determinant()
    print("  Determinant:", d, file=out)
    if (d < 0 and change_of_basis != "to_inverse_hand"):
      print("  **************************************************************", file=out)
      print("  W A R N I N G: This change of basis operator changes the hand!", file=out)
      print("  **************************************************************", file=out)
    if(eliminate_invalid_indices):
      sel = cb_op.apply_results_in_non_integral_indices(
        miller_indices=miller_array.indices())
      toss = flex.bool(miller_array.indices().size(),sel)
      keep = ~toss
      keep_array = miller_array.select(keep)
      toss_array = miller_array.select(toss)
      print("  Mean value for kept reflections:", \
        flex.mean(keep_array.data()), file=out)
      if (len(toss_array.data()) > 0):
        print("  Mean value for invalid reflections:", \
          flex.mean(toss_array.data()), file=out)
      miller_array=keep_array
    processed_array = miller_array.change_basis(cb_op=cb_op)
    print("  Crystal symmetry after change of basis:", file=out)
    crystal.symmetry.show_summary(processed_array, prefix="    ", f=out)
    return processed_array, cb_op

  def change_symmetry(self,
      space_group_symbol=None,
      space_group_info=None,
      volume_warning_threshold=0.001,
      expand_to_p1_if_necessary=True,
      remove_systematic_absences=True,
      merge_non_unique=True,
      log=None):
    """
    Encapsulates all operations required to convert the original data to a
    different symmetry (e.g. as suggested by Xtriage).  This includes
    reindexing and adjusting the unit cell parameters if necessary, and
    expansion to P1 (for moving to lower symmetry) or merging equivalents.

    :param space_group_symbol: Python str for space group symbol (any format)
    :param space_group_info: Pre-defined sgtbx.space_group_info object
    :param volume_warning_threshold: Cutoff for relative change in unit cell
      volume beyond which a warning is issued.
    :param expand_to_p1_if_necessary: When moving to lower symmetry, expand the
      data to P1 first.
    :param remove_systematic_absences: eliminate reflections that are
      systematically absent in the new symmetry.
    :param merge_non_unique: merge reflections that are no longer symmetry-
      unique under the new symmetry.
    :param log: filehandle-like object
    :returns: Miller array in the new symmetry
    """
    if (log is None):
      log = sys.stdout
    assert [space_group_symbol,space_group_info].count(None) == 1
    miller_array = self
    symm = miller_array.crystal_symmetry()
    space_group_old = symm.space_group()
    print("Current symmetry:", file=log)
    symm.show_summary(f=log, prefix="  ")
    # change to Niggli cell
    miller_array = miller_array.niggli_cell()
    symm = miller_array.crystal_symmetry()
    print("Niggli cell symmetry:", file=log)
    symm.show_summary(f=log, prefix="  ")
    if (space_group_info is None):
      space_group_info = sgtbx.space_group_info(space_group_symbol)
    space_group_new = space_group_info.group()
    if (space_group_new.n_smx() < space_group_old.n_smx()):
      if expand_to_p1_if_necessary :
        print("Changing to lower symmetry, expanding to P1 first", file=log)
        miller_array = miller_array.expand_to_p1()
      else :
        warnings.warn("This operation will result in incomplete data without "+
          "symmetry expansion!")
    unit_cell = symm.unit_cell()
    number = space_group_info.type().number()
    # reindex automatically if new lattice is P3 or P6 and angles are ~90,90,60
    # but only if the reference setting is requested!
    if (143 <= number < 195) and space_group_info.is_reference_setting():
      cb_op = sgtbx.change_of_basis_op("a,-b,-c")
      unit_cell_new = unit_cell.change_basis(cb_op)
      if (unit_cell_new.parameters()[-1] > unit_cell.parameters()[-1]):
        print("Reindexing with a,-b,-c", file=log)
        miller_array = miller_array.change_basis(cb_op)
        symm = miller_array.crystal_symmetry()
        unit_cell = symm.unit_cell()
    if (not space_group_info.group().is_compatible_unit_cell(unit_cell)):
      unit_cell_old = unit_cell
      unit_cell = space_group_info.group().average_unit_cell(unit_cell_old)
      print("Coercing unit cell into parameters compatible with %s" % \
        space_group_info, file=log)
      print("  Old cell: %s" % str(unit_cell_old.parameters()), file=log)
      print("  New cell: %s" % str(unit_cell.parameters()), file=log)
      volume_start = unit_cell_old.volume()
      volume_new = unit_cell.volume()
      volume_change_fraction = abs(volume_new - volume_start) / volume_new
      if (volume_change_fraction > volume_warning_threshold):
        warnings.warn("This operation will change the unit cell volume by "+
          "more than %.1f%%." % (volume_change_fraction*100), UserWarning)
    symm_new = crystal.symmetry(
      unit_cell=unit_cell,
      space_group_info=space_group_info)
    miller_array = miller_array.customized_copy(crystal_symmetry=symm_new)
    miller_array, cb_op = miller_array.apply_change_of_basis(
      change_of_basis="to_reference_setting",
      eliminate_invalid_indices=True,
      out=log)
    if merge_non_unique and (not miller_array.is_unique_set_under_symmetry()):
      miller_array = miller_array.merge_equivalents().array()
    if remove_systematic_absences :
      ma_old = miller_array.deep_copy()
      miller_array = miller_array.remove_systematic_absences()
      ls = ma_old.lone_set(other=miller_array)
    print("New Miller array:", file=log)
    miller_array.show_summary(f=log, prefix="  ")
    return miller_array

  def show_all_possible_systematic_absences(self, out=sys.stdout, prefix=""):
    """
    For each possible space group sharing the same basic intensity symmetry,
    show a list of possible systematically absent reflections and corresponding
    I/sigmaI.  Note that if the data have already been processed in a specific
    space group rather than the basic point group, for example P212121 instead
    of P222, all systematically absent reflections are likely to have been
    removed already.

    :returns: a systematic_absences_info object
    """
    return systematic_absences_info(self).show(out=out, prefix=prefix)

  #---------------------------------------------------------------------
  # Xtriage extensions - tested in mmtbx/scaling/tst_xtriage_twin_analyses.py
  def analyze_intensity_statistics(self, d_min=2.5,
    completeness_as_non_anomalous=None,  log=None):
    """
    Detect translational pseudosymmetry and twinning, using methods in
    Xtriage.  Returns a mmtbx.scaling.twin_analyses.twin_law_interpretation
    object.  (Requires mmtbx to be configured to be functional.)
    """
    import libtbx.load_env
    if (not libtbx.env.has_module("mmtbx")):
      raise ImportError("mmtbx is required for this functionality.")
    from mmtbx.scaling import twin_analyses as twin_analyses
    return twin_analyses.analyze_intensity_statistics(
      self=self,
      d_min=d_min,
      completeness_as_non_anomalous=completeness_as_non_anomalous,
      log=log)

  def has_twinning(self, d_min=2.5):
    """
    Convenience method for identifying twinned data.  Note that this is
    hugely inefficient if any other Xtriage analyses are planned, since it
    discards the other results.  Requires mmtbx.
    """
    return self.analyze_intensity_statistics(d_min=d_min).has_twinning()

  def average_neighbors(self, layers=1, include_origin=False,
      offset_list=None, average_with_cc=None):
    # Return array with values equal to average of values at neighboring indices
    # layers=1 means 26 or 27 neighboring points.
    # layers=0 means 6 or 7 closest points
    if average_with_cc:
      assert type(self.data()[0])==type((1+1j))  # must be complex

    sum_array=None
    sum_n_array=None
    if not offset_list:
      offset_list=get_offset_list(layers=layers,include_origin=include_origin)

    for offset in offset_list:
      offset_array=offset_indices(self,offset=offset)
      if average_with_cc:
        offset_array_matching, self_matching = offset_array.common_sets(self)
        offset_array=self.customized_copy(
           indices=offset_array_matching.indices(),
           data=offset_array_matching.data())
        self_array=self.customized_copy(
           indices=self_matching.indices(),
           data=self_matching.data())
        weight=self_array.map_correlation(offset_array)
      else:
        weight=1
      if sum_array is None:
        sum_array=offset_array*weight
      else:
        sum_array=sum_array.combine(offset_array,scale=False,
          scale_for_matches=weight)

      count_array=offset_array.customized_copy(
                      data=flex.double(offset_array.size(),abs(weight)))
      if sum_n_array is None:
        sum_n_array=count_array
      else:
        sum_n_array=sum_n_array.combine(count_array,scale=False)
    s=(sum_n_array.data() < 1.e-6)
    sum_n_array.data().set_selected(s,1)
    sum_array=sum_array.customized_copy(
      data=sum_array.data()*(1/sum_n_array.data()))

    # and extract values corresponding to original data
    sum_array_matching, self_matching = sum_array.common_sets(self)
    new_array=self.customized_copy(indices=sum_array_matching.indices(),
       data=sum_array_matching.data())
    return new_array


########################################################################
# END array class

class crystal_symmetry_is_compatible_with_symmetry_from_file:

  def __init__(self, miller_array,
         unit_cell_relative_length_tolerance=0.02,
         unit_cell_absolute_angle_tolerance=3.,
         working_point_group=None):
    self.miller_array = miller_array
    self.unit_cell_is_compatible = True
    self.space_group_is_compatible = True
    info = miller_array.info()
    if (info is None or info.crystal_symmetry_from_file is None): return
    ucf = info.crystal_symmetry_from_file.unit_cell()
    if (ucf is not None):
      uc = miller_array.unit_cell()
      if (uc is not None
          and not uc.is_similar_to(
            other=ucf,
            relative_length_tolerance=unit_cell_relative_length_tolerance,
            absolute_angle_tolerance=unit_cell_absolute_angle_tolerance)):
        self.unit_cell_is_compatible = False
    sgf = info.crystal_symmetry_from_file.space_group()
    if (sgf is not None):
      if (working_point_group is None):
        sg = miller_array.space_group()
        if (sg is not None):
          working_point_group = sg.build_derived_point_group()
      if (working_point_group is not None):
        point_group_from_file = sgf.build_derived_point_group()
        if (point_group_from_file != working_point_group):
          self.space_group_is_compatible = False

  def format_error_message(self, data_description):
    ma = self.miller_array
    what = []
    msg = ["  %s: %s" % (data_description, str(ma.info()))]
    if (not self.unit_cell_is_compatible):
      what.append("unit cell")
      msg.extend([
        "  Unit cell from file: " + str(
          ma.info().crystal_symmetry_from_file.unit_cell()),
        "    Working unit cell: " + str(ma.unit_cell())])
    if (not self.space_group_is_compatible):
      what.append("space group")
      msg.extend([
        "  Space group from file: " + str(
          ma.info().crystal_symmetry_from_file.space_group_info()),
        "    Working space group: " + str(ma.space_group_info())])
    if (len(what) != 0):
      if (len(what) == 2): what = ["crystal symmetry"]
      return "\n".join([
        "Working %s is not compatible with %s" % (what[0], what[0])
        + " from reflection file:"] + msg)
    return None


class normalised_amplitudes(object):
  """ E-values and related statistics """

  def __init__(self, miller_array, asu_contents, wilson_plot=None):
    assert miller_array.is_xray_amplitude_array()
    normalisations = miller_array.amplitude_normalisations(asu_contents,
                                                           wilson_plot)
    e = miller_array.data() / normalisations.data()
    self._array = array(
      miller_set=set(
        crystal_symmetry=miller_array.crystal_symmetry(),
        indices=miller_array.indices()).auto_anomalous(),
      data=e).set_observation_type_xray_amplitude()
    e_sq = flex.pow2(e);
    self._sum_e_sq_minus_1 = flex.sum(flex.abs(e_sq - 1))
    self._n_e_greater_than_2 = (e_sq > 4).count(True)

  def array(self):
    return self._array

  def mean_e_sq_minus_1(self):
    return self._sum_e_sq_minus_1/self._array.size()

  def percent_e_sq_gt_2(self):
    return (100.0 * self._n_e_greater_than_2)/self._array.size()


class merge_equivalents(object):
  """
  Wrapper for merging redundant observations to obtain a symmetry-unique
  array.  This also calculates some useful statistics resulting from the
  merging operation.  Normally this would not be instantiated directly, but
  instead obtained by calling array.merge_equivalents(...).
  """
  def __init__(self, miller_array, algorithm="gaussian",
               incompatible_flags_replacement=None,
               use_internal_variance=True):
    """
    :param miller_array: a non-unique array of experimental data
    :param algorithm: merging method (options are "gaussian" or "shelx")
    """
    assert algorithm in ["gaussian", "shelx"]
    self._r_linear = None
    self._r_square = None
    sigmas = None
    self._r_int = self._r_merge = self._r_meas = self._r_pim = None
    self._inconsistent_equivalents = None
    self.n_incompatible_flags = None
    data_type_str = miller_array.data().__class__.__name__
    merge_type = {
      "bool": ext.merge_equivalents_exact_bool,
      "int": ext.merge_equivalents_exact_int,
      "complex_double": ext.merge_equivalents_complex,
      "hendrickson_lattman": ext.merge_equivalents_hl,
    }.get(data_type_str, None)
    if (merge_type is not None):
      asu_array = miller_array.map_to_asu()
      perm = asu_array.sort_permutation(by_value="packed_indices")
      try :
        if data_type_str in ("bool", "int"):
          merge_ext = merge_type(
            asu_array.indices().select(perm),
            asu_array.data().select(perm),
            incompatible_flags_replacement=incompatible_flags_replacement)
        else:
          merge_ext = merge_type(
            asu_array.indices().select(perm),
            asu_array.data().select(perm))
      except RuntimeError as e :
        if ("merge_equivalents_exact: incompatible" in str(e)):
          raise Sorry(str(e) + " (mismatch between Friedel mates)")
        raise
      sigmas = None
      del asu_array
      if hasattr(merge_ext, "n_incompatible_flags"):
        self.n_incompatible_flags = merge_ext.n_incompatible_flags
    elif (isinstance(miller_array.data(), flex.double) ):
      asu_set = set.map_to_asu(miller_array)
      perm = asu_set.sort_permutation(by_value="packed_indices")
      if (miller_array.sigmas() is not None):
        if algorithm == "gaussian":
          merge_ext = ext.merge_equivalents_obs(
            asu_set.indices().select(perm),
            miller_array.data().select(perm),
            miller_array.sigmas().select(perm),
            use_internal_variance=use_internal_variance)
        elif algorithm == "shelx":
          merge_ext = ext.merge_equivalents_shelx(
            asu_set.indices().select(perm),
            miller_array.data().select(perm),
            miller_array.sigmas().select(perm))
          self._inconsistent_equivalents = merge_ext.inconsistent_equivalents
        else:
          raise RuntimeError("Programming error (should be unreachable).")
        sigmas = merge_ext.sigmas
      else:
        merge_ext = ext.merge_equivalents_real(
          asu_set.indices().select(perm),
          miller_array.data().select(perm))
        sigmas = None
      del asu_set
      self._r_linear = merge_ext.r_linear
      self._r_square = merge_ext.r_square
      self._r_int = merge_ext.r_int
      self._r_merge = merge_ext.r_merge
      self._r_meas = merge_ext.r_meas
      self._r_pim = merge_ext.r_pim
    elif (isinstance(miller_array.data(), flex.std_string) ):
      asu_array = miller_array.map_to_asu()
      perm = asu_array.sort_permutation(by_value="packed_indices")
      merge_ext = ext.merge_equivalents_string(
        asu_array.indices().select(perm),
        miller_array.data().select(perm))
      del asu_array
    else:
      raise RuntimeError(
        "cctbx.miller.merge_equivalents: unsupported array type:\n"
        "  data: %s\n"
        "  sigmas: %s" % (
          repr(miller_array.data()), repr(miller_array.sigmas())))
    self._array = array(
      miller_set=set(
        crystal_symmetry=miller_array,
        indices=merge_ext.indices,
        anomalous_flag=miller_array.anomalous_flag()),
      data=merge_ext.data,
      sigmas=sigmas).set_observation_type(miller_array)
    self._redundancies = merge_ext.redundancies

  def array(self):
    """
    Return the merged Miller array.
    """
    return self._array

  def redundancies(self):
    """
    Return an array representing the redundancy or multiplicity of each
    reflection in the merged array.
    """
    return self._array.array(data=self._redundancies)

  def r_linear(self):
    "R-linear = sum(abs(data - mean(data))) / sum(abs(data))"
    if (self._r_linear is None): return None
    return self._array.array(data=self._r_linear)

  def r_square(self):
    "R-square = sum((data - mean(data))**2) / sum(data**2)"
    if (self._r_square is None): return None
    return self._array.array(data=self._r_square)

  def r_int(self):
    return self._r_int

  def r_merge(self):
    """
    Standard (but flawed) metric of dataset internal consistency.

    .. math::
       R_{merge} = \\dfrac{\\sum_{hkl}{\\sum_{i}{|I_{i}(hkl) - \\left \\langle I_{i}(hkl) \\right \\rangle|}}}{\\sum_{hkl}{\\sum_{i}{I_{i}(hkl)}}}
    """
    return self._r_merge

  def r_meas(self):
    """
    Alternate metric of dataset internal consistency.  Explained in detail in
    Diederichs K & Karplus PA (1997) Nature Structural Biology 4:269-275.

    .. math::
       R_{meas} = \\dfrac{\\sum_{hkl}{ {\\left \\{ N(hkl) / [N(hkl) - 1] \\right \\} }^{1/2} \\times \\sum_{i}{|I_{i}(hkl) - \\left \\langle I_{i}(hkl) \\right \\rangle|}}}{\\sum_{hkl}{\\sum_{i}{I_{i}(hkl)}}}
    """
    return self._r_meas

  def r_pim(self):
    """
    Alternate metric of dataset internal consistency or quality.  Explained in
    detail in Weiss MS (2001) J Appl Cryst 34:130-135.

    .. math::
       R_{meas} = \\dfrac{\\sum_{hkl}{ {\\left \\{ 1 / [N(hkl) - 1] \\right \\} }^{1/2} \\times \\sum_{i}{|I_{i}(hkl) - \\left \\langle I_{i}(hkl) \\right \\rangle|}}}{\\sum_{hkl}{\\sum_{i}{I_{i}(hkl)}}}
    """
    return self._r_pim

  def inconsistent_equivalents(self):
    if self._inconsistent_equivalents != None:
      return self._inconsistent_equivalents
    return 0

  def r_sigma(self):
    return flex.sum(self.array().sigmas()) / flex.sum(self.array().data())

  def show_summary(self, n_bins=10, out=None, prefix=""):
    if (out is None): out = sys.stdout
    redundancies = self.redundancies().as_double()
    redundancies.setup_binner(n_bins=n_bins)
    red_mean = redundancies.mean(use_binning=True)
    selection = self.redundancies().data() > 1
    r_linear = self.r_linear()
    if (r_linear is not None):
      r_linear = r_linear.select(selection)
      r_linear.use_binning_of(redundancies)
      r_l_mean = r_linear.mean(use_binning=True)
    r_square = self.r_square()
    if (r_square is not None):
      r_square = r_square.select(selection)
      r_square.use_binning_of(redundancies)
      r_s_mean = r_square.mean(use_binning=True)
    fields = ["", "Min", "Max", "Mean"]
    if (r_linear is None): fields.append("")
    else: fields.append("R-linear")
    if (r_square is None): fields.append("")
    else: fields.append("R-square")
    lines = [fields]
    max_lengths = [len(field) for field in lines[0]]
    for i_bin in red_mean.binner.range_all():
      fields = [red_mean.binner.bin_legend(i_bin=i_bin, show_counts=False)]
      sel = red_mean.binner.selection(i_bin)
      r = self.redundancies().select(sel).data()
      if (r.size() == 0):
        fields.extend(["", ""])
      else:
        fields.append("%d" % flex.min(r))
        fields.append("%d" % flex.max(r))
      if (red_mean.data[i_bin] is None):
        fields.append("")
      else:
        fields.append("%.3f" % red_mean.data[i_bin])
      if (r_linear is None or r_l_mean.data[i_bin] is None):
        fields.append("")
      else:
        fields.append("%.4f" % r_l_mean.data[i_bin])
      if (r_square is None or r_s_mean.data[i_bin] is None):
        fields.append("")
      else:
        fields.append("%.4f" % r_s_mean.data[i_bin])
      lines.append(fields)
      max_lengths = [max(max_len,len(field))
        for max_len,field in zip(max_lengths, fields)]
    if (r_linear is not None):
      print(prefix+self.r_linear.__doc__, file=out)
    if (r_square is not None):
      print(prefix+self.r_square.__doc__, file=out)
    if (r_linear is not None or r_square is not None):
      print(prefix+"In these sums single measurements are excluded.", file=out)
    n = flex.sum(flex.int(max_lengths[1:4]))+4
    fmt = "%%%ds  %%%ds  %%%ds  %%%ds" % tuple(
      [max_lengths[0], n] + max_lengths[4:])
    fields = ["", "Redundancy"+" "*((n-10+1)//2)]
    for r in [r_linear, r_square]:
      if (r is None): fields.append("")
      else: fields.append("Mean  ")
    print(prefix + (fmt % tuple(fields)).rstrip(), file=out)
    fmt = "%%%ds  %%%ds  %%%ds  %%%ds  %%%ds  %%%ds" % tuple(max_lengths)
    for fields in lines:
      print(prefix + (fmt % tuple(fields)).rstrip(), file=out)

class fft_map(maptbx.crystal_gridding):
  """
  Container for an FFT from reciprocal space (complex double) into real space.
  Normally this is obtained by calling array.fft_map(...), not instantiated
  directly outside this module.  If the input array is anomalous, the
  resulting map will be a flex.complex_double (with grid accessor), otherwise
  it will be a flex.double.
  """
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
    map = maptbx.structure_factors.to_map(
      space_group=self.space_group(),
      anomalous_flag=self.anomalous_flag(),
      miller_indices=fourier_coefficients.indices(),
      structure_factors=fourier_coefficients.data(),
      n_real=self.n_real(),
      map_grid=flex.grid(n_complex),
      conjugate_flag=True)
    if (f_000 is not None):
      assert map.complex_map()[0] == 0j
      map.complex_map()[0] = complex(f_000)
    self._real_map_accessed = False
    if (not self.anomalous_flag()):
      self._real_map = rfft.backward(map.complex_map())
    else:
      self._complex_map = cfft.backward(map.complex_map())

  def anomalous_flag(self):
    return self._anomalous_flag

  def real_map(self, direct_access=True):
    """
    Extract the real component of the FFT'd map.

    :returns: a flex.double object with grid accessor.
    """
    if (not self.anomalous_flag()):
      assert ((self._real_map.is_padded()) or (not direct_access))
      if (direct_access):
        self._real_map_accessed = True
      return self._real_map
    else:
      return flex.real(self._complex_map)

  def as_map_manager(self, in_place=True, wrapping=True):
    '''
     Create a map_manager object from real_map_unpadded version of this map
    '''
    map_data=self.real_map_unpadded(in_place=in_place)
    from iotbx.map_manager import map_manager
    return map_manager(map_data=map_data,
      unit_cell_crystal_symmetry=self.crystal_symmetry(),
      unit_cell_grid=map_data.all(),
      wrapping=wrapping)

  def real_map_unpadded(self, in_place=True):
    """
    Extract the real component of the FFT'd map, removing any padding required
    for the FFT grid.

    :returns: a flex.double object with grid accessor.
    """
    if (in_place):
      assert (not self._real_map_accessed)
    result = self.real_map(direct_access=False)
    if (not result.is_padded()): return result
    elif (in_place):
      maptbx.unpad_in_place(map=result)
      return result
    else :
      return maptbx.copy(result, flex.grid(result.focus()))

  def complex_map(self):
    assert self.anomalous_flag()
    return self._complex_map

  def statistics(self):
    return maptbx.statistics(self.real_map(direct_access=False))

  def apply_scaling(self, scale):
    if (not self.anomalous_flag()):
      self._real_map *= scale
    else:
      self._complex_map *= scale
    return self

  def apply_fourier_scaling(self):
    return self.apply_scaling(scale=1/matrix.col(self.n_real()).product())

  def apply_volume_scaling(self):
    """
    Volume-scale the map values in place.
    """
    return self.apply_scaling(scale=1/self.unit_cell().volume())

  def apply_sigma_scaling(self):
    """
    Sigma-scale the map values in place.
    """
    statistics = self.statistics()
    if (statistics.sigma() == 0):
      return self
    scale = 1 / statistics.sigma()
    if (self.anomalous_flag()):
      scale = complex(scale)
    return self.apply_scaling(scale=scale)

  def peak_search(self, parameters=None, verify_symmetry=True):
    return self.tags().peak_search(
      parameters=parameters,
      map=self.real_map(direct_access=False),
      verify_symmetry=verify_symmetry)

  def as_ccp4_map(self,
                   file_name,
                   gridding_first=None,
                   gridding_last=None,
                   labels=["Values outside boundaries are wrapped inside",
                           "fft_map from Phenix"]):
    """
    Write the real component of the map to a CCP4-format file.
    """
    from iotbx import mrcfile
    map_data = self.real_map(direct_access=False)
    if (gridding_first is None) and (gridding_last is None):
      labels.append("wrapping_outside_cell")
    if gridding_first is None :
      gridding_first = (0,0,0)
    if gridding_last is None :
      gridding_last = tuple(self.n_real())
      # only write out the exact unit cell, without padding
      if (gridding_first == (0,0,0)):
        gridding_last = tuple([ (n-1) for n in gridding_last ])
    assert (len(labels) <= 10)
    mrcfile.write_ccp4_map(file_name=file_name,
      unit_cell=self.unit_cell(),
      space_group=self.space_group(),
      gridding_first=gridding_first,
      gridding_last=gridding_last,
      map_data=map_data,
      labels=flex.std_string(labels))

  def as_dsn6_map(self,
      file_name,
      gridding_first=None,
      gridding_last=None):
    """
    Write the real component of the map to a DSN6-format file.
    """
    from iotbx import dsn6
    map_data = self.real_map(direct_access=False)
    if gridding_first is None :
      gridding_first = (0,0,0)
    if gridding_last is None :
      gridding_last = tuple(self.n_real())
      # only write out the exact unit cell, without padding
      if (gridding_first == (0,0,0)):
        gridding_last = tuple([ (n-1) for n in gridding_last ])
    dsn6.write_dsn6_map(file_name=file_name,
      unit_cell=self.unit_cell(),
      gridding_first=gridding_first,
      gridding_last=gridding_last,
      map_data=map_data)

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
  i_patt = array(
    i_patt, data=i_patt.data() * flex.complex_double(i_patt.data().size(), 1))
  if (f_000 is not None):
    f_000 = f_000 * f_000
  return fft_map(crystal_gridding, i_patt, f_000)

def structure_factor_box_from_map(crystal_symmetry, map=None, n_real=None,
                                  anomalous_flag=False, include_000=False,
                                  f_000=None, d_min=None):
  #assert [map, n_real].count(None) in [0,2]
  if(map    is not None): assert n_real is None
  if(n_real is not None): assert map is None
  if(n_real is None):
    n_real = map.focus()
  max_index = [(i-1)//2 for i in n_real]
  complete_set = build_set(
    crystal_symmetry = crystal_symmetry,
    anomalous_flag   = anomalous_flag,
    max_index        = max_index)
  if(d_min is not None):
    complete_set = complete_set.resolution_filter(d_min=d_min)
  if(include_000 or f_000 is not None):
    indices = complete_set.indices()
    indices.append((0,0,0))
    complete_set = complete_set.customized_copy(indices = indices)
  if(map is None): return complete_set
  else:
    return complete_set.structure_factors_from_map(
      map            = map,
      use_scale      = True,
      anomalous_flag = anomalous_flag,
      use_sg         = False)

# this is tested as part of phenix.merging_statistics (note that the exact
# values are not reproducible)
def compute_cc_one_half(unmerged, n_trials=1, return_n_refl=False):
  """
  Implementation of array.cc_one_half, assuming that the reflections are
  already in the ASU.  Because the implementation uses random numbers, the
  function has the option to calculate the mean over multiple trials.
  """
  cc_all = []
  unmerged = unmerged.select(unmerged.sigmas() > 0)
  for x in range(n_trials):
    # this will obviously not be very random, but it's close enough for
    # the purpose of sampling different outcomes
    seed = 0
    if (n_trials > 1):
      seed = int(random.random()*10000)
    data_1 = data_2 = None
    split_datasets = split_unmerged(
      unmerged_indices=unmerged.indices(),
      unmerged_data=unmerged.data(),
      unmerged_sigmas=unmerged.sigmas(),
      seed=seed)
    data_1 = split_datasets.data_1
    data_2 = split_datasets.data_2
    cc = flex.linear_correlation(data_1, data_2).coefficient()
    cc_all.append(cc)
  cc_one_half = sum(cc_all) / n_trials
  if return_n_refl:
    return cc_one_half, data_1.size()
  return cc_one_half

class systematic_absences_info(object):
  """
  Container for information about possible systematically absent reflections in
  the array, trying both the current space group and all intensity-equivalent
  groups (i.e. all possible screw axis combinations).  This object would
  normally be instantiated directly from a Miller array, but is self-contained
  to enable saving as part of Xtriage results.

  :param obs: X-ray intensity (preferred) or amplitude array
  """
  def __init__(self, obs, was_filtered=None):
    self.was_filtered = was_filtered
    self.input_amplitudes = False
    assert (obs.sigmas() is not None)
    obs = obs.sort("packed_indices")
    if obs.is_xray_amplitude_array():
      obs = obs.f_as_f_sq()
      self.input_amplitudes = True
    assert obs.is_xray_intensity_array()
    if (not obs.is_unique_set_under_symmetry()):
      obs = obs.merge_equivalents().array()
    if obs.anomalous_flag():
      obs = obs.average_bijvoet_mates()
    self.space_group_info = obs.space_group_info()
    assert (self.space_group_info is not None)
    all_groups = self.space_group_info.reflection_intensity_equivalent_groups()
    point_group = self.space_group_info.group().build_derived_point_group()
    complete_sel = obs.customized_copy(
      space_group_info=point_group.info()).complete_set()
    self.space_group_symbols_and_selections = []
    self.n_possible_max = 0
    self.n_found_max = 0
    for group in all_groups :
      absent_sel = group.is_sys_absent(obs.indices()).iselection()
      all_possible = group.is_sys_absent(complete_sel.indices()).iselection()
      n_possible = len(all_possible)
      if (n_possible > self.n_possible_max):
        self.n_possible_max = n_possible
      if (len(absent_sel) > 0) : # systematic absences found
        if (len(absent_sel) > self.n_found_max):
          self.n_found_max = len(absent_sel)
        absences = obs.select(absent_sel)
        self.space_group_symbols_and_selections.append((group.info(), absences))
      elif (len(all_possible) == 0) : # no possible absences in this SG
        self.space_group_symbols_and_selections.append((group.info(), False))
      else : # absences possible, but not present in array
        self.space_group_symbols_and_selections.append((group.info(), None))

  # FIXME there must be a cleaner way to display this...
  def show(self, out=sys.stdout, prefix=""):
    """
    For each possible space group, show a list of possible systematically
    absent reflections and corresponding I/sigmaI.
    """
    if (self.n_possible_max == 0):
      print("No systematic absences possible in any intensity-equivalent groups.", file=out)
      return self
    if (self.input_amplitudes):
      print("""\
Please note that the input data were amplitudes, which means that weaker
reflections may have been modified by French-Wilson treatment or discarded
altogether, and the original intensities will not be recovered.  For best
results, use intensities as input.
""", file=out)
    for group_info, absences in self.space_group_symbols_and_selections :
      group_note = ""
      if (str(group_info) == str(self.space_group_info)):
        group_note = " (input space group)"
      if (absences == False):
        print(prefix+"%s%s: no systematic absences possible" % \
          (group_info, group_note), file=out)
      elif (absences is None):
        print(prefix+"%s%s: no absences found" % \
          (group_info, group_note), file=out)
      else :
        print(prefix+"%s%s:" % (group_info, group_note), file=out)
        for i_hkl, hkl in enumerate(absences.indices()):
          intensity = absences.data()[i_hkl]
          sigma = absences.sigmas()[i_hkl]
          indices_fmt = "(%4d, %4d, %4d)" % hkl
          if (sigma == 0):
            print(prefix+"  %s: i/sigi = undefined" % indices_fmt, file=out)
          else :
            print(prefix+"  %s: i/sigi = %6.1f" % (indices_fmt,
              intensity/sigma), file=out)
    return self
