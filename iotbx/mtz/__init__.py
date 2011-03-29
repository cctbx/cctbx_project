import cctbx.array_family.flex

import boost.python
ext = boost.python.import_ext("iotbx_mtz_ext")
from iotbx_mtz_ext import *
import iotbx_mtz_ext as ext

from iotbx.mtz import extract_from_symmetry_lib
from cctbx import xray
import cctbx.xray.observation_types
from cctbx import miller
import cctbx.crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
from libtbx.str_utils import show_string, overwrite_at, contains_one_of
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import warnings
import sys, os

expected_cmtz_struct_sizes = (
  (56, 84, 176, 500, 12336, 12488), # Linux 32-bit
  (64, 88, 184, 504, 12336, 12528), # Tru64
  #
  # 2010/04/19, addition of spg_confidence to struct SYMGRP in mtzdata.h
  (56, 84, 176, 500, 12340, 12492), # 32-bit
  (64, 88, 184, 504, 12340, 12536), # 64-bit
  #
  # 2010/06/30, addition of colsource...grpposn to struct MTZCOL
  # and xml...n_unknown_headers to struct MTZ in mtzdata.h
  (136, 84, 176, 500, 12340, 12504), # 32-bit
  (144, 88, 184, 504, 12340, 12560), # 64-bit
)
if (tuple(ext.cmtz_struct_sizes()) not in expected_cmtz_struct_sizes):
  warnings.warn("""Unexpected iotbx.mtz.cmtz_struct_sizes(): %s

os.name: %s
sys.platform: %s
boost.python.platform_info: %s

The iotbx.mtz module makes certain assumptions about the C structs in
the CMtz library. This warning appears if the sizes of the C structs
are not in a table of expected sizes. It may appear on platforms
where the iotbx.mtz module was not tested before, or if the CMtz
library structs have changed. If you see this warning, please send
the complete output to cctbx@cci.lbl.gov to give us the opportunity
to investigate.
Thank you!
""" % (
  str(tuple(ext.cmtz_struct_sizes())),
  os.name,
  sys.platform,
  boost.python.platform_info))

column_type_legend_source = "ccp4/doc/mtzformat.doc"
column_type_legend = {
  "H": "index h,k,l",
  "J": "intensity",
  "F": "amplitude",
  "D": "anomalous difference",
  "Q": "standard deviation",
  "G": "F(+) or F(-)",
  "L": "standard deviation",
  "K": "I(+) or I(-)",
  "M": "standard deviation",
  "E": "normalized amplitude",
  "P": "phase angle in degrees",
  "W": "weight (of some sort)",
  "A": "phase probability coefficients (Hendrickson/Lattman)",
  "B": "BATCH number",
  "Y": "M/ISYM, packed partial/reject flag and symmetry number",
  "I": "integer",
  "R": "real",
}

def default_column_types(miller_array):
  result = None
  if (miller_array.is_complex_array() and miller_array.sigmas() is None):
    assert not miller_array.is_xray_intensity_array()
    if (miller_array.anomalous_flag()):
      result = "GP"
    else:
      result = "FP"
  elif (miller_array.is_xray_intensity_array()):
    if (miller_array.anomalous_flag()):
      result = "K"
      if (miller_array.sigmas() is not None):
        result += "M"
    else:
      result = "J"
      if (miller_array.sigmas() is not None):
        result += "Q"
  elif (miller_array.is_xray_amplitude_array()
        or (miller_array.is_real_array()
            and miller_array.sigmas() is not None)):
    if (miller_array.anomalous_flag()):
      result = "G"
      if (miller_array.sigmas() is not None):
        result += "L"
    else:
      result = "F"
      if (miller_array.sigmas() is not None):
        result += "Q"
  elif ((   miller_array.is_bool_array()
         or miller_array.is_integer_array())
        and miller_array.sigmas() is None):
    result = "I"
  elif (miller_array.is_real_array()):
    if (miller_array.anomalous_flag()):
      result = "G"
    else:
      result = "F"
  elif (miller_array.is_hendrickson_lattman_array()):
    result = "AAAA"
  return result

anomalous_label_patterns = {
  "+": ("+", "PLUS"),
  "-": ("-", "MINU")
}

def is_anomalous_label(sign, label):
  return contains_one_of(label.upper(), anomalous_label_patterns[sign])

def are_anomalous_labels(sign, labels):
  for label in labels:
    if (not is_anomalous_label(sign, label)): return False
  return True

class label_decorator(__builtins__["object"]):

  def __init__(self,
        anomalous_plus_suffix="(+)",
        anomalous_minus_suffix="(-)",
        sigmas_prefix="SIG",
        sigmas_suffix="",
        phases_prefix="PHI",
        phases_suffix="",
        hendrickson_lattman_suffix_list=["A","B","C","D"]):
    assert len(hendrickson_lattman_suffix_list) == 4
    adopt_init_args(self, locals())

  def anomalous(self, root_label, sign=None):
    assert sign in (None, "+", "-")
    if (sign is None): return root_label
    if (sign == "+"): return root_label + self.anomalous_plus_suffix
    return root_label + self.anomalous_minus_suffix

  def sigmas(self, root_label, anomalous_sign=None):
    return self.anomalous(
      self.sigmas_prefix + root_label + self.sigmas_suffix,
      anomalous_sign)

  def phases(self, root_label, anomalous_sign=None):
    return self.anomalous(
      self.phases_prefix + root_label + self.phases_suffix,
      anomalous_sign)

  def hendrickson_lattman(self, root_label, i_coeff, anomalous_sign=None):
    assert 0 <= i_coeff < 4
    return self.anomalous(
      root_label + self.hendrickson_lattman_suffix_list[i_coeff],
      anomalous_sign)

def format_min_max(func, values):
  if (len(values) == 0): return "None"
  value = func(values)
  result = "%.2f" % value
  if (len(result) > 12): result = "%12.4E" % value
  return result

show_column_data_format_keywords = [
  "human_readable",
  "machine_readable",
  "spreadsheet"]

def tidy_show_column_data_format_keyword(input):
  if (input is None): return show_column_data_format_keywords[0]
  input = input.lower()
  for k in show_column_data_format_keywords:
    if (input.startswith(k[0])): return k
  raise Sorry(
      "Column data format keyword not recognized: %s\n" % show_string(input)
    + "  Valid keywords are: %s" % ", ".join(show_column_data_format_keywords))

class _(boost.python.injector, ext.object):

  def space_group_info(self):
    return sgtbx.space_group_info(group=self.space_group())

  def set_space_group_info(self, space_group_info, symbol=None):
    if (symbol is None):
      symbol = extract_from_symmetry_lib.ccp4_symbol(
        space_group_info=space_group_info,
        lib_name="symop.lib")
      if (symbol is None):
        symbol = "No.%d" % space_group_info.type().number()
    group = space_group_info.group()
    self.set_space_group_name(name=symbol)
    self.set_space_group_number(number=space_group_info.type().number())
    self.set_point_group_name(name=group.point_group_type())
    self.set_lattice_centring_type(
      symbol=group.conventional_centring_type_symbol())
    if (self.lattice_centring_type() == "\0"):
      self.set_lattice_centring_type(symbol="?")
    self.set_space_group(space_group=space_group_info.group())
    return self

  def set_hkl_base(self, unit_cell):
    assert self.n_crystals() == 0
    return self.add_crystal(
      name="HKL_base",
      project_name="HKL_base",
      unit_cell=unit_cell).set_id(id=0).add_dataset(
        name="HKL_base",
        wavelength=0).set_id(id=0)

  def n_columns(self):
    result = 0
    for crystal in self.crystals():
      for dataset in crystal.datasets():
        result += dataset.n_columns()
    return result

  def columns(self):
    for crystal in self.crystals():
      for dataset in crystal.datasets():
        for column in dataset.columns():
          yield column

  def column_labels(self):
    return [column.label() for column in self.columns()]

  def column_types(self):
    return [column.type() for column in self.columns()]

  def show_summary(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    p = prefix
    print >> out, p+"Title:", self.title()
    print >> out, p+"Space group symbol from file:", self.space_group_name()
    print >> out, p+"Space group number from file:", self.space_group_number()
    self.space_group_info().show_summary(
      f=out, prefix=p+"Space group from matrices: ")
    print >> out, p+"Point group symbol from file:", self.point_group_name()
    if (self.n_batches() > 0):
      print >> out, p+"Number of batches:", self.n_batches()
    print >> out, p+"Number of crystals:", self.n_crystals()
    print >> out, p+"Number of Miller indices:", self.n_reflections()
    if (self.n_crystals() > 0 and self.n_reflections() > 0):
      print >> out, p+"Resolution range: %.6g %.6g" % self.max_min_resolution()
    print >> out, p+"History:"
    for line in self.history():
      print >> out, p+" ", line.rstrip()
    for i_crystal,crystal in enumerate(self.crystals()):
      print >> out, p+"Crystal %d:" % (i_crystal+1)
      print >> out, p+"  Name:", crystal.name()
      print >> out, p+"  Project:", crystal.project_name()
      print >> out, p+"  Id:", crystal.id()
      crystal.unit_cell().show_parameters(f=out, prefix=p+"  Unit cell: ")
      print >> out, p+"  Number of datasets:", crystal.n_datasets()
      for i_dataset,dataset in enumerate(crystal.datasets()):
        print >> out, p+"  Dataset %d:" % (i_dataset+1)
        print >> out, p+"    Name:", dataset.name()
        print >> out, p+"    Id:", dataset.id()
        print >> out, p+"    Wavelength: %.6g" % dataset.wavelength()
        print >> out, p+"    Number of columns:", dataset.n_columns()
        if (dataset.n_columns() > 0):
          fields_list = [[
            "label", "#valid", "%valid", "min", "max", "type", ""]]
          max_field_lengths = [len(field) for field in fields_list[0]]
          max_field_lengths[-2] = 0
          for i_column,column in enumerate(dataset.columns()):
            valid_values = column.extract_valid_values()
            fields = [
              column.label(),
              "%d" % valid_values.size(),
              "%.2f%%" %(100.*valid_values.size()/max(1,self.n_reflections())),
              format_min_max(flex.min, valid_values),
              format_min_max(flex.max, valid_values),
              column.type()+":",
              column_type_legend.get(
                column.type(), "*** UNDEFINED column type ***")]
            fields_list.append(fields)
            for i,field in enumerate(fields):
              max_field_lengths[i] = max(max_field_lengths[i], len(field))
          format = "    %%-%ds %%%ds %%%ds %%%ds %%%ds %%%ds %%s" % tuple(
            max_field_lengths[:6])
          for fields in fields_list:
            print >> out, p+(format % tuple(fields)).rstrip()
    return self

  def _show_column_data_preparation(self):
    miller_indices = self.extract_miller_indices()
    labels = []
    pairs = []
    for column in self.columns():
      if (column.label() in ["H", "K", "L"]): continue
      labels.append(column.label())
      pairs.append((column.extract_values(), column.selection_valid()))
    return miller_indices, labels, pairs

  def show_column_data_human_readable(self, out=None):
    if (out is None): out = sys.stdout
    miller_indices, labels, pairs = self._show_column_data_preparation()
    if (miller_indices.size() == 0):
      h_width = 1
    else:
      h_width = len("-%d" % int(flex.max(flex.abs(
        miller_indices.as_vec3_double().as_double()))+.5))
    h_format = " ".join(["%%%dd" % h_width]*3)
    h_blank = " "*(h_width*3+2)
    def show_labels():
      print >> out, h_blank,
      for i,label in enumerate(labels):
        if (i and i % 4 == 0):
          print >> out
          print >> out, h_blank,
        print >> out, "%15s" % label,
      print >> out
      print >> out
    n_data_lines = 0
    print >> out, "Column data:"
    print >> out, "-"*79
    show_labels()
    for iref,h in enumerate(miller_indices):
      if (n_data_lines > 20):
        print >> out
        show_labels()
        n_data_lines = 0
      h_str = h_format % h
      print >> out, h_str,
      for i,(data,selection) in enumerate(pairs):
        if (i and i % 4 == 0):
          print >> out
          n_data_lines += 1
          print >> out, h_blank,
        if (selection[iref]): print >> out, "%15.6g" % data[iref],
        else:                 print >> out, "%15s" % "None",
      print >> out
      n_data_lines += 1
    print >> out, "-"*79
    return self

  def show_column_data_machine_readable(self, out=None):
    if (out is None): out = sys.stdout
    miller_indices, labels, pairs = self._show_column_data_preparation()
    print >> out, "Machine readable colum data:"
    print >> out, "Number of columns:", len(labels)
    print >> out, "Column labels (one per line):"
    for label in labels:
      print >> out, label
    print >> out, "Number of Miller indices:", miller_indices.size()
    print >> out, "Column data (HKL followed by data):"
    for iref,h in enumerate(miller_indices):
      print >> out, " ".join([str(i) for i in h])
      for i,(data,selection) in enumerate(pairs):
        if (selection[iref]): print >> out, "%.7g" % data[iref]
        else:                 print >> out, "None"
    print >> out, "End of column data."
    return self

  def show_column_data_spreadsheet(self, out=None):
    if (out is None): out = sys.stdout
    miller_indices, labels, pairs = self._show_column_data_preparation()
    print >> out, ",".join([label.replace(",","_")
      for label in ["H","K","L"]+labels])
    for iref,h in enumerate(miller_indices):
      row = [str(i) for i in h]
      for i,(data,selection) in enumerate(pairs):
        if (selection[iref]): row.append("%.7g" % data[iref])
        else:                 row.append("")
      print >> out, ",".join(row)
    return self

  def show_column_data(self, out=None, format="human_readable"):
    assert format in show_column_data_format_keywords
    return getattr(self, "show_column_data_"+format)(out=out)

  def change_basis_in_place(self,
        cb_op,
        new_space_group_info=None,
        assert_is_compatible_unit_cell=False):
    assert len(column_type_legend) == 17 # programmer alert
      # force update if column_type_legend is changed
    for column_type in self.column_types():
      if (column_type == "P"):
        raise RuntimeError(
          "In-place transformation of phase angles not implemented.")
      if (column_type == "A"):
        raise RuntimeError(
          "In-place transformation of Hendrickson-Lattman coefficients"
          " not implemented.")
    self.replace_miller_indices(cb_op.apply(self.extract_miller_indices()))
    if (new_space_group_info is None):
      new_space_group_info = self.space_group_info().change_basis(cb_op)
    self.set_space_group_info(space_group_info=new_space_group_info)
    for crystal in self.crystals():
      crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=crystal.unit_cell().change_basis(cb_op=cb_op),
        space_group_info=new_space_group_info,
        assert_is_compatible_unit_cell=assert_is_compatible_unit_cell)
      crystal.set_unit_cell_parameters(
        crystal_symmetry.unit_cell().parameters())
    # transform & symmetrize per-batch unit cell to support Scala 6.0 (NKS)
    # re-zero the U-matrix so nobody tries to use it downstream
    for batch in self.batches():
      batch_uc = uctbx.unit_cell(list(batch.cell()))
      batch_crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=batch_uc.change_basis(cb_op=cb_op),
        space_group_info=new_space_group_info,
        assert_is_compatible_unit_cell=assert_is_compatible_unit_cell)
      batch.set_cell(flex.float(
        batch_crystal_symmetry.unit_cell().parameters()))
      batch.set_umat(flex.float(
        (0,0,0,0,0,0,0,0,0)))

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    other_symmetry = crystal_symmetry
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="ccp4_mtz")
    result = []
    for crystal in self.crystals():
      crystal_symmetry_from_file = cctbx.crystal.symmetry(
        unit_cell=crystal.unit_cell(),
        space_group_info=self.space_group_info())
      crystal_symmetry = crystal_symmetry_from_file.join_symmetry(
        other_symmetry=other_symmetry,
        force=force_symmetry)
      for dataset in crystal.datasets():
        column_groups = self.group_columns(
          crystal_symmetry_from_file=crystal_symmetry_from_file,
          crystal_symmetry=crystal_symmetry,
          base_array_info=base_array_info,
          dataset=dataset)
        for column_group in column_groups:
          if (merge_equivalents
              and isinstance(column_group.data(), flex.double)
              and isinstance(column_group.sigmas(), flex.double)
              and column_group.sigmas().size() != 0
              and flex.min(column_group.sigmas()) > 0):
            merged_column_group = column_group.merge_equivalents().array()
            if (merged_column_group.indices().size()
                != column_group.indices().size()):
              merged_column_group.set_info(
                column_group.info().customized_copy(merged=True))
              column_group = merged_column_group
          result.append(column_group)
    return result

  def group_columns(self,
        crystal_symmetry_from_file,
        crystal_symmetry,
        base_array_info,
        dataset,
        strict=True):
    known_mtz_column_types = "".join(column_type_legend.keys())
    assert len(known_mtz_column_types) == 17 # safety guard
    all_columns = dataset.columns()
    all_column_labels = dataset.column_labels()
    all_column_types = mend_non_conforming_anomalous_column_types(
      dataset.column_types(), all_column_labels)
    groups = []
    i_column = -1
    while 1:
      i_column += 1
      if (i_column == len(all_columns)): break
      column = all_columns[i_column]
      if (strict and column.type() not in known_mtz_column_types):
        raise RuntimeError(
          'Unknown MTZ column type: "%s" (column label: "%s")' % (
            column.type(), column.label()))
      t0 = all_column_types[i_column]
      if (t0 == "H"): continue # skip h,k,l
      l0 = all_column_labels[i_column]
      remaining_types = all_column_types[i_column:]
      labels = None
      group = None
      observation_type = None
      if (t0 in "BYI"): # integer columns
        if (len(remaining_types) > 1
            and remaining_types[1] == t0
            and is_anomalous_label("+", l0)
            and is_anomalous_label("-", all_column_labels[i_column+1])):
          labels = all_column_labels[i_column:i_column+2]
          i_column += 1
          group = self.extract_integers_anomalous(*labels)
        else:
          labels = [l0]
          group = self.extract_integers(column_label=l0)
      elif (t0 in "R"): # general real column
        labels = [l0]
        group = self.extract_reals(column_label=l0)
      elif (t0 in "A"): # Hendrickson-Lattman coefficients
        if (all_column_types[i_column:i_column+4] == "AAAA"):
          if (len(remaining_types) >= 8
              and remaining_types[4:8] == "AAAA"
              and are_anomalous_labels("+",
                    all_column_labels[i_column:i_column+4])
              and are_anomalous_labels("-",
                    all_column_labels[i_column+4:i_column+8])):
            labels = all_column_labels[i_column:i_column+8]
            i_column += 7
            group = self.extract_hendrickson_lattman_anomalous(*labels)
          else:
            labels = all_column_labels[i_column:i_column+4]
            if (    are_anomalous_labels("+",
                      all_column_labels[i_column:i_column+2])
                and are_anomalous_labels("-",
                      all_column_labels[i_column+2:i_column+4])):
              group = self.extract_hendrickson_lattman_anomalous_ab_only(
                *labels)
            else:
              group = self.extract_hendrickson_lattman(*labels)
            i_column += 3
        elif (all_column_types[i_column:i_column+2] == "AA"):
          labels = all_column_labels[i_column:i_column+2]
          i_column += 1
          group = self.extract_hendrickson_lattman_ab_only(*labels)
        else:
          raise RuntimeError(
            'Invalid MTZ column combination'
            ' (incomplete Hendrickson-Lattman array),'
            + ' column labels: ' + ", ".join(['"%s"' % l
              for l in all_column_labels[i_column:i_column+4]])
            + ' column types: ' + ", ".join(['"%s"' % t
              for t in all_column_types[i_column:i_column+4]]))
      elif (remaining_types[:4] == "FQDQ"):
        labels = all_column_labels[i_column:i_column+4]
        i_column += 3
        group = self.extract_delta_anomalous(*labels)
        observation_type = xray.observation_types.reconstructed_amplitude()
      elif (t0 in "JFED"):
        # "J": "intensity"
        # "F": "amplitude"
        # "E": "normalized amplitude"
        # "D": "anomalous difference"
        # "Q": "standard deviation"
        # "P": "phase angle in degrees"
        labels = [l0]
        if (    i_column+1 < len(all_column_types)
            and all_column_types[i_column+1] in "QP"):
          labels = all_column_labels[i_column:i_column+2]
          i_column += 1
          if (all_column_types[i_column] == "Q"):
            group = self.extract_observations(*labels)
          else:
            if (t0 == "J"):
              raise RuntimeError(
                'Invalid MTZ column combination (intensity + phase angle),'
                + ' column labels: "%s" + "%s"' % tuple(labels))
            group = self.extract_complex(*labels)
        else:
          group = self.extract_reals(column_label=l0)
      elif (t0 in "GK"):
        # "G": "F(+) or F(-)"
        # "L": "standard deviation"
        # "P": "phase angle in degrees"
        # "K": "I(+) or I(-)"
        # "M": "standard deviation"
        perm = None
        if (remaining_types[:4] in ("GLGL", "GPGP", "KMKM")):
          perm = [0,1,2,3]
        elif (remaining_types[:4] in ("GGLL", "GGPP", "KKMM")):
          perm = [0,2,1,3]
        elif (remaining_types[:2] in ("GG", "KK")):
          perm = [0,1]
        else:
          raise RuntimeError('Invalid MTZ column combination:'
            + ' incomplete anomalous data for column label: "%s"' % l0)
        labels = [all_column_labels[i_column+i] for i in perm]
        i_column += len(perm)-1
        if (len(perm) == 2):
          group = self.extract_reals_anomalous(*labels)
        elif ("P" in remaining_types[:4]):
          group = self.extract_complex_anomalous(*labels)
        else:
          group = self.extract_observations_anomalous(*labels)
      else:
        labels = [l0]
        group = self.extract_reals(l0)
      groups.append(column_group(
        crystal_symmetry_from_file=crystal_symmetry_from_file,
        crystal_symmetry=crystal_symmetry,
        base_array_info=base_array_info,
        primary_column_type=t0,
        labels=labels,
        group=group,
        observation_type=observation_type))
    return groups

def column_group(
      crystal_symmetry_from_file,
      crystal_symmetry,
      base_array_info,
      primary_column_type,
      labels,
      group,
      observation_type):
  assert group.data.size() == group.indices.size()
  sigmas = getattr(group, "sigmas", None)
  if (sigmas is not None): assert sigmas.size() == group.indices.size()
  if (group.anomalous_flag):
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=group.indices,
      anomalous_flag=True)
    if (miller_set.n_bijvoet_pairs() == 0):
      # account for non-sensical files generated via
      # ccp4i "import merged data" tab with default parameters
      miller_set = miller.set(
        crystal_symmetry=crystal_symmetry,
        indices=group.indices,
        anomalous_flag=False)
  else:
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=group.indices).auto_anomalous(min_fraction_bijvoet_pairs=2/3.)
  result = (miller.array(
    miller_set=miller_set,
    data=group.data,
    sigmas=sigmas)
    .set_info(base_array_info.customized_copy(
      labels=labels,
      crystal_symmetry_from_file=crystal_symmetry_from_file)))
  if (observation_type is not None):
    result.set_observation_type(observation_type)
  elif (not result.is_complex_array()):
    if   (primary_column_type in "FG"):
      result.set_observation_type_xray_amplitude()
    elif (primary_column_type in "JK"):
      result.set_observation_type_xray_intensity()
  return result

def mend_non_conforming_anomalous_column_types(all_types, all_labels):

  replacements = {
    "JQJQ": "KMKM",
    "FQFQ": "GLGL",
    "JJQQ": "KKMM",
    "FFQQ": "GGLL",
    "FPFP": "GPGP",
    "FFPP": "GGPP",
  }

  def are_anomalous_labels():
    if (group_types[1] == "Q"):
      permutation = ((0,1),(2,3))
    else:
      permutation = ((0,2),(1,3))
    for offsets,sign in zip(permutation, ("+", "-")):
      for offs in offsets:
        if (not is_anomalous_label(sign, all_labels[i_group_start + offs])):
          return False
    return True

  def find_group():
    for group_types in replacements.keys():
      i_group_start = all_types_x.find(group_types)
      if (i_group_start >= 0):
        return group_types, i_group_start
    return None, None

  all_types = "".join(all_types)
  if (0): # for debugging only
    all_types = all_types.replace("GLGL", "FQFQ")
    all_types = all_types.replace("KMKM", "JQJQ")
  all_types_x = all_types
  while 1:
    group_types, i_group_start = find_group()
    if (group_types is None):
      break
    if (are_anomalous_labels()):
      replacement = replacements[group_types]
      all_types = overwrite_at(all_types, i_group_start, replacement)
    else:
      replacement = "X"
    all_types_x = overwrite_at(all_types_x, i_group_start, replacement)
  return all_types

class _(boost.python.injector, ext.crystal):

  def crystal_symmetry(self):
    return cctbx.crystal.symmetry(
      unit_cell=self.unit_cell(),
      space_group_info=self.mtz_object().space_group_info())

  def miller_set(self, anomalous_flag=None):
    return miller.set(
      crystal_symmetry=self.crystal_symmetry(),
      indices=self.mtz_object().extract_miller_indices(),
      anomalous_flag=anomalous_flag)

class _(boost.python.injector, ext.dataset):

  def column_labels(self):
    return [column.label() for column in self.columns()]

  def column_types(self):
    return [column.type() for column in self.columns()]

  def initialize_hkl_columns(self):
    if (self.mtz_object().n_columns() == 0):
      for label in "HKL":
        self.add_column(label=label, type="H")

  def _add_observations(self, data_label, sigmas_label, column_types,
                              indices, data, sigmas):
    mtz_reflection_indices = self.add_column(
      label=data_label,
      type=column_types[0]).set_reals(
        miller_indices=indices,
        data=data)
    if (sigmas is not None):
      self.add_column(
        label=sigmas_label,
        type=column_types[1]).set_reals(
          mtz_reflection_indices=mtz_reflection_indices,
          data=sigmas)

  def _add_complex(self, amplitudes_label, phases_label, column_types,
                         indices, data):
    mtz_reflection_indices = self.add_column(
      label=amplitudes_label,
      type=column_types[0]).set_reals(
        miller_indices=indices,
        data=flex.abs(data))
    self.add_column(
      label=phases_label,
      type=column_types[1]).set_reals(
        mtz_reflection_indices=mtz_reflection_indices,
        data=flex.arg(data, True))

  def add_miller_array(self,
        miller_array,
        column_root_label,
        column_types=None,
        label_decorator=None):
    assert column_types is None or isinstance(column_types, str)
    if (label_decorator is None):
      label_decorator = globals()["label_decorator"]()
    default_col_types = default_column_types(miller_array=miller_array)
    if (default_col_types is None):
      raise RuntimeError(
        "Conversion of given type of miller_array to MTZ format"
        " is not supported.")
    if (column_types is None):
      column_types = default_col_types
    elif (len(column_types) != len(default_col_types)):
      raise RuntimeError(
        "Invalid MTZ column_types for the given miller_array.")
    self.initialize_hkl_columns()
    if (not miller_array.anomalous_flag()):
      if (default_col_types in ["FQ", "JQ"]):
        self._add_observations(
          data_label=column_root_label,
          sigmas_label=label_decorator.sigmas(column_root_label),
          column_types=column_types,
          indices=miller_array.indices(),
          data=miller_array.data(),
          sigmas=miller_array.sigmas())
      elif (default_col_types == "FP"):
        self._add_complex(
          amplitudes_label=column_root_label,
          phases_label=label_decorator.phases(column_root_label),
          column_types=column_types,
          indices=miller_array.indices(),
          data=miller_array.data())
      elif (default_col_types in ["F", "J"]):
        self.add_column(
          label=column_root_label,
          type=column_types).set_reals(
            miller_indices=miller_array.indices(),
            data=miller_array.data())
      elif (default_col_types == "I"):
        self.add_column(
          label=column_root_label,
          type=column_types).set_reals(
            miller_indices=miller_array.indices(),
            data=miller_array.data().as_double())
      elif (default_col_types == "AAAA"):
        mtz_reflection_indices = self.add_column(
          label=label_decorator.hendrickson_lattman(column_root_label, 0),
          type=column_types[0]).set_reals(
            miller_indices=miller_array.indices(),
            data=miller_array.data().slice(0))
        for i in xrange(1,4):
          self.add_column(
            label=label_decorator.hendrickson_lattman(column_root_label, i),
            type=column_types[i]).set_reals(
              mtz_reflection_indices=mtz_reflection_indices,
              data=miller_array.data().slice(i))
      else:
        raise RuntimeError("Fatal programming error.")
    else:
      asu, matches = miller_array.match_bijvoet_mates()
      for anomalous_sign in ("+","-"):
        sel = matches.pairs_hemisphere_selection(anomalous_sign)
        sel.extend(matches.singles_hemisphere_selection(anomalous_sign))
        if (anomalous_sign == "+"):
          indices = asu.indices().select(sel)
        else:
          indices = -asu.indices().select(sel)
        data = asu.data().select(sel)
        if (default_col_types in ["GL", "KM"]):
          self._add_observations(
            data_label=label_decorator.anomalous(
              column_root_label, anomalous_sign),
            sigmas_label=label_decorator.sigmas(
              column_root_label, anomalous_sign),
            column_types=column_types,
            indices=indices,
            data=data,
            sigmas=asu.sigmas().select(sel))
        elif (default_col_types == "GP"):
          self._add_complex(
            amplitudes_label=label_decorator.anomalous(
              column_root_label, anomalous_sign),
            phases_label=label_decorator.phases(
              column_root_label, anomalous_sign),
            column_types=column_types,
            indices=indices,
            data=data)
        elif (default_col_types in ["G", "K"]):
          self.add_column(
            label=label_decorator.anomalous(column_root_label, anomalous_sign),
            type=column_types).set_reals(
              miller_indices=indices,
              data=data)
        elif (default_col_types == "I"):
          self.add_column(
            label=label_decorator.anomalous(column_root_label, anomalous_sign),
            type=column_types).set_reals(
              miller_indices=indices,
              data=data.as_double())
        elif (default_col_types == "AAAA"):
          mtz_reflection_indices = self.add_column(
            label=label_decorator.hendrickson_lattman(
              column_root_label, 0, anomalous_sign),
            type=column_types[0]).set_reals(
              miller_indices=indices,
              data=data.slice(0))
          for i in xrange(1,4):
            self.add_column(
              label=label_decorator.hendrickson_lattman(
                column_root_label, i, anomalous_sign),
              type=column_types[i]).set_reals(
                mtz_reflection_indices=mtz_reflection_indices,
                data=data.slice(i))
        else:
          raise RuntimeError("Fatal programming error.")
    return self

def miller_array_as_mtz_dataset(self,
      column_root_label,
      column_types=None,
      label_decorator=None,
      title=None,
      crystal_name="crystal",
      project_name="project",
      dataset_name="dataset",
      wavelength=1.0):
  if (title is None):
    title = str(self.info())
  if (title is None):
    title = "cctbx.miller.array"
  unit_cell = self.unit_cell()
  if (unit_cell is None):
    unit_cell = uctbx.unit_cell((1,1,1,90,90,90))
  space_group_info = self.space_group_info()
  if (space_group_info is None):
    space_group_info = sgtbx.space_group_info(symbol="P 1")
  mtz_object = object() \
    .set_title(title=title) \
    .set_space_group_info(space_group_info=space_group_info)
  mtz_object.set_hkl_base(unit_cell=unit_cell)
  return mtz_object.add_crystal(
    name=crystal_name,
    project_name=project_name,
    unit_cell=unit_cell).add_dataset(
      name=dataset_name,
      wavelength=wavelength).add_miller_array(
        miller_array=self,
        column_root_label=column_root_label,
        column_types=column_types,
        label_decorator=label_decorator)

class _(boost.python.injector, ext.batch):

  def show(self, out=None):
    if (out is None): out = sys.stdout
    print >> out, "batch number:", self.num()
    print >> out, "batch title:", self.title().rstrip()
    print >> out, "names of the three axes:", list(self.gonlab())
    print >> out, "type of orientation block:", self.iortyp()
    print >> out, "refinement flags for cell:", list(self.lbcell())
    print >> out, "number of phixyz used (0, 1, or 2):", self.misflg()
    print >> out, "reciprocal axis closest to rotation axis:", self.jumpax()
    print >> out, "crystal number:", self.ncryst()
    print >> out, "mosaicity model: 0 = isotropic, 1 = anisotropic:", \
      self.lcrflg()
    print >> out, "type of data: 2D (1), 3D (2), or Laue (3):", self.ldtype()
    print >> out, "goniostat scan axis number:", self.jsaxs()
    print >> out, "number of batch scales & Bfactors (0 if unset):", \
      self.nbscal()
    print >> out, "number of goniostat axes:", self.ngonax()
    print >> out, "flag for type of beam info:", self.lbmflg()
    print >> out, "  0: for alambd, delamb; 1: also delcor, divhd, divvd"
    print >> out, "number of detectors (current maximum 2):", self.ndet()
    print >> out, "dataset id:", self.nbsetid()
    print >> out, "cell dimensions:", list(self.cell())
    print >> out, "orientation matrix U:", list(self.umat())
    print >> out, "  in Fortranic order, i.e. U(1,1), U(2,1) ..."
    print >> out, "missetting angles at beginning and end of oscillation:", \
      list(self.phixyz())
    print >> out, "mosaicity:", list(self.crydat())
    print >> out, "datum values of goniostat axes:", list(self.datum())
    print >> out, "start of phi relative to datum:", self.phistt()
    print >> out, "end of phi relative to datum:", self.phiend()
    print >> out, "rotation axis in lab frame:", list(self.scanax())
    print >> out, "start time:", self.time1()
    print >> out, "stop time:", self.time2()
    print >> out, "batch scale:", self.bscale()
    print >> out, "batch temperature factor:", self.bbfac()
    print >> out, "sd bscale:", self.sdbscale()
    print >> out, "sd bbfac:", self.sdbfac()
    print >> out, "phi range:", self.phirange()
    print >> out, 'vectors ("Cambridge" laboratory axes) defining' \
      ' ngonax goniostat axes:'
    print >> out, "  vector 1:", list(self.e1())
    print >> out, "  vector 2:", list(self.e2())
    print >> out, "  vector 3:", list(self.e3())
    print >> out, "idealised source vector:", list(self.source())
    print >> out, "source vector:", list(self.so())
    print >> out, "wavelength (A):", self.alambd()
    print >> out, "dispersion (deltalambda / lambda):", self.delamb()
    print >> out, "correlated component:", self.delcor()
    print >> out, "horizontal beam divergence:", self.divhd()
    print >> out, "vertical beam divergence:", self.divvd()
    print >> out, "xtal to detector distance:", list(self.dx())
    print >> out, "detector tilt angle:", list(self.theta())
    print >> out, "min & max values of detector coords (pixels):", \
      list(self.detlm())
