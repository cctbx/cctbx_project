# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import iotbx.gui_tools
from iotbx import file_reader
from libtbx.utils import Sorry
import os
import six
from six.moves import range
from libtbx.phil import parse
import textwrap


def space_group_as_str(space_group):
  from cctbx import sgtbx
  if space_group is None :
    return ""
  elif isinstance(space_group, sgtbx.space_group_info):
    return str(space_group)
  else :
    sg_info = sgtbx.space_group_info(group=space_group)
    return str(sg_info)

def unit_cell_as_str(unit_cell, separator=" "):
  assert isinstance(separator, str) and separator != ""
  if unit_cell is not None :
    format = separator.join([ "%g" for i in range(6) ])
    return format % unit_cell.parameters()
  else :
    return ""

class reflections_handler(iotbx.gui_tools.manager):
  file_type = "hkl"
  file_type_label = "Reflections"
  def __init__(self,
                allowed_param_names=None,
                allowed_multiple_params=None,
                debug=False,
                minimum_data_score=4,
                prefer_amplitudes=False):
    iotbx.gui_tools.manager.__init__(self,
      allowed_param_names=allowed_param_names,
      allowed_multiple_params=allowed_multiple_params,
      debug=debug)
    self.minimum_data_score = minimum_data_score
    self.prefer_amplitudes = prefer_amplitudes

  def get_miller_array(self, labels, file_name=None, file_param_name=None):
    hkl_file = self.get_file(file_name=file_name,
      file_param_name=file_param_name)
    if (hkl_file is None):
      return None
    for array in hkl_file.file_server.miller_arrays :
      array_label_string = array.info().label_string()
      array_labels = array.info().labels
      if ((array_label_string == labels) or
          (array_labels == labels) or
          (",".join(array_labels) == labels)):
        return array
    return None

  def get_wavelength(self, *args, **kwds):
    array = self.get_miller_array(*args, **kwds)
    if (array is not None):
      info = array.info()
      if (info is not None):
        return info.wavelength
    return None

  def check_symmetry(self, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    hkl_server = hkl_file.file_server
    ma = hkl_server.miller_arrays[0]
    if (ma.space_group() is None) or (ma.unit_cell() is None):
      raise Sorry("Incomplete symmetry for this reflections file.  Please "+
        "use a format that includes both space group and unit cell, such as "+
        "an MTZ or merged Scalepack file.")
    return True

  def get_file_type_label(self, file_name=None, input_file=None):
    if (input_file is not None):
      return input_file.file_object.file_type()

  def get_all_labels(self, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    if hkl_file is None :
      labels = []
      for (file_name, hkl_file) in self.input_files():
        miller_arrays = hkl_file.file_server.miller_arrays
        labels.extend([ a.info().label_string() for a in miller_arrays ])
    else :
      miller_arrays = hkl_file.file_server.miller_arrays
      labels = [ a.info().label_string() for a in miller_arrays ]
    return labels

  def get_rfree_labels(self, *args, **kwds):
    prefer_neutron = kwds.pop("neutron", None)
    hkl_file = self.get_file(*args, **kwds)
    labels = []
    if hkl_file is not None :
      hkl_server = hkl_file.file_server
      arrays_and_flags = hkl_server.get_r_free_flags(
        file_name                = None,
        label                    = None,
        test_flag_value          = None,
        disable_suitability_test = False,
        parameter_scope          = "",
        return_all_valid_arrays  = True,
        minimum_score            = 1)
      if (prefer_neutron is not None):
        xray_arrays = []
        neutron_arrays = []
        other_arrays = []
        for array_and_flag_value in arrays_and_flags :
          array = array_and_flag_value[0]
          label = array.info().label_string().lower()
          if ("neutron" in label):
            neutron_arrays.append(array_and_flag_value)
          elif ("xray" in label):
            xray_arrays.append(array_and_flag_value)
          else :
            other_arrays.append(array_and_flag_value)
        if (prefer_neutron) and (len(neutron_arrays) > 0):
          arrays_and_flags = neutron_arrays + other_arrays
        elif (not prefer_neutron) and (len(xray_arrays) > 0):
          arrays_and_flags = xray_arrays + other_arrays
      if len(arrays_and_flags) > 0 :
        labels = [ a.info().label_string() for a, fl in arrays_and_flags ]
    return labels

  def has_rfree(self, *args, **kwds):
    return len(self.get_rfree_labels(*args, **kwds)) > 0

  def has_neutron_rfree(self, *args, **kwds):
    kwds['neutron'] = True
    labels = self.get_rfree_labels(*args, **kwds)
    return ([ "neutron" in l.lower() for l in labels ].count(True) > 0)

  def get_rfree_flag_value(self, array_name, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    flag = None
    if hkl_file is not None :
      hkl_server = hkl_file.file_server
      arrays_and_flags = hkl_server.get_r_free_flags(
        file_name                = None,
        label                    = None,
        test_flag_value          = None,
        disable_suitability_test = False,
        parameter_scope          = "",
        return_all_valid_arrays  = True,
        minimum_score            = 1)
      for miller_array, array_flag in arrays_and_flags :
        label = miller_array.info().label_string()
        if label == array_name :
          flag = array_flag
          break
    return flag

  def get_experimental_phase_labels(self, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    labels = []
    if hkl_file is not None :
      hkl_server = hkl_file.file_server
      miller_arrays = hkl_server.get_experimental_phases(
        file_name               = hkl_file.file_name,
        labels                  = None,
        ignore_all_zeros        = True,
        parameter_scope         = "",
        return_all_valid_arrays = True,
        minimum_score           = 1)
      labels = [ array.info().label_string() for array in miller_arrays ]
    return labels

  def has_phases(self, *args, **kwds):
    return len(self.get_experimental_phase_labels(*args, **kwds)) > 0

  def get_phase_arrays(self, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    labels = []
    if hkl_file is not None :
      hkl_server = hkl_file.file_server
      miller_arrays = hkl_server.get_phases_deg(
        file_name=hkl_file.file_name,
        labels=None,
        convert_to_phases_if_necessary=False,
        original_phase_units=None,
        parameter_scope=None,
        parameter_name=None,
        return_all_valid_arrays=True,
        minimum_score=1)
      return miller_arrays
    return []

  def get_phase_deg_labels(self, *args, **kwds):
    miller_arrays = self.get_phase_arrays(*args, **kwds)
    labels = []
    for array in miller_arrays :
      labels_str = array.info().label_string()
      if labels_str.startswith("FOM") : continue
      labels.append(labels_str)
    return labels

  def get_phase_column_labels(self, *args, **kwds):
    labels = []
    miller_arrays = self.get_phase_arrays(*args, **kwds)
    for array in miller_arrays :
      for label in array.info().labels :
        if label.upper().startswith("PH"):
          labels.append(label)
          break
    return labels

  def get_data_arrays(self, *args, **kwds):
    prefer_neutron = kwds.pop("neutron", None)
    hkl_file = self.get_file(*args, **kwds)
    if hkl_file is not None :
      hkl_server = hkl_file.file_server
      miller_arrays = hkl_server.get_xray_data(
        file_name               = None,
        labels                  = None,
        ignore_all_zeros        = False,
        parameter_scope         = "",
        return_all_valid_arrays = True,
        minimum_score           = self.minimum_data_score)
      if self.prefer_amplitudes :
        f_arrays = []
        other_arrays = []
        for array in miller_arrays :
          if array.is_xray_amplitude_array():
            f_arrays.append(array)
          else :
            other_arrays.append(array)
        miller_arrays = f_arrays + other_arrays
      elif (prefer_neutron is not None):
        xray_arrays = []
        neutron_arrays = []
        other_arrays = []
        for array in miller_arrays :
          label = array.info().label_string().lower()
          if ("neutron" in label):
            neutron_arrays.append(array)
          elif ("xray" in label):
            xray_arrays.append(array)
          else :
            other_arrays.append(array)
        if (prefer_neutron) and (len(neutron_arrays) > 0):
          miller_arrays = neutron_arrays + other_arrays
        elif (not prefer_neutron) and (len(xray_arrays) > 0):
          miller_arrays = xray_arrays + other_arrays
      return miller_arrays
    return []

  def get_data_labels(self, *args, **kwds):
    return extract_labels(self.get_data_arrays(*args, **kwds))

  def get_amplitude_arrays(self, *args, **kwds):
    allow_conversion = kwds.pop('allow_conversion', False)
    hkl_file = self.get_file(*args, **kwds)
    if hkl_file is not None :
      hkl_server = hkl_file.file_server
      miller_arrays = hkl_server.get_amplitudes(
        file_name               = None,
        labels                  = None,
        convert_to_amplitudes_if_necessary = allow_conversion,
        parameter_scope         = "",
        parameter_name          = "",
        return_all_valid_arrays = True,
        minimum_score           = 1,
        strict                  = True)
      return miller_arrays
    return []

  def get_amplitude_labels(self, *args, **kwds):
    return extract_labels(self.get_amplitude_arrays(*args, **kwds))

  def get_amplitude_column_labels(self, *args, **kwds):
    miller_arrays = self.get_amplitude_arrays(*args, **kwds)
    labels = []
    for array in miller_arrays :
      # XXX what about anomalous data?
      labels.extend(array.info().labels[0])
    return labels

  def get_intensity_arrays(self, *args, **kwds):
    miller_arrays = self.get_data_arrays(*args, **kwds)
    i_arrays = []
    for array in miller_arrays :
      if array.is_xray_intensity_array():
        i_arrays.append(array)
    return i_arrays

  def get_intensity_labels(self, *args, **kwds):
    return extract_labels(self.get_intensity_arrays(*args, **kwds))

  # space-separated, column labels only
  def get_data_labels_for_wizard(self, *args, **kwds):
    miller_arrays = self.get_data_arrays(*args, **kwds)
    labels = [ " ".join(array.info().labels) for array in miller_arrays ]
    return labels

  def has_data(self, *args, **kwds):
    return len(self.get_data_labels(*args, **kwds)) > 0

  def has_neutron_data(self, *args, **kwds):
    kwds['neutron'] = True
    labels = self.get_data_labels(*args, **kwds)
    return ([ "neutron" in l.lower() for l in labels ].count(True) > 0)

  def get_anomalous_data_labels(self, *args, **kwds):
    kwds = dict(kwds) # XXX gross...
    allow_dano = kwds.pop("allow_reconstructed_amplitudes", True)
    hkl_file = self.get_file(*args, **kwds)
    labels = []
    if hkl_file is not None :
      hkl_server = hkl_file.file_server
      miller_arrays = hkl_server.get_xray_data(
        file_name               = None,
        labels                  = None,
        ignore_all_zeros        = False,
        parameter_scope         = "",
        return_all_valid_arrays = True,
        minimum_score           = self.minimum_data_score)
      for array in miller_arrays :
        if (array.anomalous_flag()):
          if ((array.is_xray_reconstructed_amplitude_array()) and
              (not allow_dano)):
            continue
          labels.append(array.info().label_string())
    return labels

  def has_anomalous_data(self, *args, **kwds):
    return (len(self.get_anomalous_data_labels(*args, **kwds)) > 0)

  def get_phaser_map_fc_labels(self, *args, **kwds):
    labels = self.get_fmodel_labels(*args, **kwds)
    first_column_only = kwds.pop('first_column_only', False)
    hkl_file = self.get_file(*args, **kwds)
    if (hkl_file is not None):
      for miller_array in hkl_file.file_server.miller_arrays :
        if (miller_array.is_complex_array()):
          labels_str = miller_array.info().label_string()
          if labels_str.startswith("FWT"):
            labels.append(miller_array.info().labels[0])
    return labels

  def get_fmodel_labels(self, *args, **kwds):
    first_column_only = kwds.pop('first_column_only', False)
    hkl_file = self.get_file(*args, **kwds)
    labels_list = []
    if (hkl_file is not None):
      for miller_array in hkl_file.file_server.miller_arrays :
        if (miller_array.is_complex_array()):
          labels = miller_array.info().label_string()
          if (labels.startswith("F-model") or
              labels.upper().startswith("FMODEL") or
              labels.upper().startswith("FC")):
            if (first_column_only):
              labels_list.append(miller_array.info().labels[0])
            else :
              labels_list.append(labels)
    return labels_list

  def get_map_coeff_labels(self, *args, **kwds):
    # FIXME this is just gross...
    kwds_basic = {}
    kwds_maps = {}
    for kwd in dict(kwds):
      if (kwd in ["file_name", "file_param_name"]):
        kwds_basic[kwd] = kwds[kwd]
      else :
        kwds_maps[kwd] = kwds[kwd]
    hkl_file = self.get_file(*args, **kwds_basic)
    if hkl_file is not None :
      return get_map_coeff_labels(hkl_file.file_server, **kwds_maps)
    return []

  def get_map_coeff_labels_for_build(self, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    if hkl_file is not None :
      return get_map_coeffs_for_build(hkl_file.file_server)
    return []

  def get_map_coeff_labels_for_fft(self, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    labels_list = []
    if hkl_file is not None :
      all_labels = get_map_coeff_labels(hkl_file.file_server,
        keep_array_labels=True)
      for labels in all_labels :
        if isinstance(labels, str):
          labels_list.append(labels)
        else :
          labels_list.append(" ".join(labels))
    return labels_list

  def get_two_fofc_map_labels(self, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    labels_list = []
    if (hkl_file is not None):
      for array in hkl_file.file_server.miller_arrays :
        labels = array.info().label_string()
        if (labels.startswith("FWT") or labels.startswith("2FOFC")):
          labels_list.append(labels)
    return labels_list

  def get_fofc_map_labels(self, *args, **kwds):
    hkl_file = self.get_file(*args, **kwds)
    labels_list = []
    if (hkl_file is not None):
      for array in hkl_file.file_server.miller_arrays :
        labels = array.info().label_string()
        if (labels.startswith("DELFWT") or labels.startswith("FOFC")):
          labels_list.append(labels)
    return labels_list

  def get_amplitude_column_labels(self, *args, **kwds):
    miller_arrays = self.get_amplitude_arrays(*args, **kwds)
    labels = []
    for array in miller_arrays :
      # XXX what about anomalous data?
      labels.extend(array.info().labels[0])
    return labels

  def d_max_min(self, file_name=None, file_param_name=None,
      array_name=None, array_names=None):
    from iotbx.reflection_file_editor import get_best_resolution
    miller_arrays = []
    if (file_name is None) and (file_param_name is None):
      for phil_name, file_name in six.iteritems(self._param_files):
        input_file = self.get_file(file_name)
        if (input_file is not None):
          miller_arrays.extend(input_file.file_server.miller_arrays)
    else :
      hkl_file = self.get_file(file_name, file_param_name)
      if (hkl_file is None):
        return (None, None)
      for miller_array in hkl_file.file_server.miller_arrays :
        label = miller_array.info().label_string()
        if (array_name is not None):
          if (label == array_name):
            miller_arrays = [miller_array]
            break
        elif (array_names is not None):
          if (label in array_names):
            miller_arrays.append(miller_array)
        else :
          miller_arrays.append(miller_array)
    (d_max, d_min) = get_best_resolution(miller_arrays)
    return (d_max, d_min)

  def get_resolution_range(self, *args, **kwds):
    (d_max, d_min) = self.d_max_min(*args, **kwds)
    if (d_max is None):
      return ""
    else :
      return "(%.3f - %.3f)" % (d_max, d_min)

  def get_resolution_limits(self, *args, **kwds):
    (d_max, d_min) = self.d_max_min(*args, **kwds)
    (d_max_str, d_min_str) = ("", "")
    if (d_max is not None):
      d_max_str = "(%.3f)" % d_max
    if (d_min is not None):
      d_min_str = "(%.3f)" % d_min
    return (d_max_str, d_min_str)

  def space_group(self, *args, **kwds):
    symm = self.crystal_symmetry(*args, **kwds)
    if symm is not None :
      return symm.space_group()
    return None

  def space_group_as_str(self, *args, **kwds):
    space_group = self.space_group(*args, **kwds)
    return space_group_as_str(space_group)

  def unit_cell(self, *args, **kwds):
    symm = self.crystal_symmetry(*args, **kwds)
    if symm is not None :
      return symm.unit_cell()
    return None

  def unit_cell_as_str(self, file_name=None, file_param_name=None,
      separator=" "):
    unit_cell = self.unit_cell(file_name, file_param_name)
    return unit_cell_as_str(unit_cell, separator)

  def crystal_symmetry(self, file_name=None, file_param_name=None):
    final_symm = None
    if file_name is None and file_param_name is None :
      for param_name, file_name in six.iteritems(self._param_files):
        input_file = self.get_file(file_name)
        miller_arrays = input_file.file_server.miller_arrays
        for array in miller_arrays :
          symm = array.crystal_symmetry()
          if symm is not None :
            final_symm = symm
            break
    else :
      hkl_file = self.get_file(file_name, file_param_name)
      if hkl_file is not None :
        miller_arrays = hkl_file.file_server.miller_arrays
        if len(miller_arrays) > 0 :
          final_symm = miller_arrays[0].crystal_symmetry()
    return final_symm

  def check_symmetry_consistency(self):
    pass

def extract_labels(miller_arrays):
  return [ array.info().label_string() for array in miller_arrays ]

def get_fp_fpp_from_sasaki(guess_ha,wavelength):
  from cctbx.eltbx import sasaki
  from decimal import Decimal, ROUND_HALF_UP
  if wavelength is None or guess_ha is None :
    return None, None
  try:
    table = sasaki.table(guess_ha)
  except Exception as e :
    return None, None
  fp_fdp = table.at_angstrom(wavelength)
  f_prime=fp_fdp.fp()
  f_double_prime=fp_fdp.fdp()
  fp = Decimal(f_prime).quantize(Decimal('0.01'), ROUND_HALF_UP)
  fpp = Decimal(f_double_prime).quantize(Decimal('0.01'), ROUND_HALF_UP)
  return fp, fpp

def get_high_resolution(server):
  d_min = None #999.99
  for miller_array in server.miller_arrays :
    try :
      array_min = miller_array.d_min()
      if d_min is None or array_min < d_min :
        d_min = array_min
    except Exception :
      pass
  return d_min

def get_miller_array_symmetry(miller_array):
  from cctbx import sgtbx
  symm = miller_array.crystal_symmetry()
  if symm is None :
    return (None, None)
  unit_cell = symm.unit_cell()
  if (unit_cell is not None):
    uc = unit_cell_as_str(unit_cell)
  else :
    uc = None
  space_group = symm.space_group()
  if (space_group is not None):
    sg = str(sgtbx.space_group_info(group=space_group))
  else :
    sg = None
  return (sg, uc)

def get_array_description(miller_array):
  from iotbx import reflection_file_utils
  info = miller_array.info()
  labels = info.label_string()
  if labels in ["PHI", "PHIB", "PHIM", "PHIC"] :
    return "Phases"
  if labels in ["FOM", "FOMM"] :
    return "Weights"
  if ((miller_array.is_integer_array() or miller_array.is_bool_array()) and
      reflection_file_utils.looks_like_r_free_flags_info(info)):
    return "R-free flag"
  methods_and_meanings = [ ("is_complex_array", "Map coeffs"),
                           ("is_xray_amplitude_array", "Amplitude"),
                           ("is_xray_reconstructed_amplitude_array", "Amplitude"),
                           ("is_xray_intensity_array", "Intensity"),
                           ("is_hendrickson_lattman_array", "HL coeffs"),
                           ("is_bool_array", "Boolean"),
                           ("is_integer_array", "Integer"),
                           ("is_real_array", "Floating-point"),
                           ("is_string_array", "String")
                           ]
  for method, desc in methods_and_meanings :
    test = getattr(miller_array, method)
    if test():
      return desc
  # resort to the name of the python data type if not matching any of the above
  return type(miller_array.data()).__name__


def get_mtz_label_prefix(input_file=None, file_name=None):
  if input_file is None :
    input_file = file_reader.any_file(file_name)
  assert (input_file.file_type == "hkl" and
          input_file.file_object.file_type() == "ccp4_mtz")
  file_content = input_file.file_object.file_content()
  last_crystal = file_content.crystals()[-1]
  crystal = last_crystal.name()
  dataset = last_crystal.datasets()[-1].name()
  return "/%s/%s" % (crystal, dataset)

def extract_map_coeffs(miller_arrays, f_lab, phi_lab, fom_lab):
  f_array = None
  phi_array = None
  fom_array = None
  #print f_lab, phi_lab, fom_lab
  for miller_array in miller_arrays :
    labels = miller_array.info().label_string()
    if (labels == f_lab):
      f_array = miller_array
    elif (labels == phi_lab):
      phi_array = miller_array
    elif (labels == fom_lab):
      fom_array = miller_array
  return (f_array, phi_array, fom_array)

def map_coeffs_from_mtz_file(mtz_file, f_label="FP", phi_label="PHIM",
    fom_label="FOMM"):
  if not os.path.isfile(mtz_file):
    raise Sorry(
      "No map coefficients are available for conversion.")
  mtz_in = file_reader.any_file(mtz_file)
  mtz_in.assert_file_type("hkl")
  miller_arrays = mtz_in.file_server.miller_arrays
  (f_array, phi_array, fom_array) = extract_map_coeffs(miller_arrays,
    f_label, phi_label, fom_label)
  if (f_array is not None) and (f_array.is_complex_array()):
    map_coeffs = f_array
  else :
    if (f_array is None) or (phi_array is None):
      raise Sorry("One or more of the columns %s and %s was not found." %
        (f_label, phi_label))
    if fom_array is not None :
      weighted_f = f_array * fom_array
    else :
      weighted_f = f_array
    map_coeffs = weighted_f.phase_transfer(phi_array, deg=True)
  return map_coeffs

def extract_phenix_refine_map_coeffs(mtz_file, limit_arrays=None):
  assert (limit_arrays is None) or (isinstance(limit_arrays, list))
  if not os.path.isfile(mtz_file):
    raise Sorry("No map coefficients are available for conversion.")
  mtz_in = file_reader.any_file(mtz_file)
  mtz_in.assert_file_type("hkl")
  miller_arrays = mtz_in.file_server.miller_arrays
  assert len(miller_arrays) > 0
  map_names = {"2FOFCWT" : "2mFo-DFc",
               "FOFCWT" : "mFo-DFc",
               "2FOFCWT_no_fill" : "2mFo-DFc_no_fill",
               "FOFCWT_no_fill" : "mFo-DFc_no_fill"}
  output_arrays = []
  for miller_array in miller_arrays :
    if miller_array.is_complex_array():
      labels = miller_array.info().label_string()
      if labels.startswith("F-model"):
        continue
      if (limit_arrays is not None) and (not labels in limit_arrays):
        continue
      f_label = miller_array.info().labels[0]
      map_name = map_names.get(f_label)
      if map_name is None :
        map_name = f_label
      output_arrays.append((miller_array, map_name))
  return output_arrays

def get_map_coeff_labels(server,
    build_only=False,
    include_fom=True,
    exclude_anomalous=False,
    exclude_fmodel=False,
    keep_array_labels=False):
  all_labels = []
  phi_labels = []
  fom_labels = []
  for miller_array in server.miller_arrays :
    label = miller_array.info().label_string()
    if label.startswith("FOM") and include_fom :
      fom_labels.append(label)
  phase_arrays = server.get_phases_deg(None, None, False, None, None, None,
                                       True, 3)
  for miller_array in phase_arrays :
    labels = miller_array.info().label_string()
    if miller_array.is_hendrickson_lattman_array():
      continue
    elif miller_array.is_complex_array():
      if ((labels.startswith("F-model") or labels.startswith("FMODEL") or
           labels.startswith("FC")) and (exclude_fmodel)):
        continue
      elif (labels.upper().startswith("ANOM")) and (exclude_anomalous):
        continue
      if build_only :
        if (not labels[0:4] in ["FOFC", "DELF", "ANOM"]):
          # list these first
          if (labels in ["FWT,PHWT", "2FOFCWT,PH2FOFCWT"]):
            all_labels.insert(0, labels)
          else :
            all_labels.append(labels)
      else :
        all_labels.append(labels)
    elif miller_array.info().labels[0].startswith("PHI"):
      phi_labels.append(labels)
  amp_arrays = server.get_amplitudes(
    file_name               = None,
    labels                  = None,
    convert_to_amplitudes_if_necessary = False,
    parameter_scope         = "",
    parameter_name          = "",
    return_all_valid_arrays = True,
    minimum_score           = 2,
    strict                  = True)
  if keep_array_labels :
    for miller_array in amp_arrays :
      data_label = miller_array.info().label_string()
      for phase_label in phi_labels :
        hybrid_label = [data_label, phase_label]
        if len(fom_labels) > 0 and include_fom :
          for fom in fom_labels :
            hybrid_label_ = hybrid_label + [ fom ]
            all_labels.append(hybrid_label_)
        else :
          all_labels.append(hybrid_label)
  else :
    for miller_array in amp_arrays :
      f_label = miller_array.info().labels[0]
      if f_label[0] == "F" and f_label != "FC" :
        for phase_label in phi_labels :
          hybrid_label = "%s,%s" % (f_label, phase_label)
          if len(fom_labels) > 0 and include_fom :
            for fom in fom_labels :
              final_label = hybrid_label + ",%s" % fom
              all_labels.append(final_label)
          else :
            all_labels.append(hybrid_label)
  return all_labels

def get_map_coeffs_for_build(server):
  return get_map_coeff_labels(server, build_only=True)

def format_map_coeffs_for_resolve(f_label, phi_label, fom_label):
  return "FP=%s PHIB=%s FOM=%s" % (f_label, phi_label, fom_label)

def decode_resolve_map_coeffs(labels):
  fields = labels.strip().split()
  f_label = None
  phi_label = None
  fom_label = None
  for field in fields :
    resolve_label, array_label = field.split("=")
    if resolve_label == "FP" :
      f_label = array_label
    elif resolve_label == "PHIB" :
      phi_label = array_label
    elif resolve_label == "FOM" :
      fom_label = array_label
  return (f_label, phi_label, fom_label)



class ArrayInfo:
  """
  Extract information from a list of miller_array objects and format it for printing as a table
  To be called in a loop like:

  for i,array in enumerate(arrays):
    arrayinfo = ArrayInfo(array)
    info_fmt, headerstr, infostr = arrayinfo.get_selected_info_columns_from_phil()
    if i==0:
      print(headerstr)
    print(infostr)

  """
  def __init__(self, millarr, wrap_labels=0):
    from crys3d.hklviewer import display2
    from cctbx.array_family import flex
    from scitbx import graphics_utils
    from libtbx.math_utils import roundoff
    import math
    nan = float("nan")
    self.wrap_labels = wrap_labels
    if millarr.space_group() is None :
      self.spginf = "?"
    else:
      self.spginf = millarr.space_group_info().symbol_and_number()
    if millarr.unit_cell() is None:
      self.ucell = (nan, nan, nan, nan, nan, nan)
    else:
      self.ucell = millarr.unit_cell().parameters()
    self.ucellinf = "({:.6g}Å, {:.6g}Å, {:.6g}Å, {:.6g}°, {:.6g}°, {:.6g}°)".format(*self.ucell)
    data = millarr.deep_copy().data()
    self.maxdata = self.mindata = self.maxsigmas = self.minsigmas = nan
    self.minmaxdata = (nan, nan)
    self.minmaxsigs = (nan, nan)
    self.data_sigdata_max = nan
    self.data_sigdata = nan
    self.desc = ""
    self.arrsize = data.size()
    if not isinstance(data, flex.std_string):
      if isinstance(data, flex.hendrickson_lattman):
        data = graphics_utils.NoNansHL( data )
        # estimate minmax values of HL coefficients as a simple sum
        if self.arrsize:
          self.maxdata = max([e[0]+e[1]+e[2]+e[3] for e in data ])
          self.mindata = min([e[0]+e[1]+e[2]+e[3] for e in data ])
          self.arrsize = len([42 for e in millarr.data() if not math.isnan(e[0]+e[1]+e[2]+e[3])])
      elif isinstance(data, flex.vec3_double) or isinstance(data, flex.vec2_double):
        # XDS produces 2D or 3D arrays in some of its files
        if self.arrsize:
          self.maxdata = max([max(e) for e in data ])
          self.mindata = min([min(e) for e in data ])
      else:
        # Get a list of bools with True whenever there is a nan
        selection = ~ graphics_utils.IsNans( flex.abs( millarr.data()).as_double() )
        # count elements that are not nan values
        self.arrsize = millarr.data().select(selection).size()
        if (isinstance(data, flex.int)):
          data = flex.double([e for e in data if e!= display2.inanval])
        if millarr.is_complex_array():
          data = flex.abs(data)
        i=0
        while math.isnan(data[i]):
          i += 1 # go on until we find a data[i] that isn't NaN
        data = graphics_utils.NoNansArray( data, data[i] ) # assuming data[0] isn't NaN
        self.maxdata = flex.max( data )
        self.mindata = flex.min( data )
      if millarr.sigmas() is not None:
        data = millarr.sigmas().deep_copy()
        i=0
        while math.isnan(data[i]):
          i += 1 # go on until we find a data[i] that isn't NaN
        data = graphics_utils.NoNansArray( data, data[i] )
        self.maxsigmas = flex.max( data )
        self.minsigmas = flex.min( data )
        # Inspired by Diederichs doi:10.1107/S0907444910014836 I/SigI_asymptotic
        data_over_sigdata = millarr.data()/millarr.sigmas()
        self.data_sigdata = flex.sum(data_over_sigdata)/len(data_over_sigdata)
        self.data_sigdata_max = flex.max( data_over_sigdata)
      self.minmaxdata = (self.mindata, self.maxdata)
      self.minmaxsigs = (self.minsigmas, self.maxsigmas)
    self.labels = self.desc = self.wavelength = ""
    if millarr.info():
      self.labels = millarr.info().labels
      self.desc = get_array_description(millarr)
      self.wavelength = "{:.6g}".format(millarr.info().wavelength) if millarr.info().wavelength is not None else float("nan")
    self.span = "(?,?,?), (?,?,?)"
    self.dmin = nan
    self.dmax = nan
    if millarr.unit_cell() is not None:
      self.span = str(millarr.index_span().min()) + ", "+ str(millarr.index_span().max())
      self.dmin = millarr.d_max_min()[1]
      self.dmax = millarr.d_max_min()[0]
    self.dminmax = roundoff((self.dmin,self.dmax))
    self.issymunique = "?"
    self.isanomalous = "?"
    self.n_sys_abs = 0
    self.n_bijvoet = self.n_singletons = 0
    self.ano_mean_diff = nan
    self.ano_completeness = nan
    self.data_compl_infty = nan
    self.data_completeness = nan
    self.n_centric = nan
    # computations below done as in cctbx.miller.set.show_comprehensive_summary()
    if self.spginf != "?":
      self.issymunique = millarr.is_unique_set_under_symmetry()
      self.isanomalous = millarr.anomalous_flag()
      sys_absent_flags = millarr.sys_absent_flags().data()
      self.n_sys_abs = sys_absent_flags.count(True)
      if (self.n_sys_abs != 0):
        millarr = millarr.select(selection=~sys_absent_flags)
      self.n_centric = millarr.centric_flags().data().count(True)
    if not math.isnan(self.ucell[0]):
      if (self.spginf != "?"
          and millarr.indices().size() > 0
          and self.issymunique):
        millarr.setup_binner(n_bins=1)
        completeness_d_max_d_min = millarr.completeness(use_binning=True)
        binner = completeness_d_max_d_min.binner
        assert binner.counts_given()[0] == 0
        assert binner.counts_given()[2] == 0
        n_obs = binner.counts_given()[1]
        n_complete = binner.counts_complete()[1]
        if (n_complete != 0 and self.dmax != self.dmin):
            self.data_completeness = n_obs/n_complete
        n_complete += binner.counts_complete()[0]
        if (n_complete != 0):
            self.data_compl_infty = n_obs/n_complete
        if (self.isanomalous) and (millarr.is_xray_intensity_array() or
          millarr.is_xray_amplitude_array()):
            self.ano_completeness = millarr.anomalous_completeness()
    if (self.spginf != "?" and self.isanomalous and self.issymunique):
      asu, matches = millarr.match_bijvoet_mates()
      self.n_bijvoet = matches.pairs().size()
      self.n_singletons = matches.n_singles() - self.n_centric
      if millarr.is_real_array():
        self.ano_mean_diff = millarr.anomalous_signal()
    # break long label into list of shorter strings
    self.labelstr = ",".join(self.labels)
    if self.wrap_labels > 0:
      tlabels = textwrap.wrap(self.labelstr, width= self.wrap_labels)
      nlabl = len(tlabels)
      self.labelsformat = "{0[0]:>%d} " %(1+self.wrap_labels)
      if len(tlabels)>1:
        for i in range((len(tlabels)-1)):
          self.labelsformat += "\n{0[%d]:>%d} "%(i+1, self.wrap_labels+1)
      blanks = self.wrap_labels-5
    else:
      self.labelsformat = "{:>16} "
      if len(self.labelstr)>15:
        self.labelsformat = "{}\n                 "
      blanks= 10

    self.info_format_dict = {
      # the keys here must be verbatim copies of names of phil attributes in arrayinfo_phil_str below
      "labels":            (" %s" %self.caption_dict["labels"][0] + " "*blanks,   self.labelstr,         "{}",                   self.labelsformat),
      "description":       ("       %s      "%self.caption_dict["description"][0],self.desc,             "{}",                   "{:>16} "),
      "wavelength":        ("   %s   "%self.caption_dict["wavelength"][0],        self.wavelength,       "{}",                   "{:>8} "),
      "n_reflections":     ("  %s  " %self.caption_dict["n_reflections"][0],      self.arrsize,          "{}",                   "{:>8} "),
      "span":              (" "*15 + self.caption_dict["span"][0] + " "*14,       self.span,             "{}",                   "{:>32} "),
      "minmax_data":       ("     %s       " %self.caption_dict["minmax_data"][0],self.minmaxdata,       "{0[0]:.6}, {0[1]:.6}", "{0[0]:>11.5}, {0[1]:>11.5}"),
      "minmax_sigmas":     ("     %s     " %self.caption_dict["minmax_sigmas"][0],self.minmaxsigs,       "{0[0]:.6}, {0[1]:.6}", "{0[0]:>11.5}, {0[1]:>11.5}"),
      "data_sigdata":      (" %s" %self.caption_dict["data_sigdata"][0],          self.data_sigdata,     "{:.4g}",               "{:>9.4g} "),
      "data_sigdata_max":  ("%s" %self.caption_dict["data_sigdata_max"][0],       self.data_sigdata_max, "{:.4g}",               "{:>11.4g} "),
      "d_minmax":          ("  %s   " %self.caption_dict["d_minmax"][0],          self.dminmax,          "{0[0]:.6}, {0[1]:.6}", "{0[0]:>8.5}, {0[1]:>8.5}"),
      "unit_cell":         ("     %s      " %self.caption_dict["unit_cell"][0],   self.ucell,            "{0[0]:>7.5g},{0[1]:>7.5g},{0[2]:>7.5g},{0[3]:>7.5g},{0[4]:>7.5g},{0[5]:>7.5g}",
                                                                                              "{0[0]:>7.5g},{0[1]:>7.5g},{0[2]:>7.5g},{0[3]:>7.5g},{0[4]:>7.5g},{0[5]:>7.5g} "),
      "space_group":       ("   %s      " %self.caption_dict["space_group"][0],   self.spginf,           "{}",                   "{:>19} "),
      "n_centrics":        ("%s"%self.caption_dict["n_centrics"][0],              self.n_centric,        "{}",                   "{:>8} "),
      "n_sys_abs":         ("%s"%self.caption_dict["n_sys_abs"][0],               self.n_sys_abs,        "{}",                   "{:>9} "),
      "data_completeness": ("%s"%self.caption_dict["data_completeness"][0],       self.data_completeness,"{:.5g}",               "{:>10.5g} "),
      "data_compl_infty":  ("%s"%self.caption_dict["data_compl_infty"][0],        self.data_compl_infty, "{:.5g}",               "{:>9.5g} "),
      "is_anomalous":      ("%s"%self.caption_dict["is_anomalous"][0],            str(self.isanomalous), "{}",                   "{:>8} "),
      "is_symmetry_unique":("%s"%self.caption_dict["is_symmetry_unique"][0],      str(self.issymunique), "{}",                   "{:>8} "),
      "ano_completeness":  ("%s"%self.caption_dict["ano_completeness"][0],        self.ano_completeness, "{:.5g}",               "{:>11.5g} "),
      "ano_mean_diff":     ("%s"%self.caption_dict["ano_mean_diff"][0],           self.ano_mean_diff,    "{:.5g}",               "{:>8.5g} "),
      "n_bijvoet":         ("%s"%self.caption_dict["n_bijvoet"][0],               self.n_bijvoet,        "{}",                   "{:>8} "),
      "n_singletons":      ("%s"%self.caption_dict["n_singletons"][0],            self.n_singletons,     "{}",                   "{:>10} "),
    }

  # govern whether or not a property of the ArrayInfo should be returned by get_selected_info_columns_from_phil()
  arrayinfo_phil_str = """
wrap_labels = 15
  .type = int
  .short_caption = Wrap width for labels
  .help = Number of letters for wrapping long miller array labels. If less than 1 no wrapping is done
delimiter = "|"
  .type = str
  .short_caption = "column delimiter when printing table to standard output"
  .help = "column delimiter"
selected_info
  .help = "If values are set to True then tabulate respective properties of datasets in the reflection file."
{
    labels = True
      .type = bool
      .help = "Name of data array"
      .short_caption = "Labels"
    description = True
      .type = bool
      .help = "Type of data"
      .short_caption = "Type"
    wavelength = True
      .type = bool
      .help = "Recorded wavelength/Å"
      .short_caption = "λ/Å"
    n_reflections = True
      .type = bool
      .help = "Number of reflections"
      .short_caption = "#HKLs"
    span = True
      .type = bool
      .help = "Crude range of hkl indices"
      .short_caption = "Span"
    minmax_data = True
      .type = bool
      .help = "minimum, maximum values of data"
      .short_caption = "min,max data"
    minmax_sigmas = True
      .type = bool
      .help = "minimum, maximum values of sigmas"
      .short_caption = "min,max sigmas"
    data_sigdata = False
      .type = bool
      .help = "Average value of data/sigma"
      .short_caption = "DatSigDat"
    data_sigdata_max = False
      .type = bool
      .help = "maximum value of data/sigma"
      .short_caption = "MaxDatSigDat"
    d_minmax = True
      .type = bool
      .help = "d_min,d_max/Å"
      .short_caption = "d_min,d_max/Å"
    unit_cell = False
      .type = bool
      .help = "Unit cell parameters (a/Å, b/Å, c/Å, α°, β°, γ°)"
      .short_caption = "unit cell (a/Å, b/Å, c/Å, α°, β°, γ°)"
    space_group = False
      .type = bool
      .help = "Space group"
      .short_caption = "space group"
    n_centrics = False
      .type = bool
      .help = "Number of centric reflections"
      .short_caption = "#centrics"
    is_anomalous = True
      .type = bool
      .help = "Is data anomalous"
      .short_caption = "Anomalous"
    is_symmetry_unique = True
      .type = bool
      .help = "Is data symmetry unique"
      .short_caption = "Sym.uniq."
    n_sys_abs = False
      .type = bool
      .help = "Systematic absences"
      .short_caption = "#Syst.abs."
    data_completeness = True
      .type = bool
      .help = "Completeness in resolution range"
      .short_caption = "Data compl."
    data_compl_infty = False
      .type = bool
      .help = "Completeness with d_max=infinity"
      .short_caption = "Compl.inf."
    ano_completeness = False
      .type = bool
      .help = "Anomalous completeness in resolution range"
      .short_caption = "Ano.complete"
    ano_mean_diff = False
      .type = bool
      .help = "Mean anomalous difference."
      .short_caption = "Ano.dif. "
    n_bijvoet = False
      .type = bool
      .help = "Number of Bijvoet pairs"
      .short_caption = "#Bijvoets"
    n_singletons = False
      .type = bool
      .help = "Number of lone anomalous reflections"
      .short_caption = "#Singletons"
  }

  """

  philobj = parse(arrayinfo_phil_str)
  # make a dictionary of user friendly captions for HKLviewer GUI displaying list of table headers
  caption_dict = {}
  for objs in philobj.objects:
    if objs.name == "selected_info":
      for obj in objs.objects:
        caption_dict[obj.name] = (obj.short_caption, obj.help)


  def get_selected_info_columns_from_phil(self,philxtr=None):
    """
    Returns
    info_fmt: A list of selected property values from the miller_array object, together with their
      respective format strings for presenting in a users table such as in HKLviewer.
    headerstr: A formatted string of column names for printing a table to stdout.
    infostr: A formatted string of selected property values from the miller_array object for printing
      a table to stdout.

    If printing a table to stdout this can be done like:

    for i,array in enumerate(arrays):
      arrayinfo = ArrayInfo(array)
      info_fmt, headerstr, infostr = arrayinfo.get_selected_info_columns_from_phil()
      if i==0:
        print(headerstr)
      print(infostr)

    """
    info_format_tpl = []
    if not philxtr: # then use the default values in the arrayinfo_phil_str
      philxtr = parse(self.arrayinfo_phil_str).extract()
    delim = philxtr.delimiter
    for colname,selected in list(philxtr.selected_info.__dict__.items()):
      if not colname.startswith("__"):
        if selected:
          info_format_tpl.append( self.info_format_dict[colname] )
    # transpose info_format_tpl to return a list of headers, a list of values, and two lists of format strings
    info_fmt = list(zip(*info_format_tpl))
    # make a line of array info formatted for a table as well as a table header
    headerlst, infolst, dummy, fmtlst = info_fmt
    headerstr = ""
    for h in headerlst:
      headerstr += h + delim
    infostr = ""
    for i,info in enumerate(infolst):
      inf = info
      if i==0 and self.wrap_labels>0:
        inf = textwrap.wrap(info, width=self.wrap_labels)
      infostr += fmtlst[i].format(inf)  + delim
    return info_fmt, headerstr, infostr
