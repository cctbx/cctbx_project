
import iotbx.gui_tools
from iotbx import file_reader
from libtbx.utils import Sorry
import os

def space_group_as_str (space_group) :
  from cctbx import sgtbx
  if space_group is None :
    return ""
  elif isinstance(space_group, sgtbx.space_group_info) :
    return str(space_group)
  else :
    sg_info = sgtbx.space_group_info(group=space_group)
    return str(sg_info)

def unit_cell_as_str (unit_cell, separator=" ") :
  assert isinstance(separator, str) and separator != ""
  if unit_cell is not None :
    format = separator.join([ "%g" for i in range(6) ])
    return format % unit_cell.parameters()
  else :
    return ""

class reflections_handler (iotbx.gui_tools.manager) :
  file_type = "hkl"
  file_type_label = "Reflections"
  def __init__ (self,
                allowed_param_names=None,
                allowed_multiple_params=None,
                debug=False,
                minimum_data_score=4,
                prefer_amplitudes=False) :
    iotbx.gui_tools.manager.__init__(self,
      allowed_param_names=allowed_param_names,
      allowed_multiple_params=allowed_multiple_params,
      debug=debug)
    self.minimum_data_score = minimum_data_score
    self.prefer_amplitudes = prefer_amplitudes

  def check_symmetry (self, *args, **kwds) :
    hkl_file = self.get_file(*args, **kwds)
    hkl_server = hkl_file.file_server
    ma = hkl_server.miller_arrays[0]
    if (ma.space_group() is None) or (ma.unit_cell() is None) :
      raise Sorry("Incomplete symmetry for this reflections file.  Please "+
        "use a format that includes both space group and unit cell, such as "+
        "an MTZ or merged Scalepack file.")
    return True

  def get_file_type_label (self, file_name=None, input_file=None) :
    if (input_file is not None) :
      return input_file.file_object.file_type()

  def get_all_labels (self, *args, **kwds) :
    hkl_file = self.get_file(*args, **kwds)
    if hkl_file is None :
      labels = []
      for (file_name, hkl_file) in self.input_files() :
        miller_arrays = hkl_file.file_server.miller_arrays
        labels.extend([ a.info().label_string() for a in miller_arrays ])
    else :
      miller_arrays = hkl_file.file_server.miller_arrays
      labels = [ a.info().label_string() for a in miller_arrays ]
    return labels

  def get_rfree_labels (self, *args, **kwds) :
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
      if len(arrays_and_flags) > 0 :
        labels = [ a.info().label_string() for a, fl in arrays_and_flags ]
    return labels

  def has_rfree (self, *args, **kwds) :
    return len(self.get_rfree_labels(*args, **kwds)) > 0

  def get_rfree_flag_value (self, array_name, *args, **kwds) :
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

  def get_experimental_phase_labels (self, *args, **kwds) :
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

  def has_phases (self, *args, **kwds) :
    return len(self.get_experimental_phase_labels(*args, **kwds)) > 0

  def get_phase_deg_labels (self, *args, **kwds) :
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
      labels = [ array.info().label_string() for array in miller_arrays ]
    return labels

  def get_data_arrays (self, *args, **kwds) :
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
          if array.is_xray_amplitude_array() :
            f_arrays.append(array)
          else :
            other_arrays.append(array)
        miller_arrays = f_arrays + other_arrays
      return miller_arrays
    return []

  def get_data_labels (self, *args, **kwds) :
    return extract_labels(self.get_data_arrays(*args, **kwds))

  def get_amplitude_arrays (self, *args, **kwds) :
    hkl_file = self.get_file(*args, **kwds)
    if hkl_file is not None :
      hkl_server = hkl_file.file_server
      miller_arrays = hkl_server.get_amplitudes(
        file_name               = None,
        labels                  = None,
        ignore_all_zeros        = False,
        parameter_scope         = "",
        return_all_valid_arrays = True,
        minimum_score           = self.minimum_data_score)
      return miller_arrays
    return []

  def get_amplitude_labels (self, *args, **kwds) :
    return extract_labels(self.get_amplitude_arrays(*args, **kwds))

  def get_intensity_arrays (self, *args, **kwds) :
    miller_arrays = self.get_data_arrays(*args, **kwds)
    i_arrays = []
    for array in miller_arrays :
      if array.is_xray_intensity_array() :
        i_arrays.append(array)
    return i_arrays

  def get_intensity_labels (self, *args, **kwds) :
    return extract_labels(self.get_intensity_arrays(*args, **kwds))

  # space-separated, column labels only
  def get_data_labels_for_wizard (self, *args, **kwds) :
    miller_arrays = self.get_data_arrays(*args, **kwds)
    labels = [ " ".join(array.info().labels) for array in miller_arrays ]
    return labels

  def has_data (self, *args, **kwds) :
    return len(self.get_data_labels(*args, **kwds)) > 0

  def get_anomalous_data_labels (self, *args, **kwds) :
    kwds = dict(kwds) # XXX gross...
    allow_dano = kwds.get("allow_reconstructed_amplitudes", True)
    if ("allow_reconstructed_amplitudes" in kwds) :
      del kwds["allow_reconstructed_amplitudes"]
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
        if (array.anomalous_flag()) :
          if ((array.is_xray_reconstructed_amplitude_array()) and
              (not allow_dano)) :
            continue
          labels.append(array.info().label_string())
    return labels

  def has_anomalous_data (self, *args, **kwds) :
    return (len(self.get_anomalous_data_labels(*args, **kwds)) > 0)

  def get_map_coeff_labels (self, *args, **kwds) :
    hkl_file = self.get_file(*args, **kwds)
    if hkl_file is not None :
      return get_map_coeff_labels(hkl_file.file_server)
    return []

  def get_map_coeff_labels_for_build (self, *args, **kwds) :
    hkl_file = self.get_file(*args, **kwds)
    if hkl_file is not None :
      return get_map_coeffs_for_build(hkl_file.file_server)
    return []

  def get_map_coeff_labels_for_fft (self, *args, **kwds) :
    hkl_file = self.get_file(*args, **kwds)
    labels_list = []
    if hkl_file is not None :
      all_labels = get_map_coeff_labels(hkl_file.file_server,
        keep_array_labels=True)
      for labels in all_labels :
        if isinstance(labels, str) :
          labels_list.append(labels)
        else :
          labels_list.append(" ".join(labels))
    return labels_list

  def d_max_min (self, file_name=None, file_param_name=None,
      array_name=None, array_names=None) :
    from iotbx.reflection_file_editor import get_best_resolution
    miller_arrays = []
    if (file_name is None) and (file_param_name is None) :
      for phil_name, file_name in self._param_files.iteritems() :
        input_file = self.get_file(file_name)
        if (input_file is not None) :
          miller_arrays.extend(input_file.file_server.miller_arrays)
    else :
      hkl_file = self.get_file(file_name, file_param_name)
      if (hkl_file is None) :
        return (None, None)
      for miller_array in hkl_file.file_server.miller_arrays :
        label = miller_array.info().label_string()
        if (array_name is not None) :
          if (label == array_name) :
            miller_arrays = [miller_array]
            break
        elif (array_names is not None) :
          if (label in array_names) :
            miller_arrays.append(miller_array)
        else :
          miller_arrays.append(miller_array)
    (d_max, d_min) = get_best_resolution(miller_arrays)
    return (d_max, d_min)

  def get_resolution_range (self, *args, **kwds) :
    (d_max, d_min) = self.d_max_min(*args, **kwds)
    if (d_max is None) :
      return ""
    else :
      return "(%.3f - %.3f)" % (d_max, d_min)

  def get_resolution_limits (self, *args, **kwds) :
    (d_max, d_min) = self.d_max_min(*args, **kwds)
    (d_max_str, d_min_str) = ("", "")
    if (d_max is not None) :
      d_max_str = "(%.3f)" % d_max
    if (d_min is not None) :
      d_min_str = "(%.3f)" % d_min
    return (d_max_str, d_min_str)

  def space_group (self, *args, **kwds) :
    symm = self.crystal_symmetry(*args, **kwds)
    if symm is not None :
      return symm.space_group()
    return None

  def space_group_as_str (self, *args, **kwds) :
    space_group = self.space_group(*args, **kwds)
    return space_group_as_str(space_group)

  def unit_cell (self, *args, **kwds) :
    symm = self.crystal_symmetry(*args, **kwds)
    if symm is not None :
      return symm.unit_cell()
    return None

  def unit_cell_as_str (self, file_name=None, file_param_name=None,
      separator=" ") :
    unit_cell = self.unit_cell(file_name, file_param_name)
    return unit_cell_as_str(unit_cell, separator)

  def crystal_symmetry (self, file_name=None, file_param_name=None) :
    final_symm = None
    if file_name is None and file_param_name is None :
      for param_name, file_name in self._param_files.iteritems() :
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

  def check_symmetry_consistency (self) :
    pass

def extract_labels (miller_arrays) :
  return [ array.info().label_string() for array in miller_arrays ]

def get_fp_fpp_from_sasaki (guess_ha,wavelength):
  from cctbx.eltbx import sasaki
  from decimal import Decimal, ROUND_HALF_UP
  if wavelength is None or guess_ha is None :
    return None, None
  try:
    table = sasaki.table(guess_ha)
  except Exception, e :
    return None, None
  fp_fdp = table.at_angstrom(wavelength)
  f_prime=fp_fdp.fp()
  f_double_prime=fp_fdp.fdp()
  fp = Decimal(f_prime).quantize(Decimal('0.01'), ROUND_HALF_UP)
  fpp = Decimal(f_double_prime).quantize(Decimal('0.01'), ROUND_HALF_UP)
  return fp, fpp

def get_high_resolution (server) :
  d_min = None #999.99
  for miller_array in server.miller_arrays :
    try :
      array_min = miller_array.d_min()
      if d_min is None or array_min < d_min :
        d_min = array_min
    except Exception :
      pass
  return d_min

def get_miller_array_symmetry (miller_array) :
  from cctbx import sgtbx
  symm = miller_array.crystal_symmetry()
  if symm is None :
    return (None, None)
  unit_cell = symm.unit_cell()
  if (unit_cell is not None) :
    uc = unit_cell_as_str(unit_cell)
  else :
    uc = None
  space_group = symm.space_group()
  if (space_group is not None) :
    sg = str(sgtbx.space_group_info(group=space_group))
  else :
    sg = None
  return (sg, uc)

def get_array_description (miller_array) :
  from iotbx import reflection_file_utils
  info = miller_array.info()
  labels = info.label_string()
  if labels in ["PHI", "PHIB", "PHIM", "PHIC"] :
    return "Phases"
  if labels in ["FOM", "FOMM"] :
    return "Weights"
  if ((miller_array.is_integer_array() or miller_array.is_bool_array()) and
      reflection_file_utils.looks_like_r_free_flags_info(info)) :
    return "R-free flag"
  methods_and_meanings = [ ("is_complex_array", "Map coeffs"),
                           ("is_xray_amplitude_array", "Amplitude"),
                           ("is_xray_intensity_array", "Intensity"),
                           ("is_hendrickson_lattman_array", "HL coeffs"),
                           ("is_bool_array", "Boolean"),
                           ("is_integer_array", "Integer"),
                           ("is_real_array", "Floating-point"), ]
  for method, desc in methods_and_meanings :
    test = getattr(miller_array, method)
    if test() :
      return desc
  return "Unknown"

def get_mtz_label_prefix (input_file=None, file_name=None) :
  if input_file is None :
    input_file = file_reader.any_file(file_name)
  assert (input_file.file_type == "hkl" and
          input_file.file_object.file_type() == "ccp4_mtz")
  file_content = input_file.file_object.file_content()
  last_crystal = file_content.crystals()[-1]
  crystal = last_crystal.name()
  dataset = last_crystal.datasets()[-1].name()
  return "/%s/%s" % (crystal, dataset)

def extract_map_coeffs (miller_arrays, f_lab, phi_lab, fom_lab) :
  f_array = None
  phi_array = None
  fom_array = None
  #print f_lab, phi_lab, fom_lab
  for miller_array in miller_arrays :
    labels = miller_array.info().label_string()
    if (labels == f_lab) :
      f_array = miller_array
    elif (labels == phi_lab) :
      phi_array = miller_array
    elif (labels == fom_lab) :
      fom_array = miller_array
  return (f_array, phi_array, fom_array)

def map_coeffs_from_mtz_file (mtz_file, f_label="FP", phi_label="PHIM",
    fom_label="FOMM") :
  if not os.path.isfile(mtz_file) :
    raise Sorry(
      "No map coefficients are available for conversion.")
  mtz_in = file_reader.any_file(mtz_file)
  mtz_in.assert_file_type("hkl")
  miller_arrays = mtz_in.file_server.miller_arrays
  (f_array, phi_array, fom_array) = extract_map_coeffs(miller_arrays,
    f_label, phi_label, fom_label)
  if (f_array is not None) and (f_array.is_complex_array()) :
    map_coeffs = f_array
  else :
    if (f_array is None) or (phi_array is None) :
      raise Sorry("One or more of the columns %s and %s was not found." %
        (f_label, phi_label))
    if fom_array is not None :
      weighted_f = f_array * fom_array
    else :
      weighted_f = f_array
    map_coeffs = weighted_f.phase_transfer(phi_array, deg=True)
  return map_coeffs

def extract_phenix_refine_map_coeffs (mtz_file, limit_arrays=None) :
  assert (limit_arrays is None) or (isinstance(limit_arrays, list))
  if not os.path.isfile(mtz_file) :
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
    if miller_array.is_complex_array() :
      labels = miller_array.info().label_string()
      if labels.startswith("F-model") :
        continue
      if (limit_arrays is not None) and (not labels in limit_arrays) :
        continue
      f_label = miller_array.info().labels[0]
      map_name = map_names.get(f_label)
      if map_name is None :
        map_name = f_label
      output_arrays.append((miller_array, map_name))
  return output_arrays

def get_map_coeff_labels (server, build_only=False, include_fom=True,
    keep_array_labels=False) :
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
    if miller_array.is_hendrickson_lattman_array() :
      continue
    elif miller_array.is_complex_array() :
      # note: Phaser outputs FWT/DELFWT for *anomalous difference* map!
      if build_only :
        if (not labels.startswith("FOFC")) and (labels != "FWT,DELFWT") :
          all_labels.append(labels)
      else :
        all_labels.append(labels)
    elif miller_array.info().labels[0].startswith("PHI") :
      phi_labels.append(labels)
  amp_arrays = server.get_amplitudes(None, None, False, None, None, True, 4)
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

def get_map_coeffs_for_build (server) :
  return get_map_coeff_labels(server, build_only=True)

def format_map_coeffs_for_resolve (f_label, phi_label, fom_label) :
  return "FP=%s PHIB=%s FOM=%s" % (f_label, phi_label, fom_label)

def decode_resolve_map_coeffs (labels) :
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
