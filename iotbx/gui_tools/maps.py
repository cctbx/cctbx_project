
from libtbx.math_utils import ifloor, iceil
import os

MAP_TYPE_F_OBS = 1
MAP_TYPE_DIFF = 2
MAP_TYPE_ANOMALOUS = 4
MAP_TYPE_ISO_DIFF = 8
MAP_TYPE_NEUTRON = 16
MAP_TYPE_F_CALC = 32
MAP_TYPE_FILLED = 64

class server (object) :
  def __init__ (self, file_name=None, miller_arrays=()) :
    assert (file_name is not None) or (len(miller_arrays) > 0)
    self.file_name = file_name
    self.miller_arrays = miller_arrays
    if (len(self.miller_arrays) == 0) :
      from iotbx import file_reader
      f = file_reader.any_file(file_name, force_type="hkl")
      f.assert_file_type("hkl")
      self.miller_arrays = f.file_server.miller_arrays
    assert (len(self.miller_arrays) > 0)
    self.map_coeffs = []
    self.array_labels = []
    self.data_arrays = []
    self.phase_arrays = []
    self.weight_arrays = []
    self.rfree_arrays = []
    self.fcalc_arrays = []
    self._pdb_cache = {}
    from iotbx.reflection_file_utils import looks_like_r_free_flags_info
    for array in self.miller_arrays :
      labels = array.info().label_string()
      self.array_labels.append(labels)
      if (array.is_xray_amplitude_array() or array.is_xray_intensity_array()) :
        if (labels.startswith("FC")) :
          self.fcalc_arrays.append(array)
        else :
          self.data_arrays.append(array)
      elif (array.is_complex_array()) :
        if (labels.startswith("FC") or labels.startswith("FMODEL") or
            labels.startswith("F-model")) :
          self.fcalc_arrays.append(array)
        else :
          self.map_coeffs.append(array)
      elif (array.is_real_array()) :
        if (labels.startswith("PH")) :
          self.phase_arrays.append(array)
        elif (labels.startswith("FOM")) :
          self.weight_arrays.append(array)
      elif (array.is_bool_array()) :
        self.rfree_arrays.append(array)
      elif (array.is_integer_array()) :
        if (looks_like_r_free_flags_info(array.info())) :
          self.rfree_arrays.append(array)

  def is_phenix_refine_maps (self) :
    return "2FOFCWT,PH2FOFCWT" in self.array_labels

  def is_cctbx_maps (self) :
    return "2mFoDFc,P2mFoDFc" in self.array_labels

  def is_resolve_map (self) :
    return ((("FP" in self.array_labels) or ("FP,SIGFP" in self.array_labels))
          and ("PHIM" in self.array_labels) and ("FOMM" in self.array_labels))

  def is_ccp4_style_map (self) :
    return (("FWT,PHWT" in self.array_labels) or
            ("FWT,PHIFWT" in self.array_labels))

  def is_solve_map (self) :
    return ((("FP" in self.array_labels) or ("FP,SIGFP" in self.array_labels))
          and ("PHIB" in self.array_labels) and ("FOM" in self.array_labels))

  def is_anomalous_map (self) :
    for array in self.map_coeffs :
      if array.info().label_string().startswith("ANOM") :
        return True
    return False

  def contains_map_coefficients (self) :
    return (len(self.map_coeffs) > 0)

  def contains_amplitude_and_phase (self) :
    return (len(self.data_arrays) > 0) and (len(self.phase_arrays) > 0)

  def get_amplitudes_and_phases (self,
                                 f_label=None,
                                 phi_label=None,
                                 fom_label=None) :
    f_array, phi_array, fom_array = None, None, None
    for array in self.data_arrays :
      if (f_array is None) :
        labels = array.info().labels
        label_string = array.info().label_string()
        if (f_label is not None) :
          if (label_string == f_label) or (labels[0] == f_label) :
            f_array = array
        elif (labels[0] == "FP") :
          f_array = array
        elif (labels[0] == "F") :
          f_array = array
    for array in self.phase_arrays :
      if (phi_array is None) :
        labels = array.info().label_string()
        if (phi_label is not None) :
          if (labels == phi_label) :
            phi_array = array
            continue
        elif (labels == "PHIM") :
          phi_array = array
        elif (labels == "PHIB") and (not "PHIM" in self.array_labels) :
          phi_array = array
        elif (labels == "PHIC") and not (("PHIM" in self.array_labels) or
                                         ("PHIB" in self.array_labels)) :
          phi_array = array
        elif (labels == "PHI") :
          phi_array = array
    for array in self.weight_arrays :
      if (fom_array is None) :
        labels = array.info().label_string()
        if (fom_label is not None) :
          if (labels == fom_label) :
            fom_array = array
        elif (labels == "FOMM") :
          fom_array = array
        elif (labels == "FOM") :
          fom_array = array
    return (f_array, phi_array, fom_array)

  def get_phenix_maps (self, keep_column_labels=False) :
    f_maps = []
    diff_maps = []
    default_is_filled = None
    all_maps = []
    unknown_maps = []
    have_names = []
    for array in self.map_coeffs :
      labels = array.info().label_string()
      if (labels.endswith("_no_fill")) :
        default_is_filled = True
      elif (labels.endswith("_fill")) :
        default_is_filled = False
      if labels.startswith("2FOFCWT") or labels.startswith("2mFoDFc") :
        f_maps.append(array)
      elif labels.startswith("FOFCWT") or labels.startswith("mFoDFc") :
        diff_maps.append(array)
      elif labels.startswith("ANOM") :
        map_type = MAP_TYPE_ANOMALOUS
        map_name = "anomalous"
        if ("neutron" in labels) :
          map_type |= MAP_TYPE_NEUTRON
          map_name += "_neutron"
        if (keep_column_labels) :
          map_name = labels.split(",")[0]
        elif (map_name in have_names) :
          map_name += "_%d" % (have_names.count(map_name) + 1)
        have_names.append(map_name)
        all_maps.append(map_info(map_coeffs=array,
          map_type=map_type,
          map_name=map_name))
      elif (labels.startswith("F-model")) or (labels.startswith("FC")) :
        continue
      else :
        unknown_maps.append(array)
    unfilled_maps = []
    filled_maps = []
    for array in diff_maps :
      labels = array.info().label_string()
      map_type = MAP_TYPE_DIFF
      map_name = "mFo-DFc"
      if ("neutron" in labels) :
        map_type |= MAP_TYPE_NEUTRON
        map_name += "_neutron"
      if (keep_column_labels) :
        map_name = labels.split(",")[0]
      elif (map_name in have_names) :
       map_name += "_%d" % (have_names.count(map_name) + 1)
      have_names.append(map_name)
      all_maps.append(map_info(map_coeffs=array,
        map_type=map_type,
        map_name=map_name))
    for array in f_maps :
      labels = array.info().label_string()
      map_name = "2mFo-DFc"
      map_type = MAP_TYPE_F_OBS
      if ("neutron" in labels) :
        map_name += "_neutron"
        map_type |= MAP_TYPE_NEUTRON
      if labels.endswith("_no_fill") :
        map_name += "_no_fill"
      elif labels.endswith("_fill") or labels.endswith("_filled") :
        map_type |= MAP_TYPE_FILLED
      elif (default_is_filled) :
        map_type |= MAP_TYPE_FILLED
      else :
        map_name += "_no_fill"
      if (keep_column_labels) :
        map_name = labels.split(",")[0]
      elif (map_name in have_names) :
       map_name += "_%d" % (have_names.count(map_name) + 1)
      have_names.append(map_name)
      all_maps.append(map_info(map_coeffs=array,
        map_type=map_type,
        map_name=map_name))
    return all_maps

  def get_resolve_map (self) :
    map_coeffs = self._convert_amplitudes_and_phases(f_label="FP",
      phi_label="PHIM", fom_label="FOMM", weighted=True)
    from cctbx.miller import array_info
    info = array_info(labels=["FWT","PHWT"])
    return map_coeffs.set_info(info)

  def get_solve_map (self) :
    map_coeffs = self._convert_amplitudes_and_phases(f_label="FP",
      phi_label="PHIB", fom_label="FOM", weighted=True)
    from cctbx.miller import array_info
    info = array_info(labels=["FWT","PHWT"])
    return map_coeffs.set_info(info)

  def get_ccp4_maps (self) :
    f_map = diff_map = None
    for array in self.map_coeffs :
      labels = array.info().label_string()
      if (labels.startswith("FWT,")) :
        f_map = array
      elif (labels.startswith("DELFWT,")) :
        diff_map = array
    return (f_map, diff_map)

  def get_anomalous_map (self) :
    for array in self.map_coeffs :
      labels = array.info().label_string()
      if (labels.upper().startswith("ANOM")) :
        return array
    return None

  def _convert_amplitudes_and_phases (self, f_label=None, phi_label=None,
      fom_label=None, weighted=True) :
    f, phi, fom = self.get_amplitudes_and_phases(f_label, phi_label, fom_label)
    if (f is None) or (phi is None) :
      raise RuntimeError(("Couldn't find amplitude or phase arrays in %s.\n"+
        "File contents:\n%s") % (self.file_name, "\n".join(self.array_labels)))
    if (f.anomalous_flag()) :
      f = f.merge_bijvoet_mates()
    if (weighted) and (fom is not None) :
      f = f * fom
    map_coeffs = f.phase_transfer(phi, deg=True) # XXX is deg always True?
    return map_coeffs

  def convert_any_map (self, f_label, phi_label, fom_label, **kwds) :
    map_coeffs = self._convert_amplitudes_and_phases(f_label, phi_label,
      fom_label)
    return self._write_ccp4_map(map_coeffs, **kwds)

  def convert_resolve_map (self, **kwds) :
    map_coeffs = self.get_resolve_map()
    kwds['map_coeffs'] = map_coeffs
    kwds['simple_file_name'] = True
    return self._write_ccp4_map(**kwds)

  def convert_solve_map (self, **kwds) :
    map_coeffs = self.get_solve_map()
    kwds['map_coeffs'] = map_coeffs
    kwds['simple_file_name'] = True
    return self._write_ccp4_map(**kwds)

  def convert_ccp4_map (self, **kwds) :
    map_coeffs_list = list(self.get_ccp4_maps())
    files = []
    #kwds['simple_file_name'] = True
    for map_coeffs in map_coeffs_list :
      if (map_coeffs is None) : continue
      kwds['map_coeffs'] = map_coeffs
      output_file = self._write_ccp4_map(**kwds)
      files.append(os.path.abspath(output_file))
    return files

  def convert_phenix_maps (self, file_base, **kwds) :
    all_maps = self.get_phenix_maps()
    files = []
    for map in all_maps :
      file_name = "%s_%s.ccp4" % (file_base, map.map_name)
      kwds['output_file'] = file_name
      self._write_ccp4_map(map.map_coeffs, **kwds)
      files.append(os.path.abspath(file_name))
    return files

  def get_pdb_file (self, file_name) :
    pdb_in, mtime = self._pdb_cache.get(file_name, (None, 0))
    if (pdb_in is None) or (os.path.getmtime(file_name) > mtime) :
      from iotbx import file_reader
      f = file_reader.any_file(file_name, force_type="pdb")
      f.assert_file_type("pdb")
      pdb_in = f.file_object
      self._pdb_cache[file_name] = (pdb_in, os.path.getmtime(file_name))
    return pdb_in

  def _write_ccp4_map (self,
                       map_coeffs,
                       pdb_file=None,
                       output_file=None,
                       simple_file_name=False,
                       resolution_factor=0.25,
                       sigma_scaling=True,
                       force=False) :
    sites_cart = None
    if (pdb_file is not None) :
      pdb_in = self.get_pdb_file(pdb_file)
      sites_cart = pdb_in.atoms().extract_xyz()
    if (output_file is None) :
      if (simple_file_name) :
        output_file = os.path.splitext(self.file_name)[0] + ".ccp4"
      else :
        lab1 = map_coeffs.info().label_string().split(",")[0]
        output_file = os.path.splitext(self.file_name)[0] + "_%s.ccp4" % lab1
    if (not force) and (os.path.exists(output_file) and
        os.path.getmtime(output_file) < os.path.getmtime(self.file_name)) :
      return output_file
    fft_map = map_coeffs.fft_map(resolution_factor=resolution_factor)
    if (sigma_scaling) :
      fft_map.apply_sigma_scaling()
    else :
      fft_map.apply_volume_scaling()
    map_data = fft_map.real_map()
    write_ccp4_map(
      sites_cart=sites_cart,
      unit_cell=map_coeffs.unit_cell(),
      map_data=map_data,
      n_real=fft_map.n_real(),
      file_name=output_file)
    return output_file

  def auto_open_maps (self, use_filled=True) :
    maps = []
    if (self.is_phenix_refine_maps()) or (self.is_cctbx_maps()) :
      maps = self.get_phenix_maps()
    elif (self.is_ccp4_style_map()) :
      f_map, diff_map = self.get_ccp4_maps()
      if (f_map is not None) :
        maps.append(map_info(f_map, MAP_TYPE_F_OBS, "2mFo-DFc"))
      if (diff_map is not None) :
        maps.append(map_info(diff_map, MAP_TYPE_DIFF, "mFo-DFc"))
    elif (self.is_solve_map()) :
      f_map = self.get_solve_map()
      maps.append(map_info(f_map, MAP_TYPE_F_OBS, "FP,PHIB,FOM"))
    elif (self.is_resolve_map()) :
      f_map = self.get_resolve_map()
      maps.append(map_info(f_map, MAP_TYPE_F_OBS, "FP,PHIM,FOMM"))
    elif (self.is_anomalous_map()) :
      anom_map = self.get_anomalous_map()
      maps.append(map_info(anom_map, MAP_TYPE_ANOMALOUS, "anomalous"))
    return maps

class map_info (object) :
  def __init__ (self, map_coeffs, map_type, map_name=None) :
    self.map_coeffs = map_coeffs
    self.map_type = map_type
    self.map_name = map_name

  def is_fobs_map (self) :
    return (self.map_type & MAP_TYPE_F_OBS)

  def is_difference_map (self) :
    return (self.map_type & MAP_TYPE_DIFF)

  def is_anomalous_map (self) :
    return (self.map_type & MAP_TYPE_ANOMALOUS)

  def is_neutron_map (self) :
    return (self.map_type & MAP_TYPE_NEUTRON)

  def is_filled_map (self) :
    return (self.map_type & MAP_TYPE_FILLED)

  def fft_map (self, resolution_factor) :
    map_coeffs = self.map_coeffs.unique_under_symmetry()
    map = map_coeffs.fft_map(resolution_factor=resolution_factor)
    map.apply_sigma_scaling()
    return map

  def __getattr__ (self, name) :
    return getattr(self.map_coeffs, name)

def write_ccp4_map (sites_cart, unit_cell, map_data, n_real, file_name,
    buffer=10) :
  import iotbx.ccp4_map
  from cctbx import sgtbx
  from scitbx.array_family import flex
  if sites_cart is not None :
    frac_min, frac_max = unit_cell.box_frac_around_sites(
      sites_cart=sites_cart,
      buffer=buffer)
  else :
    frac_min, frac_max = (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)
  gridding_first = tuple([ifloor(f*n) for f,n in zip(frac_min,n_real)])
  gridding_last = tuple([iceil(f*n) for f,n in zip(frac_max,n_real)])
  space_group = sgtbx.space_group_info("P1").group()
  iotbx.ccp4_map.write_ccp4_map(
    file_name=file_name,
    unit_cell=unit_cell,
    space_group=space_group,
    gridding_first=gridding_first,
    gridding_last=gridding_last,
    map_data=map_data,
    labels=flex.std_string(["iotbx.map_conversion.write_ccp4_map_box"]))

def convert_resolve_map (mtz_file, **kwds) :
  server_ = server(file_name=mtz_file)
  return server_.convert_resolve_map(**kwds)

def convert_solve_map (mtz_file, **kwds) :
  server_ = server(file_name=mtz_file)
  return server_.convert_solve_map(**kwds)

def convert_ccp4_map (mtz_file, **kwds) :
  server_ = server(file_name=mtz_file)
  return server_.convert_ccp4_map(**kwds)
