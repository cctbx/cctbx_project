
# XXX: these functions are utilities used by the Phenix GUI for quickly
# converting data and map formats on the fly.  the main use of this is
# in phenix.refine and any program that opens maps in PyMOL.
# tested as part of Phenix regression tests (since it requires files)

import mmtbx.maps
import iotbx.phil
from iotbx import file_reader
from libtbx.math_utils import ifloor, iceil
from libtbx.utils import Sorry
import os

def extract_map_coeffs (miller_arrays, f_lab, phi_lab, fom_lab) :
  f_array = None
  phi_array = None
  fom_array = None
  #print f_lab, phi_lab, fom_lab
  for miller_array in miller_arrays :
    labels = miller_array.info().label_string()
    if labels == f_lab :
      f_array = miller_array
    elif labels == phi_lab :
      phi_array = miller_array
    elif labels == fom_lab :
      fom_array = miller_array
  return (f_array, phi_array, fom_array)

def get_maps_from_fmodel (fmodel, use_filled=False) :
  map_manager = fmodel.electron_density_map(fill_missing_f_obs=use_filled,
                                            fill_mode="dfmodel",
                                            filled_f_obs_file_name="none")
  fwt_coeffs = map_manager.map_coefficients(map_type = "2mFo-DFc")
  if fwt_coeffs.anomalous_flag() :
    fwt_coeffs = fwt_coeffs.average_bijvoet_mates()
  delfwt_coeffs = map_manager.map_coefficients(map_type = "mFo-DFc")
  if delfwt_coeffs.anomalous_flag() :
    delfwt_coeffs = delfwt_coeffs.average_bijvoet_mates()
  return (fwt_coeffs, delfwt_coeffs)

def write_map_coeffs (fwt_coeffs, delfwt_coeffs, file_name) :
  mtz_dataset = fwt_coeffs.as_mtz_dataset(
    column_root_label="2FOFCWT")
  mtz_dataset.add_miller_array(
    miller_array=delfwt_coeffs,
    column_root_label="FOFCWT")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name=file_name)
  del mtz_object

# TODO: make more modular!
def write_xplor_map_file (coeffs, frac_min, frac_max, file_base) :
  fft_map = coeffs.fft_map(resolution_factor=1/3.0)
  fft_map.apply_sigma_scaling()
  n_real = fft_map.n_real()
  gridding_first=[ifloor(f*n) for f,n in zip(frac_min,n_real)]
  gridding_last=[iceil(f*n) for f,n in zip(frac_max,n_real)]
  title_lines=["REMARK map covering model + 3.0A buffer"]
  fft_map.as_xplor_map(
    file_name="%s.map" % file_base,
    title_lines=title_lines,
    gridding_first=gridding_first,
    gridding_last=gridding_last)

def xplor_maps_from_refine_mtz (pdb_file, mtz_file, file_base=None,
    limit_arrays=None, grid_resolution_factor=0.33) :
  if file_base is None :
    file_base = os.path.join(os.path.dirname(mtz_file), "refine")
  if not os.path.isfile(pdb_file) :
    raise Sorry("The PDB file '%s' does not exist.")
  if not os.path.isfile(mtz_file) :
    raise Sorry(
      "No map coefficients are available for conversion to XPLOR format.")
  mtz_in = file_reader.any_file(mtz_file)
  pdb_in = file_reader.any_file(pdb_file)
  mtz_in.assert_file_type("hkl")
  pdb_in.assert_file_type("pdb")
  miller_arrays = mtz_in.file_object.as_miller_arrays()
  assert len(miller_arrays) > 0
  xray_structure = pdb_in.file_object.xray_structure_simple()
  map_names = {"2FOFCWT" : "2mFo-DFc",
               "FOFCWT" : "mFo-DFc",
               "2FOFCWT_no_fill" : "2mFo-DFc_no_fill",
               "FOFCWT_no_fill" : "mFo-DFc_no_fill"}
  output_files = []
  for miller_array in miller_arrays :
    if miller_array.is_complex_array() :
      labels = miller_array.info().label_string()
      if limit_arrays is not None and not labels in limit_arrays :
        continue
      f_label = miller_array.info().labels[0]
      map_name = map_names.get(f_label)
      if map_name is None :
        map_name = f_label
      file_name = "%s_%s.map" % (file_base, map_name)
      if (not os.path.exists(file_name) or
          os.path.getmtime(file_name) < os.path.getmtime(mtz_file)) :
        xplor_map_from_coeffs(miller_array=miller_array,
          output_file=file_name,
          xray_structure=xray_structure,
          grid_resolution_factor=grid_resolution_factor)
        output_files.append(file_name)
  return output_files

def xplor_map_from_coeffs (miller_array, output_file, pdb_file=None,
    xray_structure=None, grid_resolution_factor=0.33) :
  map_phil = iotbx.phil.parse("""map.file_name = %s """ % output_file)
  master_phil = mmtbx.maps.map_and_map_coeff_master_params()
  params = master_phil.fetch(source=map_phil).extract().map[0]
  params.file_name = output_file
  params.region = "selection"
  params.scale = "sigma"
  params.grid_resolution_factor = grid_resolution_factor
  if xray_structure is None :
    if pdb_file is None or not os.path.isfile(pdb_file) :
      params.region = "cell"
    else :
      pdb_in = file_reader.any_file(pdb_file)
      pdb_in.assert_file_type("pdb")
      xray_structure = pdb_in.file_object.xray_structure_simple()
  mmtbx.maps.write_xplor_map_file(params=params,
    coeffs=miller_array,
    xray_structure=xray_structure)

# XXX: sorta gross, but for pre-calculated map coefficients like FWT,PHWT,
# pass the label string as the f_label argument.
def xplor_map_from_mtz (pdb_file, mtz_file, output_file=None,
    f_label="FP", phi_label="PHIM", fom_label="FOMM") :
  if output_file is None :
    output_file = "resolve.map"
  if not os.path.isfile(mtz_file) :
    raise Sorry(
      "No map coefficients are available for conversion to XPLOR format.")
  mtz_in = file_reader.any_file(mtz_file)
  mtz_in.assert_file_type("hkl")
  miller_arrays = mtz_in.file_object.as_miller_arrays()
  (f_array, phi_array, fom_array) = extract_map_coeffs(miller_arrays,
    f_label, phi_label, fom_label)
  if f_array.is_complex_array() :
    map_coeffs = f_array
  else :
    if f_array is None or phi_array is None :
      raise Sorry("One or more of the columns %s and %s was not found." %
        (f_label, phi_label))
    if fom_array is not None :
      weighted_f = f_array * fom_array
    else :
      weighted_f = f_array
    map_coeffs = weighted_f.phase_transfer(phi_array, deg=True)
  if (not os.path.isfile(output_file) or
      os.path.getmtime(output_file) < os.path.getmtime(mtz_file)) :
    xplor_map_from_coeffs(map_coeffs, output_file, pdb_file)
  return output_file

def xplor_map_from_resolve_mtz (pdb_file, mtz_file, force=False) :
  output_file = mtz_file[:-4] + ".map"
  if (force or not os.path.isfile(output_file) or
      os.path.getmtime(output_file) < os.path.getmtime(mtz_file)) :
    xplor_map_from_mtz(pdb_file=pdb_file,
      mtz_file=mtz_file,
      output_file=output_file,
      f_label="FP,SIGFP",
      phi_label="PHIM",
      fom_label="FOMM")
  return output_file

def xplor_map_from_solve_mtz (pdb_file, mtz_file, force=False) :
  output_file = mtz_file[:-4] + ".map"
  if (force or not os.path.isfile(output_file) or
      os.path.getmtime(output_file) < os.path.getmtime(mtz_file)) :
    xplor_map_from_mtz(pdb_file=pdb_file,
      mtz_file=mtz_file,
      output_file=output_file,
      f_label="FP,SIGFP",
      phi_label="PHIB",
      fom_label="FOM")
  return output_file

def get_map_coeff_labels (server, build_only=False) :
  all_labels = []
  phi_labels = []
  fom_labels = []
  for miller_array in server.miller_arrays :
    label = miller_array.info().label_string()
    if label.startswith("FOM") :
      fom_labels.append(label)
  phase_arrays = server.get_phases_deg(None, None, False, None, None, None,
                                       True, 3)
  for miller_array in phase_arrays :
    if miller_array.is_hendrickson_lattman_array() :
      continue
    elif miller_array.is_complex_array() :
      labels = miller_array.info().label_string()
      # note: Phaser outputs FWT/DELFWT for *anomalous difference* map!
      if build_only :
        if not labels.startswith("FOFC") and labels != "FWT,DELFWT" :
          all_labels.append(miller_array.info().label_string())
      else :
        all_labels.append(miller_array.info().label_string())
    elif miller_array.info().labels[0].startswith("PHI") :
      phi_labels.append(miller_array.info().label_string())
  amp_arrays = server.get_amplitudes(None, None, False, None, None, True, 4)
  for miller_array in amp_arrays :
    f_label = miller_array.info().labels[0]
    if f_label[0] == "F" and f_label != "FC" :
      for phase_label in phi_labels :
        hybrid_label = "%s,%s" % (f_label, phase_label)
        if len(fom_labels) > 0 :
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
