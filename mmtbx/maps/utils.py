
# XXX: these functions are utilities used by the Phenix GUI for quickly
# converting data and map formats on the fly.  the main use of this is
# in phenix.refine and any program that opens maps in PyMOL.
# tested as part of Phenix regression tests (since it requires files)

from __future__ import division
import mmtbx.maps
import mmtbx.utils
import iotbx.phil
import iotbx.xplor.map
import iotbx.ccp4_map
from iotbx import file_reader
from cctbx import sgtbx
from scitbx.array_family import flex
from libtbx.math_utils import ifloor, iceil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import os, sys

#-----------------------------------------------------------------------
# MAP COEFFICIENT MANIPULATION

class fast_maps_from_hkl_file (object) :
  def __init__ (self,
                file_name,
                xray_structure,
                f_label=None,
                map_out=None,
                log=sys.stdout,
                auto_run=True) :
    adopt_init_args(self, locals())
    if f_label is None :
      print >> log, "      no label for %s, will try default labels" % \
        os.path.basename(file_name)
    f_obs = None
    default_labels = ["F,SIGF","FOBS,SIGFOBS","F(+),SIGF(+),F(-),SIGF(-)"]
    all_labels = []
    data_file = file_reader.any_file(file_name, force_type="hkl")
    for miller_array in data_file.file_server.miller_arrays :
      labels = miller_array.info().label_string()
      all_labels.append(labels)
      if labels == f_label :
        f_obs = miller_array
        break
      elif f_label is None and labels in default_labels :
        f_obs = miller_array
        break
    if f_obs is None :
      raise Sorry(("Couldn't find %s in %s.  Please specify valid "+
        "column labels (possible choices: %s)") % (f_label, file_name,
          " ".join(all_labels)))
    self.f_obs = f_obs
    self.log = None
    if auto_run :
      self.run()

  def run (self) :
    f_obs = self.f_obs
    r_free_flags = f_obs.array(data=flex.bool(f_obs.data().size(),False))
    fmodel = mmtbx.utils.fmodel_simple(
      xray_structures=[self.xray_structure],
      f_obs=f_obs,
      r_free_flags=r_free_flags)
    (f_map, df_map) = get_maps_from_fmodel(fmodel, use_filled=True)
    if self.map_out is None :
      self.map_out = os.path.splitext(self.file_name)[0] + "_map_coeffs.mtz"
    write_map_coeffs(f_map, df_map, self.map_out)

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
      if limit_arrays is not None and not labels in limit_arrays :
        continue
      f_label = miller_array.info().labels[0]
      map_name = map_names.get(f_label)
      if map_name is None :
        map_name = f_label
      output_arrays.append((miller_array, map_name))
  return output_arrays

#-----------------------------------------------------------------------
# XPLOR MAP OUTPUT

# TODO: make more modular!
def write_xplor_map_file (coeffs, frac_min, frac_max, file_base) :
  fft_map = coeffs.fft_map(resolution_factor=1/3.0)
  fft_map.apply_sigma_scaling()
  n_real = fft_map.n_real()
  gridding_first=[ifloor(f*n) for f,n in zip(frac_min,n_real)]
  gridding_last=[iceil(f*n) for f,n in zip(frac_max,n_real)]
  title_lines=["REMARK map covering model + 3.0A buffer"]
  file_name = "%s.map" % file_base
  fft_map.as_xplor_map(
    file_name=file_name,
    title_lines=title_lines,
    gridding_first=gridding_first,
    gridding_last=gridding_last)
  return file_name

def write_xplor_map(sites_cart, unit_cell, map_data, n_real, file_name,
    buffer=10) :
  if sites_cart is not None :
    frac_min, frac_max = unit_cell.box_frac_around_sites(
      sites_cart=sites_cart,
      buffer=buffer)
  else :
    frac_min, frac_max = (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)
  gridding_first=[ifloor(f*n) for f,n in zip(frac_min,n_real)]
  gridding_last=[iceil(f*n) for f,n in zip(frac_max,n_real)]
  gridding = iotbx.xplor.map.gridding(n     = map_data.focus(),
                                      first = gridding_first,
                                      last  = gridding_last)
  iotbx.xplor.map.writer(
    file_name          = file_name,
    is_p1_cell         = True,
    title_lines        = [' None',],
    unit_cell          = unit_cell,
    gridding           = gridding,
    data               = map_data,
    average            = -1,
    standard_deviation = -1)

def xplor_maps_from_refine_mtz (pdb_file, mtz_file, file_base=None,
    limit_arrays=None, grid_resolution_factor=0.33) :
  if file_base is None :
    file_base = os.path.join(os.path.dirname(mtz_file), "refine")
  if not os.path.isfile(pdb_file) :
    raise Sorry("The PDB file '%s' does not exist.")
  output_arrays = extract_phenix_refine_map_coeffs(mtz_file, limit_arrays)
  pdb_in = file_reader.any_file(pdb_file)
  pdb_in.assert_file_type("pdb")
  xray_structure = pdb_in.file_object.xray_structure_simple()
  output_files = []
  for (map_coeffs, map_name) in output_arrays :
    file_name = "%s_%s.map" % (file_base, map_name)
    if (not os.path.exists(file_name) or
        os.path.getmtime(file_name) < os.path.getmtime(mtz_file)) :
      xplor_map_from_coeffs(miller_array=map_coeffs,
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
  map_coeffs = map_coeffs_from_mtz_file(mtz_file, f_label, phi_label, fom_label)
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

#-----------------------------------------------------------------------
# CCP4 MAP OUTPUT

def ccp4_maps_from_refine_mtz (mtz_file, file_base=None,
    limit_arrays=None, resolution_factor=0.33) :
  if file_base is None :
    file_base = os.path.join(os.path.dirname(mtz_file), "refine")
  output_arrays = extract_phenix_refine_map_coeffs(mtz_file)
  output_files = []
  for (map_coeffs, map_name) in output_arrays :
    file_name = "%s_%s.ccp4" % (file_base, map_name)
    if (not os.path.exists(file_name) or
        os.path.getmtime(file_name) < os.path.getmtime(mtz_file)) :
      fft_map = map_coeffs.fft_map(resolution_factor=resolution_factor)
      fft_map.apply_sigma_scaling()
      fft_map.as_ccp4_map(file_name=file_name)
      output_files.append(file_name)
  return output_files

def ccp4_map_from_mtz (mtz_file, output_file=None, f_label="FP",
    phi_label="PHIM", fom_label="FOMM", resolution_factor=1/3.0,
    force=True) :
  if output_file is None :
    output_file = os.path.splitext(mtz_file)[0] + ".ccp4"
  if (force or not os.path.isfile(output_file) or
      os.path.getmtime(output_file) < os.path.getmtime(mtz_file)) :
    map_coeffs = map_coeffs_from_mtz_file(mtz_file, f_label, phi_label,
      fom_label)
    assert map_coeffs.is_complex_array()
    fft_map = map_coeffs.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    fft_map.as_ccp4_map(file_name=output_file)
  return output_file

def ccp4_map_from_resolve_mtz (mtz_file, force=False, resolution_factor=1/3.0) :
  return ccp4_map_from_mtz(mtz_file=mtz_file,
    f_label="FP,SIGFP",
    phi_label="PHIM",
    fom_label="FOMM",
    resolution_factor=resolution_factor,
    force=force)

def ccp4_map_from_solve_mtz (mtz_file, force=False, resolution_factor=1/3.0) :
  return ccp4_map_from_mtz(mtz_file=mtz_file,
    f_label="FP,SIGFP",
    phi_label="PHIB",
    fom_label="FOM",
    resolution_factor=resolution_factor,
    force=force)

def write_ccp4_map (sites_cart, unit_cell, map_data, n_real, file_name,
    buffer=10) :
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
    labels=flex.std_string(["mmtbx.utils.write_ccp4_map_box"]))

#-----------------------------------------------------------------------
# LABEL HANDLING

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
        if not labels.startswith("FOFC") and labels != "FWT,DELFWT" :
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
            hybrid_label.append(fom)
            all_labels.append(hybrid_label)
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
