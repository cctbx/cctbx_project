
# XXX: these functions are utilities used by the Phenix GUI for quickly
# converting data and map formats on the fly.  the main use of this is
# in phenix.refine and any program that opens maps in PyMOL.
# tested as part of Phenix regression tests (since it requires files)

from __future__ import division
import libtbx.phil
from libtbx.math_utils import ifloor, iceil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import os
import re
import sys

#-----------------------------------------------------------------------
# MAP COEFFICIENT MANIPULATION

class fast_maps_from_hkl_file (object) :
  def __init__ (self,
                file_name,
                xray_structure,
                f_label=None,
                map_out=None,
                log=sys.stdout,
                auto_run=True,
                quiet=False,
                ) :
    adopt_init_args(self, locals())
    if f_label is None and not quiet:
      print >> log, "      no label for %s, will try default labels" % \
        os.path.basename(file_name)
    f_obs = None
    fallback_f_obs = []
    default_labels = ["F,SIGF","FOBS,SIGFOBS","F(+),SIGF(+),F(-),SIGF(-)",
      "FOBS_X"]
    all_labels = []
    from iotbx import file_reader
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
      elif miller_array.is_xray_amplitude_array() :
        fallback_f_obs.append(miller_array)
    if f_obs is None :
      if (len(fallback_f_obs) == 1) and (f_label is None) :
        f_obs = fallback_f_obs[0]
      else :
        raise Sorry(("Couldn't find %s in %s.  Please specify valid "+
          "column labels (possible choices: %s)") % (f_label, file_name,
            " ".join(all_labels)))
    self.f_obs = f_obs
    self.log = None
    if auto_run :
      self.run()

  def get_maps_from_fmodel(self):
    import mmtbx.utils
    from scitbx.array_family import flex
    f_obs = self.f_obs
    r_free_flags = f_obs.array(data=flex.bool(f_obs.data().size(),False))
    fmodel = mmtbx.utils.fmodel_simple(
      xray_structures=[self.xray_structure],
      f_obs=f_obs,
      r_free_flags=r_free_flags)
    (f_map, df_map) = get_maps_from_fmodel(fmodel, use_filled=True)
    return f_map, df_map

  def run (self) :
    (f_map, df_map) = self.get_maps_from_fmodel()
    if self.map_out is None :
      self.map_out = os.path.splitext(self.file_name)[0] + "_map_coeffs.mtz"
    write_map_coeffs(f_map, df_map, self.map_out)

def get_maps_from_fmodel (fmodel, use_filled=False) :
  map_manager = fmodel.electron_density_map(fill_missing_f_obs=use_filled,
                                            fill_mode="dfmodel")
  fwt_coeffs = map_manager.map_coefficients(map_type = "2mFo-DFc")
  if fwt_coeffs.anomalous_flag() :
    fwt_coeffs = fwt_coeffs.average_bijvoet_mates()
  delfwt_coeffs = map_manager.map_coefficients(map_type = "mFo-DFc")
  if delfwt_coeffs.anomalous_flag() :
    delfwt_coeffs = delfwt_coeffs.average_bijvoet_mates()
  return (fwt_coeffs, delfwt_coeffs)

def write_map_coeffs (fwt_coeffs, delfwt_coeffs, file_name) :
  import iotbx.mtz
  decorator = iotbx.mtz.label_decorator(phases_prefix="PH")
  mtz_dataset = fwt_coeffs.as_mtz_dataset(
    column_root_label="2FOFCWT",
    label_decorator=decorator)
  mtz_dataset.add_miller_array(
    miller_array=delfwt_coeffs,
    column_root_label="FOFCWT",
    label_decorator=decorator)
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name=file_name)
  del mtz_object

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
  import iotbx.xplor.map
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
  from iotbx import file_reader
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
  import mmtbx.maps
  from iotbx import file_reader
  map_phil = libtbx.phil.parse("""map.file_name = %s """ % output_file)
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

def convert_map_coefficients (map_coefficients,
                              mtz_file,
                              pdb_file=None,
                              grid_resolution_factor=0.33) :
  assert os.path.isfile(mtz_file)
  map_files = []
  from iotbx import file_reader
  mtz_in = file_reader.any_file(mtz_file)
  mtz_in.assert_file_type("hkl")
  xray_structure = None
  if (pdb_file is not None) :
    pdb_in = file_reader.any_file(pdb_file)
    pdb_in.assert_file_type("pdb")
    xray_structure = pdb_in.file_object.xray_structure_simple()
  for map in map_coefficients :
    array_label = map.mtz_label_amplitudes + "," + map.mtz_label_phases
    map_array = None
    for miller_array in mtz_in.file_server.miller_arrays :
      if (miller_array.info().label_string() == array_label) :
        map_array = miller_array
        break
    if (miller_array is None) :
      print "Can't find %s" % array_label
      continue
    if mtz_file.endswith("_map_coeffs.mtz") :
      base_file = re.sub("_map_coeffs.mtz", "", mtz_file)
    else :
      base_file, ext = os.path.splitext(mtz_file)
    output_file = base_file + "_%s.ccp4" % map.map_type
    ccp4_map_from_coeffs(
      miller_array=map_array,
      output_file=map_file,
      xray_structure=xray_structure,
      grid_resolution_factor=grid_resolution_factor)
    map_files.append((map_file, map.map_type))
  return map_files

def ccp4_map_from_coeffs (miller_array, output_file, pdb_file=None,
    xray_structure=None, grid_resolution_factor=0.33) :
  assert miller_array.is_complex_array()
  import mmtbx.maps
  from iotbx import file_reader
  map_phil = libtbx.phil.parse("""map.file_name = %s """ % output_file)
  if (xray_structure is None) :
    if (pdb_file is not None) and (os.path.isfile(pdb_file)) :
      pdb_in = file_reader.any_file(pdb_file)
      assert (pdb_in.file_type == "pdb")
      xray_structure = pdb_in.file_object.xray_structure_simple()
  sites_cart = None
  if (xray_structure is not None) :
    sites_cart = xray_structure.sites_cart()
  fft_map = miller_array.fft_map(resolution_factor=grid_resolution_factor)
  fft_map.apply_sigma_scaling()
  write_ccp4_map(
    sites_cart=sites_cart,
    unit_cell=miller_array.unit_cell(),
    map_data=fft_map.real_map(),
    n_real=fft_map.n_real(),
    file_name=output_file)

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
    labels=flex.std_string(["mmtbx.utils.write_ccp4_map_box"]))

class write_ccp4_maps_wrapper (object) :
  def __init__ (self, pdb_hierarchy, map_coeffs, output_files,
      resolution_factor) :
    adopt_init_args(self, locals())

  def run (self) :
    sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    for map_coeffs, file_name in zip(self.map_coeffs, self.output_files) :
      if (map_coeffs is None) :
        continue
      fft_map = map_coeffs.fft_map(resolution_factor=self.resolution_factor)
      write_ccp4_map(sites_cart,
        unit_cell=map_coeffs.unit_cell(),
        map_data=fft_map.real_map(),
        n_real=fft_map.n_real(),
        file_name=file_name)

# XXX backwards compatibility
# TODO remove these ASAP once GUI is thoroughly tested
def extract_map_coeffs (*args, **kwds) :
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.extract_map_coeffs(*args, **kwds)

def map_coeffs_from_mtz_file (*args, **kwds) :
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.map_coeffs_from_mtz_file(*args, **kwds)

def extract_phenix_refine_map_coeffs (*args, **kwds) :
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.extract_phenix_refine_map_coeffs(*args,
    **kwds)

def get_map_coeff_labels (*args, **kwds) :
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.get_map_coeff_labels(*args, **kwds)

def get_map_coeffs_for_build (server) :
  return get_map_coeff_labels(server, build_only=True)

def format_map_coeffs_for_resolve (*args, **kwds) :
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.format_map_coeffs_for_resolve(*args,
    **kwds)

def decode_resolve_map_coeffs (*args, **kwds) :
  import iotbx.gui_tools.reflections
  return iotbx.gui_tools.reflections.decode_resolve_map_coeffs(*args, **kwds)
