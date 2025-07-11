"""See phenix.maps"""
# LIBTBX_SET_DISPATCHER_NAME phenix.mtz2map
# LIBTBX_SET_DISPATCHER_NAME phenix.fft
# LIBTBX_SET_DISPATCHER_NAME mmtbx.fft

# TODO: remove R-free set from map coefficients?

from __future__ import absolute_import, division, print_function
from mmtbx.maps import utils
import iotbx.map_tools
import iotbx.phil
from iotbx import file_reader
import iotbx.pdb
from libtbx import runtime_utils
import libtbx.phil
from libtbx.utils import Sorry, Usage
import sys, os
import iotbx.phil

master_phil = iotbx.phil.parse("""
mtz_file = None
  .type = path
  .short_caption = MTZ file
  .style = bold file_type:hkl OnChange:extract_map_coeffs_for_fft
  .help = MTZ file containing map coefficients
pdb_file = None
  .type = path
  .multiple = True
  .short_caption = Model file
  .help = Model file around which to draw the map.  If not supplied, the map \
    will fill the unit cell.
labels = None
  .type = strings
  .multiple = True
  .help = Map column labels.  Common examples are "2FOFCWT,PH2FOFCWT" and \
    "FP,SIGFP" "PHIM" "FOMM".  If left blank, all maps present in the input \
    file will be used.
buffer = 5.0
  .type = float
  .short_caption = Region padding
  .help = Extra padding (in Angstroms) around the selected region (or unit \
    cell)
selection = None
  .type = str
  .input_size = 400
  .short_caption = Atom selection
  .help = Atom selection around which to draw map (plus buffer).  If left \
    blank, the entire model file will be used.
d_min = None
  .type = float
  .help = High-resolution cutoff
d_max = None
  .type = float
  .help = Low-resolution cutoff
grid_resolution_factor = 0.25
  .type = float
  .help = Grid spacing (multiplied by the high-resolution limit)
gridding = None
  .type = ints
  .help = Gridding
scale = *sigma volume
  .type = choice(multi=False)
  .expert_level = 1
  .short_caption = Map scaling
  .help = Scaling method for map values
output {
  directory = None
    .type = path
    .short_caption = Output directory
    .help = Output directory (defaults to current)
    .style = output_dir bold
  prefix = None
    .type = str
    .input_size = 400
    .short_caption = Output file prefix
    .help = Output file prefix (defaults to the MTZ file name base)
  include scope libtbx.phil.interface.tracking_params
  format = xplor *ccp4 dsn6
    .type = choice
    .caption = XPLOR CCP4 DSN6
  extension = *Auto ccp4 xplor map dsn6
    .type = choice
}
r_free_flags {
  remove = False
    .type = bool
    .short_caption = Remove R-free set from map
  file_name = None
    .type = path
    .short_caption = R-free flags
    .style = OnChange:extract_r_free_flags_for_fft
  label = None
    .type = str
    .input_size = 120
    .short_caption = R-free label
    .style = renderer:draw_rfree_label_widget
  test_flag_value = None
    .type = int
}
include_fmodel = False
  .type = bool
  .short_caption = Include F(model) if present
show_maps = False
  .type = bool
""", process_includes=True)
master_params = master_phil

def find_array(miller_arrays, labels):
  for array in miller_arrays :
    if array.info().label_string() == labels :
      return array
  return None

def run(args, log=sys.stdout, run_in_current_working_directory=False):
  import iotbx.phil # FIXME this should not be necessary!
  pdb_file = None
  mtz_file = None
  input_phil = []
  if len(args) == 0 :
    print("Parameter syntax:", file=log)
    master_phil.show(out=log, prefix="  ")
    raise Usage("phenix.mtz2map [mtz_file] [pdb_file] [param_file] " +
      "[--show_maps]")
  input_objects = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="pdb_file",
    reflection_file_def="mtz_file")
  params = input_objects.work.extract()

  if params.mtz_file is None:
    raise Sorry("Please specify an MTZ file containing map coefficients.")
  if (not run_in_current_working_directory):
    if params.output.directory is None :
      params.output.directory = os.getcwd()
    if (not os.path.exists(params.output.directory)):
      os.makedirs(params.output.directory)
    elif (not os.path.isdir(params.output.directory)):
      raise Sorry("The specified output path '%s' is not a directory!" %
        params.output.directory)
    output_dir = params.output.directory
  else :
    output_dir = os.getcwd()
  if params.output.prefix is None :
    params.output.prefix = os.path.splitext(
      os.path.basename(params.mtz_file))[0]
  mtz_file = file_reader.any_file(params.mtz_file,
    force_type="hkl",
    raise_sorry_if_errors=True)
  all_labels=[]
  if params.show_maps or len(params.labels) == 0 :
    map_labels = utils.get_map_coeff_labels(mtz_file.file_server,
      exclude_anomalous=False,
      exclude_fmodel=not params.include_fmodel,
      keep_array_labels=True)
    if (not params.include_fmodel):
      all_labels = utils.get_map_coeff_labels(mtz_file.file_server,
        exclude_anomalous=False,
        exclude_fmodel=False,
        keep_array_labels=True)
    else :
      all_labels = map_labels
    if len(all_labels) > 0 :
      print("Available map coefficients in this MTZ file:", file=log)
      for labels in all_labels :
        if isinstance(labels, str):
          labels_list = [labels]
        else :
          labels_list = labels
        if (not labels in map_labels):
          extra = " (skipping, add include_fmodel=True to include)"
        else :
          extra = ""
          params.labels.append(labels_list)
        print("  %s%s" % (labels_list, extra), file=log)
    else :
      raise Sorry("No map coefficients found in this MTZ file.")
    if params.show_maps : return False
  pdb_file = None
  if len(params.pdb_file) > 0:
    pdb_file = iotbx.pdb.input(params.pdb_file[0])
  miller_arrays = mtz_file.file_object.as_miller_arrays()
  r_free_flags = None
  if params.r_free_flags.remove :
    if params.r_free_flags.file_name is not None :
      rfree_file = file_reader.any_file(params.r_free_flags.file_name,
        force_type="hkl")
    else :
      rfree_file = mtz_file
    raw_array, flag_value = rfree_file.file_server.get_r_free_flags(
      file_name=rfree_file.file_name,
      label=params.r_free_flags.label,
      test_flag_value=params.r_free_flags.test_flag_value,
      disable_suitability_test=False,
      parameter_scope="r_free_flags")
    r_free_flags = raw_array.array(
      data=raw_array.data()==flag_value).map_to_asu().average_bijvoet_mates()
  sites_cart = None
  if pdb_file is not None :
    pdb_hierarchy = pdb_file.construct_hierarchy()
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
    if params.selection is not None :
      selection_cache = pdb_hierarchy.atom_selection_cache()
      selection = selection_cache.selection(params.selection)
      sites_cart = sites_cart.select(selection)
      if (len(sites_cart) == 0):
        raise Sorry("No atoms found matching the specified selection.")
  else :
    print("No model input - will output map(s) in unit cell.", file=log)
  file_info = []
  suffixes = []
  for i, map_labels in enumerate(params.labels):
    map_coeffs = None
    if (len(map_labels) == 1):
      map_coeffs = find_array(miller_arrays, map_labels[0])
      if (map_coeffs is None):
        all_labels = utils.get_map_coeff_labels(mtz_file.file_server,
          keep_array_labels=True,
          exclude_anomalous=False,
          exclude_fmodel=not params.include_fmodel)
        labels_out = []
        if len(all_labels) > 0 :
          for labels in all_labels :
            if isinstance(labels, str):
              labels = [labels]
            labels_out.append("  " + " ".join(labels))
          raise Sorry(("No map coefficients found with labels %s.  Possible "+
            "choices are:\n%s") % (map_labels[0], "\n".join(labels_out)))
        else :
          raise Sorry(("No map coefficients found with labels %s; this file "+
            "does not appear to contain any other map coefficients!") %
            map_labels[0])
    else :
      if len(map_labels) == 2 :
        map_labels.append(None)
      (f, phi, fom) = utils.extract_map_coeffs(miller_arrays=miller_arrays,
        f_lab=map_labels[0],
        phi_lab=map_labels[1],
        fom_lab=map_labels[2])
      assert f.is_xray_amplitude_array()
      assert phi.is_real_array()
      assert (fom is None) or (fom.is_real_array())
      map_coeffs = iotbx.map_tools.combine_f_phi_and_fom(f=f, phi=phi, fom=fom)
    print("Processing map: %s" % " ".join(map_labels), file=log)
    assert map_coeffs.is_complex_array()
    map_coeffs = map_coeffs.map_to_asu().average_bijvoet_mates()
    map_coeffs = map_coeffs.resolution_filter(d_min=params.d_min,
      d_max=params.d_max)
    if (r_free_flags is not None):
      map_coeffs, flags = map_coeffs.common_sets(other=r_free_flags)
      print("  removing %d R-free flagged reflections" % \
        flags.data().count(True), file=log)
      map_coeffs = map_coeffs.select(~(flags.data()))
    from cctbx import maptbx

    if params.gridding:
      from cctbx.maptbx import crystal_gridding
      cg=crystal_gridding(
        unit_cell=map_coeffs.crystal_symmetry().unit_cell(),
        space_group_info=
           map_coeffs.crystal_symmetry().space_group_info(),
        pre_determined_n_real=params.gridding)
    else:
      cg=None
    map = map_coeffs.fft_map(resolution_factor=params.grid_resolution_factor,
      symmetry_flags=maptbx.use_space_group_symmetry,
      crystal_gridding=cg)
    if params.scale == "sigma" :
      print("  applying sigma-scaling", file=log)
      map.apply_sigma_scaling()
    elif params.scale == "volume" :
      print("  applying volume-scaling", file=log)
      map.apply_volume_scaling()
    suffix = None
    if map_labels == ["FP,SIGFP", "PHIM", "FOMM"] :
      suffix = ""
    elif map_labels[0].startswith("2FOFCWT"):
      if map_labels[0].endswith("no_fill"):
        suffix = "_2mFo-DFc_no_fill"
      else :
        suffix = "_2mFo-DFc"
    elif map_labels[0].startswith("FOFCWT"):
      suffix = "_mFo-DFc"
    elif map_labels[0] == "FWT,PHWT" : # refmac
      suffix = "_2mFo-DFc"
    elif map_labels[0] == "DELFWT,PHDELWT" : # refmac
      suffix = "_mFo-DFc"
    elif map_labels[0].startswith("ANOM"):
      suffix = "_anom"
    elif (map_labels[0].startswith("LLG")):
      suffix = "_llg"
    elif (map_labels[0].startswith("FMODEL") or
          map_labels[0].startswith("F-model")):
      suffix = "_fmodel"
    elif (map_labels[0].startswith("FC")):
      suffix = "_fcalc"
    else :
      suffix = "_%d" % (i+1)
    if ("_no_fill" in map_labels[0]):
      suffix += "_no_fill"
    elif ("_fill" in map_labels[0]):
      suffix += "_filled"
    # check for duplicate suffixes, append a number if necessary
    if (suffix in suffixes):
      suffixes.append(suffix)
      n = suffixes.count(suffix)
      suffix += "_%d" % n
    else :
      suffixes.append(suffix)
    format = params.output.format
    if params.output.extension == "Auto" :
      if format == "ccp4" :
        extension = "ccp4"
      elif format == "dsn6" :
        extension = "omap"
      else :
        extension = "xplor"
    else :
      extension = params.output.extension
      if format == "xplor" and not extension in ["xplor", "map"] :
        raise Sorry("%s is not an appropriate extension for Xplor maps." %
          extension)
      elif format == "ccp4" and not extension in ["ccp4", "map"] :
        raise Sorry("%s is not an appropriate extension for CCP4 maps." %
          extension)
    map_file_name = os.path.join(output_dir,
      params.output.prefix + suffix + "." + extension)
    if format == "xplor" :
      utils.write_xplor_map(
        sites_cart=sites_cart,
        unit_cell=map_coeffs.unit_cell(),
        map_data=map.real_map(),
        n_real=map.n_real(),
        file_name=map_file_name,
        buffer=params.buffer)
      file_info.append((map_file_name, "XPLOR map"))
    elif (format == "dsn6"):
      if (sites_cart is not None):
        import iotbx.map_tools
        iotbx.map_tools.write_dsn6_map(
          sites_cart=sites_cart,
          unit_cell=map_coeffs.unit_cell(),
          map_data=map.real_map(),
          n_real=map.n_real(),
          file_name=map_file_name,
          buffer=params.buffer)
      else :
        map.as_dsn6_map(file_name=map_file_name)
    else :
      if sites_cart is not None :
        import iotbx.map_tools
        iotbx.map_tools.write_ccp4_map(
          sites_cart=sites_cart,
          unit_cell=map_coeffs.unit_cell(),
          map_data=map.real_map(),
          n_real=map.n_real(),
          file_name=map_file_name,
          buffer=params.buffer)
      else :
        map.as_ccp4_map(file_name=map_file_name)
      file_info.append((map_file_name, "CCP4 map"))
    print("  wrote %s" % map_file_name, file=log)
  return file_info

def finish_job(result):
  return (result, []) # XXX result is already a file name/desc. list

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.makedirs(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args),
      log=sys.stdout,
      run_in_current_working_directory=True)

def validate_params(params):
  if params.mtz_file is None:
    raise Sorry("No MTZ file was provided.")

if __name__ == "__main__" :
  run(sys.argv[1:])

