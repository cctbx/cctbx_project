"""Create maps from PDB and MTZ files"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.maps

import mmtbx.maps
from scitbx.array_family import flex
import os, sys, random
import iotbx.pdb
from libtbx.utils import Sorry
from libtbx import runtime_utils
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from cctbx import crystal
from iotbx import extract_xtal_data

random.seed(0)
flex.set_random_seed(0)


legend = """
phenix.maps: a command line tool to compute various maps and save them in most
             of known formats.

How to run the command line version:

  1. Run phenix.maps without any arguments: just type phenix.maps in the command
     line and hit Enter. This will creare a parameter file called maps.params,
     which can be renamed if desired.

  2. Edit maps.params file to specify input/output file names, data labels and
     the desired maps. It is possible to request as many maps as desired. By
     default, the file maps.params specifies 5 maps to be created: 2mFo-DFc,
     2mFo-DFc with missing Fobs filled with DFcalc, mFo-DFc and anomalous
     difference maps will be output in MTZ format, and one 2mFo-DFc map will be
     output in CCP4 format.
     NOTE: the anomalous difference map will only be created if the input
     reflection data file contains Bijvoet maps (F+/F- or I+/I-).

  3. Run this command to compute requested maps: phenix.maps maps.params

Alternately, you may specify input files (and additional parameters) directly
on the command line:

  % phenix.maps model.pdb data.mtz

and it will automatically generate the default maps as described above.

Important Facts:

  - phenix.maps is available in PHENIX GUI.

  - The scope of parameters 'map_coefficients' defines the map that will be
    output as Fourier map coefficients. The scope of parameters 'map' defines
    the maps that will be output as CCP4 or X-plor format.

  - To create several maps: duplicate either 'map_coefficients' or 'map' or both
    scopes of parameters as many times as many maps is desired. Then edit each
    of them to define the maps.

  - A map is defined by specifying a map type using 'map_type' keyword available
    within each scope of parameters: 'map_coefficients' or 'map'. The general
    supported format for 'map_type' is: [p][m]Fo+[q][D]Fc[_filled]. For
    example: 2Fo-Fc, 2mFobs-DFcalc, 3Fobs-2Fmodel, Fo-Fc, mfobs-Dfcalc, anom,
    llg.  The 'map_type' parser will automatically recognize which map is
    requested.

  - The program creates as many files with CCP4 or X-plor formatted maps as
    is requested, and it creates only one MTZ formatted file with
    all Fourier map coefficients in it.

  - The CCP4 or X-plor formatted maps can be computed in the entire unit cell
    or around selected atoms only.

  - Twinning (if detected) will be accounted for automatically. This can be
    disabled by using "skip_twin_detection=True" keyword.

  - All arrays used in map calculation, for example: Fobs, Fmodel, Fcalc, Fmask,
    m, D, etc., can be output into a CNS or MTZ formatted reflection file.

  - For those who likes to experiment: bulk solvent correction and anisotropic
    scaling can be turned off, the data can be filtered by sigma and resolution.

  - For some map types certain 'map_coefficients' or 'map' scope parameters may
    not be applicable. For example, for "map_type=anomalous" the keywords
    "fill_missing_f_obs" and some other are not applicable.

  - For LLG map calculation, if you specify the wavelength any existing heavy
    atoms (P or heavier) will be modeled as anomalous scatterers using the
    theoretical values of f' and f''.
"""

default_params = """\
maps {
  map_coefficients {
    map_type = 2mFo-DFc
    format = *mtz phs
    mtz_label_amplitudes = 2FOFCWT
    mtz_label_phases = PH2FOFCWT
    fill_missing_f_obs = False
  }
  map_coefficients {
    map_type = 2mFo-DFc
    format = *mtz phs
    mtz_label_amplitudes = 2FOFCWT_fill
    mtz_label_phases = PH2FOFCWT_fill
    fill_missing_f_obs = True
  }
  map_coefficients {
    map_type = mFo-DFc
    format = *mtz phs
    mtz_label_amplitudes = FOFCWT
    mtz_label_phases = PHFOFCWT
    fill_missing_f_obs = False
  }
  map_coefficients {
    map_type = anomalous
    format = *mtz phs
    mtz_label_amplitudes = ANOM
    mtz_label_phases = PHANOM
  }
  map {
    map_type = 2mFo-DFc
    fill_missing_f_obs = False
    grid_resolution_factor = 1/4.
  }
}
"""

def analyze_input_params(params):
  # Analyze map_coefficients
  mcp = params.maps.map_coefficients
  i = 0
  while (i < len(mcp)):
    mcp_ = mcp[i]
    if (mcp_.map_type is None):
      del mcp[i]
      continue
    if(mmtbx.map_names(mcp_.map_type).anomalous):
      mcp_.fill_missing_f_obs = False
      mcp_.acentrics_scale = 2.0
      mcp_.centrics_pre_scale = 1.0
      mcp_.sharpening = False
      mcp_.sharpening_b_factor = None
    i += 1
  # Analyze maps
  mp = params.maps.map
  i = 0
  while (i < len(mp)):
    mp_ = mp[i]
    if (mp_.map_type is None):
      del mp[i]
      continue
    if(mmtbx.map_names(mp_.map_type).anomalous):
      mp_.fill_missing_f_obs = False
      mp_.acentrics_scale = 2.0
      mp_.centrics_pre_scale = 1.0
      mp_.sharpening = False
      mp_.sharpening_b_factor = None
    i += 1

def run(args, log = sys.stdout, use_output_directory=True,
    suppress_fmodel_output=False):
  print(legend, file=log)
  print("-"*79, file=log)
  master_params = mmtbx.maps.maps_including_IO_master_params()
  if(len(args)==0 or (len(args)==1 and args[0]=="NO_PARAMETER_FILE")):
    if(not (len(args)==1 and args[0]=="NO_PARAMETER_FILE")):
      parameter_file_name = "maps.params"
      print("Creating parameter file '%s' in the following directory:\n%s"%(
        parameter_file_name, os.path.abspath('.')), file=log)
      if(os.path.isfile(parameter_file_name)):
        msg="File '%s' exists already. Re-name it or move and run the command again."
        raise Sorry(msg%parameter_file_name)
      pfo = open(parameter_file_name, "w")
    else:
      pfo = log
      print("\nAll phenix.maps parameters::\n", file=pfo)
    master_params = master_params.fetch(iotbx.phil.parse(default_params))
    master_params.show(out = pfo, prefix = " ", expert_level=1)
    return
  processed_args = mmtbx.utils.process_command_line_args(
    args=args,
    log=log,
    master_params=master_params)
  working_phil = processed_args.params
  params = working_phil.extract()
  fmodel_data_file_format = params.maps.output.fmodel_data_file_format
  if (len(params.maps.map_coefficients) == 0) and (len(params.maps.map) == 0):
    print("No map input specified - using default map types", file=log)
    working_phil = master_params.fetch(sources=[working_phil,
        iotbx.phil.parse(default_params)])
    params = working_phil.extract()
  # XXX BUG - the extra fetch will always set fmodel_data_file_format to
  # mtz; this is probaby a low-level phil problem
  if (fmodel_data_file_format is None) or (suppress_fmodel_output):
    params.maps.output.fmodel_data_file_format = None
  analyze_input_params(params=params)
  have_phil_file_input = len(processed_args.phil_file_names) > 0
  if (len(processed_args.pdb_file_names) > 1):
    raise Sorry("Only one model file is allowed as input.")
  if ((params.maps.input.pdb_file_name is None) and
      (len(processed_args.pdb_file_names) == 1)):
    params.maps.input.pdb_file_name = processed_args.pdb_file_names[0]
  if(not os.path.isfile(str(params.maps.input.pdb_file_name))):
    raise Sorry(
      "model file is not given: maps.input.pdb_file_name=%s is not a file"%\
      str(params.maps.input.pdb_file_name))
  if ((params.maps.input.reflection_data.file_name is None) and
      (params.maps.input.reflection_data.r_free_flags.file_name is None) and
      (len(processed_args.reflection_file_names) == 1)):
    params.maps.input.reflection_data.file_name = \
      processed_args.reflection_file_names[0]
  print("FORMAT:", params.maps.output.fmodel_data_file_format, file=log)
  working_phil = master_params.format(python_object=params)
  print("-"*79, file=log)
  print("\nParameters to compute maps::\n", file=log)
  working_phil.show(out = log, prefix=" ")
  pdb_inp = iotbx.pdb.input(file_name = params.maps.input.pdb_file_name)
  # get all crystal symmetries
  cs_from_coordinate_files = [pdb_inp.crystal_symmetry_from_cryst1()]
  cs_from_reflection_files = []
  for rfn in [params.maps.input.reflection_data.file_name,
             params.maps.input.reflection_data.r_free_flags.file_name]:
    if(os.path.isfile(str(rfn))):
      try:
        cs_from_reflection_files.append(crystal_symmetry_from_any.extract_from(rfn))
      except KeyboardInterrupt: raise
      except RuntimeError: pass
  crystal_symmetry = None
  try :
    crystal_symmetry = crystal.select_crystal_symmetry(
      from_coordinate_files=cs_from_coordinate_files,
      from_reflection_files=cs_from_reflection_files)
  except AssertionError as e :
    if ("No unit cell and symmetry information supplied" in str(e)):
      raise Sorry("Missing or incomplete symmetry information.  This program "+
        "will only work with reflection file formats that contain both "+
        "unit cell and space group records, such as MTZ files.")
  #
  reflection_files = []
  reflection_file_names = []
  for rfn in [params.maps.input.reflection_data.file_name,
             params.maps.input.reflection_data.r_free_flags.file_name]:
    if(os.path.isfile(str(rfn))) and (not rfn in reflection_file_names):
      reflection_files.append(reflection_file_reader.any_reflection_file(
        file_name = rfn, ensure_read_access = False))
      reflection_file_names.append(rfn)
  reflection_file_server = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry   = True,
    reflection_files = reflection_files, #[],
    err              = log)
  #
  reflection_data_master_params = extract_xtal_data.data_and_flags_master_params(
    master_scope_name="reflection_data")
  reflection_data_input_params = processed_args.params.get(
    "maps.input.reflection_data")
  reflection_data_params = reflection_data_master_params.fetch(
    reflection_data_input_params).extract().reflection_data
  #
  determine_data_and_flags_result = extract_xtal_data.run(
    reflection_file_server = reflection_file_server,
    parameters             = reflection_data_params,
    keep_going             = True)
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=None
  print("-"*79, file=log)
  print("\nInput model file:", params.maps.input.pdb_file_name, file=log)
  pdb_hierarchy = pdb_inp.construct_hierarchy(set_atom_i_seq=True)
  atom_selection_manager = pdb_hierarchy.atom_selection_cache()
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry = crystal_symmetry)
  # apply omit selection
  if(params.maps.omit.selection is not None):
    omit_selection = atom_selection_manager.selection(
      string = params.maps.omit.selection)
    keep_selection = ~omit_selection
    xray_structure = xray_structure.select(selection = keep_selection)
    pdb_hierarchy = pdb_hierarchy.select(keep_selection)
    atom_selection_manager = pdb_hierarchy.atom_selection_cache()
  #
  mmtbx.utils.setup_scattering_dictionaries(
    scattering_table = params.maps.scattering_table,
    xray_structure   = xray_structure,
    d_min            = f_obs.d_min(),
    log              = log)
  if (params.maps.wavelength is not None):
    if (params.maps.scattering_table == "neutron"):
      raise Sorry("Wavelength parameter not supported when the neutron "+
        "scattering table is used.")
    xray_structure.set_inelastic_form_factors(
      photon=params.maps.wavelength,
      table="sasaki")
  xray_structure.show_summary(f = log, prefix="  ")
  print("-"*79, file=log)
  print("Bulk solvent correction and anisotropic scaling:", file=log)
  fmodel = mmtbx.utils.fmodel_simple(
    xray_structures         = [xray_structure],
    scattering_table        = params.maps.scattering_table,
    f_obs                   = f_obs,
    r_free_flags            = r_free_flags,
    outliers_rejection      = params.maps.input.reflection_data.outliers_rejection,
    skip_twin_detection     = params.maps.skip_twin_detection,
    bulk_solvent_correction = params.maps.bulk_solvent_correction,
    anisotropic_scaling     = params.maps.anisotropic_scaling)
  fmodel_info = fmodel.info()
  fmodel_info.show_rfactors_targets_scales_overall(out = log)
  print("-"*79, file=log)
  print("Compute maps.", file=log)
  # XXX if run from the Phenix GUI, the output directory parameter is actually
  # one level up from the current directory, and use_output_directory=False
  if (params.maps.output.directory is not None) and (use_output_directory):
    assert os.path.isdir(params.maps.output.directory)
    output_dir = params.maps.output.directory
  else :
    output_dir = os.getcwd()
  if params.maps.output.prefix is not None :
    file_name_base = os.path.join(output_dir,
      os.path.basename(params.maps.output.prefix))
  else :
    file_name_base = params.maps.input.pdb_file_name
    if(file_name_base.count(".")>0):
      file_name_base = file_name_base[:file_name_base.index(".")]
  xplor_maps = mmtbx.maps.compute_xplor_maps(
    fmodel                 = fmodel,
    params                 = params.maps.map,
    atom_selection_manager = atom_selection_manager,
    file_name_prefix       = None,
    file_name_base         = file_name_base,
    pdb_hierarchy          = pdb_hierarchy)
  cmo = mmtbx.maps.compute_map_coefficients(
    fmodel = fmodel,
    params = params.maps.map_coefficients,
    pdb_hierarchy = pdb_hierarchy,
    log = log)
  map_coeff_file_name = file_name_base+"_map_coeffs.mtz"
  r_free_flags_output = None
  if (params.maps.output.include_r_free_flags):
    r_free_flags_output = fmodel.r_free_flags().average_bijvoet_mates()
  write_mtz_file_result = cmo.write_mtz_file(file_name = map_coeff_file_name,
    r_free_flags=r_free_flags_output)
  if(params.maps.output.fmodel_data_file_format is not None):
    fmodel_file_name = file_name_base + "_fmodel." + \
      params.maps.output.fmodel_data_file_format
    print("Writing fmodel arrays (Fobs, Fcalc, m, ...) to %s file."%\
      fmodel_file_name, file=log)
    fmodel_file_object = open(fmodel_file_name,"w")
    fmodel.export(out = fmodel_file_object, format =
      params.maps.output.fmodel_data_file_format)
    fmodel_file_object.close()
  print("All done.", file=log)
  if (write_mtz_file_result):
    print("Map coefficients: %s" % map_coeff_file_name, file=log)
  for file_name in xplor_maps :
    print("CCP4 or XPLOR map: %s" % file_name, file=log)
  print("-"*79, file=log)
  return (map_coeff_file_name, xplor_maps)

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args),
      log=sys.stdout,
      use_output_directory=False,
      suppress_fmodel_output=True) # XXX bug fix

def validate_params(params, callback=None):
  if params.maps.input.pdb_file_name is None :
    raise Sorry("No model file defined.")
  elif params.maps.input.reflection_data.file_name is None :
    raise Sorry("No reflection file defined.")
  elif params.maps.input.reflection_data.labels is None :
    raise Sorry("No labels chosen for reflection data.")
  elif len(params.maps.map) == 0 and len(params.maps.map_coefficients) == 0 :
    raise Sorry("You have not requested any maps for output.")
  elif ((params.maps.output.directory is not None) and
        (not os.path.isdir(params.maps.output.directory))):
    raise Sorry(("The output directory %s does not exist; please choose a "+
      "valid directory, or leave this parameter blank.") %
      params.maps.output.directory)
  if (params.maps.wavelength is not None):
    if (params.maps.scattering_table == "neutron"):
      raise Sorry("Wavelength parameter not supported when the neutron "+
        "scattering table is used.")
  validate_map_params(params.maps)
  # TODO double-check this - can we get None by accident in GUI?
  #for map_coeffs in params.maps.map_coefficients :
  #  if (map_coeffs.map_type is None):
  #    raise Sorry("One or more map coefficients is missing a map type "+
  #      "definition.")
  return True

def validate_map_params(params):
  from mmtbx import map_names
  labels = []
  for map_coeffs in params.map_coefficients :
    if (map_coeffs.map_type is not None):
      try :
        decode_map = map_names(map_coeffs.map_type)
      except RuntimeError as e :
        raise Sorry(str(e))
      f = map_coeffs.mtz_label_amplitudes
      phi = map_coeffs.mtz_label_phases
      if (f in labels) or (phi in labels):
        raise Sorry(("The map coefficients with MTZ labels %s,%s duplicates at"+
          " least one previously defined label.  You may output multiple sets "+
          "of coefficients with the same map type, but the column labels must "+
          "be unique.") % (f, phi))
      elif (None in [f, phi]):
        raise Sorry("Please specify both MTZ column labels for map_type '%s'."%
          map_coeffs.map_type)
      labels.extend([f,phi])
  if (hasattr(params, "map")):
    for map in params.map :
      if (map.grid_resolution_factor > 0.5):
        # XXX can't we enforce this in phil?
        raise Sorry("The grid resolution factor for CCP4 and X-PLOR maps "+
          "must be 0.5 or less.")
  return True

def finish_job(results):
  (mtz_file, map_files) = results
  output_files = []
  if mtz_file is not None and os.path.isfile(mtz_file):
    output_files.append((mtz_file, "MTZ file"))
  for map_file in map_files :
    if os.path.isfile(map_file):
      output_files.append((map_file, "XPLOR map"))
  return (output_files, [])

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

