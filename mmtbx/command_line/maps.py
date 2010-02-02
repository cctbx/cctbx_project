# LIBTBX_SET_DISPATCHER_NAME phenix.maps

import mmtbx.maps
import os, sys
import iotbx.pdb
from libtbx.utils import Sorry
from libtbx import runtime_utils
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any

legend = """
phenix.maps: a command line tool to compute various maps.

How to run:

  1. Run phenix.maps without any arguments: just type phenix.maps in the command
     line and hit Enter. This will creare a parameter file called maps.params,
     which can be renamed if desired.

  2. Edit maps.params file to specify input/output files and the desired maps.
     It is possible to request as many maps as desired.

  3. Run this command to compute requested maps: phenix.maps maps.params

Remarks:

  - The scope of parameters 'map_coefficients' defines the map that will be
    output as Fourier map coefficients. The scope of parameters 'map' defines
    the map that will be output as X-plor formatted map.

  - To create several maps: duplicate either 'map_coefficients' or 'map' or both
    scopes of parameters as many times as many maps is desired. Then edit each
    of them to define the maps.

  - A map is defined by specifying a map type using 'map_type' keyword available
    within each scope of parameters: 'map_coefficients' or 'map'. The general
    supported format for 'map_type' is: [p][m]Fo+[q][D]Fc[kick][filled]. For
    example: 2Fo-Fc, 2mFobs-DFcalc, 3Fobs-2Fmodel, Fo-Fc, mfobs-Dfcalc, anom.
    The 'map_type' parser will automatically recognize which map is requested.

  - The program creates as many files with X-plor formatted maps as many X-plor
    formatted maps is requested, and it creates only one MTZ formatted file with
    all Fourier map coefficients in it.

  - The X-plor formatted map can be computed in the entire unit cell or around
    selected atoms only.

  - Kick maps and missing Fobs filling is done (if requested) as described in
    Adams et al. (2010). Acta Cryst. D66, 213-221.

  - Twinning (if detected) will be accounted for automatically.

  - All arrays used in map calculation, for example: Fobs, Fmodel, Fcalc, Fmask,
    m, D, etc., can be output into a CNS or MTZ formatted reflection file.
"""

def get_atom_selection_manager(pdb_inp):
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  return pdb_hierarchy.atom_selection_cache()

def run(args, log = sys.stdout):
  print >> log, legend
  print >> log, "-"*79
  if(len(args) == 0):
    parameter_file_name = "maps.params"
    print >> log, "Creating parameter file %s in the following directory:\n%s"%(
      parameter_file_name, os.path.abspath('.'))
    pfo = open(parameter_file_name, "w")
    mmtbx.maps.maps_including_IO_master_params().show(out = pfo)
    pfo.close()
    print >> log, "-"*79
    return
  processed_args = mmtbx.utils.process_command_line_args(args = args, log = log,
    master_params = mmtbx.maps.maps_including_IO_master_params())
  print >> log, "-"*79
  print >> log, "\nParameters to compute maps::\n"
  processed_args.params.show(out = log, prefix=" ")
  params = processed_args.params.extract()
  if(not os.path.isfile(str(params.maps.input.pdb_file_name))):
    raise Sorry(
      "PDB file is not given: maps.input.pdb_file_name=%s is not a file"%\
      str(params.maps.input.pdb_file_name))
  print >> log, "-"*79
  print >> log, "\nInput PDB file:", params.maps.input.pdb_file_name
  pdb_inp = iotbx.pdb.input(file_name = params.maps.input.pdb_file_name)
  cryst1 = pdb_inp.crystal_symmetry_from_cryst1()
  if(cryst1 is None):
    raise Sorry("CRYST1 record in input PDB file is incomplete or missing.")
  else:
    if([cryst1.unit_cell(), cryst1.space_group_info()].count(None) != 0):
      raise Sorry("CRYST1 record in input PDB file is incomplete or missing.")
  xray_structure = pdb_inp.xray_structure_simple()
  xray_structure.show_summary(f = log, prefix="  ")
  print >> log, "-"*79
  crystal_symmetries = []
  crystal_symmetries.append(xray_structure.crystal_symmetry())
  reflection_files = []
  for rfn in [params.maps.input.reflection_data.file_name,
             params.maps.input.reflection_data.r_free_flags.file_name]:
    if(os.path.isfile(str(rfn))):
      reflection_files.append(reflection_file_reader.any_reflection_file(
        file_name = rfn, ensure_read_access = False))
      try:
        crystal_symmetries.append(crystal_symmetry_from_any.extract_from(rfn))
      except KeyboardInterrupt: raise
      except RuntimeError: pass
  if(len(crystal_symmetries)>1):
    cs0 = crystal_symmetries[0]
    for cs in crystal_symmetries:
     if(cs.unit_cell() is not None):
       if(not cs0.is_similar_symmetry(cs)):
         raise Sorry("Crystal symmetry mismatch between different files.")
  reflection_file_server = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetries[0],
    force_symmetry   = True,
    reflection_files = reflection_files,
    err              = log)
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server  = reflection_file_server,
    parameters              = params.maps.input.reflection_data,
    data_parameter_scope    = "maps.input.reflection_data",
    flags_parameter_scope   = "maps.input.reflection_data.r_free_flags",
    data_description        = "Reflection data",
    keep_going              = True,
    log                     = log)
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=None
  print >> log, "-"*79
  print >> log, "Bulk solvent correction and anisotropic scaling:"
  fmodel = mmtbx.utils.fmodel_simple(
    xray_structures = [xray_structure],
    f_obs           = f_obs,
    r_free_flags    = r_free_flags)
  fmodel_info = fmodel.info()
  fmodel_info.show_rfactors_targets_scales_overall(out = log)
  print >> log, "-"*79
  print >> log, "Compute maps."
  atom_selection_manager = get_atom_selection_manager(pdb_inp = pdb_inp)
  if params.maps.output.directory is not None :
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
    file_name_prefix       = params.maps.output.prefix,
    file_name_base         = file_name_base)
  cmo = mmtbx.maps.compute_map_coefficients(
    fmodel = fmodel,
    params = params.maps.map_coefficients)
  map_coeff_file_name = file_name_base+"_map_coeffs.mtz"
  #if(params.maps.output.prefix is not None and len(params.maps.output.prefix)>0):
  #  map_coeff_file_name = params.maps.output.prefix + "_" + map_coeff_file_name
  cmo.write_mtz_file(file_name = map_coeff_file_name)
  if(params.maps.output.fmodel_data_file_format is not None):
    fmodel_file_name = file_name_base + "_fmodel." + \
      params.maps.output.fmodel_data_file_format
    print >> log, "Writing fmodel arrays (Fobs, Fcalc, m, ...) to %s file."%\
      fmodel_file_name
    fmodel_file_object = open(fmodel_file_name,"w")
    fmodel.export(out = fmodel_file_object, format =
      params.maps.output.fmodel_data_file_format)
    fmodel_file_object.close()
  print >> log, "All done."
  print >> log, "-"*79
  return (map_coeff_file_name, xplor_maps)

class launcher (runtime_utils.simple_target) :
  def __call__ (self) :
    os.chdir(self.output_dir)
    return run(args=list(self.args), log=sys.stdout)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
