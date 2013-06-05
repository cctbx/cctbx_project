from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.real_space_refine

from mmtbx.refinement import real_space
import sys, os
from iotbx import file_reader
from libtbx.utils import Sorry
import iotbx.phil
from mmtbx import monomer_library
from iotbx import reflection_file_reader
from libtbx import group_args
from cctbx import miller
import mmtbx.refinement.real_space
import random
from scitbx.array_family import flex
from libtbx.utils import user_plus_sys_time
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
import mmtbx.refinement.real_space.driver
import mmtbx.utils
import mmtbx.secondary_structure

if (1):
  random.seed(0)
  flex.set_random_seed(0)
  
# XXX TO BE MOVED OR RESTRUCTURED
def get_processed_pdb_object(rama_potential, log,
      raw_records=None, pdb_file_name=None):
  rama_potential = "emsley"
  assert [raw_records, pdb_file_name].count(None) == 1
  master_params = iotbx.phil.parse(
    input_string=mmtbx.monomer_library.pdb_interpretation.master_params_str,
    process_includes=True).extract()
  if(rama_potential is not None):
    master_params.peptide_link.ramachandran_restraints=True
    master_params.peptide_link.rama_potential=rama_potential
  master_params.nonbonded_weight=200
  mon_lib_srv = mmtbx.monomer_library.server.server()
  result = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = monomer_library.server.server(),
    ener_lib                 = monomer_library.server.ener_lib(),
    params                   = master_params,
    file_name                = pdb_file_name,
    raw_records              = raw_records,
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = log)
  return result

def get_geometry_restraints_manager(processed_pdb_file, xray_structure, log):
  sctr_keys=xray_structure.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  #
  sec_str = mmtbx.secondary_structure.process_structure(
    params             = None,
    processed_pdb_file = processed_pdb_file,
    tmp_dir            = os.getcwd(),
    log                = log,
    assume_hydrogens_all_missing=(not has_hd))
  sec_str.initialize(log=log)
  build_proxies = sec_str.create_hbond_proxies(
    log          = log,
    hbond_params = None)
  hbond_params = build_proxies.proxies
  hbond_params=None
  #
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    plain_pairs_radius           = 5,
    hydrogen_bond_proxies        = hbond_params,
    assume_hydrogens_all_missing = not has_hd)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)
  restraints_manager.crystal_symmetry = xray_structure.crystal_symmetry()
  return restraints_manager  
  
def process_inputs(args, log=None):
  got_pdb = False
  got_hkl = False
  got_map = False
  processed_pdb_file = None
  xray_structure     = None
  ccp4_map           = None
  map_coefficients   = None
  pdb_file_name      = None
  for arg in args :
    if(os.path.isfile(arg)):
      fro = file_reader.any_file(arg)
      if(iotbx.pdb.is_pdb_file(arg)):
        broadcast(m="Processing input PDB file:", log=log)
        processed_pdb_file = get_processed_pdb_object(pdb_file_name=fro.file_name,
          rama_potential=None, log = log)
        xray_structure = processed_pdb_file.xray_structure(show_summary = True)
        if(xray_structure is None):
          raise Sorry("Cannot not extract xray_structure.")
        pdb_file_name = fro.file_name
        got_pdb = True
      elif(fro.file_type == "ccp4_map"):
        broadcast(m="Processing input CCP4 map file: %s"%fro.file_name, log=log)
        ccp4_map = iotbx.ccp4_map.map_reader(file_name=fro.file_name)
        ccp4_map.show_summary(prefix="  ")
        got_map = True
      else:
        reflection_file = reflection_file_reader.any_reflection_file(
          file_name=arg, ensure_read_access=False)
        if(reflection_file.file_type() is not None):
          map_coefficients = reflection_file_utils.extract_miller_array_from_file(
            file_name = arg,
            #label     = params.label,
            type      = "complex",
            log       = log)
          got_hkl = True
    else:
      try: arg_phil = iotbx.phil.parse(arg)
      except Exception: raise Sorry("Bad parameter: %s"%arg)
      cmdline_phil.append(arg_phil)
  if(not got_pdb): raise Sorry("PDB file is needed.")
  if([got_map, got_hkl].count(True) == 2):
    raise Sorry("Only one, map or map coefficients, file is needed.")
  if([got_map, got_hkl].count(True) == 0):
    raise Sorry("Map or map coefficients file is needed.")
  #
  return group_args(
    processed_pdb_file = processed_pdb_file,
    xray_structure     = xray_structure,
    ccp4_map           = ccp4_map,
    map_coefficients   = map_coefficients,
    pdb_file_name      = pdb_file_name)

def validate_inputs(inputs):
  if([inputs.ccp4_map, inputs.map_coefficients].count(None) != 1):
    raise Sorry("Map or map coefficents must be provided.")
  if(inputs.map_coefficients is not None and not
     inputs.xray_structure.crystal_symmetry().is_similar_symmetry(
     inputs.map_coefficients.crystal_symmetry())):
    raise Sorry("Crystal symmetry mismatch: PDB model and map coefficients.")

def extract_target_map_data_and_crystal_gridding(inputs, resolution_factor):
  if(inputs.ccp4_map is not None):
    miller_array = mmtbx.utils.structure_factors_from_map(
      map_data          = inputs.ccp4_map.data.as_double(),
      unit_cell_lengths = inputs.ccp4_map.unit_cell_parameters[:3],
      n_real            = inputs.ccp4_map.unit_cell_grid,
      crystal_symmetry  = inputs.xray_structure.crystal_symmetry())
    target_map_data = inputs.ccp4_map.data.as_double()
    # XXX Dirty work-around to keep going, ask Ralf/Sacha for a better solution
    # XXX How to infer resolution from 3D map?
    a,b,c = inputs.ccp4_map.unit_cell_parameters[:3]
    nx,ny,nz = inputs.ccp4_map.unit_cell_grid
    d1,d2,d3 = a/nx/resolution_factor,b/ny/resolution_factor,c/nz/resolution_factor
    d_min_guess_from_map = min(d1,d2,d3)
    complete_set = miller.build_set(
      crystal_symmetry = inputs.xray_structure.crystal_symmetry(),
      anomalous_flag   = False,
      d_min            = d_min_guess_from_map)
    miller_array = complete_set.structure_factors_from_map(
      map            = inputs.ccp4_map.data.as_double(),
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = True)
    fft_map = miller_array.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    target_map_data = fft_map.real_map_unpadded()
  else:
    miller_array = inputs.map_coefficients
    fft_map = inputs.map_coefficients.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    target_map_data = fft_map.real_map_unpadded()
  return group_args(
    data             = target_map_data,
    miller_array     = miller_array,
    crystal_gridding = fft_map)

def broadcast(m, log):
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def format_usage_message(log):
  print >> log, "-"*79
  msg = """\
phenix.real_space_refine: tool for extensive real-space refinement of atomic
                          coordinates against provided map

Usage:
  phenix.real_space_refine model.pdb ccp4_formatted_map.map
  or
  phenix.real_space_refine model.pdb map.mtz
  or
  phenix.real_space_refine model.pdb map.mtz label=['2FOFCWT', 'PH2FOFCWT']

Feedback:
  PAfonine@lbl.gov or phenixbb@phenix-online.org
"""
  print >> log, msg
  print >> log, "-"*79

def run(args, log = None, resolution_factor=1./4):
  timer = user_plus_sys_time()
  if(log is None): log = sys.stdout
  format_usage_message(log = log)
  if(len(args)==0): return
  inputs = process_inputs(args = args, log = log)
  validate_inputs(inputs = inputs)
  broadcast(m="Creating geometry restraints:", log=log)
  geometry_restraints_manager = get_geometry_restraints_manager(
    processed_pdb_file = inputs.processed_pdb_file,
    xray_structure     = inputs.xray_structure,
    log                = log)
  target_map = extract_target_map_data_and_crystal_gridding(inputs=inputs,
    resolution_factor=resolution_factor)
  pdb_hierarchy = inputs.processed_pdb_file.all_chain_proxies.pdb_hierarchy
  time_startup = timer.elapsed()
  broadcast(m="Refinement start:", log=log)
  xray_structure_refined = mmtbx.refinement.real_space.driver.run(
    target_map                  = target_map,
    pdb_hierarchy               = pdb_hierarchy,
    xray_structure              = inputs.xray_structure,
    geometry_restraints_manager = geometry_restraints_manager,
    max_iterations = 100,
    macro_cycles   = 3)
  pdb_hierarchy.adopt_xray_structure(xray_structure_refined)
  pdb_hierarchy.write_pdb_file(file_name=inputs.pdb_file_name[:-4]+"_real_space_refined.pdb",
    crystal_symmetry = xray_structure_refined.crystal_symmetry())
  broadcast(m="Almost done... Run time infromation:", log=log)
  mmtbx.refinement.real_space.driver.show_time(
    external=[["  time_startup      : %6.3f", time_startup]])

if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  run(sys.argv[1:])
  print "Total time: %8.3f" % timer.elapsed()
  print "All done."
