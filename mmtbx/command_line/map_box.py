from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_box

import mmtbx.utils
from mmtbx.refinement import print_statistics
import iotbx.pdb
import libtbx.phil
from libtbx.utils import Sorry
import os, sys
from iotbx import reflection_file_utils
from cctbx import maptbx

master_phil = libtbx.phil.parse("""
  pdb_file = None
    .type = path
  map_coefficients_file = None
    .type = path
  label = None
    .type = str
  ccp4_map_file = None
    .type = str
  selection = all
    .type = str
  selection_radius = 3.0
    .type = float
  box_cushion = 3.0
    .type = float
  resolution_factor = 1./4
    .type = float
  output_format = *xplor *mtz *ccp4
    .type=choice(multi=True)
  output_file_name_prefix=None
    .type = str
  density_select = False
    .type = bool
    .help = Select boundaries based on where density is located.
  density_select_threshold = 0.05
    .type = float
    .help = Choose region where density is this fraction of maximum or greater
  ncs_file = None
    .type = path
    .help = NCS file (to be offset based on origin shift)

""")

def run(args, log=None):
  h = "phenix.map_box: extract box with model and map around selected atoms"
  if(log is None): log = sys.stdout
  print_statistics.make_header(h, out=log)
  default_message="""\

%s.

Usage:
  phenix.map_box model.pdb map_coefficients.mtz selection="chain A and resseq 1:10"

or

  phenix.map_box map.ccp4 density_select=True

Parameters:"""%h
  if(len(args) == 0):
    print default_message
    master_phil.show(prefix="  ")
    return
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_phil)
  params = inputs.params.extract()
  # PDB file
  if(len(inputs.pdb_file_names)!=1 and not params.density_select):
    raise Sorry("PDB file is needed unless density_select is set.")
  print_statistics.make_sub_header("pdb model", out=log)
  if len(inputs.pdb_file_names)>0:
    pdb_inp = iotbx.pdb.input(file_name=inputs.pdb_file_names[0])
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
  else:
    pdb_hierarchy=None
  # Map or map coefficients
  map_coeff = None
  if(inputs.ccp4_map is None):
    if(len(inputs.reflection_file_names)!=1):
      raise Sorry("Map or map coefficients file is needed.")
    map_coeff = reflection_file_utils.extract_miller_array_from_file(
      file_name = inputs.reflection_file_names[0],
      label     = params.label,
      type      = "complex",
      log       = log)
    fft_map = map_coeff.fft_map(resolution_factor=params.resolution_factor)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    map_or_map_coeffs_prefix=os.path.basename(
       inputs.reflection_file_names[0][:-4])
  else:
    print_statistics.make_sub_header("CCP4 map", out=log)
    ccp4_map = inputs.ccp4_map
    ccp4_map.show_summary(prefix="  ")
    map_data = ccp4_map.map_data()
    if inputs.ccp4_map_file_name.endswith(".ccp4"):
      map_or_map_coeffs_prefix=os.path.basename(
       inputs.ccp4_map_file_name[:-5])
    else:
      map_or_map_coeffs_prefix=os.path.basename(
       inputs.ccp4_map_file_name[:-4])
  #
  if len(inputs.pdb_file_names)>0:
    output_prefix=os.path.basename(inputs.pdb_file_names[0])[:-4]
  else:
    output_prefix=map_or_map_coeffs_prefix

  if pdb_hierarchy:
    xray_structure = pdb_hierarchy.extract_xray_structure(
      crystal_symmetry=inputs.crystal_symmetry)
    xray_structure.show_summary(f=log)
    #
    print_statistics.make_sub_header("atom selection", out=log)
    print >> log, "Selection string: selection='%s'"%params.selection
    selection = pdb_hierarchy.atom_selection_cache().selection(
      string = params.selection)
    print >> log, "  selects %d atoms from total %d atoms."%(selection.count(True),
      selection.size())
    #
    if params.density_select:
      print_statistics.make_sub_header(
      "extracting box around selected density and writing output files", out=log)
    else:
     print_statistics.make_sub_header(
      "extracting box around selected atoms and writing output files", out=log)
    #
    sites_cart_all = xray_structure.sites_cart()
    sites_cart = sites_cart_all.select(selection)
    selection = xray_structure.selection_within(
      radius    = params.selection_radius,
      selection = selection)
  else:
    xray_structure=None
    selection=None

  #
  box = mmtbx.utils.extract_box_around_model_and_map(
    xray_structure   = xray_structure,
    pdb_hierarchy    = pdb_hierarchy,
    map_data         = map_data.as_double(),
    box_cushion      = params.box_cushion,
    crystal_symmetry = inputs.crystal_symmetry,
    selection        = selection,
    density_select   = params.density_select,
    threshold        = params.density_select_threshold)

  if box.initial_shift:
    print >>log,"\nInitial coordinate shift will be (%.1f,%.1f,%.1f)\n" %(
      -box.initial_shift[0],-box.initial_shift[1],-box.initial_shift[2])
  if box.total_shift:
    print >>log,"Final coordinate shift: (%.1f,%.1f,%.1f)" %(
      -box.total_shift[0],-box.total_shift[1],-box.total_shift[2])

  print >>log,"Final cell dimensions: (%.1f,%.1f,%.1f)\n" %(
      box.box_crystal_symmetry.unit_cell().parameters()[:3])

  if box.pdb_outside_box_msg:
    print >> log, box.pdb_outside_box_msg

  if(params.output_file_name_prefix is None):
    file_name = "%s_box.pdb"%output_prefix
  else: file_name = "%s.pdb"%params.output_file_name_prefix
  if pdb_hierarchy is not None:
    box.write_pdb_file(file_name=file_name)

  if("ccp4" in params.output_format):
    if(params.output_file_name_prefix is None):
      file_name = "%s_box.ccp4"%output_prefix
    else: file_name = "%s.ccp4"%params.output_file_name_prefix
    print >> log, "writing map to CCP4 formatted file:   %s"%file_name
    box.write_ccp4_map(file_name=file_name)
  if("xplor" in params.output_format):
    if(params.output_file_name_prefix is None):
      file_name = "%s_box.xplor"%output_prefix
    else: file_name = "%s.xplor"%params.output_file_name_prefix
    print >> log, "writing map to X-plor formatted file: %s"%file_name
    box.write_xplor_map(file_name=file_name)
  if("mtz" in params.output_format):
    if(params.output_file_name_prefix is None):
      file_name = "%s_box.mtz"%output_prefix
    else: file_name = "%s.mtz"%params.output_file_name_prefix
    print >> log, "writing map coefficients to MTZ file: %s"%file_name
    if(map_coeff is not None): d_min = map_coeff.d_min()
    else:
      d_min = maptbx.d_min_from_map(map_data=box.map_box,
        unit_cell=box.xray_structure_box.unit_cell())
    box.map_coefficients(d_min=d_min,
      resolution_factor=params.resolution_factor, file_name=file_name)

  if params.ncs_file:

    if(params.output_file_name_prefix is None):
      output_ncs_file = "%s_box.ncs_spec"%output_prefix
    else: output_ncs_file = "%s.ncs_spec"%params.output_file_name_prefix
    print >>log,"\nOffsetting NCS in %s and writing to %s" %(
       params.ncs_file,output_ncs_file)
    from mmtbx.ncs.ncs import ncs
    ncs_object=ncs()
    ncs_object.read_ncs(params.ncs_file,log=log)
    ncs_object.display_all(log=log)
    if not ncs_object or ncs_object.max_operators()<1:
      print >>log,"Skipping...no NCS available"
    elif box.total_shift:
      print >>log,"Shifting origin for NCS"
      from scitbx.math import  matrix
      ncs_object=ncs_object.coordinate_offset(
       coordinate_offset=-1.0*matrix.col(box.total_shift))
      ncs_object.display_all(log=log)
    ncs_object.format_all_for_group_specification(
       file_name=output_ncs_file)
    box.ncs_object=ncs_object
  else:
    box.ncs_object=None
  print >> log
  return box

def select_box(params,map_data=None):
  # choose minimum and maximum points for a box based on where values are
  #  highest in the map
  frac_min=(0.21884765624999997, 0.4059863281250001, 0.066611328125)
  frac_max=(0.5963134765625, 0.754462890625, 0.9808056640625)

  return frac_min,frac_max

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
