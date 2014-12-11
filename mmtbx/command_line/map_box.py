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
""")

def run(args, log=None):
  h = "phenix.map_box: extract box with model and map around selected atoms"
  if(log is None): log = sys.stdout
  print_statistics.make_header(h, out=log)
  default_message="""\

%s.

Usage:
  phenix.map_box model.pdb map_coefficients.mtz selection="chain A and resseq 1:10"

Parameters:"""%h
  if(len(args) == 0):
    print default_message
    master_phil.show(prefix="  ")
    return
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_phil)
  params = inputs.params.extract()
  # PDB file
  if(len(inputs.pdb_file_names)!=1):
    raise Sorry("PDB file is needed.")
  print_statistics.make_sub_header("pdb model", out=log)
  pdb_inp = iotbx.pdb.input(file_name=inputs.pdb_file_names[0])
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
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
  else:
    print_statistics.make_sub_header("CCP4 map", out=log)
    ccp4_map = inputs.ccp4_map
    ccp4_map.show_summary(prefix="  ")
    map_data = ccp4_map.data
  #
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
  print_statistics.make_sub_header(
    "extracting box around selected atoms and writing output files", out=log)
  #
  sites_cart_all = xray_structure.sites_cart()
  sites_cart = sites_cart_all.select(selection)
  selection = xray_structure.selection_within(
    radius    = params.selection_radius,
    selection = selection)
  #
  box = mmtbx.utils.extract_box_around_model_and_map(
    xray_structure   = xray_structure,
    pdb_hierarchy    = pdb_hierarchy,
    map_data         = map_data.as_double(),
    box_cushion      = params.box_cushion,
    selection        = selection)
  output_prefix=os.path.basename(inputs.pdb_file_names[0])[:-4]
  if(params.output_file_name_prefix is None):
    file_name = "%s_box.pdb"%output_prefix
  else: file_name = "%s.pdb"%params.output_file_name_prefix
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
  print >> log

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
