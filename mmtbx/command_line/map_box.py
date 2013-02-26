from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_box

import mmtbx.utils
from mmtbx.refinement import print_statistics
import iotbx.pdb
import libtbx.phil
from libtbx.utils import Sorry
import os, sys
from iotbx import file_reader
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils

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

def match_labels(target, label):
  label = label.lower()
  for s in target:
    s = s.lower()
    if(s.count(label)): return True
  return False

def run (args, log=None):
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
  cmdline_phil = []
  got_pdb = False
  got_hkl = False
  got_map = False
  for arg in args :
    if(os.path.isfile(arg)):
      fro = file_reader.any_file(arg)
      if(iotbx.pdb.is_pdb_file(arg)):
        pdb_phil = libtbx.phil.parse("pdb_file=%s" % os.path.abspath(arg))
        cmdline_phil.append(pdb_phil)
        got_pdb = True
      elif(fro.file_type == "ccp4_map"):
        map_phil = libtbx.phil.parse("ccp4_map_file=%s" % os.path.abspath(arg))
        cmdline_phil.append(map_phil)
        got_map = True
      else:
        reflection_file = reflection_file_reader.any_reflection_file(
          file_name=arg, ensure_read_access=False)
        if(reflection_file.file_type() is not None):
          reflection_phil = libtbx.phil.parse(
            "map_coefficients_file=%s" % os.path.abspath(arg))
          cmdline_phil.append(reflection_phil)
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
  print_statistics.make_sub_header("parameters", out=log)
  working_phil, unused = master_phil.fetch(sources=cmdline_phil,
    track_unused_definitions=True)
  if(len(unused)>0):
    for u in unused:
      print str(u)
    raise Sorry("Unused parameters: see above.")
  working_phil.show(out = log)
  params = working_phil.extract()
  if(params.label is None and params.ccp4_map_file is None):
    raise Sorry("'label' has to be defined. Example: label=2mFo-DFc")
  if(params.selection is None):
    raise Sorry("'selection' has to be defined. Example: selection='chain A and resseq 1:10' ")
  #
  print_statistics.make_sub_header("pdb model", out=log)
  pdb_inp = iotbx.pdb.input(file_name=params.pdb_file)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry=pdb_inp.crystal_symmetry())
  xray_structure.show_summary(f=log)
  #
  map_coeff = None
  if(params.map_coefficients_file is not None):
    map_coeff = reflection_file_utils.extract_miller_array_from_file(
      file_name = params.map_coefficients_file,
      label     = params.label,
      type      = "complex",
      log       = log)
  #
  if(got_map):
    print_statistics.make_sub_header("CCP4 map", out=log)
    af = file_reader.any_file(params.ccp4_map_file)
    af.try_as_ccp4_map()
    ccp4_map_object = af.file_object
    ccp4_map_object.show_summary()
    map_data = ccp4_map_object.data
  #
  if(got_hkl):
    fft_map = map_coeff.fft_map(resolution_factor=params.resolution_factor)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
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
  output_prefix=os.path.basename(params.pdb_file)[:-4]
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
      d_min = mmtbx.utils.structure_factors_from_map(
        map_data          = map_data.as_double(),
        unit_cell_lengths = ccp4_map_object.unit_cell_parameters[:3],
        n_real            = ccp4_map_object.unit_cell_grid,
        crystal_symmetry  = xray_structure.crystal_symmetry(),
        resolution_factor = params.resolution_factor).d_min()
    box.map_coefficients(d_min=d_min,
      resolution_factor=params.resolution_factor, file_name=file_name)
  print >> log

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
