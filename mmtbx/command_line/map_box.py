from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_box

import mmtbx.utils
from mmtbx.refinement import print_statistics
import iotbx.pdb
import libtbx.phil
from libtbx.utils import Sorry
import os, sys
from iotbx import reflection_file_utils
from iotbx.file_reader import any_file
from cctbx import maptbx

master_phil = libtbx.phil.parse("""
  include scope libtbx.phil.interface.tracking_params
  pdb_file = None
    .type = path
    .short_caption = Model file
    .style = file_type:pdb bold input_file
  map_coefficients_file = None
    .type = path
    .style = file_type:hkl bold input_file process_hkl child:map_labels:label
  label = None
    .type = str
    .short_caption = Map labels
    .style = renderer:draw_map_arrays_widget parent:file_name:map_coefficients_file
  ccp4_map_file = None
    .type = str
  selection = all
    .type = str
    .input_size = 400
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
  get_half_height_width = True
    .type = bool
    .help = Use 4 times half-width at half-height as estimate of max size
  ncs_file = None
    .type = path
    .help = NCS file (to be offset based on origin shift)
  mask_atoms=False
    .type=bool
    .help = Set map values to 0 outside molecular mask
  mask_atoms_atom_radius = 3.
    .type=float
  value_outside_atoms = None
    .type = str
    .help = Set to 'mean' to make average value same inside and outside mask.
  soft_mask=False
    .type=bool
    .help = Use Gaussian mask in mask_atoms
  soft_mask_radius=3
    .type=float
    .help = Gaussian mask smoothing radius
  lower_bounds = None
    .type = ints
    .help = Lower bounds for cut out box
  upper_bounds = None
    .type = ints
    .help = Upper bounds for cut out box
  keep_map_size = False
    .type=bool
    .help = Keep original map gridding (do not cut anything out). \
            Use to apply soft_mask and/or mask_atoms keeping same map size.
  restrict_map_size = False
    .type=bool
    .help = Do not go outside original map boundaries

  gui
    .help = "GUI-specific parameter required for output directory"
  {
    output_dir = None
    .type = path
    .style = output_dir
  }
""", process_includes=True)

master_params = master_phil

def run(args, crystal_symmetry=None,
     pdb_hierarchy=None,
     map_data=None,
     lower_bounds=None,
     upper_bounds=None,
     write_output_files=True,
     log=None):
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
  if(len(args) == 0 and not pdb_hierarchy):
    print default_message
    master_phil.show(prefix="  ")
    return
  inputs = mmtbx.utils.process_command_line_args(args = args,
    cmd_cs=crystal_symmetry,
    master_params = master_phil)
  params = inputs.params.extract()
  master_phil.format(python_object=params).show(out=log)

  # Overwrite params with parameters in call if available
  if lower_bounds:
     params.lower_bounds=lower_bounds
  if upper_bounds:
     params.upper_bounds=upper_bounds

  # PDB file
  if params.pdb_file and not inputs.pdb_file_names and not pdb_hierarchy:
    inputs.pdb_file_names=[params.pdb_file]
  if(len(inputs.pdb_file_names)!=1 and not params.density_select and not
    pdb_hierarchy and not params.keep_map_size and not params.upper_bounds):
    raise Sorry("PDB file is needed unless density_select or keep_map_size or bounds are set .")
  if (len(inputs.pdb_file_names)!=1 and not pdb_hierarchy and \
       (params.mask_atoms or params.soft_mask )):
    raise Sorry("PDB file is needed for mask_atoms or soft_mask")
  if (params.density_select and params.keep_map_size):
    raise Sorry("Cannot set both density_select and keep_map_size")
  if (params.density_select and params.upper_bounds):
    raise Sorry("Cannot set both density_select and bounds")
  if (params.keep_map_size and params.upper_bounds):
    raise Sorry("Cannot set both keep_map_size and bounds")
  if (params.upper_bounds and not params.lower_bounds):
    raise Sorry("Please set lower_bounds if you set upper_bounds")
  print_statistics.make_sub_header("pdb model", out=log)
  if len(inputs.pdb_file_names)>0:
    pdb_inp = iotbx.pdb.input(file_name=inputs.pdb_file_names[0])
    pdb_hierarchy = pdb_inp.construct_hierarchy()
  if pdb_hierarchy:
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
  else:
    pdb_hierarchy=None
  # Map or map coefficients
  map_coeff = None
  if (not map_data):
    # read first mtz file
    if ( (len(inputs.reflection_file_names) > 0) or
         (params.map_coefficients_file is not None) ):
      # file in phil takes precedent
      if (params.map_coefficients_file is not None):
        if (len(inputs.reflection_file_names) == 0):
          inputs.reflection_file_names.append(params.map_coefficients_file)
        else:
          inputs.reflection_file_names[0] = params.map_coefficients_file
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
    # or read CCP4 map
    elif ( (inputs.ccp4_map is not None) or
           (params.ccp4_map_file is not None) ):
      if (params.ccp4_map_file is not None):
        af = any_file(params.ccp4_map_file)
        if (af.file_type == 'ccp4_map'):
          inputs.ccp4_map = af.file_content
          inputs.ccp4_map_file_name = params.ccp4_map_file
      print_statistics.make_sub_header("CCP4 map", out=log)
      ccp4_map = inputs.ccp4_map
      ccp4_map.show_summary(prefix="  ",out=log)
      map_data = ccp4_map.data#map_data()
      if inputs.ccp4_map_file_name.endswith(".ccp4"):
        map_or_map_coeffs_prefix=os.path.basename(
          inputs.ccp4_map_file_name[:-5])
      else:
        map_or_map_coeffs_prefix=os.path.basename(
          inputs.ccp4_map_file_name[:-4])
  else: # have map_data
    map_or_map_coeffs_prefix=None

  # final check that map_data exists
  if(map_data is None):
    raise Sorry("Map or map coefficients file is needed.")

  if len(inputs.pdb_file_names)>0:
    output_prefix=os.path.basename(inputs.pdb_file_names[0])[:-4]
  else:
    output_prefix=map_or_map_coeffs_prefix

  if not pdb_hierarchy: # get an empty hierarchy
    from cctbx.array_family import flex
    pdb_hierarchy=iotbx.pdb.input(
      source_info='',lines=flex.split_lines('')).construct_hierarchy()
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry=inputs.crystal_symmetry)
  xray_structure.show_summary(f=log)
  #
  selection = pdb_hierarchy.atom_selection_cache().selection(
    string = params.selection)
  if selection.size():
    print_statistics.make_sub_header("atom selection", out=log)
    print >> log, "Selection string: selection='%s'"%params.selection
    print >> log, \
        "  selects %d atoms from total %d atoms."%(selection.count(True),
        selection.size())
  sites_cart_all = xray_structure.sites_cart()
  sites_cart = sites_cart_all.select(selection)
  selection = xray_structure.selection_within(
    radius    = params.selection_radius,
    selection = selection)
#
  if params.density_select:
    print_statistics.make_sub_header(
    "Extracting box around selected density and writing output files", out=log)
  else:
   print_statistics.make_sub_header(
    "Extracting box around selected atoms and writing output files", out=log)
  #
  if params.value_outside_atoms=='mean':
    print >>log,"\nValue outside atoms mask will be set to mean inside mask"
  if params.get_half_height_width and params.density_select:
    print >>log,"\nHalf width at half height will be used to id boundaries"
  if params.soft_mask and sites_cart_all.size()>0:
    print >>log,"\nSoft mask will be applied to model-based mask"
  if params.keep_map_size:
    print >>log,"\nEntire map will be kept (not cutting out region)"
  if params.restrict_map_size:
    print >>log,"\nOutput map will be within input map"
  if params.lower_bounds and params.upper_bounds:
    print >>log,"Bounds for cut out map are (%s,%s,%s) to (%s,%s,%s)" %(
     tuple(list(params.lower_bounds)+list(params.upper_bounds)))

  box = mmtbx.utils.extract_box_around_model_and_map(
    xray_structure   = xray_structure,
    map_data         = map_data.as_double(),
    box_cushion      = params.box_cushion,
    selection        = selection,
    density_select   = params.density_select,
    threshold        = params.density_select_threshold,
    get_half_height_width = params.get_half_height_width,
    mask_atoms       = params.mask_atoms,
    soft_mask        = params.soft_mask,
    soft_mask_radius = params.soft_mask_radius,
    mask_atoms_atom_radius = params.mask_atoms_atom_radius,
    value_outside_atoms = params.value_outside_atoms,
    keep_map_size         = params.keep_map_size,
    restrict_map_size         = params.restrict_map_size,
    lower_bounds          = params.lower_bounds,
    upper_bounds          = params.upper_bounds,
    )
  if box.shift_cart:
    print >>log,"Final coordinate shift: (%.1f,%.1f,%.1f)" %(
      tuple(box.shift_cart))

  print >>log,"Final cell dimensions: (%.1f,%.1f,%.1f)\n" %(
      box.box_crystal_symmetry.unit_cell().parameters()[:3])

  if box.pdb_outside_box_msg:
    print >> log, box.pdb_outside_box_msg


  ph_box = pdb_hierarchy.select(selection)
  ph_box.adopt_xray_structure(box.xray_structure_box)

  box.hierarchy=ph_box

  if (inputs and inputs.crystal_symmetry and inputs.ccp4_map and
    inputs.crystal_symmetry.unit_cell().parameters() and
     inputs.ccp4_map.unit_cell_parameters  ) and (
       inputs.crystal_symmetry.unit_cell().parameters() !=
       inputs.ccp4_map.unit_cell_parameters):
    print >>log,"\nNOTE: Mismatch of unit cell parameters from CCP4 map:"
    print >>log,"Unit cell from CCP4 map 'unit cell parameters': "+\
      "%.1f, %.1f, %.1f, %.1f, %.1f, %.1f)" %tuple(
        inputs.ccp4_map.unit_cell_parameters)
    print >>log,"Unit cell from CCP4 map 'map grid':             "+\
      "%.1f, %.1f, %.1f, %.1f, %.1f, %.1f)" %tuple(
       inputs.crystal_symmetry.unit_cell().parameters())
    print >>log,"\nInterpreting this as the 'unit cell parameters' was "+\
      "original map \ndimension and 'map grid' is the "+\
      "portion actually in the map that was supplied here.\n"
    box.unit_cell_parameters_from_ccp4_map=inputs.ccp4_map.unit_cell_parameters
    box.unit_cell_parameters_deduced_from_map_grid=\
       inputs.crystal_symmetry.unit_cell().parameters()

  else:
    box.unit_cell_parameters_from_ccp4_map=None
    box.unit_cell_parameters_deduced_from_map_grid=None

  if not write_output_files:
    return box

  if(params.output_file_name_prefix is None):
    file_name = "%s_box.pdb"%output_prefix
  else: file_name = "%s.pdb"%params.output_file_name_prefix
  ph_box = pdb_hierarchy.select(selection)
  ph_box.adopt_xray_structure(box.xray_structure_box)
  ph_box.write_pdb_file(file_name=file_name, crystal_symmetry =
    box.xray_structure_box.crystal_symmetry())

  if("ccp4" in params.output_format):
    if(params.output_file_name_prefix is None):
      file_name = "%s_box.ccp4"%output_prefix
    else: file_name = "%s.ccp4"%params.output_file_name_prefix
    print >> log, "Writing map to CCP4 formatted file:   %s"%file_name
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
    #ncs_object.display_all(log=log)
    print >>log,"Total of %s operators read" %(ncs_object.max_operators())
    if not ncs_object or ncs_object.max_operators()<1:
      print >>log,"Skipping...no NCS available"
    elif box.shift_cart:
      from scitbx.math import  matrix
      print >>log,"Shifting NCS operators "+\
        "based on coordinate shift of (%7.1f,%7.1f,%7.1f)" %(
        tuple(box.shift_cart))
      ncs_object=ncs_object.coordinate_offset(
       coordinate_offset=matrix.col(box.shift_cart))
      #ncs_object.display_all(log=log)
      print >>log,"Total of %s operators" %(ncs_object.max_operators())
    ncs_object.format_all_for_group_specification(
       file_name=output_ncs_file)
    box.ncs_object=ncs_object
  else:
    box.ncs_object=None
  print >> log
  return box

# =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
from wxGUI2 import utils

def validate_params(params):
  return True

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    result = run(args=self.args, log=sys.stdout)
    return result

# =============================================================================

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
