from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_box

import mmtbx.utils
import mmtbx.model
from mmtbx.refinement import print_statistics
import iotbx.pdb
import libtbx.phil
from libtbx.utils import Sorry
import os, sys
from iotbx import reflection_file_utils
from iotbx.file_reader import any_file
from cctbx import maptbx
from scitbx.matrix import col

master_phil = libtbx.phil.parse("""
  include scope libtbx.phil.interface.tracking_params
  pdb_file = None
    .type = path
    .help = Model file
    .short_caption = Model file
    .style = file_type:pdb bold input_file
  map_coefficients_file = None
    .type = path
    .help = Input map coefficients file (alternative to ccp4 map file)
    .short_caption = Map coefficients
    .style = file_type:hkl bold input_file process_hkl child:map_labels:label
  label = None
    .type = str
    .short_caption = Map labels
    .help = Labels for map coefficients file
    .style = renderer:draw_map_arrays_widget parent:file_name:map_coefficients_file
  ccp4_map_file = None
    .help = Input map file (CCP4/mrc format).
    .short_caption = Input map file
    .type = str
  selection = all
    .type = str
    .help = Atom selection to be applied to input PDB file
    .short_caption = Atom selection
    .input_size = 400
  selection_radius = 3.0
    .type = float
    .help = Atoms within selection_radius of a selected atom model will be\
             kept as part of the selection.
    .short_caption = Selection radius
  box_cushion = 3.0
    .type = float
    .help = If mask_atoms is False, a box of density will be cut out around\
            the input model (after selections are applied to the model). \
            The size of the box of density will be box_cushion bigger than \
            the model.
    .short_caption = Box cushion

  mask_atoms=False
    .type=bool
    .help = Set map values to 0 outside molecular mask
    .short_caption = Mask atoms

  mask_atoms_atom_radius = 3.
    .type=float
     .help = Radius for masking around atoms
     .short_caption = Mask atoms atom radius

  value_outside_atoms = None
    .type = str
    .help = Set to 'mean' to make average value same inside and outside mask.
    .short_caption = Value outside atoms
  resolution_factor = 1./4
    .type = float
    .help = Resolution factor for calculation of map coefficients
    .short_caption = Resolution factor
  resolution = None
    .type = float
    .help = Resolution for calculation of output map coefficients. Default is \
            based on the gridding of the map (and may be higher-resolution than\
            you want).
    .short_caption = Resolution
  output_format = *xplor *mtz *ccp4
    .type=choice(multi=True)
    .help = Output format(s) for boxed map.
    .short_caption = Output format

  output_file_name_prefix=None
    .type = str
    .help = Prefix for output file names. Default is name of the pdb file \
            without the ".pdb" suffix.
    .short_caption = Output file name prefix

  density_select = False
    .type = bool
    .help = Select boundaries based on where density is located.
    .short_caption = Density select

  density_select_threshold = 0.05
    .type = float
    .help = Choose region where density is this fraction of maximum or greater
    .short_caption = density_select threshold

  get_half_height_width = True
    .type = bool
    .help = Use 4 times half-width at half-height as estimate of max size
    .short_caption = Use half-height width

  symmetry = None
    .type = str
    .help = Optional symmetry (e.g., D7, I, C2) to be used if extract_unique\
            is set.  Alternative to symmetry_file.
    .short_caption = Symmetry
  symmetry_file = None
    .type = path
    .help = Symmetry file to be offset based on origin shift.\
            Symmetry or symmetry_file required if extract_unique=True.  \
            May be a \
            Phenix .ncs_spec file or BIOMTR records or a resolve ncs file.
    .short_caption = Symmetry file

  sequence_file = None
    .type = path
    .help = Sequence file (any standard format). Can be unique part or \
            all copies.  Either sequence file or \
            molecular mass required if extract_unique is True.
    .short_caption = Sequence file

  molecular_mass = None
    .type = float
    .help = Molecular mass of object in map in Da (i.e., 33000 for 33 Kd).\
              Used in identification \
            of unique part of map. Either a sequence file or molecular mass\
            is required if extract_unique is True.
    .short_caption = Molecular mass

  solvent_content = None
    .type = float
    .help = Optional fraction of volume of map that is empty.  \
            Used in identification \
            of unique part of map.
    .short_caption = Solvent content

  extract_unique = False
    .type = bool
    .help = Extract unique part of map. Requires symmetry_file or symmetry and\
            either sequence file or molecular mass to be supplied.
    .short_caption = Extract unique


  soft_mask=False
    .type=bool
    .help = Use Gaussian mask in mask_atoms
    .short_caption = Soft mask

  soft_mask_radius=3
    .type=float
    .help = Gaussian mask smoothing radius
    .short_caption = Soft mask radius

  lower_bounds = None
    .type = ints
    .help = Lower bounds for cut out box. You can specify them directly.
    .short_caption = Lower bounds

  upper_bounds = None
    .type = ints
    .help = Upper bounds for cut out box.  You can specify them directly.
    .short_caption = Upper bounds

  keep_map_size = False
    .type=bool
    .help = Keep original map gridding (do not cut anything out). \
            Use to apply soft_mask and/or mask_atoms keeping same map size.
    .short_caption = Keep map size

  keep_origin = True
    .type=bool
    .help = If True, write out map, map_coefficients, and model \
            with origin in original location.  \
            If false, shifted origin to (0,0,0).  \
            NOTE: The unit_cell for all output will always be the \
              map_box unit cell.  Only the origin is kept/shifted.\
    .short_caption = Keep origin

  restrict_map_size = False
    .type=bool
    .help = Do not go outside original map boundaries
    .short_caption = Restrict map size

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
     ncs_object=None,
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
    pdb_hierarchy and not params.keep_map_size and not params.upper_bounds
     and not params.extract_unique):
    raise Sorry("PDB file is needed unless extract_unique, "+
      "density_select, keep_map_size \nor bounds are set .")
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
  if (params.extract_unique and not params.resolution):
    raise Sorry("Please set resolution for extract_unique")
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
      if not crystal_symmetry: crystal_symmetry=map_coeff.crystal_symmetry()
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
      if not crystal_symmetry: crystal_symmetry=ccp4_map.crystal_symmetry()
      map_data = ccp4_map.data#map_data()
      if inputs.ccp4_map_file_name.endswith(".ccp4"):
        map_or_map_coeffs_prefix=os.path.basename(
          inputs.ccp4_map_file_name[:-5])
      else:
        map_or_map_coeffs_prefix=os.path.basename(
          inputs.ccp4_map_file_name[:-4])
  else: # have map_data
    map_or_map_coeffs_prefix=None

  if crystal_symmetry and not inputs.crystal_symmetry:
    inputs.crystal_symmetry=crystal_symmetry

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

  if not ncs_object:
    from mmtbx.ncs.ncs import ncs
    ncs_object=ncs()
    if params.symmetry_file:
      ncs_object.read_ncs(params.symmetry_file,log=log)
      print >>log,"Total of %s operators read" %(ncs_object.max_operators())
  if not ncs_object or ncs_object.max_operators()<1:
      print >>log,"No symmetry available"
  if ncs_object:
    n_ops=max(1,ncs_object.max_operators())
  else:
    n_ops=1

  # Get sequence if extract_unique is set
  if params.extract_unique:
    if params.sequence_file and not params.molecular_mass:
      if n_ops > 1: # get unique part of sequence and multiply
        remove_duplicates=True
      else:
        remove_duplicates=False
      from iotbx.bioinformatics import get_sequences
      sequence=n_ops * (" ".join(get_sequences(file_name=params.sequence_file,
        remove_duplicates=remove_duplicates)))
      # get molecular mass from sequence
      from iotbx.bioinformatics import text_from_chains_matching_chain_type
      n_protein=len(text_from_chains_matching_chain_type(
        text=sequence,chain_type='PROTEIN'))
      n_rna=len(text_from_chains_matching_chain_type(
        text=sequence,chain_type='RNA'))
      n_dna=len(text_from_chains_matching_chain_type(
        text=sequence,chain_type='DNA'))
      params.molecular_mass=n_protein*110+(n_rna+n_dna)*330
    elif not params.molecular_mass:
      raise Sorry("Need a sequence file or molecular mass for extract_unique")
  else:
    molecular_mass=None
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
    restrict_map_size     = params.restrict_map_size,
    lower_bounds          = params.lower_bounds,
    upper_bounds          = params.upper_bounds,
    extract_unique        = params.extract_unique,
    solvent_content       = params.solvent_content,
    molecular_mass        = params.molecular_mass,
    resolution            = params.resolution,
    ncs_object            = ncs_object,
    symmetry              = params.symmetry,

    )

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

  # ncs_object is original
  #  box.ncs_object is shifted by shift_cart

  print >>log,"Box cell dimensions: (%.2f, %.2f, %.2f) A" %(
      box.box_crystal_symmetry.unit_cell().parameters()[:3])
  if params.keep_origin:
     print>>log,"Box origin is at grid position of : (%d, %d, %d) " %(
        tuple(box.origin_shift_grid_units(reverse=True)))
     print>>log,"Box origin is at coordinates: (%.2f, %.2f, %.2f) A" %(
       tuple(-col(box.shift_cart)))

  if box.pdb_outside_box_msg:
    print >> log, box.pdb_outside_box_msg

  # NOTE: box object is always shifted to place origin at (0,0,0)

  # For output files ONLY:
  #   keep_origin==False leave origin at (0,0,0)
   #  keep_origin==True: we shift everything back to where it was,

  if (not params.keep_origin):
    if box.shift_cart:
      print >>log,\
        "Final coordinate shift for output files: (%.2f,%.2f,%.2f) A" %(
        tuple(box.shift_cart))
    else:
      print >>log,"\nOutput files are in same location as original: origin "+\
        "is at (0,0,0)"
  else: # keep_origin
    print >>log,"\nOutput files are in same location as original, just cut out."
    print >>log,"Note that output maps are only valid in the cut out region.\n"

  if params.keep_origin:
    ph_box_original_location = ph_box.deep_copy()
    sites_cart = box.shift_sites_cart_back(
      box.xray_structure_box.sites_cart())
    xrs_offset = ph_box_original_location.extract_xray_structure(
        crystal_symmetry=box.xray_structure_box.crystal_symmetry()
          ).replace_sites_cart(new_sites = sites_cart)
    ph_box_original_location.adopt_xray_structure(xrs_offset)
    box.hierarchy_original_location=ph_box_original_location
  else:
    box.hierarchy_original_location=None

  if write_output_files:
    # Write PDB file
    if ph_box.overall_counts().n_residues>0:

      if(params.output_file_name_prefix is None):
        file_name = "%s_box.pdb"%output_prefix
      else: file_name = "%s.pdb"%params.output_file_name_prefix

      if params.keep_origin:  # Keeping origin
        print >> log, "Writing boxed PDB with box unit cell and in "+\
          "original\n    position to:   %s"%(
          file_name)
        ph_box_original_location.write_pdb_file(file_name=file_name,
          crystal_symmetry = box.xray_structure_box.crystal_symmetry())

      else: # write box PDB in box cell
        print >> log, "Writing shifted boxed PDB to file:   %s"%file_name
        ph_box.write_pdb_file(file_name=file_name, crystal_symmetry =
          box.xray_structure_box.crystal_symmetry())

    # Write NCS file if NCS
    if ncs_object and ncs_object.max_operators()>0:
      if(params.output_file_name_prefix is None):
        output_symmetry_file = "%s_box.ncs_spec"%output_prefix
      else:
        output_symmetry_file = "%s.ncs_spec"%params.output_file_name_prefix

      if params.keep_origin:
        if params.symmetry_file:
          print >>log,"\nDuplicating symmetry in %s and writing to %s" %(
            params.symmetry_file,output_symmetry_file)
        else:
          print >>log,"\nWriting symmetry to %s" %(output_symmetry_file)
        ncs_object.format_all_for_group_specification(
          file_name=output_symmetry_file)

      else:
        print >>log,"\nOffsetting symmetry in %s and writing to %s" %(
          params.symmetry_file,output_symmetry_file)
        box.ncs_object.format_all_for_group_specification(
          file_name=output_symmetry_file)

    # Write ccp4 map.  Shift back to original location if keep_origin=True
    if("ccp4" in params.output_format):
     if(params.output_file_name_prefix is None):
       file_name = "%s_box.ccp4"%output_prefix
     else: file_name = "%s.ccp4"%params.output_file_name_prefix

     if params.keep_origin:
       print >> log, "Writing boxed map with box unit_cell and "+\
          "original\n    position to CCP4 formatted file:   %s"%file_name
     else:
       print >> log, "Writing box map shifted to (0,0,0) to CCP4 "+\
          "formatted file:   %s"%file_name

     box.write_ccp4_map(file_name=file_name,
       shift_back=params.keep_origin)

    # Write xplor map.  Shift back if keep_origin=True
    if("xplor" in params.output_format):
     if(params.output_file_name_prefix is None):
       file_name = "%s_box.xplor"%output_prefix
     else: file_name = "%s.xplor"%params.output_file_name_prefix
     if params.keep_origin:
       print >> log, "Writing boxed map with box unit_cell and original "+\
         "position\n    to X-plor formatted file: %s"%file_name
     else:
       print >> log, "Writing box_map shifted to (0,0,0) to X-plor "+\
         "formatted file: %s"%file_name
     box.write_xplor_map(file_name=file_name,
         shift_back=params.keep_origin)

    # Write mtz map coeffs.  Shift back if keep_origin=True
    if("mtz" in params.output_format):
     if(params.output_file_name_prefix is None):
       file_name = "%s_box.mtz"%output_prefix
     else: file_name = "%s.mtz"%params.output_file_name_prefix

     if params.keep_origin:
       print >> log, "Writing map coefficients with box_map unit_cell"+\
         " but position matching\n   "+\
         " original position to MTZ file: %s"%file_name
     else:
       print >> log, "Writing box_map coefficients shifted to (0,0,0) "+\
          "to MTZ file: %s"%file_name
     if(map_coeff is not None):
       d_min = map_coeff.d_min()
     elif params.resolution is not None:
       d_min = params.resolution
     else:
       d_min = maptbx.d_min_from_map(map_data=box.map_box,
         unit_cell=box.xray_structure_box.unit_cell())
     box.map_coefficients(d_min=d_min,
       resolution_factor=params.resolution_factor, file_name=file_name,
       shift_back=params.keep_origin)

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
