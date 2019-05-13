from __future__ import division
from __future__ import print_function
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
    .help = Optional model file used to define region to be cut out
    .short_caption = Model file (optional)
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
    .short_caption = Atom selection (optional)
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
  map_scale_factor = None
    .type = float
    .help = Scale factor to apply to map
    .short_caption = Map scale factor
  scale_max = 99999
    .type = float
    .help = Maximum value of amplitudes for output mtz file. If None, apply\
             volume scaling
    .short_caption = Scale max
  resolution = None
    .type = float
    .help = Resolution for calculation of output map coefficients. Default is \
            based on the gridding of the map (and may be higher-resolution than\
            you want).
    .short_caption = Resolution
  output_format = xplor *mtz *ccp4
    .type=choice(multi=True)
    .help = Output format(s) for boxed map. Note that mtz format is only\
            available if keep_origin=False or keep_map_size=True. (These \
            are the cases where the map is cut down to size and placed \
            at the origin or there is a full unit cell of data.)
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
    .help = Use 4 times half-width at half-height as estimate of max size \
              in density_select
    .short_caption = Use half-height width

  symmetry = None
    .type = str
    .help = Optional symmetry (e.g., D7, I, C2) to be used if extract_unique\
            is set.  Alternative to symmetry_file.  To find symmetry \
            automatically specify symmetry=ALL.
    .short_caption = Symmetry
  symmetry_file = None
    .type = path
    .help = Symmetry file.\
            Symmetry or symmetry_file required if extract_unique=True.  \
            May be a \
            Phenix .ncs_spec file or BIOMTR records or a resolve ncs file.
    .short_caption = Symmetry file

  sequence_file = None
    .type = path
    .help = Sequence file (any standard format). Can be unique part or \
            all copies.  Used in identification of unique part of map.
    .short_caption = Sequence file (optional)

  molecular_mass = None
    .type = float
    .help = Molecular mass of object in map in Da (i.e., 33000 for 33 Kd).\
              Used in identification \
            of unique part of map.
    .short_caption = Molecular mass (optional)

  solvent_content = None
    .type = float
    .help = Optional fraction of volume of map that is empty.  \
            Used in identification \
            of unique part of map.
    .short_caption = Solvent content

  extract_unique = False
    .type = bool
    .help = Extract unique part of map. Requires either symmetry_file or \
            symmetry and\
            either sequence file or molecular mass to be supplied. If chain \
            type is not protein it should be set as well.
    .short_caption = Extract unique

  box_buffer = 5
    .type = int
    .help = Padding around unique region
    .short_caption = Padding around unique region

  soft_mask_extract_unique = True
    .type = bool
    .help = Create soft mask at edges of extract_unique box (feather map into \
               edge of box). Uses resolution as mask_radius
      .short_caption = Soft mask in extract unique

  mask_expand_ratio = 1
      .type = int
      .help = Mask expansion relative to resolution for extract_unique
      .short_caption = Mask expand ratio

  regions_to_keep = None
    .type = int
    .short_caption = Regions to keep
    .help = You can specify a limit to the number of regions to keep\
             when generating the asymmetric unit of density.

  keep_low_density = True
    .type = bool
    .help = Get remainder (weak density) with extract_unique.
    .short_caption = Get remainder


  chain_type = None *PROTEIN DNA RNA
    .type = choice
    .help = Chain type. Only used if extract_unique is set. Has minor effect \
            in setting thresholds for identification of molecular region.\
            Use None if there is a mixture.
    .short_caption = Chain type

  soft_mask = False
    .type=bool
    .help = Use Gaussian mask in mask_atoms
    .short_caption = Soft mask

  soft_mask_radius=3
    .type=float
    .help = Gaussian mask smoothing radius
    .short_caption = Soft mask radius

  lower_bounds = None
    .type = ints
    .help = Lower bounds for cut out box. You can specify them directly.\
            NOTE: lower and upper bounds refer to grid points after shifting \
            the map to place the origin at (0,0,0).
    .short_caption = Lower bounds

  upper_bounds = None
    .type = ints
    .help = Upper bounds for cut out box.  You can specify them directly.\
            NOTE: lower and upper bounds refer to grid points after shifting \
            the map to place the origin at (0,0,0).
    .short_caption = Upper bounds

  keep_map_size = False
    .type=bool
    .help = Keep original map gridding (do not cut anything out). \
            Use to apply soft_mask and/or mask_atoms keeping same map size.
    .short_caption = Keep map size

  keep_origin = True
    .type=bool
    .help = Write out map, map_coefficients, and model \
            with origin in original location.  \
            If False, shift the origin to (0,0,0).  \
            NOTE: to cut out a part of a map, shift the origin to (0,0,0),\
               and make a new small map use keep_origin=False\
               keep_input_unit_cell_and_grid=False
    .short_caption = Keep origin

  keep_input_unit_cell_and_grid = True
     .type = bool
     .help = Keep the input unit_cell dimensions and unit_cell_grid. \
             If False, use the dimensions and grid of the cut out box as the \
              unit cell map_box dimensions and grid.\
            NOTE: to cut out a part of a map, shift the origin to (0,0,0),\
               and make a new small map set keep_origin=False and \
               keep_input_unit_cell_and_grid=False
     .short_caption = Keep input unit cell and grid

  output_unit_cell = None
     .type = floats
     .help = You can specify the unit cell for your map with 3 numbers. \
              This should normally\
             not be necessary. It can be used to fix a map that has the \
             wrong unit cell.
     .short_caption = Output unit cell
     .expert_level = 3

  output_unit_cell_grid = None
    .type = ints
    .help = You can specify the grid (3 integers) corresponding to the \
              output unit cell. \
              This can be used to specify the full grid for the unit cell. \
              if output_unit_cell is not specified, new unit cell parameters\
              will be generated to maintain the grid spacing.
    .short_caption = Output unit cell grid
    .expert_level = 3

  output_origin_grid_units = None
    .type = ints
    .help = You can specify the origin of your output map.  Normally you \
           should use keep_origin=True or False to specify your origin \
           but if you want to move it to a specific grid point you can do that.\
    .short_caption = Output origin
    .expert_level = 3

  output_origin_match_this_file = None
    .type = path
    .help = As output_origin_grid_units, but use origin from this file
    .short_caption = File with origin info

  output_external_origin = None
    .type = floats
    .help = Write ORIGIN record to map file (this is an external origin \
            used to specify relationship to external files such as model files.\
            Three floating point numbers (A).
    .short_caption = output external origin

  restrict_map_size = False
    .type=bool
    .help = Do not go outside original map boundaries
    .short_caption = Restrict map size

  ignore_symmetry_conflicts = False
    .type=bool
    .help = Ignore unit cell from model if it conflicts with the map.
    .short_caption = Ignore symmetry conflicts

  output_ccp4_map_mean = None
    .type = float
    .help = Choose mean and SD of output CCP4 map
    .short_caption = Mean of output CCP4 map

  output_ccp4_map_sd = None
    .type = float
    .help = Choose mean and SD of output CCP4 map
    .short_caption = SD of output CCP4 map

  output_map_labels = None
    .type = str
    .multiple = True
    .help = Add this label to output map
    .short_caption = Add label

  gui
    .help = "GUI-specific parameter required for output directory"
  {
    output_dir = None
    .type = path
    .style = output_dir
  }
""", process_includes=True)

master_params = master_phil

def remove_element(text_list,element=None):
    new_text_list=[]
    for x in text_list:
      if x != element:
        new_text_list.append(x)
    return new_text_list

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
    print(default_message)
    master_phil.show(prefix="  ")
    return

  # Process inputs ignoring symmetry conflicts just to get the value of
  #   ignore_symmetry_conflicts...

  inputs = mmtbx.utils.process_command_line_args(args = args,
      cmd_cs=crystal_symmetry,
      master_params = master_phil,
      suppress_symmetry_related_errors=True)
  params = inputs.params.extract()

  # Now process inputs for real and write a nice error message if necessary.
  try:
    inputs = mmtbx.utils.process_command_line_args(args = args,
      cmd_cs=crystal_symmetry,
      master_params = master_phil,
      suppress_symmetry_related_errors=params.ignore_symmetry_conflicts)
  except Exception as e:
    if str(e).find("symmetry mismatch ")>1:
      raise Sorry(str(e)+"\nTry 'ignore_symmetry_conflicts=True'")
    else:
      raise e

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
  if (params.extract_unique):
    if (not params.resolution):
      raise Sorry("Please set resolution for extract_unique")
    if (not params.symmetry) and (not params.symmetry_file) and \
        (not ncs_object):
      raise Sorry(
        "Please supply a symmetry file or symmetry for extract_unique (you "+
       "\ncan try symmetry=ALL if you do not know your symmetry or "+
        "symmetry=C1 if \nthere is none)")
      from mmtbx.ncs.ncs import ncs
      ncs_object=ncs()
      ncs_object.set_unit_ncs()

  if params.keep_input_unit_cell_and_grid and (
      (params.output_unit_cell_grid is not None ) or
      (params.output_unit_cell is not None ) ):
    raise Sorry("If you set keep_input_unit_cell_and_grid then you cannot "+\
       "set \noutput_unit_cell_grid or output_unit_cell")

  if (write_output_files) and ("mtz" in params.output_format) and (
       (params.keep_origin) and (not params.keep_map_size)):
    print("\nNOTE: Skipping write of mtz file as keep_origin=True and \n"+\
       "keep_map_size is False\n")
    params.output_format=remove_element(params.output_format,element='mtz')

  if (write_output_files) and ("mtz" in params.output_format) and (
       (params.extract_unique)):
    print("\nNOTE: Skipping write of mtz file as extract_unique=True\n")
    params.output_format=remove_element(params.output_format,element='mtz')


  if params.output_origin_match_this_file:

    af = any_file(params.output_origin_match_this_file)
    if (af.file_type == 'ccp4_map'):
      origin=af.file_content.data.origin()
      params.output_origin_grid_units=origin
      print("Origin of (%s,%s,%s) taken from %s" %(
         origin[0],origin[1],origin[2],params.output_origin_match_this_file))
    if not params.ccp4_map_file:
      raise Sorry(
       "Need to specify ccp4_map_file=xxx if you use "+
           "output_origin_match_this_file=xxxx")
  if params.output_origin_grid_units is not None and params.keep_origin:
    params.keep_origin=False
    print("Setting keep_origin=False as output_origin_grid_units is set")
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
  input_unit_cell_grid=None
  input_unit_cell=None
  input_map_labels=None
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
      map_data = ccp4_map.data #map_data()
      input_unit_cell_grid=ccp4_map.unit_cell_grid
      input_unit_cell=ccp4_map.unit_cell_parameters
      input_map_labels=ccp4_map.get_labels()

      if inputs.ccp4_map_file_name.endswith(".ccp4"):
        map_or_map_coeffs_prefix=os.path.basename(
          inputs.ccp4_map_file_name[:-5])
      else:
        map_or_map_coeffs_prefix=os.path.basename(
          inputs.ccp4_map_file_name[:-4])
  else: # have map_data
    map_or_map_coeffs_prefix=None

  if params.map_scale_factor:
    print("Applying scale factor of %s to map data on read-in" %(
       params.map_scale_factor))
    map_data=map_data*params.map_scale_factor

  if params.output_origin_grid_units is not None:
    origin_to_match=tuple(params.output_origin_grid_units)
  else:
    origin_to_match=None

  if origin_to_match:
    sc=[]
    for x,o,a in zip(crystal_symmetry.unit_cell().parameters()[:3],
        origin_to_match,
        map_data.all()):
      sc.append(-x*o/a)
    shift_cart_for_origin_to_match=tuple(sc)
  else:
    origin_to_match=None
    shift_cart_for_origin_to_match=None



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
  if not params.selection: params.selection="all"
  selection = pdb_hierarchy.atom_selection_cache().selection(
    string = params.selection)
  if selection.size():
    print_statistics.make_sub_header("atom selection", out=log)
    print("Selection string: selection='%s'"%params.selection, file=log)
    print("  selects %d atoms from total %d atoms."%(selection.count(True),
        selection.size()), file=log)
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
      print("Total of %s operators read" %(ncs_object.max_operators()), file=log)
  if not ncs_object or ncs_object.max_operators()<1:
      print("No symmetry available", file=log)
  if ncs_object:
    n_ops=max(1,ncs_object.max_operators())
  else:
    n_ops=1

  # Get sequence if extract_unique is set
  sequence=None
  if params.extract_unique:
    if params.sequence_file:
      if n_ops > 1: # get unique part of sequence
        remove_duplicates=True
      else:
        remove_duplicates=False
      from iotbx.bioinformatics import get_sequences
      sequence=(" ".join(get_sequences(file_name=params.sequence_file,
        remove_duplicates=remove_duplicates)))

    if params.chain_type in ['None',None]: params.chain_type=None
    if sequence and not params.molecular_mass:
      # get molecular mass from sequence
      from iotbx.bioinformatics import text_from_chains_matching_chain_type
      if params.chain_type in [None,'PROTEIN']:
        n_protein=len(text_from_chains_matching_chain_type(
          text=sequence,chain_type='PROTEIN'))
      else:
        n_protein=0
      if params.chain_type in [None,'RNA']:
        n_rna=len(text_from_chains_matching_chain_type(
          text=sequence,chain_type='RNA'))
      else:
        n_rna=0
      if params.chain_type in [None,'DNA']:
        n_dna=len(text_from_chains_matching_chain_type(
         text=sequence,chain_type='DNA'))
      else:
        n_dna=0
      params.molecular_mass=n_ops*(n_protein*110+(n_rna+n_dna)*330)
      print("\nEstimate of molecular mass is %.0f " %(params.molecular_mass), file=log)
  if params.density_select:
    print_statistics.make_sub_header(
    "Extracting box around selected density and writing output files", out=log)
  else:
   print_statistics.make_sub_header(
    "Extracting box around selected atoms and writing output files", out=log)
  #
  if params.value_outside_atoms=='mean':
    print("\nValue outside atoms mask will be set to mean inside mask", file=log)
  if params.get_half_height_width and params.density_select:
    print("\nHalf width at half height will be used to id boundaries", file=log)
  if params.soft_mask and sites_cart_all.size()>0:
    print("\nSoft mask will be applied to model-based mask", file=log)
  if params.keep_map_size:
    print("\nEntire map will be kept (not cutting out region)", file=log)
  if params.restrict_map_size:
    print("\nOutput map will be within input map", file=log)
  if params.lower_bounds and params.upper_bounds:
    print("Bounds for cut out map are (%s,%s,%s) to (%s,%s,%s)" %(
     tuple(list(params.lower_bounds)+list(params.upper_bounds))), file=log)

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
    regions_to_keep       = params.regions_to_keep,
    box_buffer            = params.box_buffer,
    soft_mask_extract_unique = params.soft_mask_extract_unique,
    mask_expand_ratio = params.mask_expand_ratio,
    keep_low_density      = params.keep_low_density,
    chain_type            = params.chain_type,
    sequence              = sequence,
    solvent_content       = params.solvent_content,
    molecular_mass        = params.molecular_mass,
    resolution            = params.resolution,
    ncs_object            = ncs_object,
    symmetry              = params.symmetry,

    )

  ph_box = pdb_hierarchy.select(selection)
  ph_box.adopt_xray_structure(box.xray_structure_box)
  box.hierarchy=ph_box

  if (inputs and
    inputs.crystal_symmetry and inputs.ccp4_map and
    inputs.crystal_symmetry.unit_cell().parameters() and
     inputs.ccp4_map.unit_cell_parameters  ) and (
       inputs.crystal_symmetry.unit_cell().parameters() !=
       inputs.ccp4_map.unit_cell_parameters):
    print("\nNOTE: Input CCP4 map is only part of unit cell:", file=log)
    print("Full unit cell ('unit cell parameters'): "+\
      "(%.1f, %.1f, %.1f, %.1f, %.1f, %.1f) A" %tuple(
        inputs.ccp4_map.unit_cell_parameters), file=log)
    print("Size of CCP4 map 'map unit cell':        "+\
      "(%.1f, %.1f, %.1f, %.1f, %.1f, %.1f) A" %tuple(
       inputs.crystal_symmetry.unit_cell().parameters()), file=log)
    print("Full unit cell as grid units: (%s, %s, %s)" %(
      inputs.ccp4_map.unit_cell_grid), file=log)
    print("Map unit cell as grid units:  (%s, %s, %s)" %(
      map_data.all()), file=log)

    box.unit_cell_parameters_from_ccp4_map=inputs.ccp4_map.unit_cell_parameters
    box.unit_cell_parameters_deduced_from_map_grid=\
       inputs.crystal_symmetry.unit_cell().parameters()

  else:
    box.unit_cell_parameters_from_ccp4_map=None
    box.unit_cell_parameters_deduced_from_map_grid=None



  if box.pdb_outside_box_msg:
    print(box.pdb_outside_box_msg, file=log)

  # NOTE: box object is always shifted to place origin at (0,0,0)

  # NOTE ON ORIGIN SHIFTS:  The shifts are applied locally here. The box
  #  object is not affected and always has the origin at (0,0,0)
  #  output_box is copy of box with shift_cart corresponding to the output
  #   files. Normally this is the same as the original shift_cart. However
  #   if user has specified a new output origin it will differ.

  # For output files ONLY:
  #  keep_origin==False leave origin at (0,0,0)
  #  keep_origin==True: we shift everything back to where it was,
  #  output_origin_grid_units=10,10,10: output origin is at (10,10,10)

  # ncs_object is original
  #  box.ncs_object is shifted by shift_cart
  #  output_box.ncs_object is shifted back by -new shift_cart

  # Additional note on output unit_cell and grid_units.
  # The ccp4-style output map can specify the unit cell and grid units
  #  corresponding to that cell.  This can be separate from the origin and
  #  number of grid points in the map as written.  If specified, write these
  #  out to the output ccp4 map and also use this unit cell for writing
  #  any output PDB files

  from copy import deepcopy
  output_box=deepcopy(box)  # won't use box below here except to return it


  print("\nBox cell dimensions: (%.2f, %.2f, %.2f) A" %(
      box.box_crystal_symmetry.unit_cell().parameters()[:3]), file=log)
  if box.shift_cart:
     print("Working origin moved from grid position of"+\
        ": (%d, %d, %d) to (0,0,0) " %(
        tuple(box.origin_shift_grid_units(reverse=True))), file=log)
     print("Working origin moved from  coordinates of:"+\
        " (%.2f, %.2f, %.2f) A to (0,0,0)\n" %(
         tuple(-col(box.shift_cart))), file=log)

  if (params.keep_origin):
    print("\nRestoring original position for output files", file=log)
    print("Origin will be at grid position of"+\
        ": (%d, %d, %d) " %(
        tuple(box.origin_shift_grid_units(reverse=True))), file=log)
    print("\nOutput files will be in same location as original", end=' ', file=log)
    if not params.keep_map_size:
      print("just cut out.", file=log)
    else:
      print("keeping entire map", file=log)
    print("Note that output maps are only valid in the cut out region.\n", file=log)

  else:
    if origin_to_match:
      output_box.shift_cart=shift_cart_for_origin_to_match
      if params.output_origin_grid_units:
        print("Output map origin to be shifted to match target", file=log)
      print("Placing origin at grid point (%s, %s, %s)" %(
        origin_to_match)+"\n"+ \
        "Final coordinate shift for output files: (%.2f,%.2f,%.2f) A\n" %(
        tuple(col(output_box.shift_cart)-col(box.shift_cart))), file=log)
    elif box.shift_cart:
      output_box.shift_cart=(0,0,0) # not shifting back
      print("Final origin will be at (0,0,0)", file=log)
      print("Final coordinate shift for output files: (%.2f,%.2f,%.2f) A\n" %(
        tuple(col(output_box.shift_cart)-col(box.shift_cart))), file=log)
    else:
      print("\nOutput files are in same location as original and origin "+\
        "is at (0,0,0)\n", file=log)

  ph_output_box_output_location = ph_box.deep_copy()
  if output_box.shift_cart:  # shift coordinates and NCS back by shift_cart
    # NOTE output_box.shift_cart could be different than box.shift_cart if
    #  there is a target position for the origin and it is not the same as the
    #  original origin.
    sites_cart = output_box.shift_sites_cart_back(
      output_box.xray_structure_box.sites_cart())
    xrs_offset = ph_output_box_output_location.extract_xray_structure(
        crystal_symmetry=output_box.xray_structure_box.crystal_symmetry()
          ).replace_sites_cart(new_sites = sites_cart)
    ph_output_box_output_location.adopt_xray_structure(xrs_offset)

    if output_box.ncs_object:
      output_box.ncs_object=output_box.ncs_object.coordinate_offset(
         tuple(-col(output_box.shift_cart)))
    shift_back=True
  else:
    shift_back=False

  if params.keep_input_unit_cell_and_grid and \
       (input_unit_cell_grid is not None) and \
       (input_unit_cell is not None):
    params.output_unit_cell=input_unit_cell
    params.output_unit_cell_grid=input_unit_cell_grid
    print("Setting output unit cell parameters and unit cell grid to"+\
      " match\ninput map file", file=log)

  if params.output_unit_cell: # Set output unit cell parameters
    from cctbx import crystal
    output_crystal_symmetry=crystal.symmetry(
      unit_cell=params.output_unit_cell, space_group="P1")
    output_unit_cell=output_crystal_symmetry.unit_cell()
    print("Output unit cell set to: %.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" %tuple(
        output_crystal_symmetry.unit_cell().parameters()), file=log)
  else:
    output_crystal_symmetry=None

  # =============  Check/set output unit cell grid and cell parameters =======
  if params.output_unit_cell_grid or output_crystal_symmetry:
    if params.output_unit_cell_grid:
      output_unit_cell_grid=params.output_unit_cell_grid
    else:
      output_unit_cell_grid=output_box.map_box.all()
    print("Output unit cell grid set to: (%s, %s, %s)" %tuple(
        output_unit_cell_grid), file=log)

    expected_output_abc=[]
    box_spacing=[]
    output_spacing=[]
    box_abc=output_box.xray_structure_box.\
         crystal_symmetry().unit_cell().parameters()[:3]
    if output_crystal_symmetry:
      output_abc=output_crystal_symmetry.unit_cell().parameters()[:3]
    else:
      output_abc=[None,None,None]
    for a_box,a_output,n_box,n_output in zip(
        box_abc,
        output_abc,
        output_box.map_box.all(),
        output_unit_cell_grid):
      expected_output_abc.append(a_box*n_output/n_box)
      box_spacing.append(a_box/n_box)
      if output_crystal_symmetry:
        output_spacing.append(a_output/n_output)
      else:
        output_spacing.append(a_box/n_box)

    if output_crystal_symmetry: # make sure it is compatible...
      r0=expected_output_abc[0]/output_abc[0]
      r1=expected_output_abc[1]/output_abc[1]
      r2=expected_output_abc[2]/output_abc[2]
      from libtbx.test_utils import approx_equal
      if not approx_equal(r0,r1,eps=0.001) or not approx_equal(r0,r2,eps=0.001):
        print("WARNING: output_unit_cell and cell_grid will "+\
          "change ratio of grid spacing.\nOld spacings: "+\
         "(%.2f, %.2f, %.2f) A " %(tuple(box_spacing))+\
        "\nNew spacings:  (%.2f, %.2f, %.2f) A \n" %(tuple(output_spacing)), file=log)
    else:
      output_abc=expected_output_abc

    from cctbx import crystal
    output_crystal_symmetry=crystal.symmetry(
          unit_cell=list(output_abc)+[90,90,90], space_group="P1")
    print("Output unit cell will be:  (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)\n"%(
       tuple(output_crystal_symmetry.unit_cell().parameters())), file=log)

  else:
    output_unit_cell_grid = output_box.map_box.all()
    output_crystal_symmetry=output_box.xray_structure_box.crystal_symmetry()
  # ==========  Done check/set output unit cell grid and cell parameters =====

  if write_output_files:
    # Write PDB file
    if ph_box.overall_counts().n_residues>0:

      if(params.output_file_name_prefix is None):
        file_name = "%s_box.pdb"%output_prefix
      else: file_name = "%s.pdb"%params.output_file_name_prefix
      ph_output_box_output_location.write_pdb_file(file_name=file_name,
          crystal_symmetry = output_crystal_symmetry)
      print("Writing boxed PDB with box unit cell to %s" %(
          file_name), file=log)

    # Write NCS file if NCS
    if output_box.ncs_object and output_box.ncs_object.max_operators()>0:
      if(params.output_file_name_prefix is None):
        output_symmetry_file = "%s_box.ncs_spec"%output_prefix
      else:
        output_symmetry_file = "%s.ncs_spec"%params.output_file_name_prefix
      output_box.ncs_object.format_all_for_group_specification(
          file_name=output_symmetry_file)

      print("\nWriting symmetry to %s" %( output_symmetry_file), file=log)

    # Write ccp4 map.
    if("ccp4" in params.output_format):
     if(params.output_file_name_prefix is None):
       file_name = "%s_box.ccp4"%output_prefix
     else: file_name = "%s.ccp4"%params.output_file_name_prefix
     from iotbx.mrcfile import create_output_labels
     if params.extract_unique:
       program_name='map_box using extract_unique'
       limitations=["extract_unique"]
     else:
       program_name='map_box'
       limitations=[]
     labels=create_output_labels(program_name=program_name,
       input_file_name=inputs.ccp4_map_file_name,
       input_labels=input_map_labels,
       limitations=limitations,
       output_labels=params.output_map_labels)
     output_box.write_ccp4_map(file_name=file_name,
       output_crystal_symmetry=output_crystal_symmetry,
       output_mean=params.output_ccp4_map_mean,
       output_sd=params.output_ccp4_map_sd,
       output_unit_cell_grid=output_unit_cell_grid,
       shift_back=shift_back,
       output_map_labels=labels,
       output_external_origin=params.output_external_origin)
     print("Writing boxed map "+\
          "to CCP4 formatted file:   %s"%file_name, file=log)

    # Write xplor map.  Shift back if keep_origin=True
    if("xplor" in params.output_format):
     if(params.output_file_name_prefix is None):
       file_name = "%s_box.xplor"%output_prefix
     else: file_name = "%s.xplor"%params.output_file_name_prefix
     output_box.write_xplor_map(file_name=file_name,
         output_crystal_symmetry=output_crystal_symmetry,
         output_unit_cell_grid=output_unit_cell_grid,
         shift_back=shift_back,)
     print("Writing boxed map "+\
         "to X-plor formatted file: %s"%file_name, file=log)

    # Write mtz map coeffs.  Shift back if keep_origin=True
    if("mtz" in params.output_format):
     if(params.output_file_name_prefix is None):
       file_name = "%s_box.mtz"%output_prefix
     else: file_name = "%s.mtz"%params.output_file_name_prefix

     print("Writing map coefficients "+\
         "to MTZ file: %s"%file_name, file=log)
     if(map_coeff is not None):
       d_min = map_coeff.d_min()
     elif params.resolution is not None:
       d_min = params.resolution
     else:
       d_min = maptbx.d_min_from_map(map_data=output_box.map_box,
         unit_cell=output_box.xray_structure_box.unit_cell())
     output_box.map_coefficients(d_min=d_min,
       scale_max=params.scale_max,
       resolution_factor=params.resolution_factor, file_name=file_name,
       shift_back=shift_back)

  print(file=log)
  return box

# =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
from wxGUI2 import utils

def validate_params(params):
  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    result = run(args=self.args, log=sys.stdout)
    return result

# =============================================================================

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
