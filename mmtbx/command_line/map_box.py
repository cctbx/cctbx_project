"""Simple cut out map around a PDB file"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.map_box

import mmtbx.utils
import mmtbx.model
from libtbx import adopt_init_args
from mmtbx.refinement import print_statistics
import libtbx.phil
from libtbx.utils import Sorry
import os, sys
from iotbx import reflection_file_utils
from cctbx import maptbx
from scitbx.matrix import col
from six.moves import zip
from iotbx.map_manager import map_manager

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
    .type = path
    .style = file_type:ccp4_map bold input_file
  mask_file_name = None
    .help = Input mask file (CCP4/mrc format).
    .short_caption = Input mask file
    .type = str
  target_ncs_au_file = None
    .help = File with model indicating which au to choose in extract_unique
    .short_caption = Input target ncs au file
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
    .help = If model is supplied, a box of density will be cut out around\
            the input model (after selections are applied to the model). \
            The size of the box of density will be box_cushion bigger than \
            the model.  Box cushion also applied if density_select is set.\
            Box cushion is also used in extract_unique.
    .short_caption = Box cushion

  mask_atoms = False
    .type = bool
    .help = Set map values to 0 outside molecular mask
    .short_caption = Mask atoms

  mask_atoms_atom_radius = 3.
    .type = float
     .help = Radius for masking around atoms
     .short_caption = Mask atoms atom radius

  write_mask_file = False
     .type = bool
     .help = Write mask file. Requires setting mask_atoms
     .short_caption = Write mask file

  set_outside_to_mean_inside = False
    .type = bool
    .help = Set value outside mask equal to mean inside
    .short_caption = Set outside to mean inside

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
    .type = choice(multi = True)
    .help = Output format(s) for boxed map. Note that mtz format is only\
            available if keep_origin = False or keep_map_size = True. (These \
            are the cases where the map is cut down to size and placed \
            at the origin or there is a full unit cell of data.)
    .short_caption = Output format

  output_file_name_prefix = None
    .type = str
    .help = Prefix for output file names. Default is name of the pdb file \
            without the ".pdb" suffix.
    .short_caption = Output file name prefix

  mask_select = False
    .type = bool
    .help = Select boundaries (min, max in x, y, z) based on auto-mask
    .short_caption = Mask select

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
            automatically specify symmetry = ALL.
    .short_caption = Symmetry
  symmetry_file = None
    .type = path
    .help = Symmetry file.\
            Symmetry or symmetry_file required if extract_unique = True.  \
            May be a \
            Phenix .ncs_spec file or BIOMTR records or a resolve ncs file.
    .short_caption = Symmetry file

  sequence_file = None
    .type = path
    .help = Sequence file (any standard format). Can be unique part or \
            all copies.  Used in identification of unique part of map \
            and in masking with mask_select
    .short_caption = Sequence file (optional)

  molecular_mass = None
    .type = float
    .help = Molecular mass of object in map in Da (i.e., 33000 for 33 Kd).\
              Used in identification \
            of unique part of map and in masking by mask_select.
    .short_caption = Molecular mass (optional)

  solvent_content = None
    .type = float
    .help = Optional fraction of volume of map that is empty.  \
            Used in identification \
            of unique part of map and in masking by mask_select
    .short_caption = Solvent content

  extract_unique = False
    .type = bool
    .help = Extract unique part of map. Requires either symmetry_file or \
            symmetry and\
            either sequence file or molecular mass to be supplied. If chain \
            type is not protein it should be set as well.
    .short_caption = Extract unique

  use_cubic_boxing = None
    .type = bool
    .help = When boxing map, create a box that is cubic (same number \
             of grid points in each direction). Can be used with \
             stay_inside_current_map to keep the box inside the current \
             map.  Note that if box is bigger than the current map and \
             wrapping is False, any points outside original \
             map will be set to zero. Can be used with \
             extract_unique, mask_select, density_select and default\
             boxing around molecule.  Incompatible with lower_bounds \
             and upper_bounds.

  stay_inside_current_map = None
    .type = bool
    .help = If True, then when applying boxing, keep the new box inside the \
             current map.  If False, allow new box to go outside the \
             current map. Values outside will be zero if wrapping=False \
             as determined from input map or the parameter wrapping, \
             and wrapped values if wrapping= True

  increase_box_cushion_and_atom_radius_for_soft_mask = True
    .type = bool
    .help = Expand cushion and atom radii by soft_mask_radius
    .short_caption = Increase box cushion and atom radius for soft mask

  soft_mask_extract_unique = True
    .type = bool
    .help = Create soft mask at edges of extract_unique box (feather map into \
               edge of box). Uses resolution as mask_radius
      .short_caption = Soft mask in extract unique

  mask_expand_ratio = 1
      .type = float
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
    .type = bool
    .help = Use Gaussian mask in mask_atoms and on outside surface of box
    .short_caption = Soft mask

  invert_mask = False
    .type = bool
    .help = Make mask with 1 outside model (only applies to mask_around_atoms)
    .short_caption = Invert mask to outside atoms

  soft_mask_radius = 3
    .type = float
    .help = Gaussian mask smoothing radius
    .short_caption = Soft mask radius

  lower_bounds = None
    .type = ints
    .help = Lower bounds for cut out box. You can specify them directly.\
            NOTE: lower and upper bounds refer to grid points after shifting \
            the map to place the origin at (0, 0, 0). To refer to absolute \
            values specify bounds_are_absolute = True.
    .short_caption = Lower bounds

  upper_bounds = None
    .type = ints
    .help = Upper bounds for cut out box.  You can specify them directly.\
            NOTE: lower and upper bounds refer to grid points after shifting \
            the map to place the origin at (0, 0, 0). To refer to absolute \
            values specify bounds_are_absolute = True.
    .short_caption = Upper bounds

  bounds_are_absolute = False
    .type = bool
    .help = Define lower and upper bounds as absolute. \
            NOTE: lower and upper bounds refer to grid points after shifting \
            the map to place the origin at (0, 0, 0). To refer to absolute \
            values specify bounds_are_absolute = True.
    .short_caption = Bounds are absolute

  keep_map_size = False
    .type = bool
    .help = Keep original map gridding (do not cut anything out). \
            Use to apply soft_mask and/or mask_atoms keeping same map size.
    .short_caption = Keep map size

  keep_origin = True
    .type = bool
    .help = Write out map, map_coefficients, and model \
            with origin in original location.  \
            If False, shift the origin to (0, 0, 0).  \
            NOTE: to cut out a part of a map, shift the origin to (0, 0, 0), \
               and make a new small map use keep_origin = False\
               keep_input_unit_cell_and_grid = False
    .short_caption = Keep origin

  keep_input_unit_cell_and_grid = True
     .type = bool
     .help = Keep the input unit_cell dimensions and unit_cell_grid. \
             If False, use the dimensions and grid of the cut out box as the \
              unit cell map_box dimensions and grid.\
            NOTE: to cut out a part of a map, shift the origin to (0, 0, 0), \
               and make a new small map set keep_origin = False and \
               keep_input_unit_cell_and_grid = False
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
           should use keep_origin = True or False to specify your origin \
           but if you want to move it to a specific grid point you can do that.\
    .short_caption = Output origin
    .expert_level = 3

  output_origin_match_this_file = None
    .type = path
    .help = As output_origin_grid_units, but use origin from this file
    .short_caption = File with origin info

  bounds_match_this_file = None
    .type = path
    .help = Take the lower and upper bounds from this map file and apply them \
             to the input map file.
    .short_caption = File with bounds to match

  output_external_origin = None
    .type = floats
    .help = Write ORIGIN record to map file (this is an external origin \
            used to specify relationship to external files such as model files.\
            Three floating point numbers (A).
    .short_caption = output external origin

  restrict_map_size = False
    .type = bool
    .help = Do not go outside original map boundaries
    .short_caption = Restrict map size

  ignore_symmetry_conflicts = False
    .type = bool
    .help = Ignore unit cell from model if it conflicts with the map.
    .short_caption = Ignore symmetry conflicts

  wrapping = False
    .type = bool
    .help = If wrapping, map wraps around at map boundaries.
    .short_caption = Wrapping

  check_wrapping = False
    .type = bool
    .help = Check that wrapping is consistent with map if it is set to True
    .short_caption = Check wrapping

  output_ccp4_map_mean = None
    .type = float
    .help = Choose mean and SD of output CCP4 map
    .short_caption = Mean of output CCP4 map

  output_ccp4_map_sd = None
    .type = float
    .help = Choose mean and SD of output CCP4 map
    .short_caption = SD of output CCP4 map

  output_map_label = None
    .type = str
    .multiple = True
    .help = Add this label to output map
    .short_caption = Add label

  remove_output_map_labels = None
    .type = bool
    .help = Remove all output map labels
    .short_caption = Remove labels

  invert_hand = False
    .type = bool
    .help = Just before writing out the map, swap the order of all sections \
             in Z.  This will change the hand of the map. Note that this\
             removes any correspondence to models (these are not inverted). \
             If you use this, be sure to apply it to all your starting maps.\
    .short_caption = Invert hand of map

  gui
    .help = "GUI-specific parameter required for output directory"
  {
    output_dir = None
    .type = path
    .style = output_dir
  }
""", process_includes = True)

master_params = master_phil

def remove_element(text_list, element = None):
    new_text_list = []
    for x in text_list:
      if x !=  element:
        new_text_list.append(x)
    return new_text_list

def get_model_from_inputs(
    model = None,
    pdb_hierarchy = None,
    file_names = None,
    crystal_symmetry = None,
    log = sys.stdout):
  print_statistics.make_sub_header("pdb model", out = log)

  if pdb_hierarchy:  # convert to model object . XXX should come in this way
    model =  pdb_hierarchy.as_model_manager(crystal_symmetry = crystal_symmetry)

  if len(file_names)>0:
    file_name = file_names[0]
    if not file_name or not os.path.isfile(file_name):
      raise Sorry("The file %s is missing" %(file_name))
    print("Reading model from %s" %(file_name), file = log)
    from iotbx.data_manager import DataManager
    dm = DataManager()
    dm.set_overwrite(True)
    dm.process_model_file(file_name)
    model = dm.get_model(file_name)
    if crystal_symmetry and not model.crystal_symmetry():
      model.set_crystal_symmetry(crystal_symmetry)

  return model

def get_map_manager_objects(
    params = None,
    inputs = None,
    crystal_symmetry = None,
    ccp4_map = None,
    mask_as_map_manager = None,
     map_data = None,  # XXX
     mask_data = None, # XXX
     log = sys.stdout):

  # Map or map coefficients
  map_coeff = None
  resolution_from_map_coeffs = None
  input_unit_cell_grid = None
  input_unit_cell = None
  input_map_labels = None
  map_or_map_coeffs_prefix = None

  if map_data and not ccp4_map:  # convert to map_manager
    # Called with map_data.  We do not know for sure if map_data is
    #  wrapped or not. Require wrapping to be set to define it.

    assert isinstance(params.wrapping, bool)
    ccp4_map = map_manager(map_data = map_data,
      unit_cell_grid = map_data.all(),
      unit_cell_crystal_symmetry = crystal_symmetry,
      wrapping = params.wrapping)
  elif (not ccp4_map):

    # read first mtz file
    if ( (len(inputs.reflection_file_names) > 0) or
         (params.map_coefficients_file is not None) ):
      # Here with MTZ input, wrapping default is True (crystallographic map)
      if not isinstance(params.wrapping, bool):
        params.wrapping = True
      # file in phil takes precedent
      if (params.map_coefficients_file is not None):
        if (len(inputs.reflection_file_names)  ==  0):
          inputs.reflection_file_names.append(params.map_coefficients_file)
        else:
          inputs.reflection_file_names[0] = params.map_coefficients_file
      map_coeff = reflection_file_utils.extract_miller_array_from_file(
        file_name = inputs.reflection_file_names[0],
        label     = params.label,
        type      = "complex",
        log       = log)
      resolution_from_map_coeffs = map_coeff.d_min()
      if not crystal_symmetry: crystal_symmetry = map_coeff.crystal_symmetry()
      fft_map = map_coeff.fft_map(resolution_factor = params.resolution_factor)
      fft_map.apply_sigma_scaling()
      map_data = fft_map.real_map_unpadded()
      map_or_map_coeffs_prefix = os.path.basename(
         inputs.reflection_file_names[0][:-4])
      # Convert map_data to map_manager object
      ccp4_map = map_manager(map_data = map_data,
        unit_cell_grid = map_data.all(),
        unit_cell_crystal_symmetry = crystal_symmetry,
        wrapping = params.wrapping)

    # or read CCP4 map
    elif ( (inputs.ccp4_map is not None) or
           (params.ccp4_map_file is not None) ):
      # Here wrapping comes from map file; no need to set default. If not
      #  specified in map labels, wrapping will be False for map file.

      if (params.ccp4_map_file is not None):
        inputs.ccp4_map = read_map_file_with_data_manager(params.ccp4_map_file)
        inputs.ccp4_map_file_name = params.ccp4_map_file
      print_statistics.make_sub_header("CCP4 map", out = log)
      ccp4_map = inputs.ccp4_map
      ccp4_map.show_summary(prefix = "  ", out = log)
      if not crystal_symmetry: crystal_symmetry = ccp4_map.crystal_symmetry()
      map_data = ccp4_map.map_data()
      input_unit_cell_grid = ccp4_map.unit_cell_grid
      input_unit_cell = ccp4_map.unit_cell().parameters()
      input_map_labels = ccp4_map.get_labels()

      if inputs.ccp4_map_file_name.endswith(".ccp4"):
        map_or_map_coeffs_prefix = os.path.basename(
          inputs.ccp4_map_file_name[:-5])
      else:
        map_or_map_coeffs_prefix = os.path.basename(
          inputs.ccp4_map_file_name[:-4])

  if mask_data and (not mask_as_map_manager):
    mask_as_map_manager = map_manager(map_data = mask_data.as_double(),
        unit_cell_grid = mask_data.all(),
        unit_cell_crystal_symmetry = crystal_symmetry,
        wrapping = params.wrapping)
  if (not mask_as_map_manager) and params.mask_file_name:
    mask_as_map_manager = read_map_file_with_data_manager(
        params.mask_file_name)

  if len(inputs.pdb_file_names)>0:
    output_prefix = os.path.basename(inputs.pdb_file_names[0])[:-4]
  else:
    output_prefix = map_or_map_coeffs_prefix

  return ccp4_map, mask_as_map_manager, \
    output_prefix, resolution_from_map_coeffs

def read_map_file_with_data_manager(file_name):
    from iotbx.data_manager import DataManager
    dm = DataManager()
    dm.set_overwrite(True)
    dm.process_real_map_file(file_name)
    return dm.get_real_map(file_name)

def check_parameters(inputs = None, params = None,
   model = None,
   ncs_object = None,
   pdb_hierarchy = None,
   log = sys.stdout):

  if(len(inputs.pdb_file_names)!= 1 and not params.density_select and not
    params.mask_select and not pdb_hierarchy and
     not model and not params.keep_map_size and not params.upper_bounds
     and not params.extract_unique and not params.bounds_match_this_file):
    raise Sorry("PDB file is needed unless extract_unique, "+
      "density_select, mask_select, keep_map_size \nor bounds are set .")
  if (len(inputs.pdb_file_names)!= 1 and not pdb_hierarchy and not model and\
       (params.mask_atoms )):
    raise Sorry("PDB file is needed for mask_atoms")
  if params.soft_mask and (not params.resolution) and \
        (len(inputs.pdb_file_names)!= 1 and not pdb_hierarchy and not model):
    raise Sorry("Need resolution for soft_mask without PDB file")
  if (params.density_select and params.extract_unique):
    raise Sorry("Cannot set both density_select and extract_unique")
  if ((params.density_select or params.mask_select) and params.keep_map_size):
    raise Sorry("Cannot set both density_select/mask_select and keep_map_size")
  if ((params.density_select or params.mask_select) and params.upper_bounds):
    raise Sorry("Cannot set both density_select/mask_select and bounds")
  if (params.keep_map_size and params.upper_bounds):
    raise Sorry("Cannot set both keep_map_size and bounds")
  if (params.upper_bounds and not params.lower_bounds):
    raise Sorry("Please set lower_bounds if you set upper_bounds")
  if (params.lower_bounds and params.use_cubic_boxing):
    raise Sorry("Cubic boxing is not compatible with "+
       "setting lower or upper bounds")
  if (params.extract_unique):
    if (not params.resolution):
      raise Sorry("Please set resolution for extract_unique")
    if (not params.symmetry) and (not params.symmetry_file) and \
        (not ncs_object):
      params.symmetry="ALL"
      print("Setting symmetry=ALL as no symmetry information supplied",
        file = log)
    if params.mask_atoms:
      raise Sorry("You cannot set mask_atoms with extract_unique")

  if params.keep_input_unit_cell_and_grid and (
      (params.output_unit_cell_grid is not None ) or
      (params.output_unit_cell is not None ) ):
    raise Sorry("If you set keep_input_unit_cell_and_grid then you cannot "+\
       "set \noutput_unit_cell_grid or output_unit_cell")

def print_default_message(log = sys.stdout):
  h = "phenix.map_box: extract box with model and map around selected atoms"
  print_statistics.make_header(h, out = log)
  default_message = """\

%s.

Usage:
  phenix.map_box model.pdb map_coefficients.mtz selection = "chain A and resseq 1:10"

or

  phenix.map_box map.ccp4 density_select = True

Parameters:"""%h

def process_inputs(args = None,
  crystal_symmetry = None,
  log = sys.stdout):

  # Process inputs ignoring symmetry conflicts

  inputs = mmtbx.utils.process_command_line_args(args = args,
      cmd_cs = crystal_symmetry,
      master_params = master_phil,
      suppress_symmetry_related_errors = True)
  params = inputs.params.extract()

  master_phil.format(python_object = params).show(out = log)

  return inputs, params

def get_origin_or_bounds_from_ccp4_file(params = None, log = sys.stdout):
  if params.output_origin_match_this_file:
    fn = params.output_origin_match_this_file
    if params.bounds_match_this_file:
      raise Sorry("Cannot match origin and bounds at same time")
  else:
    fn = params.bounds_match_this_file
  if not params.ccp4_map_file:
    raise Sorry(
     "Need to specify your input file with ccp4_map_file = xxx if you use "+
      "output_origin_match_this_file = xxxx or bounds_match_this_file = xxxx")

  ccp4_map = read_map_file_with_data_manager(fn)

  if (ccp4_map):
    origin = ccp4_map.map_data().origin()
    if params.output_origin_match_this_file:
      params.output_origin_grid_units = origin
      print("Origin of (%s, %s, %s) taken from %s" %(
         origin[0], origin[1], origin[2], fn))
    else:
      all = ccp4_map.map_data().all()
      params.lower_bounds = origin
      print("Lower bounds of (%s, %s, %s) taken from %s" %(
         params.lower_bounds[0], params.lower_bounds[1],
           params.lower_bounds[2], fn))
      params.upper_bounds = list(col(origin)+col(all)-col((1, 1, 1)))
      print("upper bounds of (%s, %s, %s) taken from %s" %(
         params.upper_bounds[0], params.upper_bounds[1],
          params.upper_bounds[2], fn))
      params.bounds_are_absolute = True
  else:
    raise Sorry("Unable to interpret %s as map file" %(fn))

  return params

def modify_params(params = None,
    inputs = None,
    model = None,
    pdb_hierarchy = None,
    write_output_files = None,
    upper_bounds = None,
    lower_bounds = None,
    wrapping = None,
    log = sys.stdout):

  #  Update wrapping if specified
  if isinstance(wrapping, bool):
    params.wrapping = wrapping

  # PDB file
  if params.pdb_file and not inputs.pdb_file_names and not pdb_hierarchy \
      and not model:
    inputs.pdb_file_names = [params.pdb_file]

  # Overwrite params with parameters in call if available
  if lower_bounds:
     params.lower_bounds = lower_bounds
  if upper_bounds:
     params.upper_bounds = upper_bounds

  if (write_output_files) and ("mtz" in params.output_format) and (
       (params.keep_origin) and (not params.keep_map_size)):
    print("\nNOTE: Skipping write of mtz file as keep_origin = True and \n"+\
       "keep_map_size is False\n", file = log)
    params.output_format = remove_element(params.output_format, element = 'mtz')

  if params.output_external_origin and (not params.keep_origin):
    raise Sorry(
      "If you specify an external origin you must set keep_origin=True")

  if (write_output_files) and ("mtz" in params.output_format) and (
       (params.extract_unique)):
    print("\nNOTE: Skipping write of mtz file as extract_unique = True\n",
      file = log)
    params.output_format = remove_element(params.output_format, element = 'mtz')


  # XXX Get origin or bounds from a ccp4 file Instead use data_manager to read

  if params.output_origin_match_this_file or params.bounds_match_this_file:
    params = get_origin_or_bounds_from_ccp4_file(params = params, log = log)


  # Check that bounds are None or tuples of three

  if params.lower_bounds and len(params.lower_bounds) != 3:
    raise Sorry("Need 3 values for lower_bounds")
  if params.upper_bounds and len(params.upper_bounds) != 3:
    raise Sorry("Need 3 values for upper_bounds")

  if params.output_origin_grid_units is not None and params.keep_origin:
    params.keep_origin = False
    print("Setting keep_origin = False as output_origin_grid_units is set",
       file = log)

  if params.soft_mask_radius is None:
    params.soft_mask_radius = params.resolution

  if (params.soft_mask and params.mask_atoms and
      params.increase_box_cushion_and_atom_radius_for_soft_mask):
    params.box_cushion+= params.mask_atoms_atom_radius
    print ("Increasing box_cushion by mask_atoms_atom_radius for soft mask",
       file = log)

  return params

def get_origin_to_match(
   params = None,
   n_real = None,
   crystal_symmetry = None):

  if params.output_origin_grid_units is not None:
    origin_to_match = tuple(params.output_origin_grid_units)
  else:
    origin_to_match = None

  if origin_to_match:
    sc = []
    for x, o, a in zip(crystal_symmetry.unit_cell().parameters()[:3],
        origin_to_match, n_real):
      sc.append(-x*o/a)
    shift_cart_for_origin_to_match = tuple(sc)
  else:
    origin_to_match = None
    shift_cart_for_origin_to_match = None

  return origin_to_match, shift_cart_for_origin_to_match

def apply_selection_to_model(params = None, model = None, log = sys.stdout):
  if not model:
     return

  if not params.selection or params.selection == "all":
    return model

  selection = model.selection(params.selection)
  model = model.select(selection)
  return model

def get_ncs_object(params = None,
    ncs_object = None,
    log = sys.stdout):
  if not ncs_object:
    from mmtbx.ncs.ncs import ncs
    ncs_object = ncs()
    if params.symmetry_file:
      ncs_object.read_ncs(params.symmetry_file, log = log)
      print("Total of %s operators read" %(ncs_object.max_operators()), file = log)
  if ncs_object.max_operators()<1:
      print("No symmetry available", file = log)

  return ncs_object

def get_sequence_and_molecular_mass(params = None, ncs_object = None,
   log = sys.stdout):

  if (not params.extract_unique) and (not params.mask_select):
    return  params, None  # no sequence

  if params.sequence_file:
    if ncs_object.max_operators()> 1: # get unique part of sequence
      remove_duplicates = True
    else:
      remove_duplicates = False
    from iotbx.bioinformatics import get_sequences
    sequence = (" ".join(get_sequences(file_name = params.sequence_file,
      remove_duplicates = remove_duplicates)))
  else:
    sequence = None

  if params.chain_type in ['None', None]: params.chain_type = None
  if sequence and not params.molecular_mass:
    # get molecular mass from sequence
    from iotbx.bioinformatics import text_from_chains_matching_chain_type
    if params.chain_type in [None, 'PROTEIN']:
      n_protein = len(text_from_chains_matching_chain_type(
        text = sequence, chain_type = 'PROTEIN'))
    else:
      n_protein = 0
    if params.chain_type in [None, 'RNA']:
      n_rna = len(text_from_chains_matching_chain_type(
        text = sequence, chain_type = 'RNA'))
    else:
      n_rna = 0
    if params.chain_type in [None, 'DNA']:
      n_dna = len(text_from_chains_matching_chain_type(
       text = sequence, chain_type = 'DNA'))
    else:
      n_dna = 0
    params.molecular_mass = ncs_object.max_operators()*(
         n_protein*110+(n_rna+n_dna)*330)
    print("\nEstimate of molecular mass is %.0f " %(
        params.molecular_mass), file = log)
  return params, sequence

def print_what_will_happen(
   params = None,
   model = None,
   log = sys.stdout):


  if params.density_select or params.mask_select:
    print_statistics.make_sub_header(
    "Extracting box around selected density and writing output files", out = log)
  else:
   print_statistics.make_sub_header(
    "Extracting box around selected atoms and writing output files", out = log)
  #
  if params.set_outside_to_mean_inside:
    print("\nValue outside atoms mask will be set to mean inside mask", file = log)
  if params.get_half_height_width and params.density_select:
    print("\nHalf width at half height will be used to id boundaries", file = log)

  if params.soft_mask and model and \
      model.get_xray_structure().sites_cart().size()>0:
    print("\nSoft mask will be applied to model-based mask", file = log)
  elif params.soft_mask:
    print ("\nSoft mask will be applied to outside of map box", file = log)
  if params.keep_map_size:
    print("\nEntire map will be kept (not cutting out region)", file = log)
  if params.restrict_map_size:
    print("\nOutput map will be within input map", file = log)
  if params.lower_bounds and params.upper_bounds:
    print("Bounds for cut out map are (%s, %s, %s) to (%s, %s, %s)" %(
     tuple(list(params.lower_bounds)+list(params.upper_bounds))), file = log)

  if params.use_cubic_boxing:
    print("Output boxed map will be a cubic box", file = log)
    if params.stay_inside_current_map:
      print("Cubic map will be inside current map")
    elif params.wrapping:
      print("Cubic map may go outside current map;"+
           " values outside will be wrapped")
    else:
      print("Cubic map may go outside current map;"+
           " values outside will be set to zero")

def print_notes(params = None,
    mam = None,
    crystal_symmetry = None,
    ccp4_map = None,
    log = sys.stdout):

  if params.mask_select and hasattr(mam, 'get_solvent_content') and \
        mam.get_solvent_content():
    print("\nSolvent content used in mask_select: %.3f " %(
      mam.get_solvent_content()), file = log)

  if (ccp4_map and
    crystal_symmetry and
    crystal_symmetry.unit_cell().parameters() and
     ccp4_map.unit_cell().parameters()  ) and (
       crystal_symmetry.unit_cell().parameters() !=
       ccp4_map.unit_cell().parameters()):
    print("\nNOTE: Input CCP4 map is only part of unit cell:", file = log)
    print("Full unit cell ('unit cell parameters'): "+\
      "(%.1f, %.1f, %.1f, %.1f, %.1f, %.1f) A" %tuple(
        ccp4_map.unit_cell().parameters()), file = log)
    print("Size of CCP4 map 'map unit cell':        "+\
      "(%.1f, %.1f, %.1f, %.1f, %.1f, %.1f) A" %tuple(
       crystal_symmetry.unit_cell().parameters()), file = log)
    print("Full unit cell as grid units: (%s, %s, %s)" %(
      tuple(ccp4_map.unit_cell_grid)), file = log)
    print("Map unit cell as grid units:  (%s, %s, %s)" %(
      tuple(ccp4_map.map_data().all())), file = log)

    if params.invert_hand:
      print("\nOutput map will be inverted hand "+
         "(swapping order of sections in Z)", file = log)

def run(args,
     ncs_object = None,  # ncs object
     model = None,  # model.manager object
     ccp4_map = None,  # map_manager object
     mask_as_map_manager = None, # map_manager object
     crystal_symmetry = None,  # XXX remove
     pdb_hierarchy = None, #XXX remove
     map_data = None,  # XXX remove
     mask_data = None, # XXX remove
     lower_bounds = None,
     upper_bounds = None,
     wrapping = None,  # Alternative way to specify wrapping
     write_output_files = True,
     log = None):

  if (log is None): log = sys.stdout


  print_default_message(log = log)
  if(len(args)  ==  0 and not pdb_hierarchy):
    master_phil.show(prefix = "  ")
    return

  # Read files with file reader and get parameters
  inputs, params = process_inputs(args = args,
    crystal_symmetry = crystal_symmetry,
    log = log)


  # Custom changes in parameters based on input files and supplied bounds
  params = modify_params(params = params,
    inputs = inputs,
    model = model,
    pdb_hierarchy = pdb_hierarchy, # XXX remove later
    write_output_files = write_output_files,
    upper_bounds = upper_bounds,
    lower_bounds = lower_bounds,
    wrapping = wrapping,
    log = log)

  # Check parameters and issue error messages if necessary
  check_parameters(inputs = inputs, params = params,
    model = model,
    ncs_object = ncs_object,
    pdb_hierarchy = pdb_hierarchy, # remove later XXX
    log = log)

  # Use inputs.crystal_symmetry (precedence there is for map)
  crystal_symmetry = inputs.crystal_symmetry

  # Get map_manager objects

  # XXX get rid of most of these as they are part of mm objects
  ccp4_map, mask_as_map_manager, \
    output_prefix, resolution_from_map_coeffs = get_map_manager_objects(
      params = params,
      inputs = inputs,
      ccp4_map = ccp4_map,
      mask_as_map_manager = mask_as_map_manager,
      crystal_symmetry = crystal_symmetry,
      map_data = map_data,  # XXX delete
      mask_data = mask_data, # XXX  delete
      log = log)

  if not ccp4_map:
    raise Sorry("Need a map for map_box")

  # Apply a scale factor to map data on read-in if requested
  if params.map_scale_factor:
    print("Applying scale factor of %s to map data on read-in" %(
       params.map_scale_factor))
    ccp4_map = ccp4_map.customized_copy(
      map_data = ccp4_map.map_data()*params.map_scale_factor)

  # Use ccp4_map crystal_symmetry if not set
  if ccp4_map and not crystal_symmetry:
    crystal_symmetry = ccp4_map.unit_cell_crystal_symmetry()


  # Get model object (replaces pdb_hierarchy)
  model = get_model_from_inputs(
    model = model,
    pdb_hierarchy = pdb_hierarchy,
    file_names = inputs.pdb_file_names,
    crystal_symmetry = crystal_symmetry,
    log = log)

  # Apply selection to model if desired
  if model:
    model = apply_selection_to_model(params = params, model = model, log = log)

  # Get target model object for extract_unique if present
  if params.target_ncs_au_file:
    target_ncs_au_model = get_model_from_inputs(
      file_names = [params.target_ncs_au_file],
      crystal_symmetry = crystal_symmetry,
      log = log)
  else:
    target_ncs_au_model = None

  # final check that map_data exists
  if(ccp4_map.map_data is None):
    raise Sorry("Map or map coefficients file is needed.")

  # Set wrapping if specified:
  if params.wrapping in [True, False]:
    ccp4_map.set_wrapping(params.wrapping)

  ncs_object = get_ncs_object(params = params,
      ncs_object = ncs_object, log = log)

  # Get sequence if extract_unique or mask_select is set
  params, sequence = get_sequence_and_molecular_mass(params = params,
    ncs_object = ncs_object,
    log = log)

  # Summarize what will happen
  print_what_will_happen(params = params,
    model = None,
    log = log)

  # Run now

  # Change map/model unit_cell_crystal_symmetry if requested
  if params.output_unit_cell and \
     tuple(params.output_unit_cell)!= tuple(ccp4_map.unit_cell().parameters()):
    ccp4_map, model = change_output_unit_cell(params = params,
       ccp4_map = ccp4_map,
       model = model)


  # Decide if we are going to box at the beginning:
  box = (model and (not params.keep_map_size))
  if (params.lower_bounds and params.upper_bounds):
    box = False
  if (params.extract_unique):
    box = False
  if (params.density_select):
    box = False
  if (params.mask_select):
    box = False

  from iotbx.map_model_manager import map_model_manager
  mam = map_model_manager(
    model = model,
    map_manager = ccp4_map,
    ncs_object = ncs_object,
    ignore_symmetry_conflicts = params.ignore_symmetry_conflicts)
  if box:
    mam.box_all_maps_around_model_and_shift_origin(
      use_cubic_boxing = params.use_cubic_boxing,
      stay_inside_current_map = params.stay_inside_current_map,
      box_cushion = params.box_cushion)
    if mam.warning_message():
      print (mam.warning_message(), file = log)

  # Map and model and ncs are boxed if requested and
  #   shifted to place origin at (0, 0, 0)
  # Now box the map if desired and shift origin if requested
  # Shift bounds if bounds_are_absolute and original origin is not zero:
  if params.bounds_are_absolute:
    params.lower_bounds = tuple([lb-o for lb, o in zip(params.lower_bounds,
      mam.map_manager().origin_shift_grid_units)])
    params.upper_bounds = tuple([lb-o for lb, o in zip(params.upper_bounds,
      mam.map_manager().origin_shift_grid_units)])
    params.bounds_are_absolute = False

  if params.lower_bounds and params.upper_bounds:  # Box it
    assert not box # should not have used boxing
    assert not params.use_cubic_boxing  # not compatible
    from cctbx.maptbx.box import with_bounds
    mam = with_bounds(mam.map_manager(), # actually a box
         params.lower_bounds,
         params.upper_bounds,
         model = mam.model(),
         log = log)
    if mam.warning_message():
      print (mam.warning_message(), file = log)

  elif params.density_select:  # Box it with density_select
    assert not box # should not have used boxing
    from cctbx.maptbx.box import around_density
    mam = around_density(mam.map_manager(), # actually a box
         box_cushion = params.box_cushion,
         threshold = params.density_select_threshold,
         get_half_height_width = params.get_half_height_width,
         model = mam.model(),
         stay_inside_current_map = params.stay_inside_current_map,
         use_cubic_boxing = params.use_cubic_boxing,
         log = log)
    if mam.warning_message():
      print (mam.warning_message(), file = log)

  elif params.mask_select:  # Box it with mask_select
    assert not box # should not have used boxing
    from cctbx.maptbx.box import around_mask
    if not mask_as_map_manager: # Generate it
      mm = mam.map_manager().deep_copy()
      mm.create_mask_around_density(
        resolution = params.resolution,
        molecular_mass = params.molecular_mass,
        sequence = sequence,
        solvent_content = params.solvent_content,
        )
      cm=mm._created_mask
      mask_as_map_manager = mm.get_mask_as_map_manager()
      if not mask_as_map_manager:
        raise Sorry("Unable to auto-generate mask")

    mam = around_mask(mam.map_manager(), # actually a box, shifted
         mask_as_map_manager = mask_as_map_manager,
         box_cushion = params.box_cushion,
         model = mam.model(),
         stay_inside_current_map = params.stay_inside_current_map,
         use_cubic_boxing = params.use_cubic_boxing,
         log = log)
    if mam.warning_message():
      print (mam.warning_message(), file = log)

  # Now mask map if requested

  if (params.extract_unique):  # mask around unique part of map and rebox
    # NOTE: actually returns box not mam XXX
    mam = apply_around_unique(mam, params = params,
       sequence = sequence,
       target_ncs_au_model = target_ncs_au_model,
       log = log)

  elif (params.mask_atoms):  # mask around atoms, optionally soft
    mam = apply_mask_around_atoms(mam, params = params, log = log)

  elif (params.soft_mask):  # apply soft mask to outside of box
    mam = apply_mask_around_edge_of_box(mam, params = params, log = log)

  if params.write_mask_file:
    if not hasattr(mam.map_manager(),'_created_mask') or not \
       mam.map_manager()._created_mask:
      raise Sorry("Cannot create mask file if no mask has been created")
    mask_map_manager = mam.map_manager()._created_mask.map_manager().deep_copy()
  else:
    mask_map_manager = None

  # Shift origin of output file if requested
  if params.output_origin_grid_units or params.output_unit_cell_grid or\
       (not params.keep_origin) or (not params.keep_input_unit_cell_and_grid):

    if params.output_origin_grid_units:
      print ("Setting origin of final map to be at %s" %(
       str(params.output_origin_grid_units)), file = log)
    elif not params.keep_origin:
      print ("Setting origin of final map to be at (0, 0, 0)", file = log)
      params.output_origin_grid_units = (0, 0, 0)

    if (not params.keep_input_unit_cell_and_grid) and (
         not params.output_unit_cell_grid):
      params.output_unit_cell_grid = mam.map_manager().map_data().all()
    if params.output_unit_cell_grid:
      print ("Setting gridding of unit cell of final map to be at %s" %(
       str(params.output_unit_cell_grid)), file = log)
    mam.map_manager().set_original_origin_and_gridding(
       original_origin = params.output_origin_grid_units,
       gridding = params.output_unit_cell_grid)
    if mam.map_manager().ncs_object():
      mam.map_manager().ncs_object().set_shift_cart(
        mam.map_manager().shift_cart())

    if mam.model():
      mam.model().shift_model_and_set_crystal_symmetry(
       shift_cart = (0,0,0),
       crystal_symmetry = mam.map_manager().crystal_symmetry())
      mam.model().set_shift_cart((0,0,0)) # next line requires shift-cart=0,0,0
      mam.model().set_unit_cell_crystal_symmetry(
        mam.map_manager().crystal_symmetry())
      mam.model().set_shift_cart(mam.map_manager().shift_cart())

  if params.wrapping in [True, False] and mam.map_manager().is_full_size():
    mam.map_manager().set_wrapping(params.wrapping)
    if params.wrapping and params.check_wrapping and (
       not mam.map_manager().is_consistent_with_wrapping()):
      print("\nWARNING: This map is not consistent with wrapping but wrapping"+
        " is set to True",file = log)

  # Adjust labels
  if params.output_map_label:  # add it
    for label in params.output_map_label:
      mam.map_manager().add_label(label)

  if params.remove_output_map_labels:
    mam.map_manager().remove_labels()

  if params.invert_hand:
    mam.map_manager().invert_hand()

  # Print out any notes about the output files
  print_notes(params = params, mam = mam,
    crystal_symmetry = crystal_symmetry,
    ccp4_map = ccp4_map,
    log = log)

  # For output files ONLY:
  #  keep_origin == False leave origin at (0, 0, 0)
  #  keep_origin == True: we shift everything back to where it was,
  #  output_origin_grid_units = 10, 10, 10: output origin is at (10, 10, 10)
  #  output_external_origin = 10,10,10; set output_external_origin value

  print("\nBox cell dimensions: (%.2f, %.2f, %.2f) A" %(
      mam.map_manager().crystal_symmetry().unit_cell().parameters()[:3]),
       file = log)

  if mam.map_manager().origin_shift_grid_units:
     print("Working origin moved from grid position of"+\
        ": (%d, %d, %d) to (0, 0, 0) " %(
        tuple( mam.map_manager().origin_shift_grid_units)),
        file = log)
     print("Working origin moved from  coordinates of:"+\
        " (%.2f, %.2f, %.2f) A to (0, 0, 0)\n" %(
        tuple([-x for x in mam.map_manager().shift_cart()])),
           file = log)

  #  For now, need to return a box object
  if mam.model():
    xrs = mam.model().get_xray_structure()
    hierarchy = mam.model().get_hierarchy()
  else:
    xrs = None
    hierarchy = None
  output_box = box_object(
      shift_cart = mam.map_manager().shift_cart(),
      xray_structure_box = xrs,
      hierarchy = hierarchy,
      ncs_object = mam.map_manager().ncs_object(),
      map_box = mam.map_manager().map_data(),
      map_data = ccp4_map.map_data(),
      map_box_half_map_list = None,
      box_crystal_symmetry = mam.map_manager().crystal_symmetry(),
      pdb_outside_box_msg = "",
      gridding_first = getattr(mam, 'gridding_first', (0, 0, 0)),
      gridding_last = getattr(mam,
          'gridding_last', mam.map_manager().map_data().all()),
      solvent_content = params.solvent_content,
      origin_shift_grid_units = [
         -x for x in mam.map_manager().origin_shift_grid_units],
      )
  if write_output_files:
    model = mam.model()
    map_manager = mam.map_manager()
    ncs_object = mam.map_manager().ncs_object()
    from iotbx.data_manager import DataManager
    dm = DataManager(datatypes = [
       'model', 'ncs_spec', 'real_map', 'miller_array'])
    dm.set_overwrite(True)

    if params.output_external_origin:
      assert (isinstance(params.output_external_origin,tuple) or \
             isinstance(params.output_external_origin,list)) and \
             len(params.output_external_origin) == 3
      map_manager.set_output_external_origin(params.output_external_origin)
      print("Set output_external_origin to %s" %(
       str(params.output_external_origin)), file = log)

    # Write PDB file
    if model:
      if(params.output_file_name_prefix is None):
        filename = "%s_box"%output_prefix
      else: filename = "%s"%params.output_file_name_prefix
      full_filename = dm.write_model_file(
        model, filename = filename, format = "pdb")
      print("Writing boxed PDB with box unit cell to %s" %(
          "%s" %full_filename), file = log)

    # Write NCS file if NCS
    if ncs_object and ncs_object.max_operators()>0:
      if(params.output_file_name_prefix is None):
        filename =  "%s_box.ncs_spec"%output_prefix
      else:
        filename =  "%s.ncs_spec"%params.output_file_name_prefix
      dm.write_ncs_spec_file(ncs_object, filename = filename)
      print("\nWriting symmetry to %s" %( filename), file = log)
      # we are writing out new location

    # Write ccp4 map.
    if("ccp4" in params.output_format):
      if(params.output_file_name_prefix is None):
        filename = "%s_box.ccp4"%output_prefix
      else:
        filename = "%s.ccp4"%params.output_file_name_prefix
      dm.write_real_map_file(
         map_manager, filename = filename)
      print("\nWriting map to %s" %( filename), file = log)

    # Write ccp4 mask.
    if("ccp4" in params.output_format and mask_map_manager):
      if(params.output_file_name_prefix is None):
        filename = "%s_mask_box.ccp4"%output_prefix
      else:
        filename = "%s_mask.ccp4"%params.output_file_name_prefix
      dm.write_real_map_file(
         mask_map_manager, filename = filename)
      print("\nWriting mask to %s" %( filename), file = log)

    # Write xplor map.  Shift back if keep_origin = True
    if("xplor" in params.output_format):
     if(params.output_file_name_prefix is None):
       file_name = "%s_box.xplor"%output_prefix
     else: file_name = "%s.xplor"%params.output_file_name_prefix
     output_box.write_xplor_map(file_name = file_name,
         output_crystal_symmetry = mam.map_manager().crystal_symmetry(),
         output_unit_cell_grid = mam.map_manager().unit_cell_grid,
         shift_back = (output_box.shift_cart !=  (0, 0, 0)) )
     print("Writing boxed map "+\
         "to X-plor formatted file: %s"%file_name, file = log)

    # Write mtz map coeffs.  Shift back if keep_origin = True
    if("mtz" in params.output_format):
     if(params.output_file_name_prefix is None):
       file_name = "%s_box.mtz"%output_prefix
     else: file_name = "%s.mtz"%params.output_file_name_prefix

     print("Writing map coefficients "+\
         "to MTZ file: %s"%file_name, file = log)
     if(resolution_from_map_coeffs is not None):
       d_min = resolution_from_map_coeffs
     elif params.resolution is not None:
       d_min = params.resolution
     else:
       d_min = maptbx.d_min_from_map(map_data = mam.map_manager().map_data(),
         unit_cell = mam.map_manager().crystal_symmetry().unit_cell())
     map_coeffs = mam.map_manager().map_as_fourier_coefficients(
       d_min = d_min)
     mtz_dataset = map_coeffs.as_mtz_dataset(column_root_label = 'F')
     mtz_object = mtz_dataset.mtz_object()
     dm.write_miller_array_file(mtz_object, filename = file_name)

  print(file = log)
  return output_box

class box_object(object):
  '''
    Temporary holder that replaces box in map_box
  '''
  def __init__(self,
      shift_cart = None,
      origin_shift_grid_units = None,
      hierarchy = None,
      xray_structure_box = None,
      ncs_object = None,
      map_box = None,  # boxed map_data
      map_data = None,  # original map_data
      map_box_half_map_list = None,
      box_crystal_symmetry = None,
      pdb_outside_box_msg = "",
      solvent_content = None,
      gridding_first = None,
      gridding_last = None,
      ):
    adopt_init_args(self, locals())
    del self.origin_shift_grid_units
    self._origin_shift_grid_units = origin_shift_grid_units

  def show_summary(self, log = sys.stdout):
     print("Box object summary", file = log)
     print("Value of shift_cart: ", self.shift_cart, file = log)
     print("Value of origin_shift_grid_units: ", self.origin_shift_grid_units(),
       file = log)
     print("Value of map_box.origin(): ", self.map_box.origin(), file = log)
     print("Value of map_box.all(): ", self.map_box.all(), file = log)
     print("Value of map_data.origin(): ", self.map_data.origin(), file = log)
     print("Value of map_data.all(): ", self.map_data.all(), file = log)
     print("Value of solvent_content: ", self.solvent_content, file = log)
     print("Value of gridding_first: ", self.gridding_first, file = log)
     print("Value of gridding_last: ", self.gridding_last, file = log)

  def shift_sites_cart_back(self, sites_cart):
    from scitbx.matrix import col
    return sites_cart-col(self.shift_cart)

  def get_solvent_content(self): # XXX need to save it from
    return self.solvent_content

  def origin_shift_grid_units(self, reverse = True):
    if reverse:
      return tuple([-x for x in self._origin_shift_grid_units])
    else:
      return self._origin_shift_grid_units

  def shift_map_back(self, map_data):
    # Shift map from map_box cell to original coordinate system
    #  Note this map only applies in the region of the map_box cell (the
    #   map may be repeated in space but only one copy is valid).
    # The dimensions of this map are the same as the box map.
    from scitbx.matrix import col
    new_origin = self.origin_shift_grid_units(reverse = True)
    new_all = list(col(self.map_box.all())+col(new_origin))
    shifted_map_data = map_data.deep_copy()
    from scitbx.array_family import flex
    shifted_map_data.resize(flex.grid(new_origin, new_all))
    return shifted_map_data

  def write_xplor_map(self, file_name = "box.xplor", shift_back = None,
      output_unit_cell_grid = None,
      output_crystal_symmetry = None):

    # write out xplor map on same grid as ccp4 map (0 to focus-1)
    from scitbx.matrix import col
    if shift_back:
      map_data = self.shift_map_back(self.map_box)
    else:
      map_data = self.map_box

    if output_unit_cell_grid is None:
     output_unit_cell_grid = map_data.all()

    if output_crystal_symmetry is None:
      output_crystal_symmetry = self.xray_structure_box.crystal_symmetry()
    import iotbx.xplor
    gridding = iotbx.xplor.map.gridding(
        n     = output_unit_cell_grid,
        first = map_data.origin(),
        last  = tuple(col(map_data.focus())-col((1, 1, 1))))

    iotbx.xplor.map.writer(
      file_name          = file_name,
      is_p1_cell         = None, # XXX temporary flag allowing any cell
      title_lines        = ['Map in box', ],
      unit_cell          = output_crystal_symmetry.unit_cell(),
      gridding           = gridding,
      data               = map_data.as_double(),
      average            = -1,
      standard_deviation = -1)

def change_output_unit_cell(params = None,
   ccp4_map = None,
   model = None):

    '''
     Change the output unit cell as requested by user

     Check to see if user is also going to set the output_unit_cell_grid.
     If so, change that first.
    '''

    # Make sure there is no origin offset because that would change
    if ccp4_map.origin_shift_grid_units!= (0, 0, 0) or \
         ccp4_map.map_data().origin()!= (0, 0, 0):
       raise Sorry("Input map cannot have an origin "+
         "shift if output_unit_cell is to be changed")

    # Does user want to simultaneously change output_unit_cell_grid
    #  If so, then they want the output pixel size to be
    #      (output_unit_cell/output_unit_cell_grid)
    #    therefore change the unit_cell_grid first, then set crystal symmetry
    if params.output_unit_cell_grid and \
      tuple(params.output_unit_cell_grid)!= tuple(ccp4_map.unit_cell_grid):
      ccp4_map.set_original_origin_and_gridding(
        gridding = params.output_unit_cell_grid)

    from cctbx import crystal
    new_symmetry = crystal.symmetry(
      tuple(params.output_unit_cell),
      ccp4_map.crystal_symmetry().space_group_number())

    ccp4_map.set_unit_cell_crystal_symmetry(new_symmetry)

    # Now set model crystal_symmetry to match so we can combine them
    if model:
      model.set_unit_cell_crystal_symmetry(
        ccp4_map.unit_cell_crystal_symmetry())
      model.set_crystal_symmetry(ccp4_map.crystal_symmetry())
    return ccp4_map, model

def apply_around_unique(mam,
      params = None,
      sequence = None,
      target_ncs_au_model = None,
      log = None):

    from cctbx.maptbx.box import around_unique
    new_mam = around_unique(
      mam.map_manager(),
      model = mam.model(),
      target_ncs_au_model = target_ncs_au_model,
      sequence = sequence,
      regions_to_keep = params.regions_to_keep,
      solvent_content = params.solvent_content,
      resolution = params.resolution,
      molecular_mass = params.molecular_mass,
      symmetry = params.symmetry,
      chain_type = params.chain_type,
      keep_low_density = params.keep_low_density,
      box_cushion = params.box_cushion,
      soft_mask = params.soft_mask_extract_unique,
      mask_expand_ratio = params.mask_expand_ratio,
      )
    return new_mam  # XXX actually it is a box not an mam

def apply_mask_around_atoms(mam, params = None, log = None):
    assert mam.model() is not None
    if (params.soft_mask and
       params.increase_box_cushion_and_atom_radius_for_soft_mask):
      # add soft_mask_radius to atom radius
      print ("Mask radius around atoms increased by soft_mask_radius", file = log)
      mask_atoms_atom_radius = \
         params.mask_atoms_atom_radius+params.soft_mask_radius
    else: # use atom radius
      mask_atoms_atom_radius = params.mask_atoms_atom_radius
    print ("Applying mask around atoms with radius of %.1f A" %(
       mask_atoms_atom_radius), file = log)
    if params.set_outside_to_mean_inside:
      print("Value outside mask will be set to mean inside", file = log)
    if params.invert_mask:
      print("Mask will be inverted (zero inside region of atoms, one outside)",
        file=log)
    mam.map_manager().create_mask_around_atoms(model = mam.model(),
        mask_atoms_atom_radius = mask_atoms_atom_radius,
        invert_mask = params.invert_mask)
    if (params.soft_mask): # make it a soft mask
      mam.map_manager().soft_mask(soft_mask_radius = params.soft_mask_radius)
      print ("Mask will be soft with radius of %.1f A" %(
         params.soft_mask_radius), file = log)
    mam.map_manager().apply_mask(
      set_outside_to_mean_inside = params.set_outside_to_mean_inside)
    return mam


def apply_mask_around_edge_of_box(mam, params = None, log = None):
    print ("Applying soft mask around edge of box with radius of %.1f A" %(
       params.soft_mask_radius), file = log)
    if params.set_outside_to_mean_inside:
      print("Value outside mask will be set to mean inside", file = log)

    mam.map_manager().create_mask_around_edges(
          boundary_radius = params.soft_mask_radius)
    mam.map_manager().soft_mask(soft_mask_radius = params.soft_mask_radius)
    mam.map_manager().apply_mask(
      set_outside_to_mean_inside = params.set_outside_to_mean_inside)
    return mam

#  =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
from wxGUI2 import utils

def validate_params(params):
  if params.write_mask_file and not params.mask_atoms:
    raise Sorry("You need to set mask_atoms for write_mask_file")
  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    result = run(args = self.args, log = sys.stdout)
    return 0

#  =============================================================================

if (__name__  ==  "__main__"):
  run(args = sys.argv[1:])
