from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.polder
import time
import sys
from cStringIO import StringIO
import mmtbx.f_model
import mmtbx.utils
import mmtbx.masks
from mmtbx import utils # No need to import it twice. Previous import is better.
from mmtbx import map_tools
from iotbx import pdb
from iotbx import ccp4_map
from iotbx import file_reader
from iotbx import reflection_file_utils
import iotbx.phil
import iotbx.pdb # Second attempt to import it. This way is better.
import iotbx.mtz
from cctbx import xray
from cctbx import maptbx
from cctbx.array_family import flex
from libtbx import phil # Second import. Previous is better.
from libtbx.utils import Sorry

# in general, imports as mmtbx.utils are preferable, because in is obvious
# in the code what is called. There are surprising number of modules called
# "utils" in cctbx. On the opposite, from libtbx.utils import Sorry is
# ok, because Sorry is rather unique, and it would be excessively long
# call otherwise.

legend = """\

phenix.polder:
Computes omit maps by excluding the bulk solvent in the area around a
selection. One example of application are ligand omit maps. Polder omit maps
can be helpful if the ligand density is weak and obscured by bulk solvent
in conventional omit maps (where the ligand is deleted from the model).

Inputs:
  - File with reflection data (Fobs or Iobs) and R-free flags. It can
    be in most of known formats and spread across multiple files;
  - label(s) selecting which reflection data arrays should be used (in case
    there are multiple choices in input file; otherwise there is no need to
    provide labels);
  - Model file;
  - Atom selection (such as ligand)

Usage examples:
  1. phenix.polder model.cif data.mtz selection="chain A and resseq 1"
  2. phenix.polder model.pdb data.hkl data_labels="FP" selection="chain A"
  3. phenix.polder a.hkl b.hkl model.pdb selection="resseq 435"

Output:
  MTZ file with map coefficients:

  - mFo-DFc_polder    : polder difference map coefficients
  - PHImFo-DFc_polder : corresponding phases
  - mFo-DFc_omit      : omit difference map coefficients
  - PHImFo-DFc_omit   : corresponding phases

Optional output:
  CCP4 files with mask data:

  - mask_all.ccp4    : mask of original model
  - mask_omit.ccp4   : mask when ligand is omitted
  - mask_polder.ccp4 : mask obtained by polder procedure
"""

master_params_str = """
model_file_name = None
  .type = str
  .multiple = False
  .help = Model file name
solvent_exclusion_mask_selection = None
  .type = str
  .help = Atoms around which bulk solvent mask is set to zero
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
r_free_flags_labels = None
  .type = str
  .help = Labels for free reflections.
sphere_radius = 5
  .type = float
  .help = Radius of sphere around atoms where solvent mask is reset to zero
high_resolution = None
  .type = float
  .help = High resolution limit
low_resolution = None
  .type = float
  .help = Low resolution limit
scattering_table = *n_gaussian wk1995 it1992 neutron electron
  .type = choice
  .help = Scattering table for structure factors calculations
resolution_factor = 0.25
  .type = float
  .help = Used to determine the grid step = resolution_factor * high resolution
output_file_name_prefix = None
  .type = str
  .help = Prefix for output filename
mask_output = False
  .type = bool
  .help = Additional output: ccp4 maps containing the solvent mask for inital \
   model (mask_all.ccp4), when ligand is omitted (mask_omit.ccp4) and the mask \
   used for polder (mask_polder.ccp4).
debug = False
  .type = bool
  .expert_level=3
  .help = Additional output: biased omit map (ligand used for mask calculation but omitted from model)
"""

# it is not a good idea to call function exactly the same name as
# imported module. I'm not sure the presence of the function is justified
# because it is literally a single call. On top of that,
# the name is too general and
# does not provide any clue about what the function is doing.
def ccp4_map(file_name, uc, sg, map_data):
  # This is duplicated import, find it in the beginning of the file
  from iotbx import ccp4_map
  ccp4_map.write_ccp4_map(
    file_name   = file_name,
    unit_cell   = uc,
    space_group = sg,
    map_data    = map_data,
    labels      = flex.std_string([""]))

def output_map(f_obs,r_free_flags, xray_structure, mask_data, filename, params):
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure = xray_structure).f_calc()
  mask = f_obs.structure_factors_from_map(map = mask_data,
    # strings should be 80 chars long if there is no strong reason to breake
    # this rule. Comment could be here.
    use_scale = True, anomalous_flag = False, use_sg = False) # is it really use_sg = false?
  fmodel = mmtbx.f_model.manager(
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    f_calc       = f_calc,
    f_mask       = mask)
  fmodel.update_all_scales()
  # Generally, it is good idea to output everything somewhere, not straight
  # to stdout.
  # Especially when it is a funciton. If it needs to output anything,
  # the "log" parameter could be supplied. If stdout output is desired,
  # sys.stdout could be passed as value of "log" parameter.
  # Then it would be
  # print >> log, "r_work=" ... etc.
  print "r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free())
  print "*"*79
  mc_diff = map_tools.electron_density_map(
    fmodel=fmodel).map_coefficients(
      map_type         = "mFo-DFc",
      isotropize       = True,
      fill_missing     = False)
  if (params.mask_output and filename != 'bias_omit' ):
    ccp4_map(file_name="mask_"+filename+".ccp4", uc=f_obs.unit_cell(),
      sg=f_obs.space_group(), map_data=mask_data)
  return mc_diff

def mask_modif(f_obs, mask_data, sites_cart, sphere_radius):
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = f_obs.crystal_symmetry().unit_cell(),
    fft_n_real = mask_data.focus(),
    fft_m_real = mask_data.all(),
    sites_cart = sites_cart,
    site_radii = flex.double(sites_cart.size(), sphere_radius))
  mask = mask_data.as_1d()
  mask.set_selected(sel, 0)
  mask.reshape(mask_data.accessor())
  return mask

def mask_from_xrs_unpadded(xray_structure, n_real):
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xray_structure,
    p1                       = True,
    shrink_truncation_radius = mask_params.shrink_truncation_radius,
    solvent_radius           = mask_params.solvent_radius,
    for_structure_factors    = True,
    n_real                   = n_real).mask_data
  maptbx.unpad_in_place(map = mask)
  return mask

# I think it would be useful to make it as convenient as possible to
# call this function from anywhere else. If it needs to return more than one
# thing, it could be rewritted as a class. In this case (of using this
# functionality from somewhere else), remark about printing to log
# becomes even more relevant.
def polder(f_obs, r_free_flags, xray_structure, pdb_hierarchy, params):
  # see above the remark about print statements.
  print "selecting atoms..."
  selection_bool = pdb_hierarchy.atom_selection_cache().selection(
    string = params.solvent_exclusion_mask_selection)
  # it would be great to ouput parameter before you actually do selection.
  # if it is bad, the user would look at traceback at this point
  # and won't have a chance to see the parameter :)
  print "Selection:", params.solvent_exclusion_mask_selection
  n_selected = selection_bool.count(True)
  if(n_selected == 0):
    raise Sorry("No atoms where selected. Check selection syntax again.")
  print "Number of atoms selected:", n_selected
  ligand_str = pdb_hierarchy.select(selection_bool).as_pdb_string()
  # I would strongly recommend against using "naked" prints at all,
  # i.e. print <variable name>.
  # It is extremely difficult to find them later when they need to be removed.
  # Put something meaningful in front of them, like
  # print "atoms selected\n", ligand_str
  print ligand_str
  # when extracting cartesian coordinates, xray_structure needs to be in P1!
  sites_cart_ligand = xray_structure.select(selection_bool).expand_to_p1(
    sites_mod_positive=True).sites_cart()
  # xray structure object without ligand/selection
  xray_structure_noligand = xray_structure.select(~selection_bool)
  print "Calculating solvent mask..."
  print "*"*79
  crystal_gridding = f_obs.crystal_gridding(
    d_min             = f_obs.d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = params.resolution_factor)
  # calculate mask using ALL atoms
  # Style inconsistency: here we have spaces around "=" in parameter list,
  # 10 lines above we don't ("sites_mod_positive=True")
  mask_data_all = mask_from_xrs_unpadded(xray_structure = xray_structure,
    n_real = crystal_gridding.n_real())
  print "R factors for unmodified input model and data:"
  mc_diff_all = output_map(f_obs, r_free_flags, xray_structure,
    mask_data_all, filename = "all", params = params)
  # biased map - for developers
  if(params.debug):
    print "R factor when ligand is used for mask calculation (biased map):"
    mc_diff_bias_omit = output_map(f_obs, r_free_flags, xray_structure_noligand,
      mask_data_all, filename = "bias_omit", params = params)
  #------
  mask_polder = mask_modif(f_obs, mask_data_all, sites_cart_ligand,
    sphere_radius = params.sphere_radius)
  print "R factor for polder map"
  mc_diff_polder = output_map(f_obs, r_free_flags, xray_structure_noligand,
    mask_polder, filename = "polder", params = params)
  # calculate mask for structure without ligand
  # Example of style inconsistency in one function call. And skipping spaces
  # wasn't enough to keep the line shorter than 80 symbols. And anyway,
  # the function call is two lines :)
  mask_data_omit = mask_from_xrs_unpadded(xray_structure=xray_structure_noligand,
    n_real = crystal_gridding.n_real())
  print "R factor when ligand is excluded for mask calculation:"
  mc_diff_omit = output_map(f_obs, r_free_flags, xray_structure_noligand,
    mask_data_omit, filename = "omit", params = params)
  mtz_dataset = mc_diff_polder.as_mtz_dataset(column_root_label="mFo-DFc_polder")
  if(params.debug):
    mtz_dataset.add_miller_array(
      miller_array = mc_diff_bias_omit,
      column_root_label = "mFo-DFc_bias_omit")
  mtz_dataset.add_miller_array(
    miller_array = mc_diff_omit,
    column_root_label = "mFo-DFc_omit")
  mtz_object = mtz_dataset.mtz_object()
  polder_file_name = "polder_map_coeffs.mtz"
  if (params.output_file_name_prefix is not None):
    polder_file_name = params.output_file_name_prefix + "_" + polder_file_name
  mtz_object.write(file_name = polder_file_name)
  print "Finished."

# parse through command line arguments
# command_name parameter is not used. I don't see any reason to have it.
def cmd_run(args, command_name):
  if(len(args)==0):
    print legend
    return
  print "phenix.polder is running..."
  # There is no reason to make a new variable and not use it. Nothing is wrong
  # with original master_params_str
  master_params = master_params_str
  # The use of process_includes here is a bit misleadin. You don't have any
  # includes in your phil.
  parsed = phil.parse(master_params_str, process_includes=True)
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = parsed)
  #inputs.params.show() #check
  params = inputs.params.extract()
  if(params.scattering_table not in ["n_gaussian","wk1995","it1992","neutron",
    "electron"]):
    # I bet you won't get here. I haven't check it, but most likely,
    # if user provides wrong
    # table here it will fail during phil parsing.
    raise Sorry("Incorrect scattering_table.")
  if(params.solvent_exclusion_mask_selection is None):
    raise Sorry("No selection for mask calculation found.")
  # check model file
  if(len(inputs.pdb_file_names) == 0):
    if(params.model_file_name is None):
      raise Sorry("No model file found.")
  elif(len(inputs.pdb_file_names) == 1):
    params.model_file_name = inputs.pdb_file_names[0]
  else:
    raise Sorry("Only one model file should be given")
  # check reflection file
  reflection_files = inputs.reflection_files
  if(len(reflection_files) == 0):
    if(params.reflection_file_name is None):
      raise Sorry("No reflection file found.")
    else:
      hkl_in = file_reader.any_file(params.reflection_file_name,
        force_type="hkl")
      hkl_in.assert_file_type("hkl")
      reflection_files = [ hkl_in.file_object ]
  # crystal symmetry
  # What is the reason to assing crystal_symmetry to None first?
  crystal_symmetry = None
  crystal_symmetry = inputs.crystal_symmetry
  if(crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  f_obs, r_free_flags = None, None
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry   = True,
    reflection_files = reflection_files,
    err              = StringIO())
  parameters = mmtbx.utils.data_and_flags_master_params().extract()
  if(params.data_labels is not None):
    parameters.labels = params.data_labels
  if(params.r_free_flags_labels is not None):
    parameters.r_free_flags.label = params.r_free_flags_labels
  # This name looks extremely like a function name due to a verb in it.
  # In reality it is an object. Consider, e.g. another name:
  # "determined_data_and_flags"
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server = rfs,
    parameters             = parameters,
    keep_going             = True,
    log                    = StringIO())
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  # This assertion is wrong. It passes if f_obs is None, and r_free_flags is
  # not None. I don't know if this situation is possible, but in this case
  # the code will fail on the next line in printing f_obs labels.
  assert [f_obs, r_free_flags].count(None)<=1
  print "Input data:"
  print "  Iobs or Fobs:", f_obs.info().labels
  if(r_free_flags is not None):
    print "  Free-R flags:", r_free_flags.info().labels
  else:
    print "  Free-R flags: Not present"
  if (params.sphere_radius < 3 or params.sphere_radius > 10):
    raise Sorry("Sphere radius out of range: must be between 3 A and 10 A")
  pdb_input = iotbx.pdb.input(file_name = params.model_file_name)
  pdb_hierarchy = pdb_input.construct_hierarchy()
  # This way of extracting xray_structure is actually wrong. It will be
  # inconsistent with pdb_hierarchy, because pdb_hierarchy gets sorted and
  # pdb_input.xray_structure_simple is not. A way to test it would be to
  # put a side chain atoms in front of main chain in a residue and supply
  # selection that would select side-chain atoms. It would actually make
  # a good test ;) I'm sure atoms selected in xray_structure
  # (in polder function) will be wrong.
  # Use pdb_hierarchy.extract_xray_structure() instead.
  xray_structure = pdb_input.xray_structure_simple()
  mmtbx.utils.setup_scattering_dictionaries(
    scattering_table = params.scattering_table,
    xray_structure   = xray_structure,
    d_min            = f_obs.d_min())
  if(f_obs is not None):
    f_obs = f_obs.resolution_filter(d_min = params.high_resolution,
      d_max = params.low_resolution)
  if(r_free_flags is not None):
    r_free_flags = r_free_flags.resolution_filter(d_min=params.high_resolution,
      d_max = params.low_resolution)
  polder(f_obs          = f_obs,
         r_free_flags   = r_free_flags,
         xray_structure = xray_structure,
         pdb_hierarchy  = pdb_hierarchy,
         params         = params)

if(__name__ == "__main__"):
  t0 = time.time()
  cmd_run(
    args         = sys.argv[1:],
    command_name = "phenix.polder")
  print "Time:", round(time.time()-t0, 2)
