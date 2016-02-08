from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.polder
import time
import sys
import os
from cStringIO import StringIO
import mmtbx.f_model
import mmtbx.utils
import mmtbx.masks
from mmtbx import utils
from mmtbx import map_tools
from iotbx import pdb
from iotbx import ccp4_map
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import iotbx.phil
import iotbx.pdb
import iotbx.mtz
from cctbx import xray
from cctbx import maptbx
from cctbx.array_family import flex
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry

master_phil_string = """
model_file_name = None
  .type = str
  .multiple = False
  .help = PDB file name
solvent_exclusion_mask_selection = None
  .type = str
  .help = atoms around which bulk solvent mask is set to zero
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
sphere_radius = 5
  .type = float
  .help = radius of sphere around atoms, where solvent mask is reset to zero
high_resolution = None
  .type = float
low_resolution = None
  .type = float
output_file_name_prefix = None
  .type = str
scattering_table = *n_gaussian wk1995 it1992 neutron electron
  .type = choice
  .help = Scattering table for structure factors calculations
resolution_factor = 0.25
  .type = float
debug = True
  .type = bool
  .expert_level=3
"""

def ccp4_map(file_name, uc, sg, map_data):
  from iotbx import ccp4_map
  ccp4_map.write_ccp4_map(
    file_name   = file_name,
    unit_cell   = uc,
    space_group = sg,
    map_data    = map_data,
    labels      = flex.std_string([""]))

def output_map(f_obs,r_free_flags, xray_structure, mask_data, filename):
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure = xray_structure).f_calc()
  mask = f_obs.structure_factors_from_map(map = mask_data,
    use_scale = True, anomalous_flag = False, use_sg = False) # is it really use_sg = false?
  fmodel = mmtbx.f_model.manager(
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    f_calc       = f_calc,
    f_mask       = mask)
  fmodel.update_all_scales()
  print "r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free())
  print "*"*79
  mc_diff = map_tools.electron_density_map(
    fmodel=fmodel).map_coefficients(
      map_type         = "mFo-DFc",
      isotropize       = True,
      fill_missing     = False)
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

def polder(f_obs, r_free_flags, xray_structure, pdb_hierarchy, params):
  print "selecting atoms..."
  selection_bool = pdb_hierarchy.atom_selection_cache().selection(
    string = params.solvent_exclusion_mask_selection)
  print "Selection:", params.solvent_exclusion_mask_selection
  n_selected = selection_bool.count(True)
  if(n_selected == 0):
    raise Sorry("No atoms where selected. Check selection syntax again.")
  print "Atoms selected:", n_selected
  ligand_str = pdb_hierarchy.select(selection_bool).as_pdb_string()
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
  mask_data_all = mask_from_xrs_unpadded(xray_structure = xray_structure,
    n_real = crystal_gridding.n_real())
  print "R factors for unmodified input model and data:"
  mc_diff_all = output_map(f_obs, r_free_flags, xray_structure,
    mask_data_all, filename = "all")
  if(params.debug):
    print "R factor when ligand is used for mask calculation:"
    mc_diff_lig_omit = output_map(f_obs, r_free_flags, xray_structure_noligand,
      mask_data_all, filename = "lig_omit")
  #------
  mask_polder = mask_modif(f_obs, mask_data_all, sites_cart_ligand,
    sphere_radius = params.sphere_radius)
  print "R factor for polder map"
  mc_diff_polder = output_map(f_obs,r_free_flags, xray_structure_noligand,
    mask_polder, filename = "polder")
  # calculate mask for structure without ligand
  mask_data_omit = mask_from_xrs_unpadded(xray_structure=xray_structure_noligand,
    n_real = crystal_gridding.n_real())
  print "R factor when ligand is excluded for mask calculation:"
  mc_diff_omit = output_map(f_obs,r_free_flags, xray_structure_noligand,
    mask_data_omit, filename = "omit")
  mtz_dataset=mc_diff_polder.as_mtz_dataset(column_root_label="mFo-DFc_polder")
  if(params.debug):
    mtz_dataset.add_miller_array(
      miller_array = mc_diff_lig_omit,
      column_root_label = "mFo-DFc_lig-omit")
  mtz_dataset.add_miller_array(
    miller_array = mc_diff_omit,
    column_root_label = "mFo-DFc_omit")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "polder_map_coeffs.mtz")
  print "Finished."

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def run(params):
  print "phenix.polder is running..."
  if(params.scattering_table not in ["n_gaussian","wk1995","it1992","neutron",
    "electron"]):
    raise Sorry("Incorrect scattering_table.")
  crystal_symmetry = None
  crystal_symmetries = []
  for f in [str(params.model_file_name), str(params.reflection_file_name)]:
    cs = crystal_symmetry_from_any.extract_from(f)
    if(cs is not None): crystal_symmetries.append(cs)
  if(len(crystal_symmetries) == 1): crystal_symmetry = crystal_symmetries[0]
  elif(len(crystal_symmetries) == 0):
    raise Sorry("No crystal symmetry found.")
  else:
    if(not crystal_symmetries[0].is_similar_symmetry(crystal_symmetries[1])):
      raise Sorry("Crystal symmetry mismatch between different files.")
    crystal_symmetry = crystal_symmetries[0]
  if(params.solvent_exclusion_mask_selection is None):
    raise Sorry("No selection for mask calculation found.")
  f_obs, r_free_flags = None, None
  if(params.reflection_file_name is not None):
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name = params.reflection_file_name, ensure_read_access = True)
    rfs = reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      reflection_files = [reflection_file])
    parameters = utils.data_and_flags_master_params().extract()
    if(params.data_labels is not None):
      parameters.labels = [processed_args.data_labels]
    determine_data_and_flags_result = utils.determine_data_and_flags(
      reflection_file_server = rfs,
      parameters             = parameters,
      keep_going             = True,
      log                    = StringIO())
    f_obs = determine_data_and_flags_result.f_obs
    r_free_flags = determine_data_and_flags_result.r_free_flags
  assert [f_obs, r_free_flags].count(None)<=1
  pdb_input = iotbx.pdb.input(file_name = params.model_file_name)
  pdb_hierarchy = pdb_input.construct_hierarchy()
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

# parses through command line arguments
def cmd_run(args, command_name):
  msg = "Tool for improvement of ligand omit map."
  print msg
  master_params = master_phil_string
  master_phil = phil.parse(master_phil_string, process_includes=True)
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="polder")
  pdb = []
  mtz = []
  phils = []
  phil_args = []
  for arg in args:
    if(os.path.isfile(arg)):
      if(iotbx.pdb.is_pdb_file(arg)):
        pdb.append(arg)
      elif(arg.lower().endswith(".mtz")):
        mtz.append(arg)
      else:
        try:
          file_phil = phil.parse(file_name=arg)
        except RuntimeError:
          pass
        else:
          phils.append(file_phil)
    else:
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  working_params = working_phil.extract()
  if(working_params.model_file_name is None):
    if(len(pdb) == 1):
      working_params.model_file_name = pdb[0]
    else:
      raise Sorry("Exactly one model file should be given.")
  if(working_params.reflection_file_name is None):
    if(len(mtz) == 1):
      working_params.reflection_file_name = mtz[0]
    else:
      raise Sorry("Exactly one mtz file should be given.")
  run(working_params)

if(__name__ == "__main__"):
  t0 = time.time()
  cmd_run(
    args         = sys.argv[1:],
    command_name = "phenix.polder")
  print "Time:", time.time()-t0
