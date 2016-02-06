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
scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
"""

#def master_params():
 # return iotbx.phil.parse(master_params_str, process_includes=False)

def ccp4_map(file_name, uc, sg, map_data):
  from iotbx import ccp4_map
#  map_data.apply_sigma_scaling()
#  target_map_data = map_data.real_map_unpadded()
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
  #return mc_diff
  ccp4_map(file_name="mask_"+filename+".ccp4", uc=f_obs.unit_cell(),
    sg=f_obs.space_group(), map_data=mask_data)
  return mc_diff
  #mtz_object = mtz_dataset.mtz_object()
  #mtz_object.write(file_name = "polder_map_coeffs.mtz")

# sets grid points around sites_cart to zero
def mask_modif(f_obs, mask_data, sites_cart, params):
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = f_obs.crystal_symmetry().unit_cell(),
    fft_n_real = mask_data.focus(),
    fft_m_real = mask_data.all(),
    sites_cart = sites_cart,
    site_radii = flex.double(sites_cart.size(), params.sphere_radius))
  mask = mask_data.as_1d()
  # count the grid points which change 1 --> 0
  #print 'number of grid points total', mask.size()
  #print 'min max mean before modification', mask.min_max_mean().as_tuple()
  number = 0
  for i in sel:
    number += mask[i]
  #print 'number of points = 1:', number
 # Reset the grid values around ligand atoms to zero
  mask.set_selected(sel, 0)
  #print 'min max mean after modification', mask.min_max_mean().as_tuple()
  mask.reshape(mask_data.accessor())
  return mask

def polder(f_obs, r_free_flags, xray_structure, pdb_hierarchy, params):
  print 'now selecting atoms...'
  #no_H_select = params.solvent_exclusion_mask_selection + ' and not element H'
  selection_bool = pdb_hierarchy.atom_selection_cache().selection(string = params.solvent_exclusion_mask_selection)
  #print no_H_select
  #selection_bool = pdb_hierarchy.atom_selection_cache().selection(string = no_H_select)
  print "Selection:", params.solvent_exclusion_mask_selection #checkprint
  selection_length = xray_structure.sites_cart().select(selection_bool).size()
  if(selection_length == 0):
    raise Sorry("No atoms where selected. Check selection syntax again.")
  print "Atoms selected:", selection_length
  ligand_str = pdb_hierarchy.select(selection_bool).as_pdb_string()
  print ligand_str
  #STOP()
  # when extracting cartesian coordinates, xray_structure needs to be in P1!
  sites_cart_ligand = xray_structure.select(selection_bool).expand_to_p1(sites_mod_positive=True).sites_cart()
  # xray structure object without ligand/selection
  xray_structure_noligand = xray_structure.select(selection_bool,negate = True)
  print 'now calculating solvent mask...'
  print "*"*79
  resolution_factor = 0.25
  crystal_gridding = f_obs.crystal_gridding(
    d_min             = f_obs.d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = resolution_factor)
  # calculate mask using ALL atoms
  mask_data_all = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xray_structure,
    p1                       = True,
    solvent_radius           = 1.11,
    shrink_truncation_radius = 0.9,
    for_structure_factors    = True,
    n_real                   = crystal_gridding.n_real()).mask_data
  maptbx.unpad_in_place(map = mask_data_all)
  print "R factors for unmodified input model and data:"
  mc_diff_all = output_map(f_obs,r_free_flags, xray_structure,
    mask_data_all, filename = "all")
  # This map is not necessary for final version, but useful for tests
  print "R factor when ligand is used for mask calculation:" #del
  mc_diff_lig_omit = output_map(f_obs,r_free_flags, xray_structure_noligand, #del
    mask_data_all, filename = "lig_omit") #del
  #------
  mask_polder = mask_modif(f_obs, mask_data_all, sites_cart_ligand, params)
  print "R factor for polder map"
  mc_diff_polder = output_map(f_obs,r_free_flags, xray_structure_noligand,
    mask_polder, filename = "polder")
  # calculate mask for structure without ligand
  mask_data_omit = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xray_structure_noligand,
    p1                       = True,
    solvent_radius           = 1.11,
    shrink_truncation_radius = 0.9,
    for_structure_factors    = True,
    n_real                   = crystal_gridding.n_real()).mask_data
  maptbx.unpad_in_place(map = mask_data_omit)
  print "R factor when ligand is excluded for mask calculation:"
  mc_diff_omit = output_map(f_obs,r_free_flags, xray_structure_noligand,
    mask_data_omit, filename = "omit")
  mtz_dataset = mc_diff_polder.as_mtz_dataset(column_root_label="mFo-DFc_polder")
  mtz_dataset.add_miller_array(
    miller_array = mc_diff_lig_omit,
    column_root_label = "mFo-DFc_lig-omit")
  mtz_dataset.add_miller_array(
    miller_array = mc_diff_omit,
    column_root_label = "mFo-DFc_omit")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "polder_map_coeffs.mtz")
  print "Finished"

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def run(params):
  print "phenix.polder is running..."
  if(params.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
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
  if (params.solvent_exclusion_mask_selection is None):            # DL
        raise Sorry("No selection for mask calculation found.")    # DL
  f_obs = None
  r_free_flags = None
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
      reflection_file_server  = rfs,
      parameters              = parameters,
      keep_going              = True,
      log                     = StringIO())
    f_obs = determine_data_and_flags_result.f_obs
    r_free_flags = determine_data_and_flags_result.r_free_flags
    # add something if no rfree
    if(r_free_flags is None):
        raise Sorry("No Rfree flags found in input file.")
      #r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
      #test_flag_value=None
    f_obs, r_free_flags = f_obs.common_sets(r_free_flags) #DL
  pdb_input = iotbx.pdb.input(file_name = params.model_file_name)
  pdb_hierarchy = pdb_input.construct_hierarchy()
  xray_structure = pdb_input.xray_structure_simple()
  if(f_obs is not None):
    f_obs = f_obs.resolution_filter(d_min = params.high_resolution,
      d_max = params.low_resolution)
    r_free_flags = r_free_flags.resolution_filter(d_min = params.high_resolution,
      d_max = params.low_resolution)
  polder(f_obs             = f_obs,
         r_free_flags      = r_free_flags,
         xray_structure    = xray_structure,
         pdb_hierarchy     = pdb_hierarchy,
         params            = params)

# parses through command line arguments
def cmd_run(args, command_name):
  msg = "Tool for improvement of ligand omit map."
  print msg
  #if len(args) != 3:
  #     raise Sorry("Sorry, wrong number of parameters.")
  master_params = master_phil_string
  master_phil = phil.parse(master_phil_string, process_includes=True)
  #print args #checkprint
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="polder")
  pdb = []
  mtz = []
  phils = []
  phil_args = []
  for arg in args:
    if os.path.isfile(arg) :
      if iotbx.pdb.is_pdb_file(arg):
        pdb.append(arg)
      elif arg.lower().endswith(".mtz"):
        mtz.append(arg)
      else:
        try :
          file_phil = phil.parse(file_name=arg)
        except RuntimeError :
          pass
        else :
          phils.append(file_phil)
    else :
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  #print phil_args #checkprint
  #print pdb, mtz #checkprint
  working_phil = master_phil.fetch(sources=phils)
  #working_phil.show()
  working_params = working_phil.extract()
  #print working_phil.format(python_object=working_params).as_str() #checkprint
  if working_params.model_file_name is None:
    if len(pdb) == 1:
      working_params.model_file_name = pdb[0]
    else:
          raise Sorry("Exactly one model file should be given.")
  if working_params.reflection_file_name is None:
    if len(mtz) == 1:
          working_params.reflection_file_name = mtz[0]
    else:
          raise Sorry("Exactly one mtz file should be given.")
  #print working_phil.format(python_object=working_params).as_str() #checkprint
  run(working_params)

if(__name__ == "__main__"):
  t0 = time.time()
  cmd_run(
    args         = sys.argv[1:],
    command_name = "phenix.polder")
  print "Time:", time.time()-t0
