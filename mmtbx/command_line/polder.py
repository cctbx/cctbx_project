from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.polder

import time
import sys
from cctbx.array_family import flex
import mmtbx.f_model
from mmtbx import utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from cStringIO import StringIO
import iotbx.phil
from iotbx import crystal_symmetry_from_any
from libtbx.utils import Sorry
import os
from iotbx import pdb
from cctbx import xray
#from mmtbx import maps
from mmtbx import map_tools
import mmtbx.masks
#
from cctbx import maptbx
from iotbx import ccp4_map
import iotbx.pdb
import iotbx.mtz

#*******************
# TO DO
#
# - error if there are no Rfree in original cif (fetch creates column with 1 everywhere)
#
# - delete checkprint
#
# - user input dmin or let program decide?
#
# - what about sigma scaling
#
# - put back again radius from params
#
#*******************


master_params_str = """\
pdb_file_name = None
  .type = str
  .multiple = False
  .help = PDB file name (model only)
no_mask_selection = None
  .type = str
  .help = atoms around which bulk solvent mask is set to zero
ligand_code = None
  .type = str
  .help = 3letter code for the ligand to be replaced by dummy atoms
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
sphere_radius = 5
  .type = float
  .help = radius of sphere around ligand atoms where solvent mask is reset to zero
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

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

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

def get_f_mask(xray_structure, miller_array, resolution_factor, sites_cart, selection,
  params, r_free_flags):
  crystal_gridding = miller_array.crystal_gridding(
    d_min             = miller_array.d_min(),    # dmin or paramas dmin? --> decide
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = resolution_factor)
  mask_data = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xray_structure,
    p1                       = True,
    solvent_radius           = 1.11,
    shrink_truncation_radius = 0.9,
    for_structure_factors    = True,
    n_real                   = crystal_gridding.n_real()).mask_data
  #print 'mask_data is padded', mask_data.is_padded() # checkprint
  # The next line is important for mask = mask_data.as_1d(); otherwise failure
  maptbx.unpad_in_place(map=mask_data)
  print 'min max mean total mask', mask_data.as_1d().min_max_mean().as_tuple()
  #print 'mask_data is padded after command', mask_data.is_padded() # checkprint
  # save solvent mask for whole structure
  ccp4_map(file_name="mask-all.ccp4", uc=miller_array.unit_cell(),
        sg=miller_array.space_group(), map_data=mask_data) # checkprint
  #-----------------
  # test: calculate mask for xray_structure when ligand is ommited --> solvent goes into cavity
  xray_structure_noligand = xray_structure.select(selection, negate = True)
  mask_data_omit = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xray_structure_noligand,
    p1                       = True,
    solvent_radius           = 1.11,
    shrink_truncation_radius = 0.9,
    for_structure_factors    = True,
    n_real                   = crystal_gridding.n_real()).mask_data
  maptbx.unpad_in_place(map=mask_data_omit)
  print 'min max mean omit mask', mask_data_omit.as_1d().min_max_mean().as_tuple()
  ccp4_map(file_name="mask-omit.ccp4", uc=miller_array.unit_cell(),
        sg=miller_array.space_group(), map_data=mask_data_omit)
#-----------------
#  mc = miller_array.structure_factors_from_map(
#    map            = mask_data,
#    use_scale      = True,
#    anomalous_flag = False,
#    use_sg         = False)
#  mtz_dataset = mc.as_mtz_dataset(column_root_label="MASK")
#  mtz_object = mtz_dataset.mtz_object()
#  mtz_object.write(file_name = "mask.mtz")

  # select grid points around atoms from params.no_mask_selection
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = miller_array.crystal_symmetry().unit_cell(),
    fft_n_real = mask_data.focus(),
    fft_m_real = mask_data.all(),
    sites_cart = sites_cart,
    site_radii = flex.double(sites_cart.size(), 5.0))
  print 'size selection', sel.size()
  #print sel[22] # this one changes from 1 to zero, 3ksw
  print mask_data.size()
  mask = mask_data.as_1d()
  # count the grid points which change 1 --> 0
  print 'number of grid points total', mask.size()
  print 'min max mean before modification', mask.min_max_mean().as_tuple()
  number = 0
  for i in sel:
    number += mask[i]
  print 'number of points = 1:', number
 # Reset the grid values around ligand atoms to zero
  mask.set_selected(sel, 0)
  #print mask[sel[22]]
  print 'min max mean after modification', mask.min_max_mean().as_tuple()
  mask.reshape(mask_data.accessor())
  #a = mask_data.as_1d().min_max_mean().mean() * mask_data.size()
  #b = mask.as_1d().min_max_mean().mean() * mask_data.size()
  #print a, b, a-b
  ccp4_map(file_name="mask-hole.ccp4", uc=miller_array.unit_cell(),
        sg=miller_array.space_group(), map_data=mask) # checkprint
  return miller_array.structure_factors_from_map(map = mask,
    use_scale = True, anomalous_flag = False, use_sg = False)
  # is it really use_sg = false

# not necessary - for comparison with "usual" omit map
def test_map(f_obs, r_free_flags, xray_structure, pdb_hierarchy, params):
  print 'now calculating test map for comparison'
  crystal_gridding = f_obs.crystal_gridding(
    d_min             = f_obs.d_min(),    # dmin or paramas dmin? --> decide
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 0.25)
  mask_data = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xray_structure,
    p1                       = True,
    solvent_radius           = 1.11,
    shrink_truncation_radius = 0.9,
    for_structure_factors    = True,
    n_real                   = crystal_gridding.n_real()).mask_data
  maptbx.unpad_in_place(map = mask_data)
  print mask_data.size()
  print 'min max mean total mask', mask_data.as_1d().min_max_mean().as_tuple()
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure = xray_structure).f_calc()
  mask = f_obs.structure_factors_from_map(map = mask_data,
    use_scale = True, anomalous_flag = False, use_sg = False)
  fmodel = mmtbx.f_model.manager(
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    f_calc       = f_calc,
    f_mask       = mask)
  fmodel.update_all_scales()
  print "r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free())
  fmodel.show(show_header=False, show_approx=False)
  mc_diff2 = map_tools.electron_density_map(
    fmodel=fmodel).map_coefficients(
      map_type         = "mFo-DFc",
      isotropize       = True,
      fill_missing     = False)
  mtz_dataset2 = mc_diff2.as_mtz_dataset(column_root_label="mFo-DFc")
  mtz_object2 = mtz_dataset2.mtz_object()
  mtz_object2.write(file_name = "no-polder.mtz")

def polder(f_obs, r_free_flags, xray_structure, pdb_hierarchy, params):
  crystal_symmetry = xray_structure.crystal_symmetry() # maybe delete
  print 'now selecting atoms...' # checkprint
  selection_bool = pdb_hierarchy.atom_selection_cache().selection(string = params.no_mask_selection)
  print xray_structure.show_summary()
  print "Selection:", params.no_mask_selection #checkprint
  print "Atoms selected:", xray_structure.sites_cart().select(selection_bool).size()   # checkprint
  # when getting coordinates, xray structure needs to be in P1
  sites_cart_ligand = xray_structure.select(selection_bool).expand_to_p1(sites_mod_positive=True).sites_cart()
  xray_structure_noligand = xray_structure.select(selection_bool,negate = True)
  print xray_structure_noligand.show_summary()

  print 'now calculating solvent mask...' # checkprint
  solvent_mask = get_f_mask(
      xray_structure    = xray_structure,
      miller_array      = f_obs,
      resolution_factor = 0.25,
      sites_cart        = sites_cart_ligand,
      selection         = selection_bool,
      params            = params,
      r_free_flags      = r_free_flags)
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure = xray_structure_noligand).f_calc()
  fmodel = mmtbx.f_model.manager(
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    f_calc       = f_calc,
    f_mask       = solvent_mask)
  fmodel.update_all_scales()
  print "r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free())
  fmodel.show(show_header=False, show_approx=False)
  mc_diff = map_tools.electron_density_map(
    fmodel=fmodel).map_coefficients(
      map_type         = "mFo-DFc",
      isotropize       = True,
      fill_missing     = False)
  mtz_dataset = mc_diff.as_mtz_dataset(column_root_label="mFo-DFc")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "polder_map_coeffs.mtz")
  test_map(f_obs, r_free_flags, xray_structure_noligand, pdb_hierarchy, params)

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def cmd_run(args, command_name):
  msg = """\

Tool for improvement of ligand omit map.

How to use:
1: Run this command: phenix.polder
2: Copy, save into a file (e.g. parameters.txt) and edit the parameters
   shown between the lines *** below. Do not include *** lines.
3: Run the command with this parameters file:
   phenix.polder parameters.txt
"""
  if(len(args) == 0):
    print msg
    print "*"*79
    master_params().show()
    print "*"*79
    return
  else:
    if(not os.path.isfile(args[0]) or len(args)>1):
      print "Parameter file is expected at input. This is not a parameter file:\n", \
        args
      print "Run phenix.polder without argumets for running instructions."
      return
    processed_args = utils.process_command_line_args(args = args,
      master_params = master_params(), log = None)
    params = processed_args.params.extract()
    if(params.pdb_file_name is None):
      assert len(processed_args.pdb_file_names)==1
      params.pdb_file_name = processed_args.pdb_file_names[0]
    run(processed_args = processed_args, params = params)

def run(processed_args, params):
  print "phenix.polder is running..."
  if(params.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
    raise Sorry("Incorrect scattering_table.")
  crystal_symmetry = None
  crystal_symmetries = []
  for f in [str(params.pdb_file_name), str(params.reflection_file_name)]:
    cs = crystal_symmetry_from_any.extract_from(f)
    if(cs is not None): crystal_symmetries.append(cs)
  if(len(crystal_symmetries) == 1): crystal_symmetry = crystal_symmetries[0]
  elif(len(crystal_symmetries) == 0):
    raise Sorry("No crystal symmetry found.")
  else:
    if(not crystal_symmetries[0].is_similar_symmetry(crystal_symmetries[1])):
      raise Sorry("Crystal symmetry mismatch between different files.")
    crystal_symmetry = crystal_symmetries[0]
  if (params.no_mask_selection is None):                       # DL
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
      r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
      test_flag_value=None
    f_obs, r_free_flags = f_obs.common_sets(r_free_flags) #DL
  pdb_input = iotbx.pdb.input(file_name = params.pdb_file_name)
  pdb_hierarchy = pdb_input.construct_hierarchy()
  xray_structure = pdb_input.xray_structure_simple()
  print xray_structure.crystal_symmetry().show_summary() #checkprint
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

if(__name__ == "__main__"):
  t0 = time.time()
  cmd_run(
    args         = sys.argv[1:],
    command_name = "phenix.polder")
  print "Time:", time.time()-t0
