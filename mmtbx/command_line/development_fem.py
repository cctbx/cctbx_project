from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.development.fem

import sys, time
import mmtbx.utils
import iotbx.phil
import mmtbx.f_model
import iotbx.pdb
from iotbx import reflection_file_utils
from cStringIO import StringIO
import mmtbx.maps
from scitbx.array_family import flex
from cctbx import adptbx
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry

master_params_str="""\
f_obs_label = None
  .type = str
r_free_flags_label = None
  .type = str
scattering_table = wk1995  it1992  *n_gaussian  neutron
  .type = choice
"""

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def run(args, command_name = "phenix.development.fem", log = sys.stdout):
  if(len(args)==0):
    print >> log, \
"""
Usage:
  phenix.development.fem model.pdb data.mtz
  phenix.development.fem data.mtz model.pdb f_obs_label=F
"""
  parsed = iotbx.phil.parse(master_params_str)
  processed_args = mmtbx.utils.process_command_line_args(args = args,
    log = log, master_params = parsed)
  params = processed_args.params.extract()
  #
  reflection_files = processed_args.reflection_files
  if(len(reflection_files) == 0):
    raise Sorry("No reflection file found.")
  crystal_symmetry = processed_args.crystal_symmetry
  if(crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  if(len(processed_args.pdb_file_names) == 0):
    raise Sorry("No PDB file found.")
  pdb_file_names = processed_args.pdb_file_names
  #
  rfs = reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    reflection_files = reflection_files)
  parameters = mmtbx.utils.data_and_flags_master_params().extract()
  if(params.f_obs_label is not None):
    parameters.labels = params.f_obs_label
  if(params.r_free_flags_label is not None):
    parameters.r_free_flags.label = params.r_free_flags_label
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server  = rfs,
    parameters              = parameters,
    keep_going              = True,
    log                     = StringIO())
  f_obs = determine_data_and_flags_result.f_obs
  number_of_reflections = f_obs.indices().size()
  #
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=None
  # convert to non-anomalous
  f_obs = f_obs.average_bijvoet_mates()
  r_free_flags = r_free_flags.average_bijvoet_mates()
  f_obs, r_free_flags = f_obs.common_sets(r_free_flags)
  #
  pdb_file_names = processed_args.pdb_file_names
  assert len(pdb_file_names)==1
  xray_structure = iotbx.pdb.input(
    file_name=pdb_file_names[0]).xray_structure_simple()
  if(params.scattering_table == "neutron"):
    print >> log, "Using scattering table: neutron"
    xray_structure.switch_to_neutron_scattering_dictionary()
  #
  fmodel = mmtbx.f_model.manager(
    f_obs = f_obs,
    r_free_flags = r_free_flags,
    xray_structure = xray_structure)
  fmodel.update_all_scales(update_f_part1_for=None)
  fmodel.show(show_approx=False)
  ### BEGIN: compute THE SIMPLEST possible 2mFo-DFc (the most original one)
  mc_orig = fmodel.electron_density_map(
    update_f_part1 = False).map_coefficients(
      map_type     = "2mFo-DFc",
      isotropize   = True,
      fill_missing = False)
  ### END: compute THE SIMPLEST possible 2mFo-DFc (the most original one)
  print >> log, "r_work: %6.4f r_free: %6.4f"%(fmodel.r_work(), fmodel.r_free())
  ### b-factor sharpen
  xrs = fmodel.xray_structure
  b_iso_min = flex.min(xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
  print "Max B subtracted from atoms and used to sharpen map:", b_iso_min
  xrs.shift_us(b_shift=-b_iso_min)
  b_iso_min = flex.min(xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1))
  assert approx_equal(b_iso_min, 0, 1.e-3)
  fmodel.update_xray_structure(
    xray_structure = xrs,
    update_f_calc = True)
  ###
  fmodel.update_all_scales(update_f_part1_for="map", map_neg_cutoff=0)
  cs=fmodel.f_obs().crystal_symmetry()
  crystal_gridding = fmodel.f_obs().crystal_gridding(
    d_min                   = fmodel.f_obs().d_min(),
    resolution_factor       = 0.25,
    grid_step               = None,
    symmetry_flags          = None,
    mandatory_factors       = None,
    max_prime               = 5,
    assert_shannon_sampling = True)
  #
  ko = mmtbx.maps.kick(
    fmodel           = fmodel.deep_copy(),
    crystal_gridding = crystal_gridding)
  mc_orig_for_fem = fmodel.electron_density_map(
    update_f_part1 = True).map_coefficients(
      map_type     = "2mFo-DFc",
      isotropize   = True,
      fill_missing = True)
  #### Compute FEM start
  fem = mmtbx.maps.fem(ko=ko, crystal_gridding=crystal_gridding,
    mc_orig=mc_orig_for_fem, fmodel=fmodel)
  #### Compute FEM end
  mtz_dataset = mc_orig.as_mtz_dataset(column_root_label="2mFoDFc")
  mtz_dataset.add_miller_array(
    miller_array=mc_orig_for_fem,
    column_root_label="2mFoDFc_FilSharp")
  mtz_dataset.add_miller_array(
    miller_array=ko.map_coefficients,
    column_root_label="KICK")
  mtz_dataset.add_miller_array(
    miller_array=fem,
    column_root_label="FEM")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "fem.mtz")


if(__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print "Time: %6.4f"%(time.time()-t0)
