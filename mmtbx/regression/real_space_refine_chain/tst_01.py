from __future__ import absolute_import, division, print_function
import time, os, random
import iotbx.pdb
import mmtbx.utils
from mmtbx import monomer_library
from scitbx.array_family import flex
import mmtbx.refinement.real_space.explode_and_refine
from mmtbx.geometry_restraints import reference
from iotbx import reflection_file_reader
import libtbx.load_env

#if (1):
#  random.seed(0)
#  flex.set_random_seed(0)

def ccp4_map(crystal_symmetry, file_name, map_data):
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
      file_name=file_name,
      unit_cell=crystal_symmetry.unit_cell(),
      space_group=crystal_symmetry.space_group(),
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))

def run(prefix="tst_00"):
  # Poor model that we want to refine so it matches the answer
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/real_space_refine_chain/poor_model.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/real_space_refine_chain/poor_map.mtz",
    test=os.path.isfile)
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  ph_poor = pdb_inp.construct_hierarchy()
  ph_poor.atoms().reset_i_seq()
  xrs_poor = pdb_inp.xray_structure_simple()
  # Initialize states accumulator
  states = mmtbx.utils.states(pdb_hierarchy=ph_poor)
  states.add(sites_cart = xrs_poor.sites_cart())
  # Compute target map
  mas = reflection_file_reader.any_reflection_file(file_name =
    mtz_file).as_miller_arrays()
  assert len(mas)==1
  fc = mas[0]

  fft_map = fc.fft_map(resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  target_map_data = fft_map.real_map_unpadded()
  ccp4_map(crystal_symmetry=fc.crystal_symmetry(), file_name="map.ccp4",
    map_data=target_map_data)
  # Build geometry restraints
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.nonbonded_weight=200
  #params.peptide_link.ramachandran_restraints=True
  #params.peptide_link.rama_potential="oldfield"
  #print dir(params)
  #STOP()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = monomer_library.server.server(),
    ener_lib                 = monomer_library.server.ener_lib(),
    file_name              = pdb_file,
    params                   = params,
    #crystal_symmetry = fc.crystal_symmetry(),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)

  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    plain_pairs_radius           = 5,
    assume_hydrogens_all_missing = True)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)

  #for a in ph_answer.atoms():
  #  print a.i_seq, a.name, a.xyz
    #STOP()

  #ref_xyz = flex.vec3_double([(14.323, 35.055, 14.635), (16.099, 12.317, 16.37)])
  #selection = flex.size_t([1,76])
  #
  #restraints_manager.geometry.adopt_reference_coordinate_restraints_in_place(
  #  reference.add_coordinate_restraints(
  #    sites_cart = ref_xyz,
  #    selection = selection,
  #    sigma = 0.1))

  # Do real-space refinement
  t0=time.time()
  ear = mmtbx.refinement.real_space.explode_and_refine.run(
    xray_structure          = xrs_poor,
    pdb_hierarchy           = ph_poor,
    map_data                = target_map_data,
    restraints_manager      = restraints_manager,
    states                  = states)
  print("Time: %6.4f"%(time.time()-t0))
  ear.pdb_hierarchy.write_pdb_file(file_name="%s_refined.pdb"%prefix)
  states.write(file_name="%s_refined_all_states.pdb"%prefix)

if (__name__ == "__main__"):
  rs = 3292014
  random.seed(rs)
  flex.set_random_seed(rs)
  run()
  print("OK")
