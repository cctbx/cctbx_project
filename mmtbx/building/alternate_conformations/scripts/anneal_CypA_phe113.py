
# XXX for development purposes only
#
# anneal Phe113 from 3k0n into partial omit map - if the protocol is working
# properly, the sidechain should usually pop into the secondary conformation.

from __future__ import division
from __future__ import print_function
from mmtbx.building import disorder
from mmtbx import building
import libtbx.load_env
import os

def exercise():
  import mmtbx.utils
  from iotbx.file_reader import any_file
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/CypA_refine_3.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/3k0n.mtz",
    test=os.path.isfile)
  pdb_in = any_file(pdb_file)
  hierarchy = pdb_in.file_object.construct_hierarchy()
  xrs = pdb_in.file_object.xray_structure_simple()
  mtz_in = any_file(mtz_file)
  f_obs = mtz_in.file_server.miller_arrays[0].f_sq_as_f()
  flags = mtz_in.file_server.miller_arrays[-1]
  flags = flags.customized_copy(data=(flags.data()==0))
  assert ((flags.data().count(True) / len(flags.data())) <= 0.1)
  f_obs, flags = f_obs.common_sets(other=flags)
  fmodel = mmtbx.utils.fmodel_simple(
    f_obs=f_obs,
    r_free_flags=flags,
    scattering_table="n_gaussian",
    xray_structures=[xrs])
  fmodel.info().show_targets()
  selection = hierarchy.atom_selection_cache().selection(
    "chain A and resseq 111:115")
  two_fofc_map, fofc_map, = disorder.get_partial_omit_map(
    fmodel=fmodel,
    selection=selection,
    map_file_name="omit.mtz",
    partial_occupancy=0.7)
  xrs.shake_sites_in_place(0.1, selection=selection)
  sites_new = building.run_real_space_annealing(
    xray_structure=xrs,
    pdb_hierarchy=hierarchy,
    selection=selection,
    target_map=fofc_map,
    d_min=f_obs.d_min(),
    target_map_rsr=two_fofc_map,
    debug=True)
  hierarchy.atoms().set_xyz(sites_new)
  hierarchy.write_pdb_file("anneal.pdb", crystal_symmetry=xrs)
  print("wrote anneal.pdb and omit.mtz")

if (__name__ == "__main__"):
  exercise()
