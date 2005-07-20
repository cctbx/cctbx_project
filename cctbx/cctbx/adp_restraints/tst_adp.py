from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import sys, math
from libtbx.test_utils import approx_equal
from iotbx import pdb
from cctbx.array_family import flex
from mmtbx.monomer_library import pdb_interpretation
import libtbx.load_env
import os
import cctbx.adp_restraints.manager



def run():
  pdb_file_name = libtbx.env.find_in_repositories(
        relative_path="regression/pdb/adp_restraints.pdb", test=os.path.isfile)
  processed_pdb_file = pdb_interpretation.process(
                               mon_lib_srv = monomer_library.server.server(),
                               ener_lib    = monomer_library.server.ener_lib(),
                               file_name   = pdb_file_name)
  xray_structure = processed_pdb_file.xray_structure()
  grm = processed_pdb_file.geometry_restraints_manager(
                                                      plain_pairs_radius = 5.0)
  pair_proxies = grm.pair_proxies(sites_cart=xray_structure.sites_cart())
  manager = cctbx.adp_restraints.manager.iso(xray_structure,
                                             grm,
                                             sphere_radius = 1.6,
                                             distance_power = 0.0,
                                             wilson_b = None,
                                             mean_power = 0.0,
                                             normalize=False)

  u_iso_restraints = grm.harmonic_restraints(
                      variables    = xray_structure.extract_u_iso_or_u_equiv(),
                      type_indices = None,
                      type_weights = 1.0)
  assert approx_equal(u_iso_restraints.residual_sum, manager.target())
  assert approx_equal(u_iso_restraints.gradients, manager.gradients())

if (__name__ == "__main__"):
    run()
    print "iso_restraints: OK"
