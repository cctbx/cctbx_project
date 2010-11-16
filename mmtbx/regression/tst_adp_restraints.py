from mmtbx import monomer_library
import mmtbx.monomer_library.server
from mmtbx.monomer_library import pdb_interpretation
import cctbx.adp_restraints
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
import os
from cctbx.array_family import flex

def run():
  pdb_file_name = libtbx.env.find_in_repositories(
        relative_path="phenix_regression/pdb/adp_restraints.pdb", test=os.path.isfile)
  processed_pdb_file = pdb_interpretation.process(
                               mon_lib_srv = monomer_library.server.server(),
                               ener_lib    = monomer_library.server.ener_lib(),
                               file_name   = pdb_file_name)
  xray_structure = processed_pdb_file.xray_structure()
  grm = processed_pdb_file.geometry_restraints_manager(
                                                      plain_pairs_radius = 5.0)
  grm.pair_proxies(sites_cart=xray_structure.sites_cart())
  class parameters:
    sphere_radius = 1.6
    distance_power = 0.0
    average_power = 0.0
    min_u_sum = 1.e-6
  sel = flex.bool(xray_structure.scatterers().size(), True)
  xray_structure.scatterers().flags_set_grad_u_iso(sel.iselection())
  for use_hd in [True, False]:
    energies_adp = cctbx.adp_restraints.energies_iso(
      geometry_restraints_manager=grm,
      xray_structure=xray_structure,
      use_hd = use_hd,
      use_u_local_only = False,
      parameters=parameters)
    u_iso_restraints = grm.harmonic_restraints(
                      variables    = xray_structure.extract_u_iso_or_u_equiv(),
                      type_indices = None,
                      type_weights = 1.0)
    assert approx_equal(u_iso_restraints.residual_sum, energies_adp.residual_sum)
    assert approx_equal(u_iso_restraints.gradients, energies_adp.gradients)

if (__name__ == "__main__"):
  run()
  print format_cpu_times()
