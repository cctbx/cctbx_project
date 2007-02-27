from cctbx.array_family import flex
import math, time, sys, os
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
import mmtbx.model
from libtbx import introspection
import libtbx.load_env
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cStringIO import StringIO

def exercise():
  # initial setup
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  cif_file = libtbx.env.find_in_repositories(
            relative_path="phenix_regression/pdb/tyr.cif", test=os.path.isfile)
  mon_lib_srv.process_cif(file_name= cif_file)
  ener_lib.process_cif(file_name= cif_file)
  pdb_file = libtbx.env.find_in_repositories(
            relative_path="phenix_regression/pdb/ygg.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  xray_structure = processed_pdb_file.xray_structure()
  aal = processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = True,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  bond_proxies_simple = geometry.pair_proxies(
                  sites_cart = xray_structure.sites_cart()).bond_proxies.simple

  mol = mmtbx.model.manager(restraints_manager     = restraints_manager,
                            restraints_manager_ini = restraints_manager,
                            xray_structure         = xray_structure,
                            atom_attributes_list   = aal)
  # get dbe manager
  mol.setup_dbe_manager(fmodel = None, show = False, dbe_params = None)
  mol.use_dbe_true()
  assert mol.dbe_manager.need_restraints
  assert mol.dbe_selection.size() == 86
  assert mol.dbe_selection.count(True) == 45
  assert mol.dbe_selection.count(False) == 41
  mol.xray_structure.shake_sites_in_place(mean_distance = 0.5)
  mol.dbe_manager.target_and_gradients(
                                  sites_cart = mol.xray_structure.sites_cart(),
                                  dbe_selection = mol.dbe_selection)
  # do fd test
  scat_types = mol.xray_structure.scatterers().extract_scattering_types()
  restraints_selection = mol.dbe_manager.restraints_selection
  delta=0.00001
  unit_cell = mol.xray_structure.unit_cell()
  derivatives = flex.vec3_double()
  counter = 0
  for i_seq, sel in enumerate(mol.dbe_selection):
    if(sel):
       d_target_d_site = [0,0,0]
       for ix in xrange(3):
         target_values = []
         for d_sign in (-1, 1):
           modified_structure = mol.xray_structure.deep_copy_scatterers()
           ms = modified_structure.scatterers()[i_seq]
           site = list(ms.site)
           site_cart = list(unit_cell.orthogonalize(site))
           site_cart[ix] += d_sign * delta
           site = unit_cell.fractionalize(site_cart)
           ms.site = site
           mol.dbe_manager.target_and_gradients(
                                  sites_cart = modified_structure.sites_cart(),
                                  dbe_selection = mol.dbe_selection)
           target_values.append(mol.dbe_manager.target)
         derivative = (target_values[1] - target_values[0]) / (2 * delta)
         d_target_d_site[ix] = derivative
       g1 = ["%8.4f"%i for i in d_target_d_site]
       g2 = ["%8.4f"%i for i in mol.dbe_manager.gradients[i_seq]]
       assert approx_equal(d_target_d_site, mol.dbe_manager.gradients[i_seq], 1.e-5)
       counter += 1

def run():
  exercise()

if (__name__ == "__main__"):
  run()
