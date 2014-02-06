from __future__ import division
from cctbx.array_family import flex
import os
import mmtbx.model
import libtbx.load_env
from libtbx import easy_pickle
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cStringIO import StringIO
from mmtbx import utils
from libtbx.utils import format_cpu_times, null_out
from libtbx.test_utils import approx_equal
import math
from cctbx import multipolar


print dir(multipolar)

def should_have_imported_proper_methods():
  assert hasattr(multipolar, 'multipolar_test' ) == True, 'multipolar has no method multipolar_test'
  assert hasattr(multipolar, 'assign_atom_types' ) == True, 'multipolar has no method assign_atom_types'

def should_return_atom_types():
  atom_number = flex.int()
  coordinates = flex.vec3_double()
  for i in range (3):
    atom_number.append(i)
    coordinates.append( (i,i,i) )
  atom_types = multipolar.assign_atom_types( atom_number, coordinates )
  #print atom_types
  for atom_type in atom_types:
    print atom_type

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
                   relative_path="phenix_regression/pdb/enk.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)

  #assert hasattr(multipolar, 'multipolar_test' ) == True, 'multipolar has no method multipolar_test'


  #print dir(geometry)
  #print geometry.pair_proxies

  xrs = processed_pdb_file.xray_structure()

  for scatterer in xrs.scatterers():
    print scatterer.scattering_type


  bond_proxies_simple = geometry.pair_proxies(sites_cart =
    xrs.sites_cart()).bond_proxies.simple

  sites_cart = xrs.sites_cart()

  bond_deltas = flex.double()
  #loop over bonds
  for proxy in bond_proxies_simple:
    #print proxy.i_seqs, proxy.distance_ideal
    i,j = proxy.i_seqs
    site_1 = sites_cart[i]
    site_2 = sites_cart[j]
    dist_model = math.sqrt((site_1[0]-site_2[0])**2+(site_1[1]-site_2[1])**2+(site_1[2]-site_2[2])**2)
    bond_deltas.append(proxy.distance_ideal-dist_model)
  #print list(flex.abs(bond_deltas))

  STOP()
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  mol.xray_structure.scattering_type_registry(table = "wk1995")


def run():
  should_have_imported_proper_methods()
  should_return_atom_types()
  exercise()

if (__name__ == "__main__"):
  run()
