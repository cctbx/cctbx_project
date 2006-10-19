from iotbx import crystal_symmetry_from_any
from iotbx.pdb import crystal_symmetry_from_pdb
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.monomer_library import server
import iotbx.pdb
from iotbx.option_parser import option_parser
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.utils import Sorry
from cctbx.xray import ext
from cctbx.xray.structure import structure as cctbx_xray_structure
from stdlib import math
from iotbx import pdb
from libtbx import adopt_init_args
from cctbx import geometry_restraints
from libtbx.test_utils import approx_equal
import sys, os, math
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx import geometry_restraints

class dbe(object):

   def __init__(self, geometry_restraints_manager,
                      atom_attributes,
                      xray_structure,
                      dbe_parameters = None):
     adopt_init_args(self, locals())
     self.sites_cart = self.xray_structure.sites_cart()
     if (self.xray_structure.scatterers().count_anisotropic() > 0):
       raise RuntimeError(
         "Handling of anisotropic scatterers not implemented.")
     self.b_iso = self.xray_structure.scatterers().extract_u_iso()*math.pi**2*8.
     self.q = self.xray_structure.scatterers().extract_occupancies()
     self.pair_proxies = self.geometry_restraints_manager.pair_proxies(
                                                  sites_cart = self.sites_cart)
     self.dbe_as_xray_structure = self.xray_structure.deep_copy_scatterers()
     self.dbe_as_xray_structure.erase_scatterers()

   def bond_table(self):
     line_counter = 0
     for proxy in self.pair_proxies.bond_proxies.simple:
       line_counter += 1
       i_seqs = proxy.i_seqs
       atom_i = self.atom_attributes[i_seqs[0]]
       atom_j = self.atom_attributes[i_seqs[1]]
       b_iso_i = self.b_iso[i_seqs[0]]
       b_iso_j = self.b_iso[i_seqs[1]]
       q_i = self.q[i_seqs[0]]
       q_j = self.q[i_seqs[1]]
       print dir(atom_i)
       sys.exit(0)
       bond = geometry_restraints.bond(sites_cart = self.sites_cart,
                                       proxy      = proxy)
       print "%6d:"%line_counter,                \
             "%6d "%i_seqs[0],                   \
              atom_i.name,                       \
              atom_i.element,                    \
              atom_i.resName,                    \
              atom_i.resSeq,                     \
              "%6.2f "%b_iso_i,                  \
              "%4.2f "%q_i,                      \
              "<< %5.3f (%5.3f) >>"%(bond.distance_model, bond.distance_ideal), \
              "%6d "%i_seqs[1],                  \
              atom_j.name,                       \
              atom_j.element,                    \
              atom_j.resName,                    \
              atom_j.resSeq,                     \
              "%6.2f "%b_iso_j,                  \
              "%4.2f "%q_j,                      \
              atom_i.line_number



def run():
    log = sys.stdout
    file_name = "enkmolad_pdb_iso.cns"
    crystal_symmetry = crystal_symmetry_from_pdb.extract_from(
                       file_name = file_name)

    mon_lib_srv = server.server()
    ener_lib = server.ener_lib()

    processed_pdb_file = pdb_interpretation.process(
           mon_lib_srv                           = mon_lib_srv,
           ener_lib                              = ener_lib,
           file_name                             = file_name,
           strict_conflict_handling              = False,
           crystal_symmetry                      = crystal_symmetry,
           force_symmetry                        = True,
           log = log)

    all_chain_proxies = monomer_library.pdb_interpretation.build_all_chain_proxies(
          mon_lib_srv = mon_lib_srv,
          ener_lib    = ener_lib,
          file_name   = file_name,
          log         = log)
    sites_cart = processed_pdb_file.xray_structure().sites_cart()

    atoms = all_chain_proxies.stage_1.atom_attributes_list

    geometry_restraints_manager = \
             processed_pdb_file.geometry_restraints_manager()

    pair_proxies = geometry_restraints_manager.pair_proxies(sites_cart = sites_cart)

    dbe_manager = dbe(geometry_restraints_manager = geometry_restraints_manager,
                      atom_attributes             = atoms,
                      xray_structure              = processed_pdb_file.xray_structure(),
                      dbe_parameters              = None)
    dbe_manager.bond_table()


    if(0):
     for proxy in pair_proxies.bond_proxies.asu:
       # NOT EXERCISED
       i_seqs = (proxy.i_seq, proxy.j_seq)
       atom_i, atom_j = atoms[i_seqs[0]], atoms[i_seqs[1]]
       atom_i.coordinates, atom_i.Ucart
       bond = geometry_restraints.bond(
         sites_cart=sites_cart,
         asu_mappings = pair_proxies.bond_proxies.asu_mappings(),
         proxy=proxy)
       bond.distance_model
       bond.distance_ideal
       print "%6d:"%line_counter,                \
             "%6d "%i_seqs[0],                   \
              atom_i.name,                       \
              atom_i.element,                    \
              atom_i.resName,                    \
              atom_i.resSeq,                     \
              "%6.2f "%atom_i.tempFactor,        \
              "%4.2f "%atom_i.occupancy,         \
              "<< %5.3f (%5.3f) >>"%(bond.distance_model, bond.distance_ideal), \
              "%6d "%i_seqs[1],                  \
              atom_j.name,                       \
              atom_j.element,                    \
              atom_j.resName,                    \
              atom_j.resSeq,                     \
              "%6.2f "%atom_j.tempFactor,        \
              "%4.2f "%atom_j.occupancy,atom_i.line_number

if (__name__ == "__main__"):
  run()
