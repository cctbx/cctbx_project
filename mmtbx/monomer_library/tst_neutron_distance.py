from __future__ import division
from cctbx.array_family import flex
import mmtbx.model
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import format_cpu_times 

pdb_str = """
CRYST1   16.242   15.015   14.726  90.00  90.00  90.00 P 1
ATOM    396  N   TYR A  28       6.525   8.629   9.726  1.00 10.00           N
ATOM    397  CA  TYR A  28       5.932   8.949   8.446  1.00 10.00           C
ATOM    398  CB  TYR A  28       7.030   9.299   7.419  1.00 10.00           C
ATOM    399  CG  TYR A  28       7.962   8.140   7.011  1.00 10.00           C
ATOM    400  CD1 TYR A  28       7.605   7.249   6.028  1.00 10.00           C
ATOM    401  CE1 TYR A  28       8.436   6.191   5.650  1.00 10.00           C
ATOM    402  CZ  TYR A  28       9.678   6.043   6.246  1.00 10.00           C
ATOM    403  OH  TYR A  28      10.519   5.000   5.846  1.00 10.00           O
ATOM    404  CE2 TYR A  28      10.067   6.942   7.247  1.00 10.00           C
ATOM    405  CD2 TYR A  28       9.212   7.966   7.631  1.00 10.00           C
ATOM    406  C   TYR A  28       5.000   7.875   7.893  1.00 10.00           C
ATOM    407  O   TYR A  28       5.007   6.705   8.336  1.00 10.00           O
ATOM      0  HA  TYR A  28       5.367   9.722   8.604  1.00 10.00           H
ATOM      0  HB2 TYR A  28       6.603   9.645   6.620  1.00 10.00           H
ATOM      0  HB3 TYR A  28       7.573  10.015   7.783  1.00 10.00           H
ATOM      0  HD1 TYR A  28       6.785   7.353   5.601  1.00 10.00           H
ATOM      0  HD2 TYR A  28       9.470   8.547   8.310  1.00 10.00           H
ATOM      0  HE1 TYR A  28       8.156   5.588   5.000  1.00 10.00           H
ATOM      0  HE2 TYR A  28      10.898   6.853   7.654  1.00 10.00           H
ATOM      0  HH  TYR A  28      11.242   5.032   6.292  1.00 10.00           H
TER
HETATM    1  O   HOH A   1      -0.354   0.000  -0.573  1.00 20.00      A    O  
HETATM    2  H1  HOH A   1       1.392  -0.000  -0.481  1.00 20.00      A    H  
HETATM    3  H2  HOH A   1      -1.038  -0.000   1.054  1.00 20.00      A    H  
END
"""

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    raw_records    = flex.std_string(pdb_str.splitlines()),
    use_neutron_distances=True,
    force_symmetry = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies = False, plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry = geometry,
    normalization = True)
  xray_structure = processed_pdb_file.xray_structure()
  ph = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  ph.write_pdb_file(file_name="input.pdb")
  mmodel = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure     = xray_structure,
    pdb_hierarchy      = ph)
  mmodel.geometry_minimization(
    bond      = True,
    nonbonded = True,
    angle     = True,
    dihedral  = True,
    chirality = True,
    planarity = True)
  ph.adopt_xray_structure(mmodel.xray_structure)
  ph.write_pdb_file(file_name="output.pdb")

def run():
  exercise()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
