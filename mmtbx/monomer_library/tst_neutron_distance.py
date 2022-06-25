from __future__ import absolute_import, division, print_function
import mmtbx.model
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import format_cpu_times
import iotbx.pdb
import iotbx.phil
from libtbx.test_utils import approx_equal
import math
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from mmtbx.refinement import geometry_minimization

pdb_str = """
CRYST1   15.360   12.804   11.094  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A  28       4.991   8.496   8.094  1.00 10.00           N
ATOM      2  CA  TYR A  28       4.269   8.897   6.892  1.00 10.00           C
ATOM      3  CB  TYR A  28       5.247   9.154   5.742  1.00 10.00           C
ATOM      4  CG  TYR A  28       6.111   7.964   5.384  1.00 10.00           C
ATOM      5  CD1 TYR A  28       5.722   7.071   4.393  1.00 10.00           C
ATOM      6  CE1 TYR A  28       6.509   5.983   4.062  1.00 10.00           C
ATOM      7  CZ  TYR A  28       7.702   5.779   4.723  1.00 10.00           C
ATOM      8  OH  TYR A  28       8.488   4.698   4.396  1.00 10.00           O
ATOM      9  CE2 TYR A  28       8.111   6.652   5.709  1.00 10.00           C
ATOM     10  CD2 TYR A  28       7.318   7.736   6.033  1.00 10.00           C
ATOM     11  C   TYR A  28       3.254   7.833   6.490  1.00 10.00           C
ATOM     12  O   TYR A  28       3.415   6.655   6.811  1.00 10.00           O
ATOM     13  HA  TYR A  28       3.000   9.699   6.300  1.00 10.00           H
ATOM     14  HB2 TYR A  28       5.686   9.708   5.036  1.00 10.00           H
ATOM     15  HB3 TYR A  28       6.560   9.804   5.984  1.00 10.00           H
ATOM     16  HD1 TYR A  28       4.141   7.504   4.239  1.00 10.00           H
ATOM     17  HD2 TYR A  28       8.391   8.181   6.844  1.00 10.00           H
ATOM     18  HE1 TYR A  28       6.770   5.906   3.690  1.00 10.00           H
ATOM     19  HE2 TYR A  28       9.769   7.010   5.242  1.00 10.00           H
ATOM     20  HH ATYR A  28       8.000   3.505   4.139  0.50 10.00           H
ATOM     21  DH BTYR A  28       7.469   4.025   3.773  0.50 10.00           D
TER
HETATM   25  O   HOH C   1      11.274   5.928   4.656  1.00 10.00           O
HETATM   26  D1  HOH C   1      12.360   4.699   5.619  1.00 10.00           D
HETATM   27  D2  HOH C   1      10.422   5.411   6.216  1.00 10.00           D
TER
END
"""

def dist(a,b):
  return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def exercise(tolerance=0.01):
  for use_neutron_distances in [True, False]:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    inp = iotbx.pdb.input(lines=pdb_str.splitlines(), source_info=None)
    params = iotbx.phil.parse(
          input_string=grand_master_phil_str, process_includes=True).extract()
    params.pdb_interpretation.use_neutron_distances = use_neutron_distances
    model = mmtbx.model.manager(model_input = inp)
    model.process(pdb_interpretation_params=params, make_restraints=True)
    ph = model.get_hierarchy()
    ph.write_pdb_file(file_name="input.pdb")

    m = geometry_minimization.run2(
      restraints_manager = model.get_restraints_manager(),
      pdb_hierarchy = ph,
      correct_special_position_tolerance=1.0,
      bond      = True,
      nonbonded = True,
      angle     = True,
      dihedral  = True,
      chirality = True,
      planarity = True)
    model.set_sites_cart_from_hierarchy()
    suffix = "X"
    if(use_neutron_distances): suffix = "N"
    ph.write_pdb_file(file_name="output_%s.pdb"%suffix)
  # check X-H distances: x-ray
  cntr = 0
  awl = iotbx.pdb.input(file_name="output_X.pdb").atoms_with_labels()
  for a1 in awl:
    n1 = a1.name.strip()
    for a2 in awl:
      n2 = a2.name.strip()
      if([n1,n2]==["CD1","HD1"] or [n1,n2]==["CD2","HD2"] or
         [n1,n2]==["CE1","HE1"] or [n1,n2]==["CE2","HE2"]):
        assert approx_equal(dist(a1.xyz, a2.xyz), 0.93, tolerance)
        cntr += 1
      if(n1=="OH" and n2 in ["HH","DH"]):
        assert approx_equal(dist(a1.xyz, a2.xyz), 0.84, tolerance)
        cntr += 1
      if([n1,n2]==["CB","HB2"] or [n1,n2]==["CB","HB3"] or
         [n1,n2]==["CA","HA"]):
        assert approx_equal(dist(a1.xyz, a2.xyz), 0.97, tolerance)
        cntr += 1
      if(n1=="O" and n2 in ["D1","D2"] and a1.resname=="HOH"):
        assert approx_equal(dist(a1.xyz, a2.xyz), 0.85, tolerance)
        cntr += 1
  assert cntr == 11, cntr
  # check X-H distances: neutron
  cntr = 0
  awl = iotbx.pdb.input(file_name="output_N.pdb").atoms_with_labels()
  for a1 in awl:
    n1 = a1.name.strip()
    for a2 in awl:
      n2 = a2.name.strip()
      if([n1,n2]==["CD1","HD1"] or [n1,n2]==["CD2","HD2"] or
         [n1,n2]==["CE1","HE1"] or [n1,n2]==["CE2","HE2"]):
        assert approx_equal(dist(a1.xyz, a2.xyz), 1.08, tolerance), """
        dist(a1.xyz, a2.xyz) : %f
        tolerance            : %f
        """ % (dist(a1.xyz, a2.xyz),
               tolerance,
          )
        cntr += 1
      if(n1=="OH" and n2 in ["HH","DH"]):
        assert approx_equal(dist(a1.xyz, a2.xyz), 0.98, tolerance)
        cntr += 1
      if([n1,n2]==["CB","HB2"] or [n1,n2]==["CB","HB3"] or
         [n1,n2]==["CA","HA"]):
        assert approx_equal(dist(a1.xyz, a2.xyz), 1.09, tolerance)
        cntr += 1
      if(n1=="O" and n2 in ["D1","D2"] and a1.resname=="HOH"):
        assert approx_equal(dist(a1.xyz, a2.xyz), 0.98, tolerance)
        cntr += 1
  assert cntr == 11, cntr

def run():
  exercise()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
