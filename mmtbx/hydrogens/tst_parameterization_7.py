from __future__ import division
from __future__ import print_function
from builtins import zip
import time
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
from cctbx import geometry_restraints
from mmtbx.hydrogens import riding
import iotbx.cif

#-----------------------------------------------------------------------------
# This test makes sure that H are corrected according to angle restraint
# NH2 group is planar, but there are no dihedral restraints
# Input model has wrong NH2 configuration --> idealize can correct it.
#-----------------------------------------------------------------------------

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()

  cif_object = iotbx.cif.reader(input_string = cif_str).model()
  for srv in [mon_lib_srv, ener_lib]:
    srv.process_cif_object(cif_object = cif_object)

  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = None,
    raw_records    = pdb_str,
    force_symmetry = True)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xray_structure = processed_pdb_file.xray_structure()

  geometry_restraints = processed_pdb_file.geometry_restraints_manager(
    show_energies = False)

  sites_cart = xray_structure.sites_cart()
  atoms = pdb_hierarchy.atoms()

  riding_h_manager = riding.manager(
    pdb_hierarchy       = pdb_hierarchy,
    geometry_restraints = geometry_restraints)
  h_para = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart = sites_cart,
    threshold  = 0.05)
  h_distances   = diagnostics.h_distances
  unk_list      = diagnostics.unk_list
  number_h_para = diagnostics.number_h_para
  type_list     = diagnostics.type_list

# number of H atoms in structure
  number_h = 0
  for h_bool in xray_structure.hd_selection():
    if h_bool: number_h += 1

  assert (number_h_para == number_h), 'Not all H atoms are parameterized'
  assert(len(unk_list) == 0), \
    'Some H atoms are parameterized with an unknown type'

  for ih in h_distances:
    # One atom is expected to be moved
    if (ih == 16):
      continue
    labels = atoms[ih].fetch_labels()
    assert (h_distances[ih] < 0.1), \
      'distance too large: %s  atom: %s (%s) residue: %s ' \
      % (h_para[ih].htype, atoms[ih].name, ih, labels.resseq.strip())

  for type1, type2 in zip(type_list, type_list_known):
    assert (type1 == type2)
    #print "'%s'," % type1,

# Several fragments from pdb with double conformations which caused crashes
# and which should now be mended
pdb_str = """\
REMARK iotbx.pdb.box_around_molecule --buffer-layer=5 "02.updated.nomin..pdb"
REMARK Date 2017-07-10 Time 12:11:01 PDT -0700 (1499713861.46 s)
CRYST1   15.938   15.073   11.550  90.00  90.00  90.00 P 1
SCALE1      0.062743  0.000000  0.000000        0.00000
SCALE2      0.000000  0.066344  0.000000        0.00000
SCALE3      0.000000  0.000000  0.086580        0.00000
HETATM    1  C1  3HA   401       7.842   9.143   5.222  1.00 83.17           C
HETATM    2  C2  3HA   401       7.240   7.911   5.459  1.00 81.93           C
HETATM    3  C3  3HA   401       8.027   6.770   5.589  1.00 82.00           C
HETATM    4  C4  3HA   401       9.411   6.864   5.486  1.00 82.52           C
HETATM    5  C5  3HA   401      10.012   8.097   5.252  1.00 83.82           C
HETATM    6  C6  3HA   401       9.226   9.238   5.118  1.00 84.24           C
HETATM    7  C7  3HA   401       5.719   7.834   5.619  1.00 81.09           C
HETATM    8  N10 3HA   401       7.463   5.583   5.793  1.00 82.31           N
HETATM    9  O11 3HA   401      10.183   5.752   5.619  1.00 82.09           O
HETATM   10  O8  3HA   401       5.197   6.868   6.174  1.00 81.05           O
HETATM   11  O9  3HA   401       5.000   8.734   5.190  1.00 78.80           O
HETATM   12  H1  3HA   401       7.315   9.910   5.134  1.00 83.17           H
HETATM   13  H11 3HA   401      10.719   5.854   6.276  1.00 82.09           H
HETATM   14  H5  3HA   401      10.938   8.151   5.134  1.00 83.82           H
HETATM   15  H6  3HA   401       9.630  10.073   5.000  1.00 84.24           H
HETATM   16 H101 3HA   401       7.572   5.175   6.550  1.00 82.31           H
HETATM   17 H102 3HA   401       7.435   5.000   5.152  1.00 82.31           H
TER
END
"""


cif_str = """
#
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
3HA 3HA "Unknown                  " ligand 17 11 .
#
data_comp_3HA
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
3HA        O8      O   OC    -1 .         -2.5304   -0.3544    2.0294
3HA        C7      C   C      0 .         -1.2760   -0.4611    2.0390
3HA        O9      O   O      0 .         -0.6571   -0.5146    3.1342
3HA        C2      C   CR6    0 .         -0.5023   -0.4939    0.7224
3HA        C1      C   CR16   0 .         -0.2734   -1.7018    0.0834
3HA        C6      C   CR16   0 .          0.4107   -1.7296   -1.1210
3HA        C5      C   CR16   0 .          0.8659   -0.5495   -1.6865
3HA        C4      C   CR6    0 .          0.6370    0.6583   -1.0475
3HA        O11     O   OH1    0 .          1.1017    1.8508   -1.6167
3HA        C3      C   CR6    0 .         -0.0471    0.6861    0.1569
3HA        N10     N   NH2    0 .         -0.3039    1.9536    0.8168
3HA        H1      H   HCR6   0 .         -0.6294   -2.6246    0.5256
3HA        H6      H   HCR6   0 .          0.5887   -2.6740   -1.6213
3HA        H5      H   HCR6   0 .          1.4009   -0.5713   -2.6283
3HA        H11     H   HOH1   0 .          2.0033    1.7449   -1.8770
3HA        H101    H   HNH2   0 .          0.0644    2.1253    1.7332
3HA        H102    H   HNH2   0 .         -0.8528    2.6556    0.3576
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
3HA  O8      C7     deloc         1.259 0.020     1.259
3HA  C7      O9     deloc         1.259 0.020     1.259
3HA  C7      C2     single        1.527 0.020     1.527
3HA  C2      C1     aromatic      1.386 0.020     1.386
3HA  C2      C3     aromatic      1.385 0.020     1.385
3HA  C1      C6     aromatic      1.385 0.020     1.385
3HA  C1      H1     single        0.930 0.020     1.080
3HA  C6      C5     aromatic      1.385 0.020     1.385
3HA  C6      H6     single        0.930 0.020     1.080
3HA  C5      C4     aromatic      1.385 0.020     1.385
3HA  C5      H5     single        0.930 0.020     1.080
3HA  C4      O11    single        1.401 0.020     1.401
3HA  C4      C3     aromatic      1.385 0.020     1.385
3HA  O11     H11    single        0.850 0.020     0.980
3HA  C3      N10    single        1.452 0.020     1.452
3HA  N10     H101   single        0.860 0.020     1.020
3HA  N10     H102   single        0.860 0.020     1.020
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
3HA  C2      C7      O9           120.00 3.000
3HA  C2      C7      O8           119.99 3.000
3HA  O9      C7      O8           120.00 3.000
3HA  C3      C2      C1           120.00 3.000
3HA  C3      C2      C7           120.00 3.000
3HA  C1      C2      C7           120.00 3.000
3HA  H1      C1      C6           120.00 3.000
3HA  H1      C1      C2           120.00 3.000
3HA  C6      C1      C2           120.00 3.000
3HA  H6      C6      C5           120.00 3.000
3HA  H6      C6      C1           120.00 3.000
3HA  C5      C6      C1           120.00 3.000
3HA  H5      C5      C4           120.00 3.000
3HA  H5      C5      C6           120.00 3.000
3HA  C4      C5      C6           120.00 3.000
3HA  C3      C4      O11          120.00 3.000
3HA  C3      C4      C5           120.00 3.000
3HA  O11     C4      C5           120.00 3.000
3HA  H11     O11     C4           109.47 3.000
3HA  N10     C3      C4           120.00 3.000
3HA  N10     C3      C2           120.00 3.000
3HA  C4      C3      C2           120.00 3.000
3HA  H102    N10     H101         120.00 3.000
3HA  H102    N10     C3           120.00 3.000
3HA  H101    N10     C3           120.00 3.000
#
loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
3HA CONST_01      C5      C6      C1      C2             0.00   0.0 0
3HA CONST_02      C5      C4      C3      C2             0.00   0.0 0
3HA CONST_03      C4      C3      C2      C1            -0.00   0.0 0
3HA CONST_04      C4      C5      C6      C1            -0.00   0.0 0
3HA CONST_05      C3      C2      C1      C6             0.00   0.0 0
3HA CONST_06      C3      C4      C5      C6             0.00   0.0 0
3HA CONST_07      C6      C1      C2      C7           179.02   0.0 0
3HA CONST_08      C4      C3      C2      C7          -179.02   0.0 0
3HA CONST_09      O11     C4      C3      C2          -179.76   0.0 0
3HA CONST_10      N10     C3      C2      C1           179.11   0.0 0
3HA CONST_11      O11     C4      C5      C6           179.76   0.0 0
3HA CONST_12      N10     C3      C4      C5          -179.11   0.0 0
3HA CONST_13      H6      C6      C1      C2          -179.93   0.0 0
3HA CONST_14      H5      C5      C6      C1           180.00   0.0 0
3HA CONST_15      H1      C1      C6      C5          -180.00   0.0 0
3HA CONST_16      H101    N10     C3      C2            61.62   0.0 0
3HA CONST_17      H102    N10     C3      C2          -118.12   0.0 0
3HA Var_01        C1      C2      C7      O8           -88.86  30.0 2
#
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
3HA plan-1  C7     0.020
3HA plan-1  C2     0.020
3HA plan-1  C1     0.020
3HA plan-1  C6     0.020
3HA plan-1  C5     0.020
3HA plan-1  C4     0.020
3HA plan-1  O11    0.020
3HA plan-1  C3     0.020
3HA plan-1  N10    0.020
3HA plan-1  H1     0.020
3HA plan-1  H6     0.020
3HA plan-1  H5     0.020
3HA plan-2  C3     0.020
3HA plan-2  N10    0.020
3HA plan-2  H101   0.020
3HA plan-2  H102   0.020
3HA plan-3  O8     0.020
3HA plan-3  C7     0.020
3HA plan-3  O9     0.020
3HA plan-3  C2     0.020
"""

type_list_known = ['flat_2neigbs', 'alg1b', 'flat_2neigbs', 'flat_2neigbs',
  'alg1a', 'alg1a']

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK. Time: %8.3f"%(time.time()-t0))
