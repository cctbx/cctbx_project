from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file

pdb_str = """\
HETATM    1  C01 LIG A   1      -1.100   0.000  -0.436  1.00 20.00      A    C
HETATM    2  C02 LIG A   1       0.428   0.000  -0.436  1.00 20.00      A    C
HETATM    3  C03 LIG A   1       1.007   0.000   0.977  1.00 20.00      A    C
HETATM    4  O04 LIG A   1       1.245  -1.090   1.560  1.00 20.00      A    O
HETATM    5  O07 LIG A   1       1.245   1.090   1.560  1.00 20.00      A    O-1
HETATM    6 H011 LIG A   1      -1.461  -0.885   0.075  1.00 20.00      A    H
HETATM    7 H012 LIG A   1      -1.461   0.885   0.075  1.00 20.00      A    H
HETATM    8 H013 LIG A   1      -1.461   0.000  -1.457  1.00 20.00      A    H
HETATM    9 H021 LIG A   1       0.778   0.881  -0.959  1.00 20.00      A    H
HETATM   10 H022 LIG A   1       0.778  -0.881  -0.959  1.00 20.00      A    H
"""

restr_string = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LIG        LIG 'Unknown                  ' ligand 10 5 .
#
data_comp_LIG
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
LIG         C01    C   CH3    0    .      -1.0996    0.0000   -0.4358
LIG         C02    C   CH2    0    .       0.4275    0.0000   -0.4358
LIG         C03    C   C      0    .       1.0067    0.0000    0.9772
LIG         O04    O   O      0    .       1.2453   -1.0903    1.5597
LIG         O07    O   OC    -1    .       1.2454    1.0903    1.5597
LIG        H011    H   HCH3   0    .      -1.4607   -0.8845    0.0748
LIG        H012    H   HCH3   0    .      -1.4607    0.8845    0.0749
LIG        H013    H   HCH3   0    .      -1.4607    0.0000   -1.4572
LIG        H021    H   HCH2   0    .       0.7784    0.8815   -0.9588
LIG        H022    H   HCH2   0    .       0.7784   -0.8815   -0.9588
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
LIG   C02     C01   single        1.527 0.020     1.527
LIG   C03     C02   single        1.527 0.020     1.527
LIG   O04     C03   deloc         1.259 0.020     1.259
LIG   O07     C03   deloc         1.259 0.020     1.259
LIG  H011     C01   single        0.970 0.020     1.090
LIG  H012     C01   single        0.970 0.020     1.090
LIG  H013     C01   single        0.970 0.020     1.090
LIG  H021     C02   single        0.970 0.020     1.090
LIG  H022     C02   single        0.970 0.020     1.090
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
LIG  H013     C01    H012         109.47 3.000
LIG  H013     C01    H011         109.47 3.000
LIG  H012     C01    H011         109.47 3.000
LIG  H013     C01     C02         109.47 3.000
LIG  H012     C01     C02         109.47 3.000
LIG  H011     C01     C02         109.47 3.000
LIG  H022     C02    H021         108.92 3.000
LIG  H022     C02     C03         108.90 3.000
LIG  H021     C02     C03         108.90 3.000
LIG  H022     C02     C01         108.90 3.000
LIG  H021     C02     C01         108.90 3.000
LIG   C03     C02     C01         112.29 3.000
LIG   O07     C03     O04         120.00 3.000
LIG   O07     C03     C02         120.00 3.000
LIG   O04     C03     C02         120.00 3.000
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
LIG Var_01         O04     C03     C02     C01          90.00  30.0 2
LIG Var_02        H011     C01     C02     C03         -60.00  30.0 3
"""

def tst_1(prefix="tst_ref_model_ligands_1"):
  with open("%s.pdb" % prefix, 'w') as f:
    f.write(pdb_str)
  with open("%s.cif" % prefix, 'w') as f:
    f.write(restr_string)
  cmd = " ".join([
      "phenix.geometry_minimization",
      "%s.pdb" % prefix,
      "%s.cif" % prefix,
      "reference_model.enabled=True",
      "reference_model.file=%s.pdb" % prefix,
      ">%s.log" % prefix])
  print(cmd)
  assert not easy_run.call(cmd)
  assert_lines_in_file(
      file_name="%s.log" % prefix,
      lines=
          """Model:              Reference:
          LIG A   1  <=====>  LIG A   1
          Total # of matched residue pairs: 1
          Total # of reference model restraints: 1""")

def run():
  tst_1()
  print("OK")

if (__name__ == "__main__"):
  run()
