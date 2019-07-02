from __future__ import absolute_import, division, print_function
import sys
from six.moves import cStringIO as StringIO

from libtbx import easy_run

preamble='apply_cif_restraints'

files = {"%s.pdb" % preamble : '''
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
HETATM    1  C01 LIG B   1       9.290  -0.103  -0.646  1.00 20.00      A    C
HETATM    2  C02 LIG B   1      10.237  -0.103  -0.646  1.00 20.00      A    C
HETATM    3  C03 LIG B   1      10.816  -0.103   0.767  1.00 20.00      A    C
HETATM    4  O04 LIG B   1      10.767  -1.086   1.427  1.00 20.00      A    O
HETATM    5  O05 LIG B   1      11.405   1.058   1.282  1.00 20.00      A    O
HETATM    6 H011 LIG B   1       9.652  -1.043  -0.246  1.00 20.00      A    H
HETATM    7 H012 LIG B   1       9.652   0.713  -0.032  1.00 20.00      A    H
HETATM    8 H013 LIG B   1       9.652   0.020  -1.660  1.00 20.00      A    H
HETATM    9 H021 LIG B   1      10.588   0.778  -1.169  1.00 20.00      A    H
HETATM   10 H022 LIG B   1      10.588  -0.985  -1.169  1.00 20.00      A    H
HETATM   11 H051 LIG B   1      11.846   0.855   2.091  1.00 20.00      A    H
''',
  '%s_01.cif' % preamble : '''
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
  ''',
  '%s_02.cif' % preamble : '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LIG        LIG 'Unknown                  ' ligand 11 5 .
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
LIG         C01    C   CH3    0    .      -1.2904   -0.1034   -0.6460
LIG         C02    C   CH2    0    .       0.2367   -0.1034   -0.6460
LIG         C03    C   C      0    .       0.8158   -0.1034    0.7671
LIG         O04    O   O      0    .       0.7666   -1.0857    1.4274
LIG         O05    O   OH1    0    .       1.4049    1.0580    1.2825
LIG        H011    H   HCH3   0    .      -1.6516   -1.0432   -0.2460
LIG        H012    H   HCH3   0    .      -1.6516    0.7129   -0.0321
LIG        H013    H   HCH3   0    .      -1.6515    0.0201   -1.6598
LIG        H021    H   HCH2   0    .       0.5875    0.7781   -1.1689
LIG        H022    H   HCH2   0    .       0.5875   -0.9849   -1.1689
LIG        H051    H   HOH1   0    .       1.8459    0.8548    2.0907
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
LIG   O04     C03   double        1.185 0.020     1.185
LIG   O05     C03   single        1.401 0.020     1.401
LIG  H011     C01   single        0.970 0.020     1.090
LIG  H012     C01   single        0.970 0.020     1.090
LIG  H013     C01   single        0.970 0.020     1.090
LIG  H021     C02   single        0.970 0.020     1.090
LIG  H022     C02   single        0.970 0.020     1.090
LIG  H051     O05   single        0.850 0.020     0.980
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
LIG   O05     C03     O04         120.00 3.000
LIG   O05     C03     C02         120.00 3.000
LIG   O04     C03     C02         120.00 3.000
LIG  H051     O05     C03         109.48 3.000
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
LIG CONST_01      H051     O05     C03     C02        -169.90   0.0 0
LIG Var_01         O04     C03     C02     C01          73.23  30.0 2
LIG Var_02        H011     C01     C02     C03         -66.94  30.0 3
#
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
LIG plan-1    C02 0.020
LIG plan-1    C03 0.020
LIG plan-1    O04 0.020
LIG plan-1    O05 0.020
LIG plan-1   H051 0.020
  ''',
  '%s.eff' % preamble : '''
pdb_interpretation {
  apply_cif_restraints {
    restraints_file_name=%s_02.cif
    residue_selection="chain B and resname LIG"
  }
}
  ''' % preamble,
}

def run():
  for name, text in files.items():
    f=open(name, 'w')
    f.write(text)
    f.close()
  cmd = 'phenix.geometry_minimization %(preamble)s.pdb %(preamble)s_01.cif' % {'preamble' : preamble}
  print(cmd)
  ero = easy_run.fully_buffered(command=cmd)
  err = StringIO()
  ero.show_stderr(out=err)
  assert err.getvalue()
  print('ok')
  cmd += ' %(preamble)s.eff' % {'preamble' : preamble}
  print(cmd)
  ero = easy_run.fully_buffered(command=cmd)
  err = StringIO()
  ero.show_stderr(out=err)
  assert not err.getvalue()
  print('ok')
  return err.getvalue()

if __name__=="__main__":
  run()#sys.argv[1])
