from __future__ import division, print_function
import os

from iotbx.cif import restraint_file_merge

scn1 = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 SCN SCN THIOCYANATEION NON-POLYMER 3 3 "."
#
data_comp_SCN
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
 SCN S S S1
 SCN C C CSP
 SCN N N N
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 SCN S C single 1.574 0.02
 SCN C N triple 1.173 0.02
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 SCN S C N 178.4 1.5
 '''

scn2 = '''
 data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
SCN        SCN 'thiocyanate              ' ligand 3 3 .
#
data_comp_SCN
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
SCN         S      S   S     -1    .      -1.0850    0.0000    0.0000
SCN         C      C   CSP    0    .       0.7250    0.0000    0.0000
SCN         N      N   NS     0    .       1.8610    0.0000    0.0000
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
SCN   S       C     single        1.810 0.020     1.810
SCN   C       N     triple        1.136 0.020     1.136
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
SCN   N       C       S           180.00 3.000
'''

def main():
  restraint_file_merge.run([scn1, scn2],
                           'restraint_file_merge_01.cif',
                           no_file_access=True,
                           )
  restraint_file_merge.run([scn2, scn1],
                           'restraint_file_merge_02.cif',
                           no_file_access=True,
                           )
  restraint_file_merge.run([scn1, scn2],
                           'restraint_file_merge_03.cif',
                           no_file_access=True,
                           no_file_but_save=True,
                           )
  restraint_file_merge.run([scn2, scn1],
                           'restraint_file_merge_04.cif',
                           no_file_access=True,
                           no_file_but_save=True,
                           )
  filenames = os.listdir('.')
  assert 'restraint_file_merge_01.cif' not in filenames
  assert 'restraint_file_merge_02.cif' not in filenames
  assert 'restraint_file_merge_03.cif' in filenames
  assert 'restraint_file_merge_04.cif' in filenames

if __name__ == '__main__':
  main()
