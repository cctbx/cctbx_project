from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run

water = '''CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
HETATM  108  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM  109  H1  HOH A   8      -5.757   5.261   7.584  1.00 22.62           H
HETATM  110  H2  HOH A   8      -6.841   5.979   7.266  1.00 22.62           H
END
'''

hoh = '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
HOH        HOH 'water                    ' ligand 3 1 .
#
data_comp_HOH
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
HOH         O      O   OH2    0    .      -0.2095    0.0000   -0.2963
HOH         H1     H   HOH2   0    .       0.7334    0.0000   -0.2963
HOH         H2     H   HOH2   0    .      -0.5239    0.0000    0.5926
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
HOH   O       H1    single        %s.943 0.020
HOH   O       H2    single        %s.943 0.020
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
HOH   H2      O       H1          180.00 3.000
'''

inp='tst_user_supplied_library.pdb'
base_cmd = 'phenix.geometry_minimization %s' % inp
cmds = [
    base_cmd,
    '%s HOH.cif' % (base_cmd),
    '%s pdb_interpretation.restraints_library.user_supplied.path=usl' % (base_cmd),
    '%s pdb_interpretation.restraints_library.user_supplied.path=usl' % (base_cmd),
    '%s pdb_interpretation.restraints_library.user_supplied.path=usl' % (base_cmd),
    '%s HOH.cif pdb_interpretation.restraints_library.user_supplied.path=usl' % (base_cmd),
    '%s HOH.cif pdb_interpretation.restraints_library.user_supplied.path=usl' % (base_cmd),
    '%s HOH.cif pdb_interpretation.restraints_library.user_supplied.path=usl' % (base_cmd),
        ]
cmds[3]+=' pdb_interpretation.restraints_library.user_supplied.action=pre'
cmds[4]+=' pdb_interpretation.restraints_library.user_supplied.action=post'
cmds[6]+=' pdb_interpretation.restraints_library.user_supplied.action=pre'
cmds[7]+=' pdb_interpretation.restraints_library.user_supplied.action=post'

results = [
  '''bond pdb=" O   HOH A   8 "
     pdb=" H2  HOH A   8 "
  ideal  model  delta    sigma   weight residual
  0.850  0.850''',
  '''bond pdb=" O   HOH A   8 "
     pdb=" H2  HOH A   8 "
  ideal  model  delta    sigma   weight residual
  1.943  1.943''',
  '''bond pdb=" O   HOH A   8 "
     pdb=" H1  HOH A   8 "
  ideal  model  delta    sigma   weight residual
  10.943 10.943''',
  '''bond pdb=" O   HOH A   8 "
     pdb=" H1  HOH A   8 "
  ideal  model  delta    sigma   weight residual
  10.943 10.943''',
  '''bond pdb=" O   HOH A   8 "
     pdb=" H2  HOH A   8 "
  ideal  model  delta    sigma   weight residual
  0.850  0.850''',
  '''bond pdb=" O   HOH A   8 "
     pdb=" H2  HOH A   8 "
  ideal  model  delta    sigma   weight residual
  5.943  5.943''',
  '''bond pdb=" O   HOH A   8 "
     pdb=" H2  HOH A   8 "
  ideal  model  delta    sigma   weight residual
  6.943  6.943''',
  '''bond pdb=" O   HOH A   8 "
     pdb=" H2  HOH A   8 "
  ideal  model  delta    sigma   weight residual
  7.943  7.943''',
  ]

def write_water(a,b):
  f=open('HOH.cif', 'w')
  f.write(hoh % (a,b))
  del f

def main():
  f=open(inp, 'w')
  f.write(water)
  del f
  if not os.path.exists('usl'): os.mkdir('usl')
  os.chdir('usl')
  write_water(10,10)
  os.chdir('..')
  for i, cmd in enumerate(cmds):
    print(i, cmd)
    write_water(i,i)
    easy_run.go(cmd)
    f=open(inp.replace('.pdb', '_minimized.geo'), 'r')
    lines=f.read()
    del f
    assert lines.find(results[i])>-1, '"%s"\n not found' % results[i]

if __name__ == '__main__':
  main()
