from __future__ import absolute_import, division, print_function
from mmtbx.regression import model_1yjp, model_1yjp_with_waters
from iotbx.data_manager import DataManager

count_1yjp = {
              'Bond | covalent geometry | restraints': 59,
              'Bond angle | covalent geometry | restraints': 79,
              'Dihedral angle | covalent geometry | restraints': 22,
              '    harmonic':  7,
              '  sinusoidal': 15,
              'Planarity | covalent geometry | restraints': 13,
              'Chirality | covalent geometry | restraints': 6,
              # 'C-Beta improper torsion angle restraints': 12,
              'Dihedral angle | C-Beta improper | restraints': 12,
              # 'Parallelity restraints': 0,
              #'User supplied restraints': 0,
              #'User supplied torsion angle restraints': 0,
              #'User supplied angle restraints': 0,
              #'Metal coordination angle restraints': 0,
              #'Disulphide bridge torsion angle restraints': 0,
              #'Side chain torsion angle restraints': 0,
              #'Secondary Structure restraints around h-bond angle restraints': 0,
              #'Bond-like restraints': 0,
              #'Metal coordination restraints': 0,
              #'Disulphide bridge angle restraints': 0,
              #'Disulphide bridge restraints': 0,
              'Nonbonded | unspecified | interactions': 990,
              #
              # 'Bond | User supplied | restraints': -1,
              }
count_1yjp_with_waters = count_1yjp.copy()
count_1yjp_with_waters['Nonbonded | unspecified | interactions'] = 1178

edits = '''
refinement.geometry_restraints.edits {
  bond {
    action = *add
    atom_selection_1 = resname HOH and resid 10 and name O
    atom_selection_2 = resname ASN and resid 2 and name ND2
    symmetry_operation = None
    distance_ideal = 2.1
    sigma = 0.02
    slack = None
  }
}
'''

def check_geo(geo_lines):
  rc = {}
  for line in geo_lines.splitlines():
    tmp = line.split(':')
    tmp = list(filter(None, tmp))
    if len(tmp)!= 2: continue
    rc[tmp[0]]=int(tmp[1])
  return rc

def check_diff(d1,d2):
  outl = 'diff\n%s\n%s\n' % (d1,d2)
  for key in d1:
    print(d1)
    print(d2)
    if d1[key]!=d2[key]:
      outl += '%s : %d %d\n' % (key, d1[key], d2[key])
  return outl

def main():
  dm = DataManager()
  dm.process_model_str('testing', model_1yjp)
  model = dm.get_model()
  rc = model.restraints_as_geo(force=True)
  rc = check_geo(rc)
  assert rc == count_1yjp, check_diff(rc, count_1yjp)

  dm = DataManager()
  dm.process_model_str('testing', model_1yjp_with_waters)
  model = dm.get_model()
  rc = model.restraints_as_geo(force=True)
  rc = check_geo(rc)
  assert rc == count_1yjp_with_waters, '%s != %s' % (count_1yjp_with_waters, rc)

  params = model.get_default_pdb_interpretation_params()
  edits_1yjp = params.geometry_restraints.edits

  edits_1yjp.bond[0].action='add'
  edits_1yjp.bond[0].atom_selection_1='resname HOH and resid 10 and name O'
  edits_1yjp.bond[0].atom_selection_2='resname ASN and resid 2 and name ND2'
  edits_1yjp.bond[0].distance_ideal=2.1
  edits_1yjp.bond[0].sigma=0.1
  model.process(pdb_interpretation_params=params,
                make_restraints=True)
  rc = model.restraints_as_geo(force=True)
  rc = check_geo(rc)
  current = count_1yjp_with_waters.copy()
  # current['User supplied restraints'] = 1
  current['Bond | User supplied | restraints'] = 1
  current['Nonbonded | unspecified | interactions']   = 1176
  assert rc == current, check_diff(rc, current)

  edits_1yjp.angle[0].action='add'
  edits_1yjp.angle[0].atom_selection_1='resname HOH and resid 10 and name O'
  edits_1yjp.angle[0].atom_selection_2='resname ASN and resid 2 and name ND2'
  edits_1yjp.angle[0].atom_selection_3='resname ASN and resid 2 and name CG'
  edits_1yjp.angle[0].angle_ideal=21.9
  edits_1yjp.angle[0].sigma=1.1
  model.process(pdb_interpretation_params=params,
                make_restraints=True)
  rc = model.restraints_as_geo(force=True)
  rc = check_geo(rc)
  current = count_1yjp_with_waters.copy()
  current['Bond | User supplied | restraints'] = 1
  current['Bond angle | User supplied | restraints'] = 1
  current['Nonbonded | unspecified | interactions']   = 1176
  assert rc == current, check_diff(rc, current)

  edits_1yjp.dihedral[0].action='add'
  edits_1yjp.dihedral[0].atom_selection_1='resname HOH and resid 10 and name O'
  edits_1yjp.dihedral[0].atom_selection_2='resname ASN and resid 2 and name ND2'
  edits_1yjp.dihedral[0].atom_selection_3='resname ASN and resid 2 and name CG'
  edits_1yjp.dihedral[0].atom_selection_4='resname ASN and resid 2 and name CB'
  edits_1yjp.dihedral[0].angle_ideal=121.9
  edits_1yjp.dihedral[0].sigma=1.12
  edits_1yjp.dihedral[0].periodicity=10
  model.process(pdb_interpretation_params=params,
                make_restraints=True)
  rc = model.restraints_as_geo(force=True)
  rc = check_geo(rc)
  current = count_1yjp_with_waters.copy()
  current['Bond | User supplied | restraints'] = 1
  current['Bond angle | User supplied | restraints'] = 1
  current['Dihedral angle | User supplied | restraints'] = 1
  #current['  sinusoidal'] = 16
  current['Nonbonded | unspecified | interactions']   = 1176
  assert rc == current, check_diff(rc, current)
  print('OK')

if __name__ == '__main__':
  main()
