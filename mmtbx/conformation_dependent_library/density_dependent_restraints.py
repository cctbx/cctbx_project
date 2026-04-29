from __future__ import absolute_import, division, print_function
import sys

ddw_master_params = '''
ddw
  .style = hidden
  .short_caption = Density dependent restraints
{
  enable = False
    .type = bool
  bond_weighting = True
    .type = bool
  angle_weighting = False
    .type = bool
}
'''


def setup_restraints():
  assert 0

def get_weight_factor_bond(normalised_map_weights, i, j, atoms=None):
  if atoms: print(i,j,normalised_map_weights[i], normalised_map_weights[j])
  if normalised_map_weights[i]<.5 and normalised_map_weights[j]<.5:
    return 2
  return 1

def update_restraints(hierarchy,
                      geometry,
                      normalised_map_weights,
                      bond_weighting=True,
                      angle_weighting=False,
                      log=None,
                      verbose=False,
                      ):
  atoms = hierarchy.atoms()
  assert len(normalised_map_weights)==len(atoms), 'normalised_map_weights array not the length of atoms'
  bond_params_table = geometry.bond_params_table
  n_bonds=0
  if bond_weighting:
    for i, bonded in enumerate(bond_params_table):
      for j in bonded:
        bond = bond_params_table.lookup(i, j)
        factor=get_weight_factor_bond(normalised_map_weights, i, j, atoms=None)
        if factor==1: continue
        weight=bond.weight
        bond.weight*=factor
        n_bonds+=1
        if verbose:
          print(' bond %s-%s %5.2f %9.2f ~> %9.2f' % (
            atoms[i].quote(),
            atoms[j].quote(),
            bond.distance_ideal,
            weight,
            bond.weight,
            ))

  n_angles=[]
  if angle_weighting:
    assert 0
    for angle_proxy in geometry.angle_proxies:
      if (angle_proxy.i_seqs[0] not in total_i_seqs and
          angle_proxy.i_seqs[1] not in total_i_seqs and
          angle_proxy.i_seqs[2] not in total_i_seqs): continue
      angle_proxy.weight*=factor
      n_angles.append(angle_proxy)
      if verbose:
        i,j,k = angle_proxy.i_seqs
        print('  angle %s-%s-%s %s %s' % (atoms[i].quote(),
                                          atoms[j].quote(),
                                          atoms[k].quote(),
                                          angle_proxy.angle_ideal,
                                          angle_proxy.weight))

  print('  Upweighted %d bonds and %d angles' % (n_bonds+1, len(n_angles)), file=log)

def random_atom_weighting(atoms):
  import random
  rc=[]
  for atom in atoms:
    if atom.parent().resname in ['HOH']:
      rc.append(random.randrange(1, 10)/10)
    elif atom.name.strip() in ['CA', 'C', 'N', 'O']:
      rc.append(.9)
    else:
      rc.append(.3)
  return rc

def run(filename):
  import iotbx
  from mmtbx.model.model import manager
  from iotbx import pdb

  pdb_inp = iotbx.pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()

  m = manager(pdb_inp)
  m.process(make_restraints=True)
  restraints_manager = m.get_restraints_manager()
  normalised_map_weights=random_atom_weighting(hierarchy.atoms())
  print(normalised_map_weights)
  update_restraints(hierarchy,
                    restraints_manager.geometry,
                    normalised_map_weights=normalised_map_weights,
                    verbose=True,
    )


if __name__=="__main__":
  run(sys.argv[1])
