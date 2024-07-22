from __future__ import absolute_import, division, print_function
import sys

ddr_master_params = '''
ddr
  .style = hidden
  .short_caption = Density dependent restraints
{
  enable = False
    .type = bool
  weighting_factor = 0.5
    .type = float
  cc_minimum = 0.
    .type = float
  map_value_minimum = 0.5
    .type = float
  bond_weighting = True
    .type = bool
  angle_weighting = True
    .type = bool
  trans_peptide = True
    .type = bool
}
'''


def setup_restraints():
  assert 0

def trim_i_seqs(atoms, ddr_i_seqs):
  total_i_seqs = []
  for i_seq in ddr_i_seqs:
    atom = atoms[i_seq]
    if atom.parent().resname in ['HOH']: continue
    total_i_seqs.append(atom.i_seq)
  return total_i_seqs

def expand_i_seqs_to_residue(atoms, ddr_i_seqs):
  total_i_seqs = []
  for i_seq in ddr_i_seqs:
    atom = atoms[i_seq]
    for at in atom.parent().atoms():
      total_i_seqs.append(at.i_seq)
  return total_i_seqs

def get_selected_i_seqs(rg, names):
  rc = []
  for atom in rg.atoms():
    if atom.name in names:
      rc.append(atom.i_seq)
  return rc

def expand_i_seqs_to_neighbouring_residue(hierarchy, ddr_i_seqs):
  ddr_set = set(ddr_i_seqs)
  previous = None
  intersection = None
  total_i_seqs = []
  for rg in hierarchy.residue_groups():
    if intersection:
      rc = get_selected_i_seqs(rg, [' CA ', ' N  ', ' H  '])
      total_i_seqs += rc
    rg_set = []
    for atom in rg.atoms(): rg_set.append(atom.i_seq)
    intersection = ddr_set.intersection(set(rg_set))
    if intersection and previous:
      total_i_seqs += rg_set
      rc = get_selected_i_seqs(previous, [' CA ', ' C  ', ' O  '])
      total_i_seqs += rc
    previous = rg
  return total_i_seqs

def update_restraints(hierarchy,
                      geometry,
                      ddr_i_seqs,
                      factor=5.,
                      bond_weighting=True,
                      angle_weighting=True,
                      all_trans_peptide=True,
                      log=None,
                      verbose=False,
                      ):
  atoms = hierarchy.atoms()
  ddr_i_seqs = trim_i_seqs(atoms, ddr_i_seqs)
  # total_i_seqs = expand_i_seqs_to_residue(atoms, ddr_i_seqs)
  total_i_seqs = expand_i_seqs_to_neighbouring_residue(hierarchy, ddr_i_seqs)
  total_i_seqs.sort()
  # for i in total_i_seqs: print(atoms[i].quote())
  bond_params_table = geometry.bond_params_table
  n_bonds=0
  if bond_weighting:
    for i, bonded in enumerate(bond_params_table):
      if i in total_i_seqs:
        for j in bonded:
          bond = bond_params_table.lookup(i, j)
          bond.weight*=factor
          n_bonds+=1
          if verbose:
            print(' bond %s-%s %s %s' % (atoms[i].quote(),
                                         atoms[j].quote(),
                                         bond.distance_ideal,
                                         bond.weight))
  n_angles=[]
  if angle_weighting:
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
  n_dihedrals = []
  if all_trans_peptide:
    for dihedral_proxy in geometry.dihedral_proxies:
      if (dihedral_proxy.i_seqs[0] not in total_i_seqs and
          dihedral_proxy.i_seqs[1] not in total_i_seqs and
          dihedral_proxy.i_seqs[2] not in total_i_seqs and
          dihedral_proxy.i_seqs[3] not in total_i_seqs): continue
      i,j,k,l = dihedral_proxy.i_seqs
      if verbose:
        print('  dihedral %s-%s-%s-%s %s %s' % (atoms[i].quote(),
                                                atoms[j].quote(),
                                                atoms[k].quote(),
                                                atoms[l].quote(),
                                                dihedral_proxy.angle_ideal,
                                                dihedral_proxy.weight))
      names = [atoms[i].name, atoms[j].name, atoms[k].name, atoms[l].name]
      if names == [' CA ', ' C  ', ' N  ', ' CA ']:
        if abs(dihedral_proxy.angle_ideal)<20.:
          dihedral_proxy.angle_ideal=180.
          n_dihedrals.append(dihedral_proxy)

  print('  Upweighted %d bonds and %d angles' % (n_bonds+1, len(n_angles)), file=log)
  print('  Adjusted %s dihedrals to trans-peptide' % len(n_dihedrals), file=log)
  # geometry.reset_internals()

def run(filename):
  pdb_inp = iotbx.pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_serial()
  update_restraints(hierarchy,
                    restraints_manager,
                    verbose=True,
    )

if __name__=="__main__":
  run(sys.argv[1])
