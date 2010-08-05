def run(args):
  assert len(args) == 2
  import iotbx.pdb
  lists_of_atoms = []
  for file_name in args:
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    pdb_inp.construct_hierarchy().only_residue() # raises if more than one
    lists_of_atoms.append(pdb_inp.atoms())
  lookup_dict = {}
  for atom_i in lists_of_atoms[0]:
    lookup_dict[atom_i.name] = atom_i
  atom_pairs = []
  for atom_j in lists_of_atoms[1]:
    atom_i = lookup_dict.get(atom_j.name)
    if (atom_i is not None):
      atom_pairs.append((atom_i, atom_j))
  assert len(atom_pairs) > 2
  from scitbx.array_family import flex
  reference_sites = flex.vec3_double()
  other_sites = flex.vec3_double()
  for pair in atom_pairs:
    reference_sites.append(pair[0].xyz)
    other_sites.append(pair[1].xyz)
  import scitbx.math.superpose
  fit = scitbx.math.superpose.least_squares_fit(
    reference_sites=reference_sites,
    other_sites=other_sites)
  rmsd = fit.other_sites_best_fit().rms_difference(reference_sites)
  print "Number of atoms first pdb: ", lists_of_atoms[0].size()
  print "Number of atoms second pdb:", lists_of_atoms[1].size()
  print "Number of superposed atoms:", reference_sites.size()
  print "RMSD: %.3f" % rmsd

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
