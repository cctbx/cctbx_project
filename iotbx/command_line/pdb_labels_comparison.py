"""Compare labels in two models"""
from __future__ import absolute_import, division, print_function
def run(args):
  if (len(args) != 2):
    from libtbx.utils import Usage
    import libtbx.load_env
    raise Usage("%s first.pdb second.pdb" % libtbx.env.dispatcher_name)
  import iotbx.pdb
  inp_hier = []
  lbls_sets = []
  for file_name in args:
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    print("%6d atoms in %s" % (pdb_inp.atoms().size(), file_name))
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    inp_hier.append((pdb_inp, pdb_hierarchy))
      # need to keep hierarchy alive to get all labels
    lbls_set = set()
    for atom in pdb_inp.atoms():
      lbls = atom.id_str()
      if (lbls in lbls_set):
        raise RuntimeError("Duplicate atom labels: %s" % lbls)
      lbls_set.add(lbls)
    lbls_sets.append(lbls_set)
  print("%6d matching atom labels" % (
    len(lbls_sets[0].intersection(lbls_sets[1]))))
  def show_missing(i, j):
    diff = lbls_sets[j].difference(lbls_sets[i])
    print("%6d missing in %s" % (len(diff), args[i]))
    for atom in inp_hier[j][0].atoms():
      lbls = atom.id_str()
      if (lbls in diff):
        print("         %s" % lbls)
  show_missing(0, 1)
  show_missing(1, 0)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
