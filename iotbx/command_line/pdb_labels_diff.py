def run(args):
  if (len(args) != 2):
    from libtbx.utils import Usage
    import libtbx.load_env
    raise Usage("%s first.pdb second.pdb" % libtbx.env.dispatcher_name)
  import iotbx.pdb
  lbls_list_pair = []
  for file_name in args:
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    lbls_list = []
    for atom in pdb_inp.atoms():
      lbls_list.append(atom.id_str())
    lbls_list_pair.append("\n".join(lbls_list)+"\n")
  a, b = lbls_list_pair
  if (a != b):
    import difflib
    sys.stdout.write(
      "".join(difflib.unified_diff(a.splitlines(1), b.splitlines(1))))

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
