import iotbx.pdb
import sys

def run(args):
  if (len(args) == 0):
    raise RuntimeError("Please specify one or more pdb file names.")
  for file_name in args:
    pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
    pdb_obj.hierarchy.overall_counts().show()
    for model in pdb_obj.hierarchy.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          print 'resid: "%s"' % rg.resid()
          for ag in rg.atom_groups():
            print '  altloc: "%s", resname: "%s"' % (ag.altloc, ag.resname)
            for atom in ag.atoms():
              print '    ', atom.name

if (__name__ == "__main__"):
  run(sys.argv[1:])
