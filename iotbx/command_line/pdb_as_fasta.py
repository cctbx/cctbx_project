import iotbx.pdb
import sys, os

def run(args):
  for arg in args:
    assert os.path.isfile(arg)
    pdb_obj = iotbx.pdb.input(file_name=arg)
    hierarchy = pdb_obj.construct_hierarchy()
    for model in hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          f = conformer.format_fasta()
          if (f is not None):
            print "\n".join(f)

if (__name__ == "__main__"):
  run(sys.argv[1:])
