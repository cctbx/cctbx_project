from mmtbx.monomer_library import pdb_interpretation
import sys

def run():
  pdb_interpretation.run(
    args=sys.argv[1:],
    strict_conflict_handling=False)

if (__name__ == "__main__"):
  run()
