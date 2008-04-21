from iotbx import pdb
from cctbx.array_family import flex
import sys

def run(args):
  assert len(args) == 1
  pdb_inp = pdb.input(file_name=args[0])
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_atoms = pdb_hierarchy.atoms()

  # add random numbers [-0.5,0.5) to coordinates
  new_xyz = pdb_atoms.extract_xyz() + flex.vec3_double(
    flex.random_double(size=pdb_atoms.size()*3)-0.5)
  pdb_atoms.set_xyz(new_xyz=new_xyz)

  # reset B-factors (min=1, max=20)
  new_b = flex.random_double(size=pdb_atoms.size(), factor=19) + 1
  pdb_atoms.set_b(new_b=new_b)

  sys.stdout.write(pdb_hierarchy.as_pdb_string(
    crystal_symmetry=pdb_inp.crystal_symmetry(),
    append_end=True))

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
