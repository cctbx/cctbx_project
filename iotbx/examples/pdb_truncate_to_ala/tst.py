import v0_getting_started
import v1_loop_over_atoms
import v2_simple
import v3_better
import v4_with_bells_and_whistles
from libtbx.path import is_same_file
import libtbx.load_env
import sys, os

def run(args):
  assert len(args) == 0
  tutorial_dir = libtbx.env.under_dist(
    module_name="iotbx",
    path="examples/pdb_truncate_to_ala",
    test=os.path.isdir)
  if ("set" not in libtbx.forward_compatibility.__builtins__):
    libtbx.forward_compatibility.__builtins__["set"] = list
  for file_name in ["crambin_pieces.pdb", "resname_mix.pdb"]:
    file_path = os.path.join(tutorial_dir, file_name)
    if (   not os.path.isfile(file_name)
        or not is_same_file(file_names=[file_path, file_name])):
      libtbx.utils.copy_file(source=file_path, target=file_name)
    for vx in [v0_getting_started,
               v1_loop_over_atoms,
               v2_simple,
               v3_better,
               v4_with_bells_and_whistles]:
      vx.run(args=[file_name])
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
