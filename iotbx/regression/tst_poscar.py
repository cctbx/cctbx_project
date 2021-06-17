from __future__ import absolute_import, division, print_function
def exercise():
  import libtbx.load_env
  import os
  op = os.path
  inp_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/poscar",
    test=op.isdir)
  if (inp_dir is None):
    print("Skipping POSCAR tests: input files not available")
    return
  import iotbx.poscar
  for node in os.listdir(inp_dir):
    if (not node.startswith("POSCAR-")): continue
    file_name = op.join(inp_dir, node)
    with open(file_name) as f:
      lines = f.read().splitlines()
    poscar = iotbx.poscar.reader(lines=lines)
    assert poscar.make_up_types_if_necessary() is poscar
    xs = poscar.xray_structure(u_iso=0.1)
    assert str(xs.space_group_info()) == "P 1"
    assert xs.scatterers().size() == len(poscar.sites)

def run(args):
  assert len(args) == 0
  exercise()
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
