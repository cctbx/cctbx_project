from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run
import libtbx.load_env

def exercise():
  if (not libtbx.env.has_module("phenix_regression")):
    print("phenix_regression not configured, skipping tst_pdbx_mmcif_tutorial.py")
    return
  mmcif_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3orl.cif",
    test=os.path.isfile)
  py_file = libtbx.env.find_in_repositories(
    relative_path="iotbx/examples/pdbx_mmcif_tutorial.py",
    test=os.path.isfile)
  command = "iotbx.python %s %s" %(py_file, mmcif_file)
  print(command)
  result = easy_run.fully_buffered(command=command).raise_if_errors()

if __name__ == "__main__":
  exercise()
  print("OK")
