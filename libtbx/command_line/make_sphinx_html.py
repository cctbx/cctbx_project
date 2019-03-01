
from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os.path as op
import shutil
import os

if (__name__ == "__main__"):
  cctbx_base = libtbx.env.find_in_repositories("cctbx_project")
  assert (cctbx_base is not None)
  base_dir = op.dirname(cctbx_base)
  dest_dir = op.join(base_dir, "cctbx_docs")
  if op.exists(dest_dir):
    shutil.rmtree(dest_dir)
  os.chdir(op.join(cctbx_base, "sphinx"))
  easy_run.call("make html")
  shutil.move("build/html", dest_dir)
