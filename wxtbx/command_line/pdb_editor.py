"""Edit PDB files graphically"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_editor
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import wxtbx.app
from wxtbx import pdb_editor
import libtbx.load_env
import os
import sys

if (__name__ == "__main__"):
  app = wxtbx.app.CCTBXApp(0)
  bmp_path = libtbx.env.find_in_repositories(
    relative_path="gui_resources/icons/custom/tools.png",
    test=os.path.isfile)
  app.SetTaskbarIcon(bmp_path, "PDB Editor")
  frame = pdb_editor.PDBTreeFrame(None, -1, "PDB Editor")
  frame.Show()
  if (len(sys.argv) > 1):
    frame.LoadPDB(sys.argv[1])
  app.MainLoop()

