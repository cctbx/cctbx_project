# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import wxtbx.app
from wxtbx import pdb_editor
import sys

if (__name__ == "__main__") :
  app = wxtbx.app.CCTBXApp(0)
  frame = pdb_editor.PDBTreeFrame(None, -1, "PDB Editor")
  frame.Show()
  if (len(sys.argv) > 1) :
    frame.LoadPDB(sys.argv[1])
  app.MainLoop()
