# $Id$

import sys
sys.path.insert(0, "../build")
import MakeMakefilesMaster

class write_makefiles(MakeMakefilesMaster.write_makefiles):

  def __init__(self):
    self.files = (
      "global/error.cpp",
      "global/bpl_utils.cpp",
      "global/boost_array_bpl.cpp",
      "uctbx/uctbx.cpp",
      "uctbx/uctbxmodule.cpp",
      "uctbx/uctbxdriver.cpp",
      "uctbx/tst.py",
    )

    self.libraries = {
      "uctbx": ("uctbx", "error"),
    }

    self.executables = {
      "uctbxdriver": ("uctbxdriver", "uctbx", "error"),
    }

    self.boost_python_modules = {
      "uctbx": ("uctbxmodule",
                "uctbx",
                "error",
                "bpl_utils", "boost_array_bpl"),
    }

if (__name__ == "__main__"):
  write_makefiles().all_targets("uctbx")
