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
      "sgtbx/utils.cpp",
      "sgtbx/math.cpp",
      "sgtbx/matrix.cpp",
      "sgtbx/lattice_tr.cpp",
      "sgtbx/groups.cpp",
      "sgtbx/miller.cpp",
      "sgtbx/change_basis.cpp",
      "sgtbx/hall_in.cpp",
      "sgtbx/properties.cpp",
      "sgtbx/tidy.cpp",
      "sgtbx/symbols.cpp",
      "sgtbx/metric.cpp",
      "sgtbx/sgtbxmodule.cpp",
      "sgtbx/sgtbxdriver.cpp",
      "sgtbx/symbols.py",
      "sgtbx/tst_symbols.py",
      "sgtbx/tst_unitcell.py",
      "sgtbx/tst1.py",
      "sgtbx/tst2.py",
      "sgtbx/tst3.py",
      "sgtbx/tst4.py",
    )

    lib = (
      "error",
      "utils",
      "math",
      "matrix",
      "lattice_tr",
      "groups",
      "miller",
      "change_basis",
      "hall_in",
      "properties",
      "tidy",
      "symbols",
      "metric",
    )

    self.libraries = {
      "sgtbx": lib,
    }

    self.executables = {
      "sgtbxdriver": ("sgtbxdriver",) + lib,
    }

    self.boost_python_modules = {
      "sgtbx":   ("sgtbxmodule",) + lib + ("bpl_utils", "boost_array_bpl"),
    }

if (__name__ == "__main__"):
  write_makefiles().all_targets("sgtbx")
