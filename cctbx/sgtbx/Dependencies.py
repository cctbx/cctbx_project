# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

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
      "sgtbx/coordinates.cpp",
      "sgtbx/sgtbxmodule.cpp",
      "sgtbx/sgtbxdriver.cpp",
      "sgtbx/symbols.py",
      "sgtbx/tst_symbols.py",
      "sgtbx/tst_unitcell.py",
      "sgtbx/tst1.py",
      "sgtbx/tst2.py",
      "sgtbx/tst3.py",
      "sgtbx/tst4.py",
      "sgtbx/tst5.py",
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
      "coordinates",
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
