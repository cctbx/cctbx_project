# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "global/error.cpp",
      "global/bpl_utils.cpp",
      "global/boost_array_bpl.cpp",
      "uctbx/uctbx.cpp",
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
      "sgtbx/type.cpp",
      "sgtbx/normalizers.cpp",
      "sgtbx/wyckoff.cpp",
      "sgtbx/bricks.cpp",
      "sgtbx/miller_asu.cpp",
      "sgtbx/seminvariant.cpp",
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
      "sgtbx/tst6.py",
      "sgtbx/tst7.py",
      "sgtbx/tst8.py",
      "sgtbx/tst9.py",
      "sgtbx/tst10.py",
      "sgtbx/tst11.py",
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
      "type",
      "normalizers",
      "wyckoff",
      "bricks",
      "miller_asu",
      "seminvariant",
    )

    self.libraries = {
      "sgtbx": lib,
    }

    self.executables = {
      "sgtbxdriver": ("sgtbxdriver",) + lib + ("uctbx",),
    }

    self.boost_python_modules = {
      "sgtbx":   ("sgtbxmodule",) + lib
               + ("bpl_utils", "boost_array_bpl", "uctbx"),
    }
