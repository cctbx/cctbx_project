# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
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
      "sgtbx/miller_ref_asu.cpp",
      "sgtbx/miller_asu.cpp",
      "sgtbx/seminvariant.cpp",
      "sgtbx/select_generators.cpp",
      "sgtbx/sgtbxmodule.cpp",
      "sgtbx/sgtbxdriver.cpp",
      "sgtbx/symbols.py",
      "sgtbx/tst.py",
      "sgtbx/tst6.py",
      "sgtbx/tst7.py",
      "sgtbx/tst9.py",
      "sgtbx/tst10.py",
      "sgtbx/tst11.py",
      "sgtbx/tst12.py",
      "sgtbx/tst13.py",
      "sgtbx/tst14.py",
      "sgtbx/tst15.py",
    )

    lib = (
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
      "miller_ref_asu",
      "miller_asu",
      "seminvariant",
      "select_generators",
    )

    self.libraries = {
      "sgtbx": lib,
    }

    self.executables = {
      "sgtbxdriver": (("sgtbxdriver",), ("sgtbx", "uctbx", "cctbx_misc")),
    }

    self.boost_python_modules = {
      "sgtbx": (("sgtbxmodule",),
        ("sgtbx", "uctbx", "cctbx_misc", "cctbx_bpl1")),
    }
