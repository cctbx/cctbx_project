# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
      "uctbx/uctbx.cpp",
      "uctbx/uctbxmodule.cpp",
      "uctbx/uctbxdriver.cpp",
      "uctbx/tst.py",
    )

    self.libraries = {
      "uctbx": ("uctbx",),
    }

    self.executables = {
      "uctbxdriver": (("uctbxdriver", "uctbx"), ("cctbx_misc",)),
    }

    self.boost_python_modules = {
      "uctbx": (("uctbxmodule", "uctbx"), ("cctbx_misc", "cctbx_bpl1")),
    }
