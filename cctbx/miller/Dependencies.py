# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
      "miller/miller_lib.cpp",
      "miller/millermodule.cpp",
      "miller/tst.py",
    )

    self.libraries = {
      "cctbx_miller": ("miller_lib",),
    }

    self.boost_python_modules = {
      "miller": (("millermodule",),
        ("cctbx_miller", "sgtbx", "uctbx", "cctbx_misc", "cctbx_bpl1")),
    }
