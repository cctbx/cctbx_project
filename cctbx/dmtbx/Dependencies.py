# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
      "dmtbx/dmtbxmodule.cpp",
      "dmtbx/dbg.py",
    )

    self.boost_python_modules = {
      "dmtbx": (("dmtbxmodule",),
                ("sgtbx", "uctbx", "cctbx_misc", "cctbx_bpl1")),
    }
