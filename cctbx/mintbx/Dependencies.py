# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
      "mintbx/mintbxmodule.cpp",
      "mintbx/tst.py",
    )

    self.boost_python_modules = {
      "mintbx": (("mintbxmodule",), ("cctbx_misc", "cctbx_bpl1")),
    }
