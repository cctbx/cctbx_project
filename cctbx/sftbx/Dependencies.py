# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
      "sftbx/sftbxmodule.cpp",
      "sftbx/tst.py",
      "sftbx/tst_basic.py",
      "sftbx/dbg.py",
    )

    self.boost_python_modules = {
      "sftbx": (("sftbxmodule",),
                ("cctbx_miller", "sgtbx", "uctbx",
                 "cctbx_misc", "cctbx_bpl1")),
    }
