# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
      "misc/error.cpp",
      "misc/bpl_utils.cpp",
      "misc/tiny_bpl.cpp",
      "misc/small_bpl.cpp",
      "misc/devmodule.cpp",
      "misc/tst.py",
    )

    self.libraries = {
      "cctbx_misc": ("error",),
      "cctbx_bpl1": ("bpl_utils", "tiny_bpl", "small_bpl"),
    }

    self.boost_python_modules = {
      "dev": (("devmodule",), ("cctbx_misc", "cctbx_bpl1")),
    }
