# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
      "fftbx/fftbxmodule.cpp",
      "fftbx/tst.py",
    )

    self.boost_python_modules = {
      "fftbx": (("fftbxmodule",), ("cctbx_bpl1",)),
    }
