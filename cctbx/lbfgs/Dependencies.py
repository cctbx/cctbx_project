# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost"

    self.files = (
      "lbfgs/lbfgs_driver.cpp",
      "lbfgs/lbfgsmodule.cpp",
      "lbfgs/tst.py",
    )

    self.executables = {
      "lbfgs_driver": (("lbfgs_driver",), ("cctbx_misc",)),
    }

    self.boost_python_modules = {
      "lbfgs": (("lbfgsmodule",), ("cctbx_misc",)),
    }
