# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "global/error.cpp",
      "global/bpl_utils.cpp",
      "global/boost_array_bpl.cpp",
      "uctbx/uctbx.cpp",
      "uctbx/uctbxmodule.cpp",
      "uctbx/uctbxdriver.cpp",
      "uctbx/tst.py",
    )

    self.libraries = {
      "uctbx": ("uctbx", "error"),
    }

    self.executables = {
      "uctbxdriver": (("uctbxdriver", "uctbx", "error"), ()),
    }

    self.boost_python_modules = {
      "uctbx": (("uctbxmodule",
                 "uctbx",
                 "error",
                 "bpl_utils", "boost_array_bpl"), ()),
    }
