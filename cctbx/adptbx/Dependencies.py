# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "global/error.cpp",
      "global/bpl_utils.cpp",
      "global/carray_bpl.cpp",
      "adptbx/adptbxmodule.cpp",
      "adptbx/tst.py",
    )

    self.boost_python_modules = {
      "adptbx": (("adptbxmodule",
                  "error",
                  "bpl_utils", "carray_bpl"), ()),
    }
