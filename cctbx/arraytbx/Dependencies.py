# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "global/bpl_utils.cpp",
      "global/carray_bpl.cpp",
      "arraytbx/shared_storagemodule.cpp",
      "arraytbx/tst.py",
    )

    self.boost_python_modules = {
      "shared_storage": (("shared_storagemodule", "bpl_utils", "carray_bpl"),
                         ()),
    }
