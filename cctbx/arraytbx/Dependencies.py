# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "global/error.cpp",
      "global/bpl_utils.cpp",
      "global/carray_bpl.cpp",
      "arraytbx/std_vectormodule.cpp",
      "arraytbx/shared_storagemodule.cpp",
      "arraytbx/tst.py",
      "arraytbx/tst_af_helpers.cpp",
      "arraytbx/tst_af_1.cpp",
      "arraytbx/tst_af_2.cpp",
    )

    self.executables = {
      "tst_af_1": (("tst_af_1",), ()),
      "tst_af_2": (("tst_af_2",), ()),
    }

    self.boost_python_modules = {
      "std_vector":
        (("std_vectormodule", "error", "bpl_utils", "carray_bpl"), ()),
      "shared_storage":
        (("shared_storagemodule", "bpl_utils", "carray_bpl"), ()),
    }

