# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost/arraytbx"

    self.files = (
      "arraytbx/sharedmodule.cpp",
      "arraytbx/tst_shared.py",
      "arraytbx/tst_af_helpers.cpp",
      "arraytbx/tst_af_1.cpp",
      "arraytbx/tst_af_2.cpp",
      "arraytbx/tst_vec3.cpp",
      "arraytbx/tst_mat3.cpp",
    )

    self.executables = {
      "tst_af_1": (("tst_af_1",), ()),
      "tst_af_2": (("tst_af_2",), ()),
      "tst_vec3": (("tst_vec3",), ()),
      "tst_mat3": (("tst_mat3",), ("cctbx_misc",)),
    }

    self.boost_python_modules = {
      "shared":
        (("sharedmodule",), ("cctbx_misc", "cctbx_bpl1")),
    }

