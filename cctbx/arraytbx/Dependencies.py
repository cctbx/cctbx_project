# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost/arraytbx"

    self.files = (
      "arraytbx/flexmodule.cpp",
      "arraytbx/flex_picklers.cpp",
      "arraytbx/flex_utilsmodule.cpp",
      "arraytbx/tst_flex.py",
      "arraytbx/tst_flex_utils.py",
      "arraytbx/tst_af_helpers.cpp",
      "arraytbx/tst_af_1.cpp",
      "arraytbx/tst_af_2.cpp",
      "arraytbx/tst_af_3.cpp",
      "arraytbx/tst_af_4.cpp",
      "arraytbx/tst_vec3.cpp",
      "arraytbx/tst_mat3.cpp",
      "arraytbx/tst_sym_mat3.cpp",
      "arraytbx/tst_all.py",
    )

    self.executables = {
      "tst_af_1": (("tst_af_1",), ()),
      "tst_af_2": (("tst_af_2",), ("cctbx_misc",)),
      "tst_af_3": (("tst_af_3",), ()),
      "tst_af_4": (("tst_af_4",), ()),
      "tst_vec3": (("tst_vec3",), ()),
      "tst_mat3": (("tst_mat3",), ("cctbx_misc",)),
      "tst_sym_mat3": (("tst_sym_mat3",), ("cctbx_misc",)),
    }

    self.boost_python_modules = {
      "flex":
        (("flexmodule", "flex_picklers"),
        ("eltbx", "cctbx_misc", "cctbx_bpl1")),
      "flex_utils":
        (("flex_utilsmodule",),
        ("cctbx_misc", "cctbx_bpl1")),
    }

