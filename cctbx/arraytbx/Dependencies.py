# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.lib_python_subdir = "cctbx_boost/arraytbx"

    self.files = (
      "arraytbx/sharedmodule.cpp",
      "arraytbx/shared_mapmodule.cpp",
      "arraytbx/shared_picklers.cpp",
      "arraytbx/flexmodule.cpp",
      "arraytbx/flex_picklers.cpp",
      "arraytbx/tst_shared.py",
      "arraytbx/tst_shared_map.py",
      "arraytbx/tst_flex.py",
      "arraytbx/tst_af_helpers.cpp",
      "arraytbx/tst_af_1.cpp",
      "arraytbx/tst_af_2.cpp",
      "arraytbx/tst_af_3.cpp",
      "arraytbx/tst_vec3.cpp",
      "arraytbx/tst_mat3.cpp",
      "arraytbx/tst_sym_mat3.cpp",
    )

    self.executables = {
      "tst_af_1": (("tst_af_1",), ()),
      "tst_af_2": (("tst_af_2",), ()),
      "tst_af_3": (("tst_af_3",), ("cctbx_misc",)),
      "tst_vec3": (("tst_vec3",), ()),
      "tst_mat3": (("tst_mat3",), ("cctbx_misc",)),
      "tst_sym_mat3": (("tst_sym_mat3",), ("cctbx_misc",)),
    }

    self.boost_python_modules = {
      "shared":
        (("sharedmodule", "shared_picklers"),
        ("eltbx", "cctbx_misc", "cctbx_bpl1")),
      "shared_map":
        (("shared_mapmodule",),
        ("cctbx_misc", "cctbx_bpl1")),
      "flex":
        (("flexmodule", "flex_picklers"),
        ("cctbx_misc", "cctbx_bpl1")),
    }

