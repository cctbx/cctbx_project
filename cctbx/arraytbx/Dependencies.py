# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "global/error.cpp",
      "global/bpl_utils.cpp",
      "global/tiny_bpl.cpp",
      "arraytbx/sharedmodule.cpp",
      "arraytbx/tst_shared.py",
      "arraytbx/tst_af_helpers.cpp",
      "arraytbx/tst_af_1.cpp",
      "arraytbx/tst_af_2.cpp",
      "arraytbx/debug_overloadmodule.cpp",
      "arraytbx/tst_debug_overload.py",
    )

    self.executables = {
      "tst_af_1": (("tst_af_1",), ()),
    }
    if (self.platform != "vc60"):
      self.executables["tst_af_2"] = (("tst_af_2",), ())

    self.boost_python_modules = {
      "shared":
        (("sharedmodule", "error", "bpl_utils", "tiny_bpl"), ()),
      "debug_overload":
        (("debug_overloadmodule",), ()),
    }

