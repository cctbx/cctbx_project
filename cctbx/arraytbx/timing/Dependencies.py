# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "arraytbx/timing/run_timing.cpp",
      "arraytbx/timing/ublas_shared.cpp",
    )

    self.executables = {
      "run_timing": (("run_timing",), ()),
      "ublas_shared": (("ublas_shared",), ()),
    }
