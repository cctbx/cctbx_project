# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "fftbx/timing/test0.cpp",
      "fftbx/timing/fftbxtimer.cpp",
      "fftbx/timing/fftwtimer.cpp",
      "fftbx/timing/time_cmd.py",
      "fftbx/timing/do_timing.py",
      "fftbx/timing/cmp_times.py",
      "fftbx/timing/tst3d.cpp",
      "fftbx/timing/time3d.py",
      "fftbx/timing/eval_times.py",
      "fftbx/timing/debug_compare.csh",
    )

    self.executables = {
      "test0": (("test0",), ()),
      "fftbxtimer": (("fftbxtimer",), ()),
    }
    if (self.platform in ("tru64_cxx", "unix_gcc", "irix_CC")):
      self.executables["fftwtimer"] = (("fftwtimer",), ("fftw", "rfftw"))
      self.executables["tst3d"] = (("tst3d",), ("fftw", "rfftw"))
