# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):
    self.files = (
      "examples/cpp/getting_started.cpp",
      "examples/cpp/convert_ccp4_symop_lib.cpp",
      "examples/python/symop.lib",
    )

    self.examples = {
      "getting_started": ("getting_started",),
      "convert_ccp4_symop_lib": ("convert_ccp4_symop_lib",),
    }
