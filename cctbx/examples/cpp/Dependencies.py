# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):
    self.files = (
      "examples/cpp/getting_started.cpp",
    )

    self.examples = {
      "getting_started": ("getting_started",),
    }
