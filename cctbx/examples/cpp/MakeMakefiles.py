# $Id$

import MakeMakefilesMaster

class write_makefiles(MakeMakefilesMaster.write_makefiles):

  def dependencies(self):
    self.files = (
      "examples/cpp/getting_started.cpp",
    )

    self.examples = {
      "getting_started": ("getting_started",),
    }
