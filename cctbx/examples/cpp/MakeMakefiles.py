# $Id$

import sys
sys.path.insert(0, "../../build")
import MakeMakefilesMaster

class write_makefiles(MakeMakefilesMaster.write_makefiles):

  def __init__(self):
    self.files = (
      "examples/cpp/getting_started.cpp",
    )

    self.examples = {
      "getting_started": ("getting_started",),
    }

if (__name__ == "__main__"):
  write_makefiles().all_targets("examples/cpp")
