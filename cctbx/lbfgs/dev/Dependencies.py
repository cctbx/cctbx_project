# $Id$

import makefile_generator

class write_makefiles(makefile_generator.write_makefiles):

  def dependencies(self):

    self.files = (
      "lbfgs/dev/lbfgs.f",
      "lbfgs/dev/sdrive.f",
      "lbfgs/dev/lbfgs.pyf",
      "lbfgs/dev/pyfort.py",
      "lbfgs/dev/build_lbfgs_extension.csh",
      "lbfgs/dev/tst_fortran_lbfgs.py",
    )
