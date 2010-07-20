import smtbx.stl.map
import smtbx.array_family
import boost.python
boost.python.import_ext("smtbx_refinement_constraints_ext")
from smtbx_refinement_constraints_ext import *

class _parameter(boost.python.injector, parameter):

  def arguments(self):
    for i in xrange(self.n_arguments):
      yield self.argument(i)
