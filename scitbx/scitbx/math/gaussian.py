import scitbx.array_family.flex # for tuple mappings

from scitbx.math.ext import gaussian_term as term
from scitbx.math.ext import gaussian_sum as sum
from scitbx.math.ext import gaussian_fit as fit

import boost.python
import sys

class _sum(boost.python.injector, sum):

  def show(self, f=None, format=None):
    if (f is None): f = sys.stdout
    if (format is None): format = "%.8g"
    for l,v in (("a:", self.array_of_a()), ("b:", self.array_of_b())):
      print >> f, l, " ".join([format % x for x in v])
    if (self.use_c()):
      print >> f, "c:", format % self.c()
    return self
