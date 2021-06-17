from __future__ import absolute_import, division, print_function
import scitbx.array_family.flex # for tuple mappings

from scitbx.math.ext import gaussian_term as term # implicit import
from scitbx.math.ext import gaussian_sum as sum
from scitbx.math.ext import gaussian_fit as fit
from scitbx.array_family import flex

import boost_adaptbx.boost.python as bp
import sys
from six.moves import zip

@bp.inject_into(sum)
class _():

  def show(self, f=None, format=None):
    if (f is None): f = sys.stdout
    if (format is None): format = "%.8g"
    for l,v in (("a:", self.array_of_a()), ("b:", self.array_of_b())):
      print(l, " ".join([format % x for x in v]), file=f)
    if (self.use_c()):
      print("c:", format % self.c(), file=f)
    return self

  def sort(self):
    perm = flex.sort_permutation(
      data=flex.abs(flex.double(self.array_of_a())),
      reverse=True)
    return sum(
      flex.select(self.array_of_a(), perm),
      flex.select(self.array_of_b(), perm),
      self.c(),
      self.use_c())

@bp.inject_into(fit)
class _():

  def show_table(self, f=None):
    if (f is None): f = sys.stdout
    print("stol  table   fitted   delta rel_del", file=f)
    for x,y,v,e in zip(self.table_x(),
                       self.table_y(),
                       self.fitted_values(),
                       self.significant_relative_errors()):
      print("%4.2f %7.4f %7.4f %7.4f %7.4f" % (x,y,v,v-y,e), file=f)

  def __getinitargs__(self):
    return (self.table_x(), self.table_y(), self.table_sigmas(), sum(self))

  def sort(self):
    return fit(self.table_x(), self.table_y(), self.table_sigmas(),
               sum.sort(self))
