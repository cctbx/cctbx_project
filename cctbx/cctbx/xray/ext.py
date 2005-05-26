import boost.python
ext = boost.python.import_ext("cctbx_xray_ext")
from cctbx_xray_ext import *

from cctbx.array_family import flex
import sys

class _scattering_dictionary(boost.python.injector, scattering_dictionary):

  def sorted_items(self, heaviest_first=True):
    items = self.dict().items()
    weights = flex.double([val.gaussian.at_stol(0) for key,val in items])
    perm = flex.sort_permutation(weights, reverse=heaviest_first)
    return flex.select(items, permutation=perm)

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    for key,val in self.dict().items():
      gn = str(val.gaussian.n_terms())
      if (val.gaussian.c() != 0):
        gn += "+c"
      print >> f, "%s%s:%s*%d" % (prefix, key, gn, val.member_indices.size()),
      prefix = ""
    print >> f

  def show(self,
        header="Number of scattering types:",
        show_weights=True,
        show_gaussians=True,
        f=None,
        prefix=""):
    if (f is None): f = sys.stdout
    items = self.sorted_items()
    if (header is not None):
      print >> f, prefix + header, len(items)
    nk = max(3,max([len(key) for key,val in items]))
    nk_fmt = "%%-%ds " % nk
    nn = max(5,len(str(max([val.member_indices.size() for key,val in items]))))
    nn_fmt = "%%%dd" % nn
    if (len(items) > 0):
      line = prefix + "  Type%s %sNumber" % (" "*(nk-3), " "*(nn-5))
      if (show_weights): line += "   Weight"
      if (show_gaussians): line += "   Gaussians"
      print >> f, line
      for key,val in items:
        line = prefix + "   " \
             + nk_fmt%key \
             + nn_fmt%val.member_indices.size() + " "
        if (show_weights):
          line += " %8.2f" % val.gaussian.at_stol(0)
        if (show_gaussians):
          line += " %7s" % str(val.gaussian.n_terms())
          if (val.gaussian.c() != 0): line += "+c"
        print >> f, line.rstrip()

  def wilson_dict(self):
    result = {}
    for key,val in self.dict().items():
      result[key] = val.member_indices.size()
    return result
