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

class _scattering_type_registry(
        boost.python.injector, scattering_type_registry):

  def sorted_type_index_pairs(self, heaviest_first=True):
    ugs = self.unique_gaussians_as_list()
    pairs = []
    weights = flex.double()
    for t,i in self.type_index_pairs_as_dict().items():
      pairs.append((t,i))
      gaussian = ugs[i]
      if (gaussian is None):
        weights.append(0)
      else:
        weights.append(gaussian.at_stol(0))
    perm = flex.sort_permutation(weights, reverse=heaviest_first)
    return flex.select(pairs, permutation=perm)

  def show_summary(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    unique_gaussians = self.unique_gaussians_as_list()
    unique_counts = list(self.unique_counts)
    for t,i in self.sorted_type_index_pairs():
      gaussian = unique_gaussians[i]
      if (gaussian is None):
        gn = "None"
      else:
        gn = str(gaussian.n_terms())
        if (gaussian.c() != 0):
          gn += "+c"
      print >> out, "%s%s:%s*%d" % (prefix, t, gn, unique_counts[i]),
      prefix = ""
    print >> out

  def show(self,
        header="Number of scattering types:",
        show_weights=True,
        show_gaussians=True,
        out=None,
        prefix=""):
    if (out is None): out = sys.stdout
    unique_gaussians = self.unique_gaussians_as_list()
    unique_counts = list(self.unique_counts)
    tips = self.sorted_type_index_pairs()
    if (header is not None):
      print >> out, prefix + header, len(tips)
    nt = max(3,max([len(t) for t,i in tips]))
    nt_fmt = "%%-%ds " % nt
    nc = max(5,len(str(max(unique_counts))))
    nc_fmt = "%%%dd" % nc
    if (len(tips) > 0):
      line = prefix + "  Type%s %sNumber" % (" "*(nt-3), " "*(nc-5))
      if (show_weights): line += "   Weight"
      if (show_gaussians): line += "   Gaussians"
      print >> out, line
      for t,i in tips:
        line = prefix + "   " \
             + nt_fmt%t \
             + nc_fmt%unique_counts[i] + " "
        gaussian = unique_gaussians[i]
        if (show_weights):
          if (gaussian is None):
            line += "     None"
          else:
            line += " %8.2f" % gaussian.at_stol(0)
        if (show_gaussians):
          if (gaussian is None):
            line += "      None"
          else:
            line += " %7s" % str(gaussian.n_terms())
            if (gaussian.c() != 0): line += "+c"
        print >> out, line.rstrip()

  def wilson_dict(self):
    result = {}
    unique_counts = list(self.unique_counts)
    for t,i in self.sorted_type_index_pairs():
      result[t] = unique_counts[i]
    return result
