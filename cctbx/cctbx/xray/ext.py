import boost.python
ext = boost.python.import_ext("cctbx_xray_ext")
from cctbx_xray_ext import *

import sys

class _scattering_dictionary(boost.python.injector, scattering_dictionary):

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
        types_only=False,
        header="Scattering types:",
        f=None,
        prefix=""):
    if (f is None): f = sys.stdout
    if (header is not None):
      print >> f, prefix + header
    gn = ""
    for key,val in self.dict().items():
      print >> f, "%s  %s: %d" % (prefix, key, val.member_indices.size()),
      if (not types_only):
        gn = str(val.gaussian.n_terms())
        if (val.gaussian.c() != 0): gn += "+c"
        print >> f, gn,
      print >> f
    print >> f
