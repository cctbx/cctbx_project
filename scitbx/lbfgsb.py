from __future__ import absolute_import, division, print_function
import scitbx.array_family.flex

import boost.python
ext = boost.python.import_ext("scitbx_lbfgsb_ext")
from scitbx_lbfgsb_ext import *

from scitbx.array_family import flex

class minimizer(ext.minimizer):

  def __init__(self, n=None,
                     m=None,
                     l=None,
                     u=None,
                     nbd=None,
                     enable_stp_init=False,
                     factr=None,
                     pgtol=None,
                     iprint=None):
    assert [l,u,nbd].count(None) in [0,3]
    assert n is not None or l is not None
    if (n is None):
      assert u.size() == l.size() and nbd.size() == l.size()
      n = l.size()
    elif (l is None):
      l = flex.double(n, 0)
      u = l
      nbd = flex.int(n, 0)
    if (m is None): m = 5
    if (factr is None): factr = 1.0e+7
    if (pgtol is None): pgtol = 1.0e-5
    if (iprint is None): iprint = -1
    ext.minimizer.__init__(self,
      n, m, l, u, nbd, enable_stp_init, factr, pgtol, iprint)
