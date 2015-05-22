from __future__ import division
from libtbx import adopt_init_args
import sys

class flags(object):

  def __init__(self,
        rotamer=None,
        reference=None,
        default=False):
    if (rotamer is None): rotamer = default
    if (reference is None): reference = default
    adopt_init_args(self, locals())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "generic_restraints_manager.flags:"
    print >> f, "  rotamer:", self.rotamer
    print >> f, "  reference:", self.reference
