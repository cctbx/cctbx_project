from __future__ import division
from libtbx import adopt_init_args
import sys

class flags(object):

  def __init__(self,
        ramachandran=None,
        rotamer=None,
        reference=None,
        den=None,
        default=False):
    if (ramachandran is None): ramachandran = default
    if (rotamer is None): rotamer = default
    if (reference is None): reference = default
    if (den is None): den = False #always off by default
    adopt_init_args(self, locals())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "generic_restraints_manager.flags:"
    print >> f, "  ramachandran:", self.ramachandran
    print >> f, "  rotamer:", self.rotamer
    print >> f, "  reference:", self.reference
    print >> f, "  den:", self.den
