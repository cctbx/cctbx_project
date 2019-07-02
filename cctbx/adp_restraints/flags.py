from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args
import sys

class flags(object):

  def __init__(self,
        adp_similarity=None,
        rigid_bond=None,
        isotropic_adp=None,
        default=False):
    if (adp_similarity is None): adp_similarity = default
    if (rigid_bond is None): rigid_bond = default
    if (isotropic_adp is None): isotropic_adp = default
    adopt_init_args(self, locals())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print("adp_restraints.manager.flags:", file=f)
    print("  adp_similarity:", self.adp_similarity, file=f)
    print("  rigid_bond:", self.rigid_bond, file=f)
    print("  isotropic_adp:", self.isotropic_adp, file=f)
