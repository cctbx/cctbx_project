from libtbx import adopt_init_args
import sys

class flags(object):

  def __init__(self,
        bond=None,
        nonbonded=None,
        angle=None,
        dihedral=None,
        reference_dihedral=None,
        chirality=None,
        planarity=None,
        bond_similarity=None,
        generic_restraints=None,
        default=False):
    if (bond is None): bond = default
    if (nonbonded is None): nonbonded = default
    if (angle is None): angle = default
    if (dihedral is None): dihedral = default
    if (reference_dihedral is None): reference_dihedral = default
    if (chirality is None): chirality = default
    if (planarity is None): planarity = default
    if (bond_similarity is None): bond_similarity = default
    if (generic_restraints is None) : generic_restraints = default
    adopt_init_args(self, locals())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "geometry_restraints.manager.flags:"
    print >> f, "  bond:", self.bond
    print >> f, "  nonbonded:", self.nonbonded
    print >> f, "  angle:", self.angle
    print >> f, "  dihedral:", self.dihedral
    print >> f, "  reference dihedral:", self.reference_dihedral
    print >> f, "  chirality:", self.chirality
    print >> f, "  planarity:", self.planarity
    print >> f, "  bond similarity:", self.bond_similarity
    print >> f, "  other (generic):", self.generic_restraints
