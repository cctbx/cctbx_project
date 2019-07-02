from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args
import sys

class flags(object):

  def __init__(self,
        bond=None,
        nonbonded=None,
        angle=None,
        dihedral=None,
        reference_coordinate=None,
        reference_dihedral=None,
        ncs_dihedral=None,
        den_restraints=None,
        chirality=None,
        planarity=None,
        parallelity=None,
        bond_similarity=None,
        ramachandran_restraints=None,
        default=False):
    if (bond is None): bond = default
    if (nonbonded is None): nonbonded = default
    if (angle is None): angle = default
    if (dihedral is None): dihedral = default
    if (reference_coordinate is None): reference_coordinate = default
    if (reference_dihedral is None): reference_dihedral = default
    if (ncs_dihedral is None): ncs_dihedral = default
    if den_restraints is None: den_restraints = default
    if (chirality is None): chirality = default
    if (planarity is None): planarity = default
    if (parallelity is None): parallelity = default
    if (bond_similarity is None): bond_similarity = default
    if ramachandran_restraints is None: ramachandran_restraints = default
    adopt_init_args(self, locals())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print("geometry_restraints.manager.flags:", file=f)
    print("  bond:", self.bond, file=f)
    print("  nonbonded:", self.nonbonded, file=f)
    print("  angle:", self.angle, file=f)
    print("  dihedral:", self.dihedral, file=f)
    print("  reference coordinate:", self.reference_coordinate, file=f)
    print("  reference dihedral:", self.reference_dihedral, file=f)
    print("  chirality:", self.chirality, file=f)
    print("  planarity:", self.planarity, file=f)
    print("  parallelity:", self.parallelity, file=f)
    print("  bond similarity:", self.bond_similarity, file=f)
    print("  ramachandran:", self.ramachandran_restraints, file=f)
    print("  DEN:", self.den_restraints, file=f)
