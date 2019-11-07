from __future__ import absolute_import, division, print_function

class SmallUnitCellVolume(Exception): pass

# cutoff=25.0 Angstrom-cubed is the smallest conceivable unit cell for any
# crystal.  Cells with smaller volume are assumed to have two parallel
# basis vectors and are rejected.  Cutoff value of 100.0 can be used
# for macromolecular work.

def unit_cell_too_small(uc,cutoff=25.):
    abc = uc.parameters()
    if abc[0]*abc[1]*abc[2]/cutoff > uc.volume() or uc.volume() < cutoff:
      raise SmallUnitCellVolume
