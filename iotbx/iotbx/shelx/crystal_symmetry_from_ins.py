"Extracts crystal symmetry from CNS input file."

from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
import sys

class read_shelx_latt:

  def __init__(self, line):
    flds = line.split()
    assert flds[0] == "LATT"
    n = int(flds[1])
    assert 1 <= abs(n) <= 7
    self.centric = (n > 0)
    self.z = "*PIRFABC"[abs(n)]

def extract_from(file):
  unit_cell = None
  space_group = None
  for line in file:
    l = line.strip()
    if (l.startswith("CELL ")):
      assert unit_cell is None
      flds = l.split()
      assert len(flds) == 8
      unit_cell = uctbx.unit_cell(" ".join(flds[2:]))
    if (l.startswith("LATT ")):
      assert space_group is None
      latt = read_shelx_latt(l)
      space_group = sgtbx.space_group()
      if (latt.centric):
        space_group.expand_inv(sgtbx.tr_vec((0,0,0)))
      space_group.expand_conventional_centring_type(latt.z)
    if (l.startswith("SYMM ")):
      assert space_group is not None
      s = sgtbx.rt_mx(l[5:])
      space_group.expand_smx(s)
  assert unit_cell is not None
  assert space_group is not None
  return crystal.symmetry(
    unit_cell=unit_cell,
    space_group=space_group)
