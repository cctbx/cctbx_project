from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
import sys

codewords = [
  "ABSC", "ACTA", "AFIX", "AFLS", "ANIS", "BASF", "BIND", "BLOC", "BOND",
  "BUMP", "CELL", "CGLS", "CHIV", "CONF", "CONN", "DAMP", "DANG", "DEFS",
  "DELU", "DFIX", "DISP", "EADP", "EGEN", "EQIV", "ESEL", "EXTI", "EXYZ",
  "FACE", "FEND", "FLAT", "FMAP", "FRAG", "FREE", "FVAR", "GRID", "HFIX",
  "HKLF", "HOPE", "HTAB", "ILSF", "INIT", "ISOR", "L.S.", "LATT", "LAUE",
  "LIST", "MERG", "MOLE", "MORE", "MOVE", "MPLA", "NCSY", "OMIT", "PART",
  "PATT", "PHAN", "PHAS", "PLAN", "PSEE", "RESI", "RTAB", "SADI", "SAME",
  "SFAC", "SHEL", "SIMU", "SIZE", "SLIM", "SPAG", "SPEC", "SPIN", "STIR",
  "SUMP", "SWAT", "SYMM", "TEMP", "TEXP", "TIME", "TITL", "TREF", "TWIN",
  "UNDO", "UNIT", "VECT", "VLEN", "VOID", "WGHT", "WPDB", "XHAB", "ZERR"]

class read_shelx_latt(object):

  def __init__(self, line):
    flds = line.split()
    assert flds[0] == "LATT"
    n = int(flds[1])
    assert 1 <= abs(n) <= 7
    self.centric = (n > 0)
    self.z = "*PIRFABC"[abs(n)]

def extract_from(file_name=None, file=None, max_characters=100000):
  assert [file_name, file].count(None) == 1
  if (file is None):
    file = open(file_name)
  unit_cell = None
  space_group = None
  codew = dict([(w+" ",0) for w in codewords])
  n_codew = 0
  n_characters = 0
  for line in file:
    if (max_characters != 0):
      n_characters += len(line)
      if (n_characters > max_characters): break
    l = line.rstrip().split("!")[0]
    if (l.startswith("CELL ")):
      assert unit_cell is None
      flds = l.split()
      assert len(flds) == 8
      unit_cell = uctbx.unit_cell(" ".join(flds[2:]))
    elif (l.startswith("LATT ")):
      assert space_group is None
      latt = read_shelx_latt(l)
      space_group = sgtbx.space_group()
      if (latt.centric):
        space_group.expand_inv(sgtbx.tr_vec((0,0,0)))
      space_group.expand_conventional_centring_type(latt.z)
    elif (l.startswith("SYMM ")):
      assert space_group is not None
      s = sgtbx.rt_mx(l[5:])
      space_group.expand_smx(s)
    else:
      w = l[:5]
      c = codew.get(w, -1)
      if (c == 0):
        codew[w] = 1
        n_codew += 1
  assert unit_cell is not None
  if (space_group is None and n_codew >= 3):
    space_group = sgtbx.space_group()
  return crystal.symmetry(
    unit_cell=unit_cell,
    space_group=space_group)
