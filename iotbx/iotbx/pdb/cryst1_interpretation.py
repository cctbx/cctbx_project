from iotbx.pdb import parser as pdb_parser
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx

def is90(angle, eps=0.01):
  return abs(angle-90) < eps

def is120(angle, eps=0.01):
  return abs(angle-120) < eps

def equiv(r,s,t, eps=0.01):
  m = (r+s+t)/3.
  return abs(r-m) < eps and abs(s-m) < eps and abs(t-m) < eps

rhombohedral = {
  "R3": "R 3",
  "H3": "R 3",
  "R32": "R 3 2",
  "H32": "R 3 2",
}

short_mono = (
  "P2",
  "P21",
  "C2",
  "A2",
  "B2",
  "I2",
)

special = {
  "A1": "Hall:  P 1 (-x,-1/2*y+1/2*z,1/2*y+1/2*z)",
  "C1211": "Hall:  C 2y (x+1/4,y+1/4,z)",
  "C21": "Hall:  C 2y (x+1/4,y+1/4,z)",
  "I1211": "Hall:  C 2y (x+1/4,y+1/4,-x+z-1/4)",
  "I21": "Hall:  C 2y (x+1/4,y+1/4,-x+z-1/4)",
  "P21212A": "Hall:  P 2 2ab (x+1/4,y+1/4,z)",
  "F422": "Hall:  I 4 2 (1/2*x+1/2*y,-1/2*x+1/2*y,z)",
  "C4212": "Hall:  P 4 2 (1/2*x+1/2*y-1/4,-1/2*x+1/2*y-1/4,z)",
}

_all = {}
for sym in rhombohedral.keys(): _all[sym] = rhombohedral
for sym in short_mono: _all[sym] = short_mono
for sym in special.keys(): _all[sym] = special

class categorize:

  def __init__(self, symbol):
    self.symbol = None
    try:
      self.symbol = symbol.strip()
      self.category = _all[self.symbol.replace(" ","").upper()]
    except:
      self.category = None

  def get_category(self):
    if (self.category == rhombohedral): return "rhombohedral"
    if (self.category == short_mono): return "short_mono"
    if (self.category == special): return "special"
    return None

  def space_group_info(self, unit_cell=None):
    if (self.category is None):
      try: return sgtbx.space_group_info(self.symbol)
      except RuntimeError: return None
    if (isinstance(unit_cell, uctbx.ext.unit_cell)):
      unit_cell = unit_cell.parameters()
    if (self.category == rhombohedral):
      if (unit_cell is None): return None
      (a, b, c, alpha, beta, gamma) = unit_cell
      if (abs(a - b) <= 0.01 and is90(alpha) and is90(beta) and is120(gamma)):
        basis_symbol = "H"
      elif (equiv(a,b,c) and equiv(alpha,beta,gamma)):
        basis_symbol = "R"
      else:
        return None
      return sgtbx.space_group_info(
        rhombohedral[self.symbol] + ":" + basis_symbol)
    if (self.category == short_mono):
      if (unit_cell is None): return None
      Z, T = self.symbol[0], self.symbol[1:]
      (a, b, c, alpha, beta, gamma) = unit_cell
      if (is90(alpha) and is90(gamma)):
        if (Z == "B"): return None
        return sgtbx.space_group_info(Z + " 1 " + T + " 1")
      if (is90(alpha) and is90(beta)):
        if (Z == "C"): return None
        return sgtbx.space_group_info(Z + " 1 1 " + T)
    if (self.category == special):
      return sgtbx.space_group_info(special[self.symbol])
    raise RuntimeError, "Programming error (should be unreachable)."

def crystal_symmetry(cryst1_record, line_number=None):
  if (isinstance(cryst1_record, str)):
    cryst1_record = pdb_parser.pdb_record(
      raw_record=cryst1_record,
      line_number=line_number)
  space_group_info = categorize(cryst1_record.sGroup).space_group_info(
    unit_cell=cryst1_record.ucparams)
  return crystal.symmetry(
    unit_cell=cryst1_record.ucparams,
    space_group_info=space_group_info)
