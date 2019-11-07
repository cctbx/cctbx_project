from __future__ import absolute_import, division, print_function
from six.moves import range
def level(space_group_info):
  """\
Determine the level of confidence in a given space group symbol. Often
data are scaled and merged in a space group that is different from the
crystal space group. For example, if the crystal space group is P212121
the data are sometimes processed in space group P222. Therefore this
function returns 0 for these two symbols. For some space groups this
kind of ambiguity is not possible, e.g. space group P1. The
corresponding return value is 2. Some symbols could be ambiguous,
but are typically not used as "default" symbols by data processing
programs, e.g. P3212. It is therefore more likely that these symbols
are chosen on purpose. The corresponding return value is 1.
This function returns None for all non-chiral space groups.
The return value for chiral space groups with unusual symbols is -1.
"""
  if (not space_group_info.group().is_chiral()):
    return None
  # Based on a table in the HKL manual.
  symbol_groups = [line.split() for line in [
    "P23 P213",
    "P432 P4132 P4232 P4332",
    "I23 I213",
    "I432 I4132",
    "F23",
    "F432 F4132",
    "R3:H",
    "R32:H",
    "P3 P31 P32",
    "P312 P3112 P3212",
    "P321 P3121 P3221",
    "P6 P61 P65 P62 P64 P63",
    "P622 P6122 P6522 P6222 P6422 P6322",
    "P4 P41 P42 P43",
    "P422 P41212 P4212 P4122 P4322 P4222 P42212 P43212",
    "I4 I41",
    "I422 I4122",
    "P222 P212121 P2221 P21212",
    "C2221 C222",
    "I222 I212121",
    "F222",
    "P121 P1211",
    "C121",
    "P1"
  ]]
  symbol = space_group_info.type().lookup_symbol().replace(" ", "")
  for symbol_group in symbol_groups:
    try: i = symbol_group.index(symbol)
    except KeyboardInterrupt: raise
    except Exception: continue
    if (len(symbol_group) == 1): return 2
    if (i < 2): return 0
    return 1
  return -1

def _test():
  from cctbx import sgtbx
  assert level(sgtbx.space_group_info("P -1")) == None
  assert level(sgtbx.space_group_info("P 1")) == 2
  assert level(sgtbx.space_group_info("P 2 2 21")) == 1
  assert level(sgtbx.space_group_info("P 2")) == 0
  assert level(sgtbx.space_group_info("Hall: C 1")) == -1
  n = 0
  for i in range(1,231):
    c = level(sgtbx.space_group_info(i))
    if (c is not None):
      assert c  in (0,1,2)
      n += 1
  assert n == 65
  print("OK")

if (__name__ == "__main__"):
  _test()
