"""Extracts crystal symmetry from CNS input file.

For example, the input

{===>} sg="P3(2)21";
{===>} a=76.4;
{===>} b=76.4;
{===>} c=180.94;
{===>} alpha=90;
{===>} beta=90;
{===>} gamma=120;

is returned as:

crystal.symmetry(
  unit_cell=(76.4, 76.4, 180.94, 90, 90, 120),
  space_group_symbol="P3(2)21")
"""

from cctbx import crystal

def extract_from(file_name=None, file=None, max_characters=1000000):
  assert [file_name, file].count(None) == 1
  if (file is None):
    file = open(file_name)
  unit_cell = [None for i in xrange(6)]
  space_group_symbol = None
  n_characters = 0
  for line in file:
    if (max_characters != 0):
      n_characters += len(line)
      if (n_characters > max_characters): break
    line = line.strip()
    if (line.startswith('{===>} sg="')):
      assert line[-2:] == '";'
      assert space_group_symbol is None, "Duplicate space group symbol."
      space_group_symbol = line[11:-2]
    else:
      i = -1
      for label in ("a", "b", "c", "alpha", "beta", "gamma"):
        i += 1
        if (line.startswith('{===>} %s=' % label)):
          assert line[-1] == ';'
          assert unit_cell[i] is None, "Duplicate unit cell parameter %s." % label
          unit_cell[i] = float(line[8+len(label):-1])
  assert unit_cell.count(None) == 0
  assert space_group_symbol is not None
  return crystal.symmetry(
    unit_cell=unit_cell,
    space_group_symbol=space_group_symbol)
