#! /usr/bin/env python

"""
Extracts crystal symmetry from CNS input file.

For example, the input

{===>} sg="P3(2)21";
{===>} a=76.4;
{===>} b=76.4;
{===>} c=180.94;
{===>} alpha=90;
{===>} beta=90;
{===>} gamma=120;

is reformatted as:

crystal.symmetry(
  unit_cell=(76.4, 76.4, 180.94, 90, 90, 120),
  space_group_symbol="P3(2)21")
"""

import sys

def run():
  assert len(sys.argv) in (1,2)
  if (len(sys.argv) == 1):
    file = sys.stdin
  else:
    file = open(sys.argv[1], "r")
  unit_cell = [None for i in xrange(6)]
  space_group_symbol = None
  for line in file:
    line = line.strip()
    if (line.startswith('{===>} sg="')):
      assert line[-2:] == '";'
      assert space_group_symbol == None, "Duplicate space group symbol."
      space_group_symbol = line[11:-2]
    else:
      i = -1
      for label in ("a", "b", "c", "alpha", "beta", "gamma"):
        i += 1
        if (line.startswith('{===>} %s=' % label)):
          assert line[-1] == ';'
          assert unit_cell[i] == None, "Duplicate unit cell parameter %s." % label
          unit_cell[i] = float(line[8+len(label):-1])
  assert space_group_symbol != None
  assert unit_cell.count(None) == 0
  print """\
crystal.symmetry(
  unit_cell=(%.6g, %.6g, %.6g, %.6g, %.6g, %.6g),
  space_group_symbol="%s")
""" % tuple(unit_cell + [space_group_symbol])

if (__name__ == "__main__"):
  run()
