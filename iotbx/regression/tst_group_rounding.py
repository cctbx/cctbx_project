from __future__ import absolute_import, division, print_function
from iotbx.pdb.hierarchy import group_rounding


def tst_1():
  """ Test cases where nothing should happen"""
  assert group_rounding([0.222222], 2) == [0.222222]
  assert group_rounding([5.222], 2) == [5.222]
  # sum less than 1
  assert group_rounding([0.222222, 0.22222], 2) == [0.222222, 0.22222]
  assert group_rounding([0.222222, 0.22222, 0.2, 0.2], 2) == [0.222222, 0.22222, 0.2, 0.2]

  # sum greater than 1
  assert group_rounding([0.5222222, 0.522222], 2) == [0.5222222, 0.522222]
  assert group_rounding([0.5222222, 0.522222, 0.2, 0.2], 2) == [0.5222222, 0.522222, 0.2, 0.2]

def tst_2():
  """ Test cases where rounding should occur"""
  # straightforward cases, no special procedure would have needed
  assert group_rounding([0.4555555, 0.5444445],3) == [0.456, 0.544]
  assert group_rounding([0.4555555, 0.5444445],2) == [0.46, 0.54]
  assert group_rounding([0.250001, 0.249999, 0.5],3) == [0.25, 0.25, 0.5]
  assert group_rounding([0.250001, 0.249999, 0.5],2) == [0.25, 0.25, 0.5]
  assert group_rounding([0.250001, 0.249999, 0.1, 0.1, 0.1, 0.1, 0.1],2) == [0.25, 0.25, 0.1, 0.1, 0.1, 0.1, 0.1]

def tst_3():
  # interesting cases, where standard rounding would exceed sum of 1
  assert group_rounding([0.545, 0.455],2) == [0.55, 0.45]
  assert group_rounding([0.125, 0.125, 0.75], 2) == [0.13, 0.12, 0.75]
  assert group_rounding([0.75, 0.125, 0.125], 2) == [0.75, 0.13, 0.12]
  assert group_rounding([0.125, 0.125, 0.5, 0.25], 2) == [0.13, 0.12, 0.5, 0.25]

  assert group_rounding([0.5005, 0.4995],3) == [0.5, 0.5]
  assert group_rounding([0.7555, 0.2445],3) == [0.756, 0.244]
  assert group_rounding([0.5555, 0.2445, 0.2],3) == [0.556, 0.244, 0.2]

if (__name__ == "__main__"):
  tst_1()
  tst_2()
  tst_3()
