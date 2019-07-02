# -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import absolute_import, division, print_function

def exercise():
  from mmtbx.ligands import xtal_screens
  s = xtal_screens.server()

  try:
    condition = s.get_condition("chrysler", "1")
    assert False
  except RuntimeError:
    pass

  condition = s.get_condition("crystal_screen", "4")
  assert condition.pH() == 8.5
  assert condition.ligands() == ["TRS", "CL", "NH4", "SO4"]

  condition = s.get_condition("crystal_screen", "E6")
  assert condition.pH() == 7.0
  assert condition.ligands() == ["IMD"]

  condition = s.get_condition("crystal_screen", "A9")
  assert condition.pH() == 5.6
  assert condition.ligands() == ["NH4", "ACT", "NA", "FLC", "PEG"]

  condition = s.get_condition("index", "5")
  assert condition.pH() == 7.5
  assert condition.ligands() == ["EPE", "NH4", "SO4"]

  condition = s.get_condition("index", "F12")
  assert condition.pH() == 7.5
  assert condition.ligands() == ["NA", "CL", "EPE", "PEG"]

  condition = s.get_condition("peg_ion", "48")
  assert condition.pH() == None
  assert condition.ligands() == ["NH4", "FLC", "PEG"]

  try:
    condition = s.get_condition("peg_ion", "49")
    assert False
  except RuntimeError:
    pass

  condition = s.get_condition("peg_ion", "G12")
  assert condition.pH() == None
  assert condition.ligands() == ["CIT", "BTB", "TME", "PEG"]

  print("OK")

if (__name__ == "__main__"):
  exercise()
