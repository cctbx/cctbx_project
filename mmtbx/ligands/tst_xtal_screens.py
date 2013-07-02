# -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division

def exercise () :
  from mmtbx.ligands import xtal_screens
  s = xtal_screens.server()

  condition = s.get_condition("Crystal Screen", "4")
  assert condition.pH() == 8.5
  assert condition.ligands() == ["TRS", "CL", "NH4", "SO4"]

  condition = s.get_condition("Crystal Screen II", "E6")
  assert condition.pH() == 7.0
  assert condition.ligands() == ["IMD"]

  condition = s.get_condition("Crystal Screen HT", "A9")
  assert condition.pH() == 5.6
  assert condition.ligands() == ["NH4", "ACT", "NA", "FLC", "PEG"]

  condition = s.get_condition("Index", "5")
  assert condition.pH() == 7.5
  assert condition.ligands() == ["EPE", "NH4", "SO4"]

  condition = s.get_condition("Index HT", "F12")
  assert condition.pH() == 7.5
  assert condition.ligands() == ["NA", "CL", "EPE", "PEG"]

  condition = s.get_condition("PEG/Ion", "48")
  assert condition.pH() == None
  assert condition.ligands() == ["MLT", "LMR", "PEG"]

  condition = s.get_condition("PEG/Ion 2", "49")
  assert condition == None

  condition = s.get_condition("PEG/Ion HT", "G12")
  assert condition.pH() == None
  assert condition.ligands() == ["CIT", "BTB", "TME", "PEG"]

  print "OK"

if (__name__ == "__main__") :
  exercise()
