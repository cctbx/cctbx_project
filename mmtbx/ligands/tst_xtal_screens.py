
from __future__ import division

def exercise () :
  from mmtbx.ligands import xtal_screens
  s = xtal_screens.server()
  condition = s.get_condition("crystal_screen", "A9")
  assert (condition.pH() == 5.6)
  assert (condition.ligands() == ['NH4', 'ACT', 'NA', 'FLC', 'PEG'])
  print "OK"

if (__name__ == "__main__") :
  exercise()
