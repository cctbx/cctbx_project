# $Id$

from tinypse import *
for key in ("Si", 13):
  entry = TinyPSE(key)
  print entry.Z()
  print entry.Symbol()
  print entry.Name()
  print entry.Weight()
