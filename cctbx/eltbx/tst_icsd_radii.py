# $Id$

from icsd_radii import *
for lbl in ("C", "O", "N", "Si", "Si4+"):
  entry = ICSD_Radius(lbl)
  print entry.Label(), entry.Radius()
