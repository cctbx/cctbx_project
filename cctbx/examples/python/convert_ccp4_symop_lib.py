# Example by Kevin Cowtan, 2001. Public Domain.
# Converted to Python by Ralf W. Grosse-Kunstleve.
# This little program reads the CCP4 symmetry library file, and interprets
# the symops. It then returns the Hall code for the spacegroup. Used to
# create a new library file to allow compatible spacegroup naming between
# CCP4 and cctbx.
#
# usage: convert_ccp4_symop_lib < symop.lib

import sys
import string
from cctbx_boost import sgtbx

while 1:
  line = sys.stdin.readline()[:-1]
  flds = string.split(line, None, 2)
  if (len(flds) == 0): break
  nspgrp = string.atoi(flds[0]) # read spacegroup number
  nsym = string.atoi(flds[1]) # read nsym
  print nspgrp, nsym, flds[2] # print it all
  SgOps = sgtbx.SpaceGroup() # now interpret the symops
  for i in xrange(nsym):
    line = sys.stdin.readline()[:-1] # get the i'th symop
    # print line
    SgOps.expandSMx( sgtbx.RTMx(line) ) # and interpret
  SgInfo = SgOps.Info()
  print SgInfo.BuildHallSymbol() # now produce the sg symbol
  print SgInfo.BuildLookupSymbol()
