from __future__ import absolute_import, division, print_function
# Example by Kevin Cowtan, 2001. Public Domain.
# Converted to Python by Ralf W. Grosse-Kunstleve.
# This little program reads the CCP4 symmetry library file, and interprets
# the symops. It then returns the Hall code for the spacegroup. Used to
# create a new library file to allow compatible spacegroup naming between
# CCP4 and cctbx.
#
# usage: convert_ccp4_symop_lib < symop.lib

from cctbx import sgtbx
import sys
from six.moves import range

def run():
  while 1:
    line = sys.stdin.readline()[:-1]
    flds = line.split(None, 2)
    if (len(flds) == 0): break
    nspgrp = int(flds[0]) # read spacegroup number
    nsym = int(flds[1]) # read nsym
    print(nspgrp, nsym, flds[2]) # print it all
    group = sgtbx.space_group() # now interpret the symops
    for i in range(nsym):
      line = sys.stdin.readline()[:-1] # get the i'th symop
      # print line
      group.expand_smx(sgtbx.rt_mx(line)) # and interpret
    info = sgtbx.space_group_info(group=group)
    print(info.type().hall_symbol()) # now produce the sg symbol
    print(info)

if (__name__ == "__main__"):
  run()
