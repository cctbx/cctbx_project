# LIBTBX_SET_DISPATCHER_NAME cxi.plotcv_parse

import sys
from scitbx.array_family import flex


tile = flex.int()
delx = flex.double()
dely = flex.double()
for line in sys.__stdin__.readlines():
  if line.find("Tile")==0:
    tokens = line.split()
    #print tokens[1], int(tokens[1][:-1]),tokens
    tile.append( int(tokens[1][:-1]) )
    delx.append( float(tokens[7]) )
    dely.append( float(tokens[9][:-1]) )
    if len(tile)==64: break

order = flex.sort_permutation(tile)
for x in xrange(64):
  #print tile[order[x]], delx[order[x]], dely[order[x]]
  pass

for x in xrange(8):
  for y in xrange(8):
    print "%5.2f,"%delx[order[x*8+y]], "%5.2f,"%dely[order[x*8+y]],
  print
print

for x in xrange(8):
  for y in xrange(8):
    print "%5.2f,"%dely[order[x*8+y]], "%5.2f,"%delx[order[x*8+y]],
  print
print

for x in xrange(8):
  for y in xrange(8):
    print "%5.2f,"%-delx[order[x*8+y]], "%5.2f,"%-dely[order[x*8+y]],
  print
print

for x in xrange(8):
  for y in xrange(8):
    print "%5.2f,"%-dely[order[x*8+y]], "%5.2f,"%-delx[order[x*8+y]],
  print
print
