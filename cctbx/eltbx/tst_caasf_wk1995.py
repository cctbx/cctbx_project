# $Id$

from caasf_wk1995 import *
sf = CAASF_WK1995("SI")
print sf.Table()
print sf.Label()
print sf.n_ab()
for i in xrange(sf.n_ab()):
  print sf.a(i), sf.b(i)
print sf.c()
start_stol = 0.
end_stol = 1.
step_stol = 0.1
i = 0
while 1:
  stol = start_stol + i * step_stol
  if (stol > end_stol + 1.e-4): break
  print "%5.3f %.6g" % (stol, sf(stol * stol))
  i = i + 1
