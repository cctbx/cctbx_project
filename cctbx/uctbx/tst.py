# $Id$

import math
import pickle
from cctbx_boost.arraytbx import flex
from cctbx_boost import uctbx

print uctbx.__version__
print uctbx.UnitCell.__converters__

print uctbx.UnitCell((12,13,14,91,92,93))
print uctbx.UnitCell([12,13,14,91,92,93])
u = uctbx.UnitCell([12,13,14,91,92,93])
print u.Q([0,0,0])

print uctbx.UnitCell(())
print uctbx.UnitCell((1,))
print uctbx.UnitCell((1,2))
print uctbx.UnitCell((1,2,3))
print uctbx.UnitCell((1,2,3,81))
print uctbx.UnitCell((1,2,3,81,82))
print uctbx.UnitCell((1,2,3,81,82,83))
print uctbx.UnitCell((2,-1,0,-1,2,0,0,0,2))
try:
  u = uctbx.UnitCell((0,0,0,0,0,0))
except RuntimeError, e:
  print e.args[0]
else:
  raise RuntimeError, 'exception expected'
try:
  u = uctbx.UnitCell((0,0,0,0,0,0,0,0,0))
except RuntimeError, e:
  print e.args[0]
else:
  raise RuntimeError, 'exception expected'
u = uctbx.UnitCell((34.698, 71.491, 54.74, 90.0, 106.549, 90.0))
print u
print u.__reduce__()
u = uctbx.UnitCell((30,40,50, 90,90,90))
print u
print u.getVolume()
u = uctbx.UnitCell((30,30,50, 90,90,120))
print u
print u.Q((0,0,0))
print u.two_stol((0,0,0))
print u.d((0,0,0))
print u.Q((3,4,5))
print u.two_stol((3,4,5))**2
print (1./u.d((3,4,5)))**2
print u.getParameters()
print u.getParameters(0)
print u.getParameters(1)
v = uctbx.UnitCell(u.getParameters(1))
print v
Xc = (7,8,9)
print Xc
Xf = u.fractionalize([7,8,9])
print Xf
print u.orthogonalize(Xf)

w = u.ChangeBasis((1,0,0, 0,1,0, 0,0,1))
print 'w=', w
w = u.ChangeBasis((1,0,0, 0,1,0, 0,0,1), 1)
print 'w=', w
w = u.ChangeBasis((0,0,12, 12,0,0, 0,12,0), 12)
print 'w=', w
w = u.ChangeBasis((1,1,0, -1,1,0, 0,0,1), 1)
print 'w=', w
v = w.ChangeBasis((1,-1,0, 1,1,0, 0,0,2), 2)
print 'v=', v
print 'u=', u

for tup in ((1,2), (1,2,3,4)):
  try:
    print u.Q(tup)
  except ValueError, e:
    print e.args[0]
  else:
    raise RuntimeError, 'exception expected'

print u.MaxMillerIndices(3)
print u.MaxMillerIndices(2)
print u.MaxMillerIndices(1)

print u.Q((0, 0, 0))
Q = u.Q((1, 2, 3))
assert Q == u.Q((1, 2, 3))
two_stol = u.two_stol((1, 2, 3))
assert two_stol == u.two_stol((1, 2, 3))
d = u.d((1, 2, 3))
assert d == u.d((1, 2, 3))
print 'Q=', Q, 'two_stol=', math.sqrt(Q), 'd=', 1./math.sqrt(Q)
print 'Q=', Q, 'two_stol=', two_stol    , 'd=', d

print 'G'
print u.getMetricalMatrix()
print 'reciprocal G'
print u.getMetricalMatrix(1)
G = u.getMetricalMatrix()
v = uctbx.UnitCell(G)
print 'u=', u
print 'v=', v

print 'FractionalizationMatrix'
print u.getFractionalizationMatrix()
print 'OrthogonalizationMatrix'
print u.getOrthogonalizationMatrix()

print u.getLongestVector2()

pstr = pickle.dumps(u)
print pstr
up = pickle.loads(pstr)
print u
print up
print u.isEqual(up)
print up.isEqual(u)
o = uctbx.UnitCell((30.1,30,50,90,90,120))
print o
print u.isEqual(o)
print o.isEqual(u)
print u.isEqual(o, 0.1)
print o.isEqual(u, 0.1)

miller_indices = flex.miller_Index(((1,2,3), (4,5,6)))
print tuple(u.Q(miller_indices))
