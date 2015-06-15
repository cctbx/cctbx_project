from __future__ import division

def distance2(a,b):
  d2 = 0
  for i in range(3):
    d2 += (a.xyz[i]-b.xyz[i])**2
  return d2

