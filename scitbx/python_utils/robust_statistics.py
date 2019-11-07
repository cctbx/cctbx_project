from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex

def percentile(x, percent):
  ## see http://cnx.rice.edu/content/m10805/latest/
  assert percent >=0.0
  assert percent <=1.0
  n = x.size()
  order = flex.sort_permutation(x)
  np = float(n+1)*percent
  ir = int(np)
  fr = np - int(np)
  tmp1 = x[ order[ir-1] ]
  tmp2 = x[ order[ir] ]
  result = tmp1 + fr*(tmp2-tmp1)
  return result

def median(x):
  result=percentile(x, 0.5)
  return( result )

def hinges(x):
  h1 = percentile( x, 0.25 )
  h2 = percentile( x, 0.75 )
  return h1, h2

def h_spread(x):
  low, high = hinges( x)
  return ( high - low )

def trimean(x):
  h1, h2 = hinges(x)
  m = median(x)
  result = (h1+h2+2.0*m)/4.0
  return result

def trimmed_mean(x,percent):
  low = percentile(x,percent/2.0)
  high = percentile(x,1.0-percent/2.0)
  tmp = x.select ( ( x > low ).iselection() )
  tmp = tmp.select( ( tmp < high ).iselection() )
  return flex.mean(tmp)
