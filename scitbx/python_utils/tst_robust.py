from scitbx.array_family import flex
from scitbx.python_utils import robust_statistics as rs

def tst_median():
  x = flex.double(range(10))
  assert rs.median(x)==4.5
  x = flex.double(range(11))
  assert rs.median(x)==5.0

def tst_hinges():
  n = 12
  x = flex.double(range(n))
  hinges = rs.hinges(x)
  assert hinges[0] == 2.25
  assert hinges[1] == 8.75

def tst_h_spread():
  n = 12
  x = flex.double(range(n))
  assert rs.h_spread(x) == 6.5

def tst_trimean():
  n = 13
  x = flex.double(range(n))
  assert rs.trimean(x) == rs.median(x)

def tst_percentile():
  n = 10
  x = flex.double(range(n))
  assert rs.percentile(x,0.50) == rs.median(x)

def tst_trimmed_mean():
  n = 10
  x = flex.double(range(n))
  assert rs.trimmed_mean(x,0.50)==4.5

def run():
  tst_median()
  tst_hinges()
  tst_trimean()
  tst_h_spread()
  tst_percentile()
  tst_trimmed_mean()
  print "OK"

run()
