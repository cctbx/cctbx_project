from boost_python_hybrid_times_ext import run_c_plus_plus
import math
import time
import sys

if (not hasattr(sys, "gettickeraccumulation")):
  print "***************************************************"
  print "WARNING: sys.gettickeraccumulation() not available."
  print "***************************************************"
  print
  def gettickeraccumulation():
    return math.ceil(time.time()*1.e6)
  sys.gettickeraccumulation = gettickeraccumulation

def factorial(n):
  fact = 1
  for i in range(2,n+1):
    fact *= i
  return fact

def power(x, n):
  pow = x
  for i in range(1,n):
    pow *= x
  return pow

def sin(x, n_terms):
  result = x
  sign = -1
  pow = 3
  for i in range(1,n_terms):
    result += power(x,pow)/(sign*factorial(pow))
    sign *= -1
    pow += 2
  return result

def run_python(n, n_terms):
  result = 0
  for j in range(n):
    for i in range(180):
      result += sin(i * 3.14159265359/180, n_terms)
  return result

class hybrid:

  def __init__(self, n_python):
    self.n_python = n_python

  def __call__(self, n, n_terms):
    result = run_python(self.n_python, n_terms)
    result += run_c_plus_plus(n-self.n_python, n_terms)
    return result

class time_per_python_tick:

  def __init__(self, worker, n, n_terms):
    time_0 = time.time()
    ticks_0 = sys.gettickeraccumulation()
    self.result = worker(n, n_terms)
    self.ticks_diff = sys.gettickeraccumulation() - ticks_0
    self.time_diff = time.time() - time_0

  def report(self, label):
    print "%-6s %6.3f %10d %10.3f" % (
      label, self.time_diff, self.ticks_diff,
      self.time_diff/self.ticks_diff*1.e6)
    return self

def end_points(n_terms):
  print "         time      ticks  time/tick"
  n = 1
  while (n <= 128):
    print "n:", n
    py = time_per_python_tick(run_python, n, n_terms).report("Python")
    cpp = time_per_python_tick(run_c_plus_plus, n, n_terms).report("C++")
    print "ratio %6.3f" % (py.time_diff / cpp.time_diff)
    assert abs(py.result - cpp.result) < 1.e-8
    print
    n *= 2

def plot(n_terms):
  for n,m in [(100,100), (1000, 100), (10000, 100), (100000, 100)]:
    print "                 time      ticks  time/tick"
    reference_result = None
    for n_python in xrange(0, m+1, m/10):
      hy = time_per_python_tick(hybrid(n_python), n, n_terms).report(
        "%6.2f%% Python" % (100.*n_python/n))
      if (reference_result is None):
        reference_result = hy.result
      else:
        assert abs(hy.result - reference_result) < 1.e-10*n
    print

def run():
  n_terms = 6 # larger values lead to integer overflows in factorial()
  end_points(n_terms)
  plot(n_terms)

if (__name__ == "__main__"):
  run()
