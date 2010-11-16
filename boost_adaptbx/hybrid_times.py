from boost_python_hybrid_times_ext import run_c_plus_plus
import sys, os

if (not hasattr(sys, "gettickeraccumulation")):
  print "***************************************************"
  print "WARNING: sys.gettickeraccumulation() not available."
  print "***************************************************"
  print
  def gettickeraccumulation():
    return 1000000
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

class hybrid(object):

  def __init__(self, n_python):
    self.n_python = n_python

  def __call__(self, n, n_terms):
    result = run_python(self.n_python, n_terms)
    result += run_c_plus_plus(n-self.n_python, n_terms)
    return result

def usr_and_sys():
  t = os.times()
  return t[0]+t[1]

class time_per_python_tick(object):

  def __init__(self, worker, n, n_terms):
    time_0 = usr_and_sys()
    ticks_0 = sys.gettickeraccumulation()
    self.result = worker(n, n_terms)
    self.ticks_diff = sys.gettickeraccumulation() - ticks_0
    self.time_diff = usr_and_sys() - time_0

  def report(self, label):
    print "%-6s %6.3f %10d %10.3f" % (
      label, self.time_diff, self.ticks_diff,
      self.time_diff/max(1,self.ticks_diff)*1.e6)
    return self

def forever(n_terms, n=200):
  sys.tracebacklimit = 0
  n_ratios = 0
  sum_py = 0
  sum_cpp = 0
  while True:
    time_0 = usr_and_sys()
    run_python(n, n_terms)
    sum_py += usr_and_sys() - time_0
    time_0 = usr_and_sys()
    run_c_plus_plus(n, n_terms)
    sum_cpp += usr_and_sys() - time_0
    n_ratios += 1
    print "mean Python: %6.3f" % (sum_py/n_ratios)
    print "mean C++:    %6.3f" % (sum_cpp/n_ratios)
    if (sum_cpp != 0):
      print "iteration %d: ratio %6.3f" % (n_ratios, sum_py/sum_cpp)
    else:
      print "iteration %d: ratio infinity" % n_ratios
    print

def end_points(n_terms):
  print "         time      ticks  time/tick"
  n = 1
  while (n <= 128):
    print "n:", n
    py = time_per_python_tick(run_python, n, n_terms).report("Python")
    cpp = time_per_python_tick(run_c_plus_plus, n, n_terms).report("C++")
    if (cpp.time_diff != 0):
      print "ratio %6.3f" % (py.time_diff / cpp.time_diff)
    else:
      print "ratio infinity"
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
  if ("--forever" in sys.argv[1:]):
    forever(n_terms)
  else:
    end_points(n_terms)
    plot(n_terms)

if (__name__ == "__main__"):
  run()
