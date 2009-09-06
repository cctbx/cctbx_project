import scitbx.math
import boost.rational
import time
import sys, os

def compare_with_boost_rational_gcd(other_gcd):
  if (other_gcd is None): return
  for a in xrange(-100,101):
    for b in xrange(-100,101):
      r = boost.rational.gcd(a, b)
      o = other_gcd(a=a, b=b)
      if (o != r):
        raise RuntimeError(str((a, b, r, o)))

def run(args):
  if (len(args) == 0):
    n = 4400
  else:
    assert len(args) == 1
    n = int(args[0])
  #
  for label in ["gcd_int_simple", "gcd_int_asm"]:
    compare_with_boost_rational_gcd(
      other_gcd=getattr(scitbx.math, label, None))
  #
  labels = ["gcd_int_boost_math", "gcd_int_simple", "gcd_int_asm"]
  impls = [getattr(scitbx.math, "time_%s" % label, None)
    for label in labels]
  for impl in impls:
    if (impl is None): continue
    impl(10)
  for label,impl in zip(labels,impls):
    if (impl is None): continue
    w0 = time.time()
    us0 = sum(os.times()[:2])
    result = impl(n)
    print "%-18s %d w=%.2f u+s=%.2f" % (
      label,
      result,
      time.time()-w0,
      sum(os.times()[:2])-us0)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
