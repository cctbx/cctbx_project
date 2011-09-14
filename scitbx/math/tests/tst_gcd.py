import scitbx.math
import boost.rational
import time
import sys, os

def compare_with_boost_rational_gcd(label):
  other_gcd = getattr(scitbx.math, label, None)
  if (other_gcd is None): return
  print "compare_with_boost_rational_gcd(%s)" % label
  samples = range(-100,101) \
          + range(-100000-10,-100000+11) \
          + range( 100000-10, 100000+11)
  for a in samples:
    for b in samples:
      r = boost.rational.gcd(a, b)
      o = other_gcd(a=a, b=b)
      if (o != r):
        raise RuntimeError(str((a, b, r, o)))

def run(args):
  if (len(args) == 0):
    n = 1000
  else:
    assert len(args) == 1
    n = int(args[0])
  #
  labels = [
    "gcd_int_boost",
    "gcd_int_simple",
    "gcd_int32_asm",
    "gcd_long_boost",
    "gcd_long_simple",
    "gcd_unsigned_long_binary",
    "gcd_long_binary",
    "gcd_int64_asm"]
  #
  for label in labels:
    compare_with_boost_rational_gcd(label=label)
  #
  impls = [getattr(scitbx.math, "time_%s" % label, None)
    for label in labels]
  for label,impl in zip(labels,impls):
    if (impl is None): continue
    w0 = time.time()
    us0 = sum(os.times()[:2])
    result = impl(n)
    print "%-24s %d w=%.2f u+s=%.2f" % (
      label,
      result,
      time.time()-w0,
      sum(os.times()[:2])-us0)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
