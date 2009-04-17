from scitbx.array_family import flex
from scitbx.lbfgs import have_lbfgs_f, fortran, raw_reference
from libtbx.utils import show_times
from libtbx.test_utils import show_diff
from libtbx import easy_run
import libtbx.load_env
import random
import sys

def exercise(lbfgs_impl, n=100, m=5, iprint=[1, 0]):
  assert n % 2 == 0
  x = flex.double(n)
  g = flex.double(n)
  diag = flex.double(n)
  diagco = 0
  eps = 1.0e-5
  xtol = 1.0e-16
  for j in xrange(0, n, 2):
    x[j] = -1.2
    x[j+1] = 1.
  size_w = n*(2*m+1)+2*m
  w = flex.double(size_w)
  iflag = 0
  icall = 0
  while 1:
    f = 0.
    for j in xrange(0, n, 2):
      t1 = 1.e0 - x[j]
      t2 = 1.e1 * (x[j+1] - x[j] * x[j])
      g[j+1] = 2.e1 * t2
      g[j] = -2.e0 * (x[j] * g[j+1] + t1)
      f = f + t1 * t1 + t2 * t2
    iflag = lbfgs_impl(
      n=n, m=m, x=x, f=f, g=g, diagco=diagco, diag=diag,
      iprint=iprint, eps=eps, xtol=xtol, w=w, iflag=iflag)
    if (iflag <= 0): break
    icall += 1
    # We allow at most 2000 evaluations of f and g
    if (icall > 2000): break

def run_and_compare_implementations(this_script, n, m, iprint):
  outputs = []
  for impl in ["fortran", "raw_reference"]:
    if (impl == "fortran" and not have_lbfgs_f):
      continue
    cmd = 'scitbx.python "%s" %s %d %d %d %d' % (
      this_script, impl, n, m, iprint[0], iprint[1])
    print cmd
    sys.stdout.flush()
    out = easy_run.fully_buffered(command=cmd)
    err = "\n".join(out.stderr_lines)
    if (len(err) != 0):
      print err
      if (err.find("== ERROR SUMMARY: 0 errors from 0 contexts") < 0):
        raise AssertionError(
          "stderr output does not appear to be valgrind output")
    out = "\n".join(out.stdout_lines)
    if (impl == "fortran"):
      out = out.replace("D-", "E-").replace("D+", "E+")
    outputs.append(out)
  if (len(outputs) != 1):
    f, r = outputs
    assert not show_diff(f, r) # may fail for other random seeds,
      # due to rounding errors

def compare_implementations():
  this_script = libtbx.env.under_dist(
    module_name="scitbx", path="lbfgs/tst_lbfgs_f.py")
  assert this_script.find('"') < 0
  rnd = random.Random(x=0)
  for iprint1 in [-1, 0, 1, 2, 3]:
    for iprint2 in [0, 1, 2, 3]:
      n = rnd.choice([100, 14, 4, 2])
      m = rnd.choice([2, 3, 4, 5, 6, 7])
      run_and_compare_implementations(
        this_script=this_script, n=n, m=m, iprint=[iprint1, iprint2])

def run(args):
  timer = show_times(time_start="now")
  if (len(args) == 5):
    exercise(
      lbfgs_impl=eval(args[0]),
      n=int(args[1]),
      m=int(args[2]),
      iprint=[int(args[3]), int(args[4])])
    return
  assert args in [[], ["--once"], ["--endless"]]
  if (not have_lbfgs_f):
    print "Skipping some tests: lbfgs.f not available."
  once = "--once" in args
  endless = "--endless" in args
  if (once or endless):
    while 1:
      if (have_lbfgs_f):
        exercise(lbfgs_impl=fortran)
      exercise(lbfgs_impl=raw_reference)
      if (not endless): break
  else:
    compare_implementations()
  timer()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
