from __future__ import division
from scitbx.array_family import flex
from scitbx.lbfgs import have_lbfgs_fem, fortran, raw_reference, raw
from libtbx.utils import show_times
from libtbx.test_utils import show_diff
from libtbx import easy_run
import libtbx.load_env
import platform
import random
import re
import sys, os

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

def run_cmd(cmd):
  print cmd
  sys.stdout.flush()
  out = easy_run.fully_buffered(command=cmd)
  err = "\n".join(out.stderr_lines)
  if (len(err) != 0):
    print err
    if (err.find("== ERROR SUMMARY: 0 errors from 0 contexts") < 0):
      raise AssertionError(
        "stderr output does not appear to be valgrind output")
  return "\n".join(out.stdout_lines)

def run_and_compare_sdrive_fem(this_script):
  sdrive_fem = libtbx.env.under_build(path="scitbx/lbfgs/sdrive_fem")
  if (not os.path.isfile(sdrive_fem)):
    return
  outputs = []
  for cmd in [sdrive_fem, 'scitbx.python "%s" fortran 100 5 1 0' % this_script]:
    outputs.append(run_cmd(cmd=cmd))
  assert not show_diff(outputs[0], outputs[1])

def truncate_floats(out):
  match_objects = re.finditer(
    "[ -][0-9]\\.[0-9][0-9][0-9]E[-+][0-9][0-9]", out)
  fragments = []
  k = 0
  for match_obj in match_objects:
    i = match_obj.start()
    j = match_obj.end()
    v = float(out[i:j])
    if (abs(v) < 1e-14):
      v = 0
    else:
      v = float(out[i:j-4])
    fmt = "%%%d.1f" % (j-4-i)
    fragments.append(out[k:i] + fmt % v)
    k = j-4
  fragments.append(out[k:])
  return "".join(fragments)

def replace_e0dd_with_edd(out):
  match_objects = re.finditer("E[-+]0[0-9][0-9]", out)
  fragments = []
  k = 0
  for match_obj in match_objects:
    j = match_obj.start() + 2
    i = out.rfind(" ", k, j)
    assert i >= 0
    if (j - i < 9):
      fragments.append(out[k:i] + " " + out[i:j])
    else:
      fragments.append(out[k:j])
    k = j + 1
  fragments.append(out[k:])
  return "".join(fragments)

def run_and_compare_implementations(this_script, n, m, iprint):
  outputs = []
  for impl in ["fortran", "raw_reference", "raw"]:
    if (impl == "fortran" and not have_lbfgs_fem):
      continue
    cmd = 'scitbx.python "%s" %s %d %d %d %d' % (
      this_script, impl, n, m, iprint[0], iprint[1])
    out = run_cmd(cmd=cmd)
    if (impl == "fortran"):
      out = out.replace("D-", "E-").replace("D+", "E+")
    out = replace_e0dd_with_edd(out=out)
    out = out.replace("E-00", "E+00")
    out = truncate_floats(out=out)
    outputs.append(out)
  assert len(outputs) >= 2
  a = outputs[0]
  for i in xrange(1, len(outputs)):
    b = outputs[i]
    assert not show_diff(a, b)

def compare_implementations():
  this_script = libtbx.env.under_dist(
    module_name="scitbx", path="lbfgs/tst_lbfgs_fem.py")
  assert this_script.find('"') < 0
  run_and_compare_sdrive_fem(this_script=this_script)
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
  if (not have_lbfgs_fem):
    print "Skipping some tests: lbfgs_fem.cpp not linked into scitbx_lbfgs_ext."
  once = "--once" in args
  endless = "--endless" in args
  if (once or endless):
    while 1:
      if (have_lbfgs_fem):
        exercise(lbfgs_impl=fortran)
      exercise(lbfgs_impl=raw_reference)
      exercise(lbfgs_impl=raw)
      if (not endless): break
  else:
    compare_implementations()
  timer()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
