from stdlib import math
import sys

def run_tests(build_dir, dist_dir, tst_list):
  import sys, os, os.path
  python_exe = sys.executable
  for tst in tst_list:
    cmd_args = ""
    if (type(tst) == type([])):
      if ("--Verbose" in sys.argv[1:]):
        cmd_args = " " + " ".join(["--Verbose"] + tst[1:])
      elif ("--Quick" in sys.argv[1:]):
        cmd_args = " " + " ".join(tst[1:])
      tst = tst[0]
    else:
      if ("--Verbose" in sys.argv[1:]):
        continue
    if (tst.startswith("$B")):
      tst_path = tst.replace("$B", build_dir)
    else:
      tst_path = tst.replace("$D", dist_dir)
    assert tst_path.find("$") < 0
    tst_path = os.path.normpath(tst_path)
    if (tst_path.endswith(".py")):
      cmd = python_exe + " " + tst_path
    else:
      cmd = tst_path
    if ("--valgrind" in sys.argv[1:]):
      cmd = "valgrind " + cmd
    elif ("--memcheck" in sys.argv[1:]):
      cmd = "valgrind --tool=memcheck " + cmd
    cmd += cmd_args
    print cmd
    sys.stdout.flush()
    os.system(cmd)
    print
    sys.stderr.flush()
    sys.stdout.flush()

def approx_equal(a1, a2, eps=1.e-6, multiplier=1.e10):
  if (hasattr(a1, "__len__")): # traverse list
    assert len(a1) == len(a2)
    for i in xrange(len(a1)):
      if (not approx_equal(a1[i], a2[i], eps, multiplier)): return False
  elif (hasattr(a1, "real")): # complex numbers
    if (not approx_equal(a1.real, a2.real, eps, multiplier)): return False
    if (not approx_equal(a1.imag, a2.imag, eps, multiplier)): return False
  else:
    d = a1 - a2
    if (abs(d) > eps):
      if (multiplier is None): return False
      am = max(a1,a2) * multiplier
      d = (am - d) - am
      if (d != 0): return False
  return True

def eps_eq(a1, a2, eps=1.e-6):
  if (hasattr(a1, "__len__")): # traverse list
    assert len(a1) == len(a2)
    for i in xrange(len(a1)):
      if (not eps_eq(a1[i], a2[i], eps)): return False
  elif (hasattr(a1, "real")): # complex numbers
    if (not eps_eq(a1.real, a2.real, eps)): return False
    if (not eps_eq(a1.imag, a2.imag, eps)): return False
  else:
    if (a1 == 0 or a2 == 0):
      if (abs(a1-a2) > eps): return False
    else:
      l1 = round(math.log(abs(a1)))
      l2 = round(math.log(abs(a2)))
      m = math.exp(-max(l1, l2))
      if (abs(a1*m-a2*m) > eps): return False
  return True

def exercise():
   assert approx_equal(1, 1)
   assert not approx_equal(1, 0)
   assert not approx_equal(1, 2)
   assert not approx_equal(1, 1+1.e-5)
   assert approx_equal(1, 1+1.e-6)
   assert not approx_equal(0, 1.e-5)
   assert approx_equal(0, 1.e-6)
   assert eps_eq(1, 1)
   assert not eps_eq(1, 0)
   assert not eps_eq(1, 2)
   assert not eps_eq(1, 1+1.e-5)
   assert eps_eq(1, 1+1.e-6)
   assert not eps_eq(0, 1.e-5)
   assert eps_eq(0, 1.e-6)

if (__name__ == "__main__"):
  exercise()
