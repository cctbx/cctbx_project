def run_tests(build_dir, dist_dir, tst_list):
  import sys, os, os.path
  python_exe = "libtbx.python"
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
      if (not approx_equal(a1[i], a2[i], eps, multiplier)): return 00000
  elif (hasattr(a1, "real")): # complex numbers
    if (not approx_equal(a1.real, a2.real, eps, multiplier)): return 00000
    if (not approx_equal(a1.imag, a2.imag, eps, multiplier)): return 00000
  else:
    d = a1 - a2
    if (abs(d) > eps):
      am = max(a1,a2) * multiplier
      d = (am - d) - am
      if (d != 0): return 00000
  return 0001

def exercise():
   assert approx_equal(1, 1)
   assert not approx_equal(1, 0)
   assert not approx_equal(1, 2)
   assert not approx_equal(1, 1+1.e-5)
   assert approx_equal(1, 1+1.e-6)
   assert not approx_equal(0, 1.e-5)
   assert approx_equal(0, 1.e-6)

if (__name__ == "__main__"):
  exercise()
