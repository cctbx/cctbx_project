def run_tests(build_dir, dist_dir, tst_list):
  import sys, os, os.path
  try: python_exe = os.environ["LIBTBX_PYTHON_EXE"]
  except: python_exe = "python"
  for tst in tst_list:
    cmd_args = ""
    if (type(tst) == type([])):
      if ("--Verbose" in sys.argv[1:]):
        cmd_args = " " + " ".join(["--Verbose"] + tst[1:])
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
    cmd += cmd_args
    print cmd
    sys.stdout.flush()
    os.system(cmd)
    print
    sys.stderr.flush()
    sys.stdout.flush()

def approx_equal(a1, a2, eps=1.e-6):
  if (hasattr(a1, "__len__")): # traverse list
    assert len(a1) == len(a2)
    for i in xrange(len(a1)):
      if (not approx_equal(a1[i], a2[i], eps)): return 0
  elif (hasattr(a1, "real")): # complex numbers
    if (not approx_equal(a1.real, a2.real, eps)): return 0
    if (not approx_equal(a1.imag, a2.imag, eps)): return 0
  else:
    if (abs(a1 - a2) > eps): return 0
  return 1
