def run_tests(build_dir, dist_dir, tst_list):
  import sys, os, os.path
  try: python_exe = os.environ["LIBTBX_PYTHON_EXE"]
  except: python_exe = "python"
  for tst in tst_list:
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
    print cmd
    sys.stdout.flush()
    os.system(cmd)
    print
    sys.stderr.flush()
    sys.stdout.flush()
