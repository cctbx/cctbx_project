import sys, os, os.path

tst_list = (
"$B/array_family/tst_af_1",
"$B/array_family/tst_af_2",
"$B/array_family/tst_af_3",
"$B/array_family/tst_af_4",
"$B/array_family/tst_af_5",
"$B/array_family/tst_vec3",
"$B/array_family/tst_mat3",
"$B/array_family/tst_sym_mat3",
"$D/array_family/boost_python/regression_test.py",
"$D/array_family/boost_python/tst_flex.py",
"$D/lbfgs/boost_python/tst_lbfgs.py",
"$D/fftpack/boost_python/tst_fftpack.py",
)

build = os.path.join(os.environ["LIBTBX_BUILD"], "scitbx")
dist = os.environ["SCITBX_DIST"]

for tst in tst_list:
  if (tst.startswith("$B")):
    tst_path = tst.replace("$B", build)
  else:
    tst_path = tst.replace("$D", dist)
  assert tst_path.find("$") < 0
  tst_path = os.path.normpath(tst_path)
  print tst_path
  sys.stdout.flush()
  if (tst_path.endswith(".py")):
    os.system("python " + tst_path)
  else:
    os.system(tst_path)
  print
  sys.stderr.flush()
  sys.stdout.flush()
