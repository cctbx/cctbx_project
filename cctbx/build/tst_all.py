import sys, os, os.path

tst_list = (
"./arraytbx/tst_af_1",
"./arraytbx/tst_af_2",
"./arraytbx/tst_af_3",
"./arraytbx/tst_vec3",
"./arraytbx/tst_mat3",
"./arraytbx/tst_sym_mat3",
"./arraytbx/tst_flex.py",
"./arraytbx/tst_flex_utils.py",
"./arraytbx/tst_shared.py",
"./arraytbx/tst_shared_map.py",
"./eltbx/tst_caasf_it1992.py",
"./eltbx/tst_caasf_wk1995.py",
"./eltbx/tst_henke.py",
"./eltbx/tst_icsd_radii.py",
"./eltbx/tst_neutron.py",
"./eltbx/tst_sasaki.py",
"./eltbx/tst_tinypse.py",
"./eltbx/tst_wavelengths.py",
"./fftbx/tst.py",
"./lbfgs/tst.py",
"./uctbx/tst.py",
"./sgtbx/tst.py",
"./miller/tst.py",
"./mintbx/tst.py",
"./adptbx/tst.py",
"./sftbx/tst_basic.py",
"./sftbx/tst.py",
"$CCTBX_DIST/cctbx/euclidean_model_matching.py",
)

cctbx_dist = os.environ["CCTBX_DIST"]

for tst in tst_list:
  tst_path = tst.replace("$CCTBX_DIST", cctbx_dist)
  print tst_path
  sys.stdout.flush()
  if (tst_path.endswith(".py")):
    os.system("python " + os.path.normpath(tst_path))
  else:
    os.system(os.path.normpath(tst_path))
  print
  sys.stderr.flush()
  sys.stdout.flush()
