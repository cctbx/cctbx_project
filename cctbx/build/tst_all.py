import os, os.path

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
"./sftbx/tst.py",
"./sftbx/tst_basic.py",
)

for tst in tst_list:
  print tst
  if (tst.endswith(".py")):
    os.system("python " + os.path.normpath(tst))
  else:
    os.system(os.path.normpath(tst))
  print
