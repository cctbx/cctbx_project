import os, shutil
import generate_operator_traits_builtin
import generate_af_std_imports
import generate_af_operator_functors
import generate_af_algebras
import generate_af_apply
import generate_af_tiny_bpl
import generate_henke_cpp
import generate_sasaki_cpp
generate_operator_traits_builtin.run()
generate_af_std_imports.run()
generate_af_operator_functors.run()
generate_af_algebras.run()
generate_af_apply.run()
generate_af_tiny_bpl.run()
generate_henke_cpp.run()
generate_sasaki_cpp.run()
array_family_include = "../include/cctbx/array_family/"
misc_src = "../misc/"
eltbx_src = "../eltbx/"
for file, dir in (
  ("operator_traits_builtin.h", array_family_include),
  ("std_imports.h", array_family_include),
  ("operator_functors.h", array_family_include),
  ("ref_algebra.h", array_family_include),
  ("tiny_algebra.h", array_family_include),
  ("small_algebra.h", array_family_include),
  ("shared_algebra.h", array_family_include),
  ("versa_algebra.h", array_family_include),
  ("tiny_plain_apply.h", array_family_include),
  ("small_plain_apply.h", array_family_include),
  ("shared_plain_apply.h", array_family_include),
  ("versa_plain_apply.h", array_family_include),
  ("tiny_apply.h", array_family_include),
  ("small_apply.h", array_family_include),
  ("shared_apply.h", array_family_include),
  ("versa_apply.h", array_family_include),
  ("ref_apply.h", array_family_include),
  ("tiny_bpl.h", array_family_include),
  ("tiny_bpl.cpp", misc_src),
  ("henke.cpp", eltbx_src),
  ("henke_tables_01_12.cpp", eltbx_src),
  ("henke_tables_13_24.cpp", eltbx_src),
  ("henke_tables_25_36.cpp", eltbx_src),
  ("henke_tables_37_48.cpp", eltbx_src),
  ("henke_tables_49_60.cpp", eltbx_src),
  ("henke_tables_61_72.cpp", eltbx_src),
  ("henke_tables_73_84.cpp", eltbx_src),
  ("henke_tables_85_92.cpp", eltbx_src),
  ("sasaki.cpp", eltbx_src),
  ("sasaki_tables_01_12.cpp", eltbx_src),
  ("sasaki_tables_13_24.cpp", eltbx_src),
  ("sasaki_tables_25_36.cpp", eltbx_src),
  ("sasaki_tables_37_48.cpp", eltbx_src),
  ("sasaki_tables_49_60.cpp", eltbx_src),
  ("sasaki_tables_61_72.cpp", eltbx_src),
  ("sasaki_tables_73_82.cpp", eltbx_src),
):
  print "Copying " + dir + file
  shutil.copy(file, dir + file)
  os.unlink(file)
