import os, shutil
import generate_operator_traits_builtin
import generate_af_std_imports
import generate_af_operator_functors
import generate_af_algebras
import generate_af_apply
import generate_flagged_value_algebra
import generate_af_tiny_bpl
generate_operator_traits_builtin.run()
generate_af_std_imports.run()
generate_af_operator_functors.run()
generate_af_algebras.run()
generate_af_apply.run()
generate_flagged_value_algebra.run()
generate_af_tiny_bpl.run()
for file, dir in (
  ("operator_traits_builtin.h", "../include/cctbx/array_family/"),
  ("std_imports.h", "../include/cctbx/array_family/"),
  ("operator_functors.h", "../include/cctbx/array_family/"),
  ("flagged_value_algebra.h", "../include/cctbx/array_family/"),
  ("ref_algebra.h", "../include/cctbx/array_family/"),
  ("tiny_algebra.h", "../include/cctbx/array_family/"),
  ("small_algebra.h", "../include/cctbx/array_family/"),
  ("shared_algebra.h", "../include/cctbx/array_family/"),
  ("versa_algebra.h", "../include/cctbx/array_family/"),
  ("tiny_plain_apply.h", "../include/cctbx/array_family/"),
  ("small_plain_apply.h", "../include/cctbx/array_family/"),
  ("shared_plain_apply.h", "../include/cctbx/array_family/"),
  ("versa_plain_apply.h", "../include/cctbx/array_family/"),
  ("tiny_apply.h", "../include/cctbx/array_family/"),
  ("small_apply.h", "../include/cctbx/array_family/"),
  ("shared_apply.h", "../include/cctbx/array_family/"),
  ("versa_apply.h", "../include/cctbx/array_family/"),
  ("tiny_bpl.h", "../include/cctbx/array_family/"),
):
  print "Copying " + dir + file
  shutil.copy(file, dir + file)
  os.unlink(file)
