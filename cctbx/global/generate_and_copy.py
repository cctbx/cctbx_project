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
array_family_include = "../include/cctbx/array_family/"
for file, dir in (
  ("operator_traits_builtin.h", array_family_include),
  ("std_imports.h", array_family_include),
  ("operator_functors.h", array_family_include),
  ("flagged_value_algebra.h", array_family_include),
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
):
  print "Copying " + dir + file
  shutil.copy(file, dir + file)
  os.unlink(file)
