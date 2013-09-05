from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = (
    "$D/refinement/constraints/tests/tst_lbfgs.py",
    "$B/refinement/constraints/tests/tst_reparametrisation",
    "$B/refinement/constraints/tests/tst_geometrical_hydrogens",
    "$B/refinement/constraints/tests/tst_special_position",
    "$D/refinement/constraints/tests/tst_reparametrisation.py",
    ["$D/refinement/constraints/tests/tst_constrained_structure.py",
     '--normal_eqns_solving_method=naive'],
    ["$D/refinement/constraints/tests/tst_constrained_structure.py",
     '--normal_eqns_solving_method=levenberg-marquardt'],
    "$D/refinement/constraints/tests/tst_rigid.py",
    "$D/refinement/constraints/tests/tst_direction.py",
    )

def run():
  build_dir = libtbx.env.under_build("smtbx")
  dist_dir = libtbx.env.dist_path("smtbx")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if __name__ == '__main__':
  run()
