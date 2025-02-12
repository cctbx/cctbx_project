from __future__ import absolute_import, division, print_function
import sys
import random
from cctbx.maptbx.tst_real_space_refinement_simple import exercise_lbfgs
from cctbx.array_family import flex
from scitbx.graph import test_cases_tardy_pdb
from libtbx.utils import null_out


def run(args):
  # Test that the real_space_refinement_simple tests run using a weight vector
  if (1):
    random.seed(0)
    flex.set_random_seed(0)
  out = null_out()
  remaining_args = []
  for arg in args:
    if (arg == "--verbose"): out = sys.stdout
    else: remaining_args.append(arg)
  test_cases = test_cases_tardy_pdb.select_test_cases(
    tags_or_indices=remaining_args)
  for test_case in test_cases:
    print("test case %d: %s" % (test_case.index, test_case.tag), file=out)
    minimized = []
    rsr_weight_vector = flex.double(len(test_case.sites),1)
    minimized.append(exercise_lbfgs(test_case, use_geo=False, out=out))
    minimized.append(exercise_lbfgs(test_case, use_geo=True,rsr_weight=rsr_weight_vector, out=out))
    m0, m1 = minimized

    assert m0.real_space_target_weight == 1
    assert all(m1.real_space_target_weight == 1)
    assert m1.f_final < m0.f_start * 0.98

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

