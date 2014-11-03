
from __future__ import division
import sys
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out

def exercise (args) :
  if 'verbose' in args:
    out=sys.stdout
  else:
    out=null_out()

  print "Running basic estimator....",
  from mmtbx.scaling.bayesian_estimator import bayesian_estimator
  run=bayesian_estimator(out=out)
  result_list=run.exercise()
  assert approx_equal(result_list,[0.975, 0.973, 0.746, 0.784],eps=1.e-3)
  print "OK"

  print "Running basic estimator with randomized data....",
  cc=run.exercise_2(out=out)
  assert approx_equal(cc,0.988,eps=1.e-3)
  print "OK"

  print "Running jacknifed estimators...",
  from mmtbx.scaling.bayesian_estimator import run_jacknife
  cc=run_jacknife(args=[],out=out)
  assert approx_equal(cc,0.889,eps=1.e-3)
  print "OK"

  print "Running group of estimators...",
  from mmtbx.scaling.bayesian_estimator import exercise_group
  cc1,cc2=exercise_group(out=out)
  assert approx_equal(cc1,0.861,eps=1.e-3)
  assert approx_equal(cc2,-0.186,eps=1.e-3)
  print "OK"



if (__name__ == "__main__") :
  exercise(sys.argv[1:])
  print "OK"
