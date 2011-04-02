import sys, os

def run(args):
  if (len(args) == 0):
    n_iterations = 1000
  elif (len(args) == 1):
    n_iterations = int(args[0])
    assert n_iterations > 0
  else:
    from libtbx.utils import Usage
    raise Usage("scitbx.show_exp_times [n_iterations]")
  evar = "LIBTBX_NO_LD_PRELOAD"
  evar_set = os.environ.has_key(evar)
  if (evar_set):
    print "%s set:" % evar
  else:
    print "%s not set:" % evar
  from scitbx.math.tests.tst_exp_functions import \
    exercise_with_random_arguments as exercise
  exercise(n_arguments=10000, n_iterations=n_iterations)
  print
  sys.stdout.flush()
  if (not evar_set):
    if (os.environ.has_key("LD_PRELOAD")):
      del os.environ["LD_PRELOAD"]
    os.environ[evar] = "1"
    from libtbx import easy_run
    easy_run.call(command="scitbx.show_exp_times %d" % n_iterations)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
