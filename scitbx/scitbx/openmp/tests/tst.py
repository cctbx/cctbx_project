from scitbx import openmp
if openmp.available:
  import boost.python
  ext = boost.python.import_ext("scitbx_openmp_tests_ext")
import libtbx.load_env
import os

def exercise_environment():
  env = openmp.environment
  assert env.dynamic in (True, False)
  assert env.nested == False
  env.dynamic = True
  assert env.dynamic == True
  env.dynamic = False
  assert env.dynamic == False
  env.nested = True
  assert env.nested == True
  for i in xrange(1,5):
    env.num_threads = i
    assert env.num_threads == i

  env.num_threads = 2
  assert ext.tst_environment() == (4,2, 4,4)

def exercise():
  if not openmp.available:
    print "Skipping OpenMP test"
  else:
    exercise_environment()
  print 'OK'

if __name__ == '__main__':
  exercise()
