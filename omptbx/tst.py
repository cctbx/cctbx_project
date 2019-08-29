from __future__ import absolute_import, division, print_function
from six.moves import range
import omptbx

def exercise():
  if (omptbx.have_omp_h): print("omptbx.have_omp_h")
  if (omptbx.have_stubs_h): print("omptbx.have_stubs_h")
  print("omtbx.omp_version:", omptbx.omp_version)
  if (omptbx.have_omp_h):
    assert not omptbx.have_stubs_h
  else:
    assert omptbx.have_stubs_h
  #
  env = omptbx.env
  #
  env.dynamic = False
  assert not env.dynamic
  env.dynamic = True
  if (omptbx.have_omp_h and omptbx.omp_version > 199819):
    assert env.dynamic
  else:
    assert not env.dynamic
  #
  env.nested = False
  assert not env.nested
  env.nested = True
  if (omptbx.have_omp_h):
    assert env.nested
  else:
    assert not env.nested
  #
  env.dynamic = False
  env.nested = False
  for i in range(1,5):
    env.num_threads = i
    if (omptbx.have_omp_h):
      assert env.num_threads == i
    else:
      assert env.num_threads == 1
  #
  env.num_threads = 2
  if (omptbx.have_omp_h):
    assert omptbx.ext.tst_environment() == (4,2, 4,4)
  else:
    assert omptbx.ext.tst_environment() == (4,1, 4,1)
  #
  print("OK")

if (__name__ == "__main__"):
  exercise()
