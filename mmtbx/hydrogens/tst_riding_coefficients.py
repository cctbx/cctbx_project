from __future__ import absolute_import, division, print_function
import time
from mmtbx_hydrogens_ext import *

# ----------------------------------------------------
# test some properties of riding_coefficients object
# in particular: deep_copy
# ----------------------------------------------------

# initialize riding_coefficients object
def set_object():
  rc1 = riding_coefficients(
    htype ='flat_2neigbs',
    ih    = 1,
    a0    = 5,
    a1    = 2,
    a2    = 3,
    a3    = 6,
    a     = 3.467,
    b     = 5.4,
    h     = 3.58,
    n     = 2,
    disth = 0.887)
  return rc1

# test getter and setter
def exercise1():
  rc1 = set_object()
  # test getter
  assert (rc1.htype == 'flat_2neigbs')
  assert (rc1.ih    == 1)
  assert (rc1.a0    == 5)
  assert (rc1.a1    == 2)
  assert (rc1.a2    == 3)
  assert (rc1.a3    == 6)
  assert (rc1.a     == 3.467)
  assert (rc1.b     == 5.4)
  assert (rc1.h     == 3.58)
  assert (rc1.n     == 2)
  assert (rc1.disth == 0.887)

  # test setter
  rc1.htype = '3neigbs'
  rc1.ih    = 10
  rc1.a0    = 11
  rc1.a1    = 12
  rc1.a2    = 13
  rc1.a3    = 14
  rc1.a     = 2.5
  rc1.b     = 3.5
  rc1.h     = 4.5
  rc1.n     = 1
  rc1.disth = 1.0

  assert (rc1.htype == '3neigbs')
  assert (rc1.ih    == 10)
  assert (rc1.a0    == 11)
  assert (rc1.a1    == 12)
  assert (rc1.a2    == 13)
  assert (rc1.a3    == 14)
  assert (rc1.a     == 2.5)
  assert (rc1.b     == 3.5)
  assert (rc1.h     == 4.5)
  assert (rc1.n     == 1)
  assert (rc1.disth == 1.0)

# shallow copy
def exercise2():
  rc1 = set_object()
  rc2 = rc1
  assert (rc1.ih == 1)
  assert (rc2.ih == 1)

  rc2.ih = 100
  assert (rc1.ih == 100)
  assert (rc2.ih == 100)

# deep copy
def exercise3():
  rc1 = set_object()
  rc2 = riding_coefficients(rc1)

  rc2.htype = '3neigbs'
  rc2.ih    = 10
  rc2.a0    = 11
  rc2.a1    = 12
  rc2.a2    = 13
  rc2.a3    = 14
  rc2.a     = 2.5
  rc2.b     = 3.5
  rc2.h     = 4.5
  rc2.n     = 1
  rc2.disth = 1.0

  assert (rc1.htype == 'flat_2neigbs')
  assert (rc1.ih    == 1)
  assert (rc1.a0    == 5)
  assert (rc1.a1    == 2)
  assert (rc1.a2    == 3)
  assert (rc1.a3    == 6)
  assert (rc1.a     == 3.467)
  assert (rc1.b     == 5.4)
  assert (rc1.h     == 3.58)
  assert (rc1.n     == 2)
  assert (rc1.disth == 0.887)

  assert (rc2.htype == '3neigbs')
  assert (rc2.ih    == 10)
  assert (rc2.a0    == 11)
  assert (rc2.a1    == 12)
  assert (rc2.a2    == 13)
  assert (rc2.a3    == 14)
  assert (rc2.a     == 2.5)
  assert (rc2.b     == 3.5)
  assert (rc2.h     == 4.5)
  assert (rc2.n     == 1)
  assert (rc2.disth == 1.0)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise1()
  exercise2()
  exercise3()
  print("OK. Time: %8.3f"%(time.time()-t0))
