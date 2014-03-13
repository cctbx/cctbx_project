from __future__ import division

import boost.python
fw_ext = boost.python.import_ext("cctbx_french_wilson_ext")
from scitbx.array_family import flex
import random

def exercise_00():
  x = flex.random_double(1000)
  y = flex.random_double(1000)
  xa = flex.double()
  ya = flex.double()
  ba = flex.bool()
  for x_, y_ in zip(x,y):
    scale1 = random.choice([1.e-6, 1.e-3, 0.1, 1, 1.e+3, 1.e+6])
    scale2 = random.choice([1.e-6, 1.e-3, 0.1, 1, 1.e+3, 1.e+6])
    b = random.choice([True, False])
    x_ = x_*scale1
    y_ = y_*scale2
    v1 = fw_ext.expectEFW(eosq=x_, sigesq=y_, centric=b)
    v2 = fw_ext.expectEsqFW(eosq=x_, sigesq=y_, centric=b)
    assert type(v1) == type(1.)
    assert type(v2) == type(1.)
    xa.append(x_)
    ya.append(y_)
    ba.append(b)
  fw_ext.is_FrenchWilson(F=xa, SIGF=ya, is_centric=ba, eps=0.001)

if (__name__ == "__main__"):
  exercise_00()
