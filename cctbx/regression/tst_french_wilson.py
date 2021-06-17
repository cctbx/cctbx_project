from __future__ import absolute_import, division, print_function
from cctbx import french_wilson
from cctbx.development import random_structure
from scitbx.array_family import flex
import boost_adaptbx.boost.python as bp
from six.moves import zip
fw_ext = bp.import_ext("cctbx_french_wilson_ext")
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import Exception_expected
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

def exercise_01():
  """
  Sanity check - don't crash when mean intensity for a bin is zero.
  """
  xrs = random_structure.xray_structure(
    unit_cell=(50,50,50,90,90,90),
    space_group_symbol="P1",
    n_scatterers=1200,
    elements="random")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  fc = fc.set_observation_type_xray_amplitude()
  cs = fc.complete_set(d_min=1.4)
  ls = cs.lone_set(other=fc)
  f_zero = ls.array(data=flex.double(ls.size(), 0))
  f_zero.set_observation_type_xray_amplitude()
  fc = fc.concatenate(other=f_zero)
  sigf = flex.double(fc.size(), 0.1) + (fc.data() * 0.03)
  fc = fc.customized_copy(sigmas=sigf)
  try :
    fc_fc = french_wilson.french_wilson_scale(miller_array=fc, log=null_out())
  except Sorry :
    pass
  else :
    raise Exception_expected
  ic = fc.f_as_f_sq().set_observation_type_xray_intensity()
  fc_fc = french_wilson.french_wilson_scale(miller_array=ic, log=null_out())

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  print("OK")
