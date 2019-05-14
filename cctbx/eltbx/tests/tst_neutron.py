from __future__ import absolute_import, division, print_function
from cctbx.eltbx import neutron
from libtbx.test_utils import approx_equal, Exception_expected

def exercise_00():
  t = neutron.neutron_news_1992_table("eu")
  assert t.label() == "Eu"
  l = t.bound_coh_scatt_length()
  assert approx_equal(l.real, 7.22)
  assert approx_equal(l.imag, -1.26)
  assert approx_equal(t.abs_cross_sect(), 4530.)
  n = 0
  for t in neutron.neutron_news_1992_table_iterator():
    n += 1
    if (n == 1):
      assert t.label() == "H"
    elif (n == 189):
      assert t.label() == "U6+"
    u = neutron.neutron_news_1992_table(t.label())
    assert u.label() == t.label()
  assert n == 189, n
  try:
    t = neutron.neutron_news_1992_table("XX")
  except ValueError as e:
    pass
  else:
    raise Exception_expected

def exercise_01():
  def strip_num_and_sign(x):
    r = ""
    for i in x:
      if(i.isalpha()): r+=i
    return r
  for t1 in neutron.neutron_news_1992_table_iterator():
    l1 = t1.label()
    l2 = strip_num_and_sign(x=l1)
    t1 = neutron.neutron_news_1992_table(l1)
    t2 = neutron.neutron_news_1992_table(l1)
    assert approx_equal(t1.abs_cross_sect(), t2.abs_cross_sect())
    assert approx_equal(t1.bound_coh_scatt_length(),t2.bound_coh_scatt_length())

def run():
  exercise_00()
  exercise_01()
  print("OK")

if (__name__ == "__main__"):
  run()
