
from __future__ import absolute_import, division, print_function
from cctbx.r_free_utils import *
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from scitbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO
from itertools import count
import random

def exercise():
  random.seed(12345)
  flex.set_random_seed(12345)
  flags_1 = assign_random_r_free_flags(n_refl=100000, fraction_free=0.05)
  assert (flags_1.count(True) == 5000)
  flags_1_ccp4 = assign_random_r_free_flags(n_refl=100000, fraction_free=0.05,
    format="ccp4")
  # XXX this is the best we can do with the current method
  assert (flags_1_ccp4.count(0) > 4000) and (flags_1_ccp4.count(0) < 6000)
  flags_1_shelx = assign_random_r_free_flags(n_refl=100000, fraction_free=0.05,
    format="shelx")
  assert (flags_1_shelx.count(-1) == 5000)
  flags_2 = assign_r_free_flags_by_shells(n_refl=100000,
    fraction_free=0.05,
    n_bins=50)
  assert (flags_2.count(True) == 5000)
  ccp4_flags = export_r_free_flags_for_ccp4(flags_1, True)
  assert (ccp4_flags.count(0) == flags_1.count(True))
  assert (flex.max(ccp4_flags) == 19)
  shelx_flags = export_r_free_flags_for_shelx(flags_1, True)
  assert ((shelx_flags==-1).all_eq((ccp4_flags==0)))
  flags_3 = assign_random_r_free_flags(n_refl=100000, fraction_free=0.025)
  assert (flags_3.count(True) == 2500)
  ccp4_flags = export_r_free_flags_for_ccp4(flags_3, True)
  assert (ccp4_flags.count(0) == flags_3.count(True))
  assert (flex.max(ccp4_flags) == 39)
  # now with an actual Miller array
  symm = crystal.symmetry(
    space_group_info=sgtbx.space_group_info("P212121"),
    unit_cell=uctbx.unit_cell((6,7,8,90,90,90)))
  set1 = miller.build_set(
    crystal_symmetry=symm,
    anomalous_flag=True,
    d_min=1.0)
  flags_4 = set1.generate_r_free_flags()
  stats = get_r_free_stats(flags_4, True)
  assert (19 <= stats[0] <= 25) # XXX is this even necessary?
  # much larger for the last few tests
  symm = crystal.symmetry(
    space_group_info=sgtbx.space_group_info("P212121"),
    unit_cell=uctbx.unit_cell((60,70,80,90,90,90)))
  set1 = miller.build_set(
    crystal_symmetry=symm,
    anomalous_flag=True,
    d_min=1.0)
  #print set1.indices().size()
  flags_5 = set1.generate_r_free_flags(fraction=0.1, max_free=None)
  flags_6 = adjust_fraction(flags_5, 0.15, log=null_out())
  frac_6 = flags_6.data().count(True) / flags_6.data().size()
  assert approx_equal(frac_6, 0.15, eps=0.001)
  flags_7 = adjust_fraction(flags_5, 0.05, log=null_out())
  frac_7 = flags_7.data().count(True) / flags_7.data().size()
  assert approx_equal(frac_7, 0.05, eps=0.001)
  n_flipped = 0
  for i_hkl, (h,k,l) in enumerate(flags_5.indices()):
    if (i_hkl % 100 == 0) and (h > 0) and (k > 0) and (l > 0):
      flag = flags_5.data()[i_hkl]
      if (not flag):
        n_flipped += 1
        flags_5.data()[i_hkl] = True
  # XXX check this for reproducibility on other systems
  #assert (n_flipped == 1559) # setting the random seed should ensure this
  try :
    flags_8 = flags_5.average_bijvoet_mates()
  except Sorry :
    pass
  else :
    raise Exception_expected
  out = StringIO()
  flags_9 = remediate_mismatches(flags_5, log=out)
  # XXX check this for reproducibility on other systems
  #assert (out.getvalue() == "  1559 reflections moved to test set\n")
  flags_10 = flags_9.average_bijvoet_mates()

if (__name__ == "__main__"):
  exercise()
  print("OK")
