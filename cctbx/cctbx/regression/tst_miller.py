from cctbx import crystal
from cctbx import miller
from cctbx import xray
from cctbx import maptbx
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import cctbx
import scitbx.math
from scitbx.python_utils import complex_math
from libtbx.test_utils import approx_equal
from cStringIO import StringIO
import pickle
import random
import math
import sys

def exercise_set():
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,2,3), (0,0,4)))
  ms = miller.set(xs, mi)
  ms = miller.set(xs, mi, False)
  ms = miller.set(xs, mi, True)
  assert ms.indices() == mi
  assert ms.anomalous_flag() == True
  mc = ms.copy()
  assert not mc is ms
  assert mc.unit_cell() is ms.unit_cell()
  assert mc.space_group_info() is ms.space_group_info()
  assert mc.indices() is ms.indices()
  assert mc.anomalous_flag() is ms.anomalous_flag()
  mc = ms.deep_copy()
  assert mc.unit_cell().is_similar_to(ms.unit_cell())
  assert mc.space_group() == ms.space_group()
  assert flex.order(mc.indices(), ms.indices()) == 0
  assert mc.anomalous_flag() == ms.anomalous_flag()
  assert tuple(ms.multiplicities().data()) == (4, 2)
  assert tuple(ms.epsilons().data()) == (1, 2)
  assert approx_equal(ms.d_star_sq().data(), [0.7211111, 0.64])
  assert approx_equal(tuple(ms.d_spacings().data()), (1.177603, 1.25))
  assert approx_equal(tuple(ms.sin_theta_over_lambda_sq().data()),
                      (0.1802778, 0.16))
  assert approx_equal(ms.d_min(), 1.177603)
  assert approx_equal(ms.resolution_range(), (1.25, 1.177603))
  p1 = ms.expand_to_p1()
  assert p1.indices().size() == 6
  b = p1.setup_binner(auto_binning=True)
  b = p1.setup_binner(reflections_per_bin=1)
  b = p1.setup_binner(n_bins=8)
  assert id(p1.binner()) == id(b)
  assert b.limits().size() == 9
  assert tuple(ms.sort().indices()) == ((0,0,4), (1,2,3))
  assert tuple(ms.sort(reverse=True).indices()) == ((1,2,3), (0,0,4))
  assert tuple(ms.sort(reverse=True)
                .sort(by_value="packed_indices").indices()) \
      == ((0, 0, 4), (1, 2, 3))
  assert tuple(ms.sort(reverse=True)
                .sort(by_value="packed_indices", reverse=True).indices()) \
      == ((1, 2, 3), (0, 0, 4))
  ms = miller.set(xs, mi, False)
  mp = ms.patterson_symmetry()
  assert str(mp.space_group_info()) == "P m m m"
  assert mp.indices() == ms.indices()
  mc = ms.complete_set()
  c = mc.completeness()
  assert c >= 1-1.e5
  assert c <= 1
  c = mc.select(~mc.all_selection()).completeness()
  assert c == 0
  ma = ms.map_to_asu()
  assert flex.order(ms.indices(), ma.indices()) == 0
  ma = ms.remove_systematic_absences()
  assert flex.order(ms.indices(), ma.indices()) == 0
  assert miller.set(xs, mi).auto_anomalous().anomalous_flag() == False
  mi.extend(flex.miller_index(((-1,-2,-3), (3,4,5), (-3,-4,-5))))
  ma = miller.set(xs, mi)
  assert ma.n_bijvoet_pairs() == 2
  assert ma.auto_anomalous().anomalous_flag() == True
  assert ma.auto_anomalous(
    min_n_bijvoet_pairs=2).anomalous_flag() == True
  assert ma.auto_anomalous(
    min_n_bijvoet_pairs=3).anomalous_flag() == False
  assert ma.auto_anomalous(
    min_fraction_bijvoet_pairs=4/5.-1.e-4).anomalous_flag() == True
  assert ma.auto_anomalous(
    min_fraction_bijvoet_pairs=4/5.+1.e-4).anomalous_flag() == False
  s = StringIO()
  mc.show_comprehensive_summary(f=s)
  assert s.getvalue() == """\
Number of Miller indices: 36
Anomalous flag: %s
Unit cell: (3, 4, 5, 90, 90, 90)
Space group: P 2 2 2 (No. 16)
Systematic absences: 0
Centric reflections: 27
Resolution range: 5 1.1776
Completeness in resolution range: 1
Completeness with d_max=infinity: 1
""" % str(False)
  c = mc.change_basis(cb_op="k,h,-l")
  assert c.unit_cell().is_similar_to(uctbx.unit_cell((4,3,5,90,90,90)))

def exercise_generate_r_free_flags():
  for anomalous_flag in [False, True]:
    miller_set = miller.build_set(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(28.174, 52.857, 68.929, 90, 90, 90),
        space_group_symbol="P 21 21 21"),
      anomalous_flag=anomalous_flag,
      d_min=8)
    for i_trial in xrange(10):
      if (i_trial == 0):
        trial_set = miller_set
      else:
        trial_set = miller_set.select(
          flex.random_permutation(size=miller_set.indices().size()))
        if (i_trial >= 5):
          trial_set = trial_set.select(
            flex.random_double(size=miller_set.indices().size()) < 0.8)
      flags = trial_set.generate_r_free_flags()
      if (not anomalous_flag):
        if (i_trial < 5):
          assert flags.indices().size() == 145
          assert flags.data().count(True) == 15
      else:
        if (i_trial < 5):
          assert flags.indices().size() == 212
        fp,fm = flags.hemispheres_acentrics()
        assert fp.data().all_eq(fm.data())
        if (i_trial < 5):
          assert fp.data().count(True) \
               + flags.select_centric().data().count(True) == 15
        flags = flags.select(~flex.bool(
          flags.indices().size(),
          flags.match_bijvoet_mates()[1].pairs_hemisphere_selection("-")))
      flags = flags.sort(by_value="resolution", reverse=True)
      isel = flags.data().iselection()
      flag_distances = isel[1:] - isel[:-1]
      assert flex.min(flag_distances) > 0
      assert flex.max(flag_distances) <= 20

def exercise_binner():
  crystal_symmetry = crystal.symmetry(
    unit_cell="14.311  57.437  20.143",
    space_group_symbol="C m c m")
  for anomalous_flag in [False, True]:
    set1 = miller.build_set(
      crystal_symmetry=crystal_symmetry,
      anomalous_flag=anomalous_flag,
      d_min=10)
    set1.setup_binner(n_bins=3)
    s = StringIO()
    set1.binner().show_summary(f=s, prefix="..")
    assert s.getvalue() == """\
..unused:         - 28.7186 [0/0]
..bin  1: 28.7186 - 14.1305 [3/3]
..bin  2: 14.1305 - 11.4473 [3/3]
..bin  3: 11.4473 - 10.0715 [2/2]
..unused: 10.0715 -         [0/0]
"""
    set2 = miller.build_set(
      crystal_symmetry=crystal_symmetry,
      anomalous_flag=anomalous_flag,
      d_min=8)
    set2.use_binning_of(set1)
    s = StringIO()
    set2.completeness(use_binning=True).show(f=s, prefix=". ")
    assert s.getvalue() == """\
. unused:         - 28.7186 [0/0]
. bin  1: 28.7186 - 14.1305 [3/3] 1.000
. bin  2: 14.1305 - 11.4473 [3/3] 1.000
. bin  3: 11.4473 - 10.0715 [2/2] 1.000
. unused: 10.0715 -         [8/8] 1.000
"""
    s = StringIO()
    set2.completeness(use_binning=True).show(show_bin_number=False, f=s)
    assert s.getvalue() == """\
        - 28.7186 [0/0]
28.7186 - 14.1305 [3/3] 1.000
14.1305 - 11.4473 [3/3] 1.000
11.4473 - 10.0715 [2/2] 1.000
10.0715 -         [8/8] 1.000
"""
    s = StringIO()
    set2.completeness(use_binning=True).show(show_d_range=False, f=s)
    assert s.getvalue() == """\
unused: [0/0]
bin  1: [3/3] 1.000
bin  2: [3/3] 1.000
bin  3: [2/2] 1.000
unused: [8/8] 1.000
"""
    s = StringIO()
    set2.completeness(use_binning=True).show(show_counts=False, f=s)
    assert s.getvalue() == """\
unused:         - 28.7186
bin  1: 28.7186 - 14.1305 1.000
bin  2: 14.1305 - 11.4473 1.000
bin  3: 11.4473 - 10.0715 1.000
unused: 10.0715 -         1.000
"""
    s = StringIO()
    set2.completeness(use_binning=True).show(
      show_bin_number=False, show_d_range=False, show_counts=False, f=s)
    assert s.getvalue() == """\

 1.000
 1.000
 1.000
 1.000
"""
  set2 = set2.customized_copy(indices=flex.miller_index())
  set2.use_binning_of(set1)
  s = StringIO()
  set2.completeness(use_binning=True).show(f=s)
  assert s.getvalue() == """\
unused:         - 28.7186 [0/0]
bin  1: 28.7186 - 14.1305 [0/3] 0.000
bin  2: 14.1305 - 11.4473 [0/3] 0.000
bin  3: 11.4473 - 10.0715 [0/2] 0.000
unused: 10.0715 -         [0/0]
"""
  try:
    pickle.dumps(set1.binner())
  except RuntimeError, e:
    assert str(e).startswith(
      "cctbx.miller.binner instances are not picklable.")
  else: raise RuntimeError("Exception expected.")
  bng = miller.binning(set1.unit_cell(), 10, set1.indices(), 0, 0)
  try:
    pickle.dumps(bng)
  except RuntimeError, e:
    assert str(e).startswith(
      "cctbx.miller.binning instances are not picklable.")
  else: raise RuntimeError("Exception expected.")

def exercise_crystal_gridding():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(95.2939, 95.2939, 98.4232, 94.3158, 115.226, 118.822),
    space_group_symbol="Hall: C 2y (x+y,-x+y+z,z)")
  f_obs = miller.build_set(crystal_symmetry, anomalous_flag=False, d_min=3.5)
  symmetry_flags = sgtbx.search_symmetry_flags(
    use_space_group_symmetry=False,
    use_space_group_ltr=0,
    use_seminvariants=True,
    use_normalizer_k2l=False,
    use_normalizer_l2n=False)
  crystal_gridding_tags = f_obs.crystal_gridding(
    symmetry_flags=symmetry_flags,
    resolution_factor=1/3.,
    mandatory_factors=(20,20,20)).tags()
  assert crystal_gridding_tags.n_real() == (100,100,100)

def exercise_array():
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  data = flex.double((1,2))
  sigmas = flex.double((0.1,0.2))
  ms = miller.set(xs, mi)
  ma = miller.array(ms)
  ma = miller.array(ms, data)
  s = StringIO()
  ma.show_array(f=s, prefix=": ")
  assert s.getvalue() == """\
: (1, -2, 3) 1.0
: (0, 0, -4) 2.0
"""
  ma = miller.array(ms, data, sigmas)
  s = StringIO()
  ma.show_array(f=s, prefix=" :")
  assert s.getvalue() == """\
 :(1, -2, 3) 1.0 0.1
 :(0, 0, -4) 2.0 0.2
"""
  ma = miller.array(ms, data, sigmas).set_info("test")
  assert ma.indices() == mi
  assert ma.data() == data
  assert ma.sigmas() == sigmas
  assert ma.info() == "test"
  assert ma.observation_type() is None
  assert ma.size() == 2
  ma.set_info(miller.array_info(
    source="file", labels=["a", "b", "c"]))
  assert str(ma.info()) == "file:a,b,c"
  assert ma.info().label_string() == "a,b,c"
  ma.set_info(miller.array_info(
    source="file", labels=["a", "b", "c"], merged=True))
  assert str(ma.info()) == "file:a,b,c,merged"
  assert ma.info().label_string() == "a,b,c,merged"
  ma.set_info("Test")
  assert ma.info() == "Test"
  ma.set_observation_type_xray_amplitude()
  assert ma.is_xray_amplitude_array()
  ma.set_observation_type_xray_intensity()
  assert ma.is_xray_intensity_array()
  ac = ma.deep_copy()
  assert flex.order(ac.data(), ma.data()) == 0
  assert flex.order(ac.sigmas(), ma.sigmas()) == 0
  assert ac.info() == "Test"
  assert ac.is_xray_intensity_array()
  aa = ac.as_amplitude_array()
  assert aa.as_amplitude_array() is aa
  ai = aa.as_intensity_array()
  assert approx_equal(ai.data(), ac.data())
  assert ai.as_intensity_array() is ai
  assert aa.eliminate_sys_absent() is aa
  aa = miller.array(
    miller_set=miller.set(
      crystal_symmetry=crystal.symmetry(
        unit_cell=aa.unit_cell(),
        space_group_symbol="F222"),
      indices=aa.indices()),
    data=aa.data())
  ae = aa.eliminate_sys_absent()
  assert ae is not aa
  assert tuple(ae.indices()) == ((0,0,-4),)
  asu = ma.map_to_asu()
  assert tuple(asu.indices()) == ((1,2,3), (0,0,4))
  mi = flex.miller_index(((1,2,3), (-1,-2,-3), (2,3,4), (-2,-3,-4), (3,4,5)))
  data = flex.double((1,2,5,3,6))
  sigmas = flex.double((0.1,0.2,0.3,0.4,0.5))
  ms = miller.set(xs, mi, anomalous_flag=True)
  ma = miller.array(ms, data, sigmas)
  assert ma.discard_sigmas().sigmas() is None
  ad = ma.anomalous_differences()
  assert tuple(ad.indices()) == ((1,2,3), (2,3,4))
  assert approx_equal(tuple(ad.data()), (-1.0, 2.0))
  assert approx_equal(tuple(ad.sigmas()), (math.sqrt(0.05), 0.5))
  for hp,hm in ((ma.hemisphere_acentrics("+"), ma.hemisphere_acentrics("-")),
                ma.hemispheres_acentrics()):
    assert tuple(hp.indices()) == ((1,2,3), (2,3,4))
    assert approx_equal(tuple(hp.data()), (1,5))
    assert approx_equal(tuple(hp.sigmas()), (0.1,0.3))
    assert tuple(hm.indices()) == ((-1,-2,-3), (-2,-3,-4))
    assert approx_equal(tuple(hm.data()), (2,3))
    assert approx_equal(tuple(hm.sigmas()), (0.2,0.4))
  assert approx_equal(ma.anomalous_signal(), 0.5063697)
  assert ma.measurability() == 1.0
  sigma_test_meas = flex.double(5,0.3)
  m_test_meas = miller.array(ms, data, sigma_test_meas )
  assert m_test_meas.measurability() == 0.5
  ms = miller.set(crystal.symmetry(), mi, anomalous_flag=True)
  ma = miller.array(ms, data, sigmas)
  ad = ma.anomalous_differences()
  assert tuple(ad.indices()) == ((1,2,3), (2,3,4))
  for hp,hm in ((ma.hemisphere_acentrics("+"), ma.hemisphere_acentrics("-")),
                ma.hemispheres_acentrics()):
    assert tuple(hp.indices()) == ((1,2,3), (2,3,4))
    assert approx_equal(tuple(hp.data()), (1,5))
    assert approx_equal(tuple(hp.sigmas()), (0.1,0.3))
    assert tuple(hm.indices()) == ((-1,-2,-3), (-2,-3,-4))
    assert approx_equal(tuple(hm.data()), (2,3))
    assert approx_equal(tuple(hm.sigmas()), (0.2,0.4))
  assert approx_equal(ma.anomalous_signal(), 0.5063697)
  assert tuple(ma.all_selection()) == (1,1,1,1,1)
  for sa in (ma.select(flex.bool((1,0,0,1,0))),
             ma.select(flex.size_t((0,3)))):
    assert tuple(sa.indices()) == ((1,2,3), (-2,-3,-4))
    assert approx_equal(tuple(sa.data()), (1,3))
    assert approx_equal(tuple(sa.sigmas()), (0.1,0.4))
  ms = miller.build_set(xs, anomalous_flag=False, d_min=1)
  ma = miller.array(ms)
  sa = ma.resolution_filter()
  assert ma.indices().size() == sa.indices().size()
  sa = ma.resolution_filter(0.5)
  assert sa.indices().size() == 0
  sa = ma.resolution_filter(d_min=2)
  assert sa.indices().size() == 10
  sa = ma.resolution_filter(d_min=2, negate=True)
  assert sa.indices().size() == 38
  ma = ma.d_spacings()
  ma = miller.array(ma, ma.data(), ma.data().deep_copy())
  assert ma.indices().size() == 48
  sa = ma.sigma_filter(0.5)
  assert sa.indices().size() == 48
  sa = ma.sigma_filter(2)
  assert sa.indices().size() == 0
  for i in (1,13,25,27,39):
    ma.sigmas()[i] /= 3
  sa = ma.sigma_filter(2)
  assert sa.indices().size() == 5
  assert approx_equal(ma.mean(False,False), 1.6460739)
  assert approx_equal(ma.mean(False,True), 1.5146784)
  ma.setup_binner(n_bins=3)
  assert approx_equal(ma.mean(True,False).data[1:-1],
    (2.228192, 1.2579831, 1.0639812))
  assert approx_equal(ma.mean(True,True).data[1:-1],
    (2.069884, 1.2587977, 1.0779636))
  assert approx_equal(ma.mean_sq(False,False), 3.3287521)
  assert approx_equal(ma.mean_sq(False,True), 2.6666536)
  assert approx_equal(ma.mean_sq(True,False).data[1:-1],
    (5.760794, 1.5889009, 1.1336907))
  assert approx_equal(ma.mean_sq(True,True).data[1:-1],
    (4.805354, 1.5916849, 1.1629777))
  assert approx_equal(ma.rms(False,False)**2, 3.3287521)
  assert approx_equal(ma.rms(False,True)**2, 2.6666536)
  assert approx_equal([x**2 for x in ma.rms(True,False).data[1:-1]],
    ma.mean_sq(True,False).data[1:-1])
  assert approx_equal([x**2 for x in ma.rms(True,True).data[1:-1]],
    ma.mean_sq(True,True).data[1:-1])
  for use_binning in (False,True):
    for use_multiplicities in (False,True):
      sa = ma.rms_filter(-1, use_binning, use_multiplicities)
      assert sa.indices().size() == 0
      sa = ma.rms_filter(100, use_binning, use_multiplicities, False)
      assert sa.indices().size() == ma.indices().size()
      sa = ma.rms_filter(-1, use_binning, use_multiplicities, negate=True)
      assert sa.indices().size() == ma.indices().size()
      sa = ma.rms_filter(100, use_binning, use_multiplicities, negate=True)
      assert sa.indices().size() == 0
      sa = ma.rms_filter(1.0, use_binning, use_multiplicities)
      assert sa.indices().size() \
          == ((36, 33), (29, 29))[use_binning][use_multiplicities]
  assert approx_equal(ma.statistical_mean(), 1.380312)
  assert approx_equal(tuple(ma.statistical_mean(True)),
                      (1.768026, 1.208446, 0.9950434))
  no = ma.remove_patterson_origin_peak()
  assert approx_equal(no.data()[0], 3.231974)
  assert approx_equal(no.data()[47], 0.004956642)
  no = ma.quasi_normalize_structure_factors(d_star_power=0)
  assert approx_equal(no.data()[0], 2.4378468)
  assert approx_equal(no.data()[47], 0.9888979)
  no = ma.quasi_normalize_structure_factors()
  assert approx_equal(no.data()[0], 2.00753806261)
  assert approx_equal(no.data()[47], 1.09976342511)
  su = ma + 3
  assert approx_equal(tuple(su.data()), tuple(ma.data() + 3))
  su = ma + ma
  assert approx_equal(tuple(su.data()), tuple(ma.data() * 2))
  assert approx_equal(tuple(su.sigmas()), tuple(ma.sigmas() * math.sqrt(2)))
  s = ma.f_as_f_sq()
  v = s.f_sq_as_f()
  assert approx_equal(tuple(ma.data()), tuple(v.data()))
  assert not approx_equal(tuple(ma.sigmas()), tuple(v.sigmas()))
  s = miller.array(ma, ma.data()).f_as_f_sq()
  v = s.f_sq_as_f()
  assert approx_equal(tuple(ma.data()), tuple(v.data()))
  assert s.sigmas() is None
  assert v.sigmas() is None
  ma = miller.array(ms)
  s = ma[:]
  assert s.data() is None
  assert s.sigmas() is None
  ma = miller.array(ms, flex.double((1,2)))
  s = ma[:]
  assert s.data().all_eq(ma.data())
  assert s.sigmas() is None
  ma = miller.array(ms, flex.double((1,2)), flex.double((3,4)))
  s = ma[:]
  assert s.data().all_eq(ma.data())
  assert s.sigmas().all_eq(ma.sigmas())
  xs = crystal.symmetry((3,4,5), "P 1 1 21")
  mi = flex.miller_index(((0,0,1), (0,0,2), (0,0,-3), (0,0,-4)))
  ms = miller.set(xs, mi)
  ma = miller.array(ms).remove_systematic_absences()
  assert tuple(ma.indices()) == ((0,0,2), (0,0,-4))
  ma = miller.array(ms).remove_systematic_absences(negate=True)
  assert tuple(ma.indices()) == ((0,0,1), (0,0,-3))
  ma = miller.array(ms, flex.double((3,4,1,-2)), flex.double((.3,.4,.1,.2)))
  sa = ma.sort(by_value="resolution")
  assert tuple(sa.indices()) == ((0,0,1), (0,0,2), (0,0,-3), (0,0,-4))
  assert approx_equal(sa.data(), (3,4,1,-2))
  assert approx_equal(sa.sigmas(), (.3,.4,.1,.2))
  sa = ma.sort(by_value="resolution", reverse=True)
  assert tuple(sa.indices()) == ((0,0,-4), (0,0,-3), (0,0,2), (0,0,1))
  assert approx_equal(sa.data(), (-2,1,4,3))
  assert approx_equal(sa.sigmas(), (.2,.1,.4,.3))
  sa = ma.sort(by_value="packed_indices", reverse=True)
  assert tuple(sa.indices()) == ((0,0,2), (0,0,1), (0,0,-3), (0,0,-4))
  sa = ma.sort(by_value="data")
  assert approx_equal(sa.data(), (4,3,1,-2))
  sa = ma.sort(by_value="data", reverse=True)
  assert approx_equal(sa.data(), (-2,1,3,4))
  sa = ma.sort(by_value="abs")
  assert approx_equal(sa.data(), (4,3,-2,1))
  sa = ma.sort(by_value="abs", reverse=True)
  assert approx_equal(sa.data(), (1,-2,3,4))
  sa = ma.sort(by_value=flex.double((3,1,4,2)))
  assert tuple(sa.indices()) == ((0,0,-3), (0,0,1), (0,0,-4), (0,0,2))
  sa = ma.sort(by_value=flex.double((3,1,4,2)), reverse=True)
  assert tuple(sa.indices()) == ((0,0,2), (0,0,-4), (0,0,1), (0,0,-3))
  aa = sa.adopt_set(ma)
  assert tuple(aa.indices()) == tuple(ma.indices())
  assert approx_equal(aa.data(), ma.data())
  assert approx_equal(aa.sigmas(), ma.sigmas())
  sa = ma.apply_scaling(target_max=10)
  assert approx_equal(flex.max(sa.data()), 10)
  assert approx_equal(flex.max(sa.sigmas()), 1)
  sa = sa.apply_scaling(factor=3)
  assert approx_equal(flex.max(sa.data()), 30)
  assert approx_equal(flex.max(sa.sigmas()), 3)
  ma = miller.array(miller.set(xs, mi, False),data,sigmas).patterson_symmetry()
  assert str(ma.space_group_info()) == "P 1 1 2/m"
  assert ma.indices() == mi
  assert ma.data() == data
  assert ma.sigmas() == sigmas
  a1 = miller.array(
    miller.set(xs, flex.miller_index(((1,-2,3), (0,0,-4)))),
    flex.double((1,2)))
  a2 = miller.array(
    miller.set(xs, flex.miller_index(((0,0,-5), (1,-2,3)))),
    flex.double((3,4)),
    flex.double((5,6)))
  c1 = a1.common_set(a2)
  assert tuple(c1.indices()) == ((1,-2,3),)
  assert tuple(c1.data()) == (1,)
  c2 = a2.common_set(a1)
  assert tuple(c2.indices()) == ((1,-2,3),)
  assert tuple(c2.data()) == (4,)
  assert tuple(c2.sigmas()) == (6,)
  l1 = a1.lone_set(a2)
  assert list(l1.indices()) == [(0,0,-4)]
  assert list(l1.data()) == [2]
  l2 = a2.lone_set(a1)
  assert list(l2.indices()) == [(0,0,-5)]
  assert list(l2.data()) == [3]
  l1, l2 = a1.lone_sets(a2)
  assert list(l1.indices()) == [(0,0,-4)]
  assert list(l2.indices()) == [(0,0,-5)]
  assert tuple(c1.adopt_set(c2).indices()) == ((1,-2,3),)
  sg = miller.array(
    miller.set(xs, flex.miller_index(((0,0,-5), (1,-2,3))), False),
    flex.double((3,4)))
  p1 = sg.expand_to_p1()
  assert p1.indices().size() == 3
  assert approx_equal(tuple(p1.data()), (3,4,4))
  assert p1.sigmas() is None
  sg = miller.array(
    miller.set(xs, flex.miller_index(((0,0,-5), (1,-2,3))), False),
    flex.double((3,4)),
    flex.double((5,6)))
  p1 = sg.expand_to_p1()
  assert p1.indices().size() == 3
  assert approx_equal(tuple(p1.data()), (3,4,4))
  assert approx_equal(tuple(p1.sigmas()), (5,6,6))
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  data = flex.double((1,2))
  a = miller.array(miller.set(xs, mi), data)
  ph = flex.double((10,20))
  b = a.phase_transfer(ph, deg=True)
  assert approx_equal(tuple(b.amplitudes().data()), a.data())
  assert approx_equal(tuple(b.phases(deg=True).data()), (10,0))
  ph = ph * math.pi/180
  b = a.phase_transfer(ph, deg=False)
  assert approx_equal(tuple(b.amplitudes().data()), a.data())
  assert approx_equal(tuple(b.phases(deg=True).data()), (10,0))
  c = a.phase_transfer(b.data())
  assert approx_equal(tuple(c.amplitudes().data()), a.data())
  assert approx_equal(tuple(c.phases(deg=True).data()), (10,0))
  a = miller.array(miller_set=a, data=flex.complex_double([1+2j,2-3j]))
  c = a.conjugate()
  assert approx_equal(a.data(), flex.conj(c.data()))
  ms = miller_set=miller.build_set(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(11,11,13,90,90,120),
      space_group_symbol="P31m"),
    anomalous_flag=True,
    d_min=3)
  ma = miller.array(
    miller_set=ms,
    data=flex.double(xrange(ms.indices().size())))
  ma.set_observation_type_xray_amplitude()
  assert ma.select_acentric().data().size() == 48
  assert ma.select_centric().data().size() == 3
  ma.setup_binner(auto_binning=True).counts_complete(
    include_acentric=False,
    include_centric=False)
  assert ma.binner().counts_complete() == [0]*10
  for i_trial in [0,1]:
    if (i_trial == 1): ma = ma.f_as_f_sq()
    assert approx_equal(ma.second_moment_of_intensities(), 1.81758415842)
    assert approx_equal(ma.wilson_ratio(), 0.742574257426)
    ma.setup_binner(auto_binning=True).counts_complete(include_centric=False)
    s = StringIO()
    ma.second_moment_of_intensities(use_binning=True).show(f=s)
    assert s.getvalue() == """\
unused:         - 13.0001 [ 0/0 ]
bin  1: 13.0001 -  5.9727 [ 7/6 ]  2.1658
bin  2:  5.9727 -  4.8197 [ 8/8 ]  1.1899
bin  3:  4.8197 -  4.2345 [ 5/4 ]  1.6278
bin  4:  4.2345 -  3.8585 [ 6/6 ]  1.2952
bin  5:  3.8585 -  3.5881 [ 4/4 ]  1.0429
bin  6:  3.5881 -  3.3805 [ 8/8 ]  1.2980
bin  7:  3.3805 -  3.2139 [ 2/2 ]  1.0234
bin  8:  3.2139 -  3.0759 [11/10]  1.2283
unused:  3.0759 -         [ 0/0 ]
"""
    s = StringIO()
    ma.wilson_ratio(use_binning=True).show(f=s)
    assert s.getvalue() == """\
unused:         - 13.0001 [ 0/0 ]
bin  1: 13.0001 -  5.9727 [ 7/6 ]  0.6007
bin  2:  5.9727 -  4.8197 [ 8/8 ]  0.9349
bin  3:  4.8197 -  4.2345 [ 5/4 ]  0.7079
bin  4:  4.2345 -  3.8585 [ 6/6 ]  0.9247
bin  5:  3.8585 -  3.5881 [ 4/4 ]  0.9892
bin  6:  3.5881 -  3.3805 [ 8/8 ]  0.9079
bin  7:  3.3805 -  3.2139 [ 2/2 ]  0.9941
bin  8:  3.2139 -  3.0759 [11/10]  0.9134
unused:  3.0759 -         [ 0/0 ]
"""
  assert ma.binner() is not None
  try:
    pickle.dumps(ma)
  except RuntimeError, e:
    assert str(e).startswith(
      "cctbx.miller.binner instances are not picklable.")
  else: raise RuntimeError("Exception expected.")
  ma.clear_binner()
  ml = pickle.loads(pickle.dumps(ma))
  sa = StringIO()
  sl = StringIO()
  ma.show_summary(f=sa).show_array(f=sa)
  ml.show_summary(f=sl).show_array(f=sl)
  assert sa.getvalue() == sl.getvalue()

def exercise_array_2(space_group_info):
  xs = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(60),
    space_group_info=space_group_info)
  for anomalous_flag in (False, True):
    st = miller.build_set(xs, anomalous_flag, d_min=1)
    for sigmas in (None, flex.double(xrange(1,st.indices().size()+1))):
      sg = miller.array(
        st,
        data=flex.double(xrange(st.indices().size())),
        sigmas=sigmas)
      p1 = sg.expand_to_p1()
      ps = miller.array(
        miller.set(xs, p1.indices(), p1.anomalous_flag()),
        p1.data(),
        p1.sigmas())
      m = ps.merge_equivalents()
      p = m.array().sort_permutation(by_value="data", reverse=True)
      assert flex.order(sg.indices(), m.array().indices().select(p)) == 0
      assert approx_equal(sg.data(), m.array().data().select(p))
      if (sigmas is not None):
        s = m.array().sigmas().select(p)
        r = m.redundancies().select(p)
        sr = s * flex.sqrt(r.as_double())
        assert approx_equal(sr, sigmas)
      #
      orig = m.array()
      if (not orig.anomalous_flag()):
        ave = orig.average_bijvoet_mates()
        assert ave.correlation(orig).coefficient() > 1-1.e-6
        if (sigmas is not None):
          assert flex.linear_correlation(
            ave.sigmas(), orig.sigmas()).coefficient() > 1-1.e-6
      else:
        if (sigmas is not None):
          # merge_equivalents uses the sigmas as weights
          orig = orig.customized_copy(
            sigmas=flex.double(orig.sigmas().size(),1))
        ave = orig.average_bijvoet_mates()
        asu, matches = orig.match_bijvoet_mates()
        vfy_indices = asu.indices().select(matches.pairs().column(0))
        vfy_data = (  asu.data().select(matches.pairs().column(0))
                    + asu.data().select(matches.pairs().column(1))) / 2
        for sign in ("+", "-"):
          sel = matches.singles(sign)
          vfy_indices.extend(asu.indices().select(sel))
          vfy_data.extend(asu.data().select(sel))
        vfy = miller.array(
          miller_set=miller.set(
            crystal_symmetry=orig,
            indices=vfy_indices,
            anomalous_flag=False),
          data=vfy_data).adopt_set(ave)
        assert vfy.correlation(ave).coefficient() > 1-1.e-6

def exercise_fft_map():
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  for anomalous_flag in (False, True):
    ms = miller.set(xs, mi, anomalous_flag=anomalous_flag)
    ma = miller.array(ms, flex.complex_double((1,2)))
    fft_map = ma.fft_map()
    assert approx_equal(fft_map.resolution_factor(), 1./3)
    assert fft_map.symmetry_flags() is None
    assert approx_equal(fft_map.max_prime(), 5)
    assert fft_map.anomalous_flag() == anomalous_flag
    assert fft_map.real_map().size() > 0
    assert not fft_map.real_map_unpadded().is_padded()
    if (anomalous_flag):
      assert fft_map.complex_map().size() > 0

def exercise_squaring_and_patterson_map(space_group_info,
                                        n_scatterers=8,
                                        d_min=2,
                                        verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=500,
    min_distance=5.,
    general_positions_only=True,
    u_iso=0.0)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  e_000 = math.sqrt(n_scatterers * structure.space_group().order_z())
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  f_calc = f_calc.sort(by_value="abs")
  f = abs(f_calc)
  assert approx_equal(f.data(), f_calc.amplitudes().data())
  assert approx_equal(f_calc.phases(deg=True).data()*math.pi/180,
                      f_calc.phases().data())
  f.setup_binner(auto_binning=True)
  e = f.quasi_normalize_structure_factors()
  grid_resolution_factor = 1/3.
  u_base = xray.calc_u_base(d_min, grid_resolution_factor)
  if (0 or verbose):
    print "u_base:", u_base
  d_star_sq = e.unit_cell().d_star_sq(e.indices())
  dw = flex.exp(d_star_sq*2*(math.pi**2)*u_base)
  eb = miller.array(miller_set=e, data=e.data()/dw)
  eb_map = eb.phase_transfer(f_calc).fft_map(
    resolution_factor=grid_resolution_factor,
    d_min=d_min,
    f_000=e_000).real_map()
  eb_map_sq = flex.pow2(eb_map)
  eb_sq = eb.structure_factors_from_map(eb_map_sq)
  mwpe = f_calc.mean_weighted_phase_error(eb_sq)
  if (0 or verbose):
    print "mean_weighted_phase_error: %.2f" % mwpe
  assert mwpe < 2
  for sharpening in (False, True):
    for origin_peak_removal in (False, True):
      patterson_map = eb.patterson_map(
        symmetry_flags=maptbx.use_space_group_symmetry,
        resolution_factor=grid_resolution_factor,
        f_000=e_000,
        sharpening=sharpening,
        origin_peak_removal=origin_peak_removal)
      grid_tags = maptbx.grid_tags(patterson_map.n_real())
      grid_tags.build(
        patterson_map.space_group_info().type(),
        maptbx.use_space_group_symmetry)
      assert grid_tags.n_grid_misses() == 0
      assert grid_tags.verify(patterson_map.real_map())

def exercise_common_set((a, b), permutation_only):
  ab = a.common_set(b)
  ba = b.common_set(a)
  if (permutation_only):
    assert ab.indices().all_eq(b.indices())
    assert ba.indices().all_eq(a.indices())
  aab,bab = a.common_sets(b)
  assert aab.indices().all_eq(bab.indices())
  assert aab.indices().all_eq(ab.indices())
  aba,bba = b.common_sets(a)
  assert aba.indices().all_eq(bba.indices())
  assert aba.indices().all_eq(ba.indices())

def exercise_array_correlation(space_group_info,
                               n_scatterers=8,
                               d_min=2,
                               verbose=0):
  arrays = []
  for i in xrange(2):
    structure = random_structure.xray_structure(
      space_group_info,
      elements=["const"]*n_scatterers)
    arrays.append(abs(structure.structure_factors(d_min=d_min-i*0.5).f_calc()))
  arrays[1] = arrays[1].select(
    flex.random_permutation(size=arrays[1].indices().size()))
  exercise_common_set((arrays[0], arrays[0].select(
      flex.random_permutation(size=arrays[0].indices().size()))),
    permutation_only=True)
  exercise_common_set(arrays, permutation_only=False)
  for anomalous_flag in [False, True]:
    if (anomalous_flag):
      arrays[0] = arrays[0].as_anomalous()
    assert approx_equal(arrays[0].correlation(arrays[0]).coefficient(), 1)
    assert approx_equal(arrays[0].correlation(arrays[1]).coefficient(),
                        arrays[1].correlation(arrays[0]).coefficient())
    arrays[0].setup_binner(auto_binning=True)
    arrays[1].use_binning_of(arrays[0])
    for corr in arrays[0].correlation(arrays[0], use_binning=True).data:
      if (corr is not None):
        assert approx_equal(corr, 1)
    corr0 = arrays[0].correlation(arrays[1], use_binning=True).data
    corr1 = arrays[1].correlation(arrays[0], use_binning=True).data
    for c0,c1 in zip(corr0,corr1):
      assert (c0 is None) == (c1 is None)
      if (c0 is None): continue
      assert approx_equal(c0, c1)

def exercise_as_hendrickson_lattman(space_group_info, n_scatterers=5, d_min=3,
                                    verbose=0):
  phase_integrator = miller.phase_integrator()
  phase_restriction = space_group_info.group().phase_restriction
  for anomalous_flag in [False, True]:
    structure = random_structure.xray_structure(
      space_group_info,
      elements=["const"]*n_scatterers,
      volume_per_atom=200,
      random_f_double_prime=False)
    f_calc = structure.structure_factors(
      d_min=d_min,
      anomalous_flag=anomalous_flag,
      algorithm="direct").f_calc()
    max_f_calc = flex.max(flex.abs(f_calc.data()))
    phase_integrals = f_calc.data() / (max_f_calc/.95)
    for h,pi_calc in zip(f_calc.indices(), phase_integrals):
      phase_info = phase_restriction(h)
      assert abs(pi_calc - phase_info.valid_structure_factor(pi_calc)) < 1.e-6
      hl = miller.as_hendrickson_lattman(
        centric_flag=phase_info.is_centric(),
        phase_integral=pi_calc,
        max_figure_of_merit=1-1.e-6)
      pi_int = phase_integrator(phase_info, hl)
      assert abs(pi_calc - pi_int) < 1.e-1

def one_random_hl(f, min_coeff=1.e-3):
  result = f * (random.random() - 0.5)
  if (result < 0):
    if (result > -min_coeff): return min_coeff
  else:
    if (result < min_coeff): return min_coeff
  return result

def generate_random_hl(miller_set, coeff_range=5):
  phase_restriction = miller_set.space_group().phase_restriction
  hl = flex.hendrickson_lattman()
  for h in miller_set.indices():
    phase_info = phase_restriction(h)
    if (phase_info.is_centric()):
      fom = max(0.01, random.random()*0.95)
      if (random.random() < 0.5): fom *= -1
      angle = phase_info.ht_angle()
      f = fom * complex(math.cos(angle), math.sin(angle))
      assert abs(f - phase_info.valid_structure_factor(f)) < 1.e-6
      hl.append(cctbx.hendrickson_lattman(
        centric_flag=True,
        phase_integral=f,
        max_figure_of_merit=1-1.e-6))
    else:
      f = 2 * coeff_range * random.random()
      hl.append([one_random_hl(f) for i in xrange(4)])
  return miller.array(miller_set=miller_set, data=hl)

def exercise_phase_integrals(space_group_info):
  crystal_symmetry = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(
      volume=250*space_group_info.group().n_ltr()),
    space_group_info=space_group_info)
  is_centric = space_group_info.group().is_centric
  for anomalous_flag in [False, True]:
    miller_set = miller.build_set(
      crystal_symmetry=crystal_symmetry,
      anomalous_flag=anomalous_flag,
      d_min=1)
    sg_hl = generate_random_hl(miller_set=miller_set)
    p1_hl = sg_hl.expand_to_p1()
    sg_phase_integrals = sg_hl.phase_integrals(n_steps=360/5)
    p1_phase_integrals = p1_hl.phase_integrals()
    p1_sg_phase_integrals = sg_phase_integrals.expand_to_p1()
    assert p1_sg_phase_integrals.indices().all_eq(p1_phase_integrals.indices())
    for h,pi_p1,pi_p1_sg in zip(p1_phase_integrals.indices(),
                                p1_phase_integrals.data(),
                                p1_sg_phase_integrals.data()):
      if (is_centric(h)):
        if (scitbx.math.phase_error(complex_math.arg(pi_p1),
                                    complex_math.arg(pi_p1_sg)) > 1.e-6):
          print "Error:", h, pi_p1, pi_p1_sg
          print "arg(pi_p1):", complex_math.arg(pi_p1)
          print "arg(pi_p1_sg):", complex_math.arg(pi_p1_sg)
          raise AssertionError
        if (not (0.5 < abs(pi_p1)/abs(pi_p1_sg) < 0.75)):
          print "Error:", h, pi_p1, pi_p1_sg
          print "abs(pi_p1):", abs(pi_p1)
          print "abs(pi_p1_sg):", abs(pi_p1_sg)
          raise AssertionError
      elif (abs(pi_p1 - pi_p1_sg) > 1.e-6):
        print "Error:", h, pi_p1, pi_p1_sg
        raise AssertionError
    #
    amplitude_array = miller.array(
      miller_set=miller_set,
      data=flex.random_double(size=miller_set.indices().size()))
    for phase_integrator_n_steps in [None, 360/5]:
      with_phases = amplitude_array.phase_transfer(
        phase_source=sg_hl,
        phase_integrator_n_steps=phase_integrator_n_steps)
      assert flex.max(  flex.abs(amplitude_array.data()
                      - abs(with_phases).data())) < 1.e-6
      assert with_phases.mean_weighted_phase_error(
        phase_source=sg_phase_integrals) < 1.e-6

def run_call_back(flags, space_group_info):
  exercise_array_2(space_group_info)
  exercise_squaring_and_patterson_map(space_group_info, verbose=flags.Verbose)
  exercise_array_correlation(space_group_info)
  exercise_as_hendrickson_lattman(space_group_info)
  exercise_phase_integrals(space_group_info)

def run():
  exercise_set()
  exercise_generate_r_free_flags()
  exercise_binner()
  exercise_array()
  exercise_crystal_gridding()
  exercise_fft_map()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
