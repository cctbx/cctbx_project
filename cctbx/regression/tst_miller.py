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
from libtbx import complex_math
from libtbx.test_utils import \
  approx_equal, not_approx_equal, show_diff, Exception_expected
from libtbx import Auto
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
  assert approx_equal(ms.d_star_cubed().data(), [0.612355017974, 0.512])
  assert approx_equal(tuple(ms.d_spacings().data()), (1.177603, 1.25))
  assert approx_equal(tuple(ms.sin_theta_over_lambda_sq().data()),
                      (0.1802778, 0.16))
  assert approx_equal(ms.d_min(), 1.177603)
  assert approx_equal(ms.d_max_min(), [1.25,1.177603])
  assert approx_equal(ms.resolution_range(), (1.25, 1.177603))
  assert approx_equal(ms.minimum_wavelength_based_on_d_min(), 2.33165395578)
  assert approx_equal(ms.minimum_wavelength_based_on_d_min(tolerance=0.1),
    2.11968541434)
  p1 = ms.expand_to_p1()
  assert p1.indices().size() == 6
  p1b, iselection = ms.expand_to_p1(return_iselection=True)
  assert p1b.indices().all_eq(p1.indices())
  assert list(iselection) == [0, 0, 0, 0, 1, 1]
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
  assert mp.min_max_indices() == ((0, 0, 3), (1, 2, 4))
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
  assert not show_diff(s.getvalue(), """\
Number of Miller indices: 36
Anomalous flag: %s
Unit cell: (3, 4, 5, 90, 90, 90)
Space group: P 2 2 2 (No. 16)
Systematic absences: 0
Centric reflections: 27
Resolution range: 5 1.1776
Completeness in resolution range: 1
Completeness with d_max=infinity: 1
""" % str(False))
  s = StringIO()
  mc.select(flex.size_t([0])).show_comprehensive_summary(f=s)
  assert not show_diff(s.getvalue(), """\
Number of Miller indices: 1
Anomalous flag: %s
Unit cell: (3, 4, 5, 90, 90, 90)
Space group: P 2 2 2 (No. 16)
Systematic absences: 0
Centric reflections: 1
Resolution range: 5 5
Completeness with d_max=infinity: 1
""" % str(False))
  #
  c = mc.change_basis(cb_op="k,h,-l")
  assert c.unit_cell().is_similar_to(uctbx.unit_cell((4,3,5,90,90,90)))
  #
  s = mc.customized_copy(space_group_info=sgtbx.space_group_info(number=19))
  assert s.sys_absent_flags().data().all_eq(flex.bool([
    True, False, True, False, True, False, False, False, False,
    False, False, False, False, True, False, True, False, False,
    False, False, False, False, False, False, False, False, False,
    False, False, False, False, False, False, False, False, False]))
  assert s.sys_absent_flags(integral_only=True).data().all_eq(False)
  r = s.reflection_intensity_symmetry()
  assert str(r.space_group_info()) == "P m m m"
  assert r.sys_absent_flags().data().all_eq(False)
  #
  s = s.customized_copy(indices=flex.miller_index([(1,2,3),(0,0,1)]))
  for deg in [False, True]:
    p = s.random_phases_compatible_with_phase_restrictions(deg=deg).data()
    assert s.space_group().is_valid_phase((0,0,1), p[1], deg)
  #
  s = miller.set(
    crystal_symmetry=xs,
    indices=flex.miller_index(((1,2,3), (0,0,0))),
    anomalous_flag=False)
  r = s.resolution_filter_selection
  assert list(r()) == [True, True]
  assert list(r(d_max=0, d_min=0)) == [True, True]
  assert list(r(d_max=1.1)) == [False, False]
  assert list(r(d_max=1.2)) == [True, False]
  assert list(r(d_min=1.1)) == [True, True]
  assert list(r(d_min=1.2)) == [False, True]
  assert list(r(d_max=1.2, d_min=1.1)) == [True, False]
  assert list(r(d_max=1.1, d_min=1.2)) == [False, False]
  #
  s = xs.miller_set(indices=flex.miller_index([(1,2,3)]), anomalous_flag=True)
  assert s.indices().size() == 1
  assert s.anomalous_flag()

def exercise_union_of_sets():
  xs = crystal.symmetry(unit_cell=(10,10,10,90,90,90), space_group_symbol="P1")
  ms = miller.set(xs, flex.miller_index())
  u = miller.union_of_sets(miller_sets=[ms])
  assert u.is_similar_symmetry(ms)
  assert len(u.indices()) == 0
  ms = miller.set(xs, flex.miller_index(((1,-2,3), (0,0,-4))))
  u = miller.union_of_sets(miller_sets=[ms])
  assert sorted(u.indices()) == [(0,0,-4), (1,-2,3)]
  u = miller.union_of_sets(miller_sets=[ms, ms])
  assert sorted(u.indices()) == [(0,0,-4), (1,-2,3)]
  ms2 = miller.set(xs, flex.miller_index(((1,-2,3), (0,0,-5))))
  u = miller.union_of_sets(miller_sets=[ms, ms2])
  assert sorted(u.indices()) == [(0,0,-5), (0,0,-4), (1,-2,3)]

def exercise_enforce_positive_amplitudes():
  from cctbx.xray import observation_types
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  data = flex.double((-1,-2))
  sigmas = flex.double((1,2))
  ms = miller.set(xs, mi)
  ma = miller.array(ms)
  ma = miller.array(ms, data=data, sigmas=sigmas).set_observation_type(
    observation_types.intensity() )
  new_ma = ma.enforce_positive_amplitudes()
  assert new_ma.data()[0]>0
  assert new_ma.data()[1]>0

def exercise_generate_r_free_flag_on_lat_sym(sg_info):
  for an_flag in [True,False]:
    full_xs = sg_info.any_compatible_crystal_symmetry(volume=50*80*100)
    low_xs = crystal.symmetry( unit_cell = full_xs.unit_cell(),
                               space_group = sgtbx.space_group("P1") )
    miller_set = miller.build_set(
      crystal_symmetry=low_xs,
      anomalous_flag=an_flag,
      d_min=3.0 )

    free_flags = miller_set.generate_r_free_flags(use_lattice_symmetry=True)
    free_flags = free_flags.select( free_flags.data() )
    fake_data_in_lat_sym = free_flags.customized_copy(
      crystal_symmetry=full_xs,
      indices = free_flags.indices(),
      data = flex.double(free_flags.indices().size(),2.0 ),
      sigmas = flex.double(free_flags.indices().size(),1.0 )
                                                     )
    fake_data_in_lat_sym = fake_data_in_lat_sym.merge_equivalents().array()
    fake_data_in_p1 = fake_data_in_lat_sym.expand_to_p1()
    assert fake_data_in_p1.indices().size()==free_flags.indices().size()
    # note that this assert will fail (unless you are lucky)
    # if max_delta is set to None

    # also check that the full symmetry doesn't do anything weird
    miller_set = miller.build_set(
      crystal_symmetry=full_xs,
      anomalous_flag=an_flag,
      d_min=3.0 )
    free_flags = miller_set.generate_r_free_flags(use_lattice_symmetry=True)

    # check the generatiuon of integer flags please
    integer_flags = miller_set.generate_r_free_flags_on_lattice_symmetry(
      fraction=None,
      max_free=None,
      return_integer_array=True,
      n_partitions=3
    )
    assert approx_equal(
      float(integer_flags.data().size())/integer_flags.data().count( 0 ),3,eps=0.1 )
    assert approx_equal(
      float(integer_flags.data().size())/integer_flags.data().count( 1 ),3,eps=0.1 )
    assert approx_equal(
      float(integer_flags.data().size())/integer_flags.data().count( 2 ),3,eps=0.1 )
    assert ( integer_flags.data().count( 3 ) == 0 )




def exercise_generate_r_free_flags(verbose=0, use_lattice_symmetry=False):
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
      flags = trial_set.generate_r_free_flags(use_lattice_symmetry=use_lattice_symmetry)
      if (i_trial == 0):
        out = StringIO()
        flags.show_r_free_flags_info(out=out, prefix="$#")
        if (verbose): sys.stdout.write(out.getvalue())
        lines = out.getvalue().splitlines()
        assert len(lines) == 13
        for line in lines:
          assert line.startswith("$#")
        accu = flags.r_free_flags_accumulation()
        assert accu.reflection_counts.size() == accu.free_fractions.size()
        if (verbose):
          print "r_free_flags_accumulation:", \
            zip(accu.reflection_counts, accu.free_fractions)
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
        flags_non_anomalous = flags.average_bijvoet_mates()
        if (i_trial < 5):
          assert flags_non_anomalous.indices().size() == 145
          assert flags_non_anomalous.data().count(True) == 15
          flags_gen = flags_non_anomalous.generate_bijvoet_mates() \
            .adopt_set(flags)
          assert flags_gen.data().all_eq(flags.data())
        flags = flags.select(~flex.bool(
          flags.indices().size(),
          flags.match_bijvoet_mates()[1].pairs_hemisphere_selection("-")))
      flags = flags.sort(by_value="resolution", reverse=True)
      isel = flags.data().iselection()
      flag_distances = isel[1:] - isel[:-1]
      assert flex.min(flag_distances) > 0
      assert flex.max(flag_distances) <= 20
      flags = trial_set.generate_r_free_flags(max_free=10, use_lattice_symmetry=use_lattice_symmetry)
      if (not anomalous_flag):
        assert flags.data().count(True) == 10
      else:
        assert 10 <= flags.data().count(True) <= 20

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
    s = StringIO()
    set1.binner().show_summary(f=s, prefix=" ", bin_range_as="stol")
    assert s.getvalue() == """\
 unused:         -  0.0174 [0/0]
 bin  1:  0.0174 -  0.0354 [3/3]
 bin  2:  0.0354 -  0.0437 [3/3]
 bin  3:  0.0437 -  0.0496 [2/2]
 unused:  0.0496 -         [0/0]
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
    set2.completeness(use_binning=True).show(
      show_bin_number=False, f=s, bin_range_as="two_theta", wavelength=0.71073)
    assert s.getvalue() == """\
        -  1.4180 [0/0]
 1.4180 -  2.8822 [3/3] 1.000
 2.8822 -  3.5579 [3/3] 1.000
 3.5579 -  4.0441 [2/2] 1.000
 4.0441 -         [8/8] 1.000
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
    s = StringIO()
    set2.completeness(use_binning=True).show(
      show_bin_range=False, show_counts=False, f=s)
    assert s.getvalue() == """\
unused:
bin  1: 1.000
bin  2: 1.000
bin  3: 1.000
unused: 1.000
"""
    s = StringIO()
    set2.completeness(use_binning=True).show(
      show_bin_number=False, bin_range_as="d_star_sq", show_counts=False, f=s)
    assert s.getvalue() == """\
        -  0.0012
 0.0012 -  0.0050 1.000
 0.0050 -  0.0076 1.000
 0.0076 -  0.0099 1.000
 0.0099 -         1.000
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
  p = pickle.dumps(set1.binner())
  l = pickle.loads(p)
  s = StringIO()
  set1.binner().show_summary(f=s)
  sl = StringIO()
  l.show_summary(f=sl)
  assert not show_diff(sl.getvalue(), s.getvalue())
  #
  expected_counts = iter([
    [0,1,0], [0,1,0], [0,1,0],
    [0,1,1,0], [0,2,0], [0,2,0], [0,2,0],
    [0,1,1,1,0], [0,2,1,0], [0,3,0], [0,3,0], [0,3,0],
    [0,1,1,1,1,0], [0,2,2,0], [0,4,0], [0,4,0], [0,4,0], [0,4,0],
    [0,1,1,1,1,1,0], [0,2,1,2,0], [0,3,2,0], [0,5,0], [0,5,0], [0,5,0], [0,5,0]
  ])
  for n in xrange(1,6):
    set2 = set1.select(flex.size_t(xrange(n)))
    for reflections_per_bin in range(1,n+1) + [n+1, n*10]:
      set2.setup_binner_counting_sorted(reflections_per_bin=reflections_per_bin)
      assert list(set2.binner().counts()) == expected_counts.next()
  set1.setup_binner_counting_sorted(
    d_max=16,
    d_min=11,
    reflections_per_bin=3)
  assert list(set1.binner().counts()) == [2,3,2,1]

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
  assert approx_equal(xs.unit_cell().parameters(),
                      ma.crystal_symmetry().unit_cell().parameters())
  assert xs.space_group().type().lookup_symbol() == \
         ma.crystal_symmetry().space_group().type().lookup_symbol()
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
  ma.set_info(ma.info().customized_copy(systematic_absences_eliminated=True))
  assert ma.info().label_string() \
      == "a,b,c,merged,systematic_absences_eliminated"
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
  assert aa.is_xray_amplitude_array()
  assert aa.as_amplitude_array() is aa
  ai = aa.as_intensity_array()
  assert ai.is_xray_intensity_array()
  assert approx_equal(ai.data(), ac.data())
  assert ai.as_intensity_array() is ai
  assert aa.eliminate_sys_absent() is aa
  assert ma.sigmas_are_sensible()
  assert ma.sigmas_are_sensible(epsilon=0.15)
  assert ma.sigmas_are_sensible(epsilon=0.15, critical_ratio=0.51)
  assert not ma.sigmas_are_sensible(epsilon=0.15, critical_ratio=0.49)
  assert not ma.sigmas_are_sensible(epsilon=0.25)
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
  a_s = aa.select_sys_absent()
  assert a_s is not aa
  assert tuple(a_s.indices()) == ((1,-2,3),)
  s = StringIO()
  aa.eliminate_sys_absent(log=s, prefix="%^")
  assert not show_diff(s.getvalue(), """\
%^Removing 1 systematic absence:
%^  Average absolute value of:
%^    Absences: 1
%^      Others: 1.41421
%^       Ratio: 0.707107

""")
  s = StringIO()
  aa.eliminate_sys_absent(integral_only=True, log=s, prefix="%^")
  assert not show_diff(s.getvalue(), """\
%^Removing 1 integral systematic absence:
%^  Average absolute value of:
%^    Absences: 1
%^      Others: 1.41421
%^       Ratio: 0.707107

""")
  aa = aa.customized_copy(data=flex.complex_double([1+1j,2-2j]))
  s = StringIO()
  aa.eliminate_sys_absent(log=s, prefix="%^")
  assert not show_diff(s.getvalue(), """\
%^Removing 1 systematic absence:
%^  Average absolute value of:
%^    Absences: 1.41421
%^      Others: 2.82843
%^       Ratio: 0.5

""")
  aa = aa.customized_copy(data=flex.int([3,7]), anomalous_flag=True) \
    .expand_to_p1().customized_copy(crystal_symmetry=aa)
  s = StringIO()
  aa.eliminate_sys_absent(log=s, prefix="%^")
  assert not show_diff(s.getvalue(), """\
%^Removing 4 systematic absences.

""")
  #
  for deg in [None, False, True]:
    asu = ma.map_to_asu(deg=deg)
    assert tuple(asu.indices()) == ((1,2,3), (0,0,4))
    if (deg is None):
      assert asu.data().all_eq(ma.data())
    else:
      assert approx_equal(asu.data(), [-1, 2])
  #
  mi = flex.miller_index(((1,2,3), (-1,-2,-3), (2,3,4), (-2,-3,-4), (3,4,5)))
  data = flex.double((1,2,5,3,6))
  sigmas = flex.double((0.1,0.2,0.3,0.4,0.5))
  ms = miller.set(xs, mi, anomalous_flag=True)
  ma = miller.array(ms, data, sigmas)
  assert ma.set().__class__.__name__ == "set"
  assert ma.set().space_group_info() is ma.space_group_info()
  assert ma.set().indices() is ma.indices()
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
  assert tuple(ma.all_selection()) == (True,True,True,True,True)
  for sa in (ma.select(flex.bool((True,False,False,True,False))),
             ma.select(flex.size_t((0,3)))):
    assert tuple(sa.indices()) == ((1,2,3), (-2,-3,-4))
    assert approx_equal(tuple(sa.data()), (1,3))
    assert approx_equal(tuple(sa.sigmas()), (0.1,0.4))
  ms = miller.build_set(xs, anomalous_flag=False, d_min=1)
  assert ms.is_unique_set_under_symmetry()
  assert ms.unique_under_symmetry() is ms
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
  #
  indices=flex.miller_index(((-3,0,0),(2,2,1),(0,3,3)))
  sa = ma.select_indices(indices=indices).sort(
    by_value="packed_indices")
  assert sa.indices().size() == 2
  assert approx_equal(sa.indices(), ((0, 3, 3),(2, 2, 1)))
  sa = ma.select_indices(indices=indices, map_indices_to_asu=True).sort(
    by_value="packed_indices")
  assert sa.indices().size() == 3
  assert approx_equal(sa.indices(), ((0, 3, 3),(2, 2, 1),(3,0,0)))
  sa = ma.select_indices(indices=indices, negate=True)
  assert sa.indices().size() == 46
  sa = ma.select_indices(indices=indices, map_indices_to_asu=True,
                         negate=True)
  assert sa.indices().size() == 45
  for index in indices:
    assert index not in sa.indices()
  ms2 = miller.set(xs, flex.miller_index(((0,0,1),(0,0,1),(0,0,2))))
  ma2 = miller.array(ms2, flex.double((4,2,3)), sigmas=flex.double((0.4,0.1,0.6)))
  indices=flex.miller_index(((0,0,1),))
  try:
    sa = ma2.select_indices(indices=indices)
  except RuntimeError, e:
    assert not show_diff(
      str(e),
      ("cctbx.miller.array.select_indices(): "
       "This method can only be used reliably on a merged array"))
  else: raise Exception_expected
  #
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
  assert approx_equal(ma.sum(False,False), 79.01154935280627)
  assert approx_equal(ma.sum(False,True), 72.704563158141255)
  assert approx_equal(ma.sum(True,False).data[1:-1],
    (46.792023266563191, 22.643695504219149, 9.5758305820238991))
  assert approx_equal(ma.sum(True,True).data[1:-1],
    (43.467564215172203, 22.658357766449797, 9.7016727390599851))
  assert approx_equal(ma.sum_sq(False,False), 159.7800992093843)
  assert approx_equal(ma.sum_sq(False,True), 127.9993737702583)
  assert approx_equal(ma.sum_sq(True,False).data[1:-1],
    (120.9766669251709, 28.600215901274829, 10.203216382938619))
  assert approx_equal(ma.sum_sq(True,True).data[1:-1],
    (100.91243638831548, 28.650328284843834, 10.466799373124193))
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
  s = StringIO()
  ma.count_and_fraction_in_bins(
    data_value_to_count=2).show(f=s)
  assert not show_diff(s.getvalue(), """\
unused:        - 5.0001 [ 0/0 ] 0 0.0000
bin  1: 5.0001 - 1.4346 [21/21] 1 0.0476
bin  2: 1.4346 - 1.1432 [18/18] 0 0.0000
bin  3: 1.1432 - 1.0000 [ 9/10] 0 0.0000
unused: 1.0000 -        [ 0/0 ] 0 0.0000
""")
  s = StringIO()
  ma.count_and_fraction_in_bins(
    data_value_to_count=2, count_not_equal=True).show(f=s)
  assert not show_diff(s.getvalue(), """\
unused:        - 5.0001 [ 0/0 ]  0 0.0000
bin  1: 5.0001 - 1.4346 [21/21] 20 0.9524
bin  2: 1.4346 - 1.1432 [18/18] 18 1.0000
bin  3: 1.1432 - 1.0000 [ 9/10]  9 1.0000
unused: 1.0000 -        [ 0/0 ]  0 0.0000
""")
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
  for f_sq_as_f_algorithm in ["xtal_3_7", "crystals"]:
    for f_as_f_sq_algorithm in ["simple", "shelxl"]:
      s = ma.f_as_f_sq(algorithm=f_as_f_sq_algorithm)
      v = s.f_sq_as_f(algorithm=f_sq_as_f_algorithm)
      assert approx_equal(ma.data(), v.data())
      assert not_approx_equal(ma.sigmas(), v.sigmas())
      s = miller.array(ma, ma.data()).f_as_f_sq(algorithm=f_as_f_sq_algorithm)
      v = s.f_sq_as_f(algorithm=f_sq_as_f_algorithm)
      assert approx_equal(ma.data(), v.data())
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
  a3 = miller.array(
    miller.set(xs, flex.miller_index(((1,0,1), (0,0,-4), (0,0,-5), (1,-2,3)))),
    flex.double((5,6,7,8)),
    flex.double((1,2,3,4)))
  m1 = a1.matching_set(other=a2, data_substitute=13)
  assert m1.indices() is a2.indices()
  assert approx_equal(m1.data(), (13, 1))
  m2 = a2.matching_set(other=a1, data_substitute=15, sigmas_substitute=17)
  assert m2.indices() is a1.indices()
  assert approx_equal(m2.data(), (4, 15))
  assert approx_equal(m2.sigmas(), (6, 17))
  m3 = a2.matching_set(other=a3, data_substitute=a3.data(),
                       sigmas_substitute=a3.sigmas())
  assert m3.indices() is a3.indices()
  assert approx_equal(m3.data(), (5,6,3,4))
  assert approx_equal(m3.sigmas(), (1,2,5,6))
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
  p1, isel = sg.expand_to_p1(return_iselection=True)
  assert list(isel) == [0,1,1]
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
  sg = sg.customized_copy(data=flex.bool([False, True]), sigmas=None)
  p1 = sg.expand_to_p1()
  assert p1.indices().size() == 3
  assert p1.data().all_eq(flex.bool([False, True, True]))
  sg = sg.customized_copy(data=flex.int([3, 5]))
  p1 = sg.expand_to_p1()
  assert p1.indices().size() == 3
  assert p1.data().all_eq(flex.int([3, 5, 5]))
  #
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  data = flex.double((1,2))
  a = miller.array(miller.set(xs, mi), data)
  ph = flex.double((10,20))
  b = a.phase_transfer(ph, deg=True)
  assert approx_equal(b.amplitudes().data(), a.data())
  assert approx_equal(b.intensities().data(), flex.pow2(a.data()))
  assert approx_equal(b.phases(deg=True).data(), (10,0))
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
  assert ma.remove_cone(fraction_percent=10).data().size() == 43
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
    d_max, d_min = ma.d_max_min()
    cl = ma.resolution_filter(d_min = 6.0,  d_max = d_max).completeness()
    ch = ma.resolution_filter(d_min = d_min,d_max = 6.0).completeness(d_max =
                                                                           6.0)
    ca = ma.resolution_filter().completeness()
    assert approx_equal(cl+ch+ca, 3.0)
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
  ml = pickle.loads(pickle.dumps(ma))
  sa = StringIO()
  sl = StringIO()
  ma.binner().show_summary(f=sa)
  ml.binner().show_summary(f=sl)
  assert not show_diff(sl.getvalue(), sa.getvalue())
  ma.clear_binner()
  assert ma.binner() is None
  ml = pickle.loads(pickle.dumps(ma))
  assert ml.binner() is None
  sa = StringIO()
  sl = StringIO()
  ma.show_summary(f=sa).show_array(f=sa)
  ml.show_summary(f=sl).show_array(f=sl)
  assert sa.getvalue() == sl.getvalue()
  #
  for i_trial in [0,1,2]:
    if (i_trial == 0):
      d = flex.random_double(size=ma.indices().size())
    elif (i_trial == 1):
      d = d < 0.5
    else:
      d = flex.int(list(flex.random_size_t(size=ma.indices().size()) % 100))
    b = ma.customized_copy(
      indices=ma.indices().concatenate(ma.indices()),
      data=d.concatenate(d))
    m = b.merge_equivalents()
    s = StringIO()
    m.show_summary(n_bins=3, out=s, prefix="@#")
    if (i_trial == 0):
      assert not show_diff(s.getvalue().replace("-0.0000", " 0.0000"), """\
@#R-linear = sum(abs(data - mean(data))) / sum(abs(data))
@#R-square = sum((data - mean(data))**2) / sum(data**2)
@#In these sums single measurements are excluded.
@#                             Redundancy       Mean      Mean
@#                           Min  Max   Mean  R-linear  R-square
@#unused:         - 13.0001
@#bin  1: 13.0001 -  4.3977    2    2  2.000    0.0000    0.0000
@#bin  2:  4.3977 -  3.5133    2    2  2.000    0.0000    0.0000
@#bin  3:  3.5133 -  3.0759    2    2  2.000    0.0000    0.0000
@#unused:  3.0759 -
""")
    else:
      assert not show_diff(s.getvalue(), """\
@#                             Redundancy
@#                           Min  Max   Mean
@#unused:         - 13.0001
@#bin  1: 13.0001 -  4.3977    2    2  2.000
@#bin  2:  4.3977 -  3.5133    2    2  2.000
@#bin  3:  3.5133 -  3.0759    2    2  2.000
@#unused:  3.0759 -
""")
  #
  ma.set_info(miller.array_info(
    source="fake",
    crystal_symmetry_from_file=crystal.symmetry(
      unit_cell=(10,10,12,90,90,120),
      space_group="P6")))
  assert ma.is_in_asu()
  s = ma.crystal_symmetry_is_compatible_with_symmetry_from_file()
  assert not s.unit_cell_is_compatible
  assert not s.space_group_is_compatible
  assert not show_diff(s.format_error_message(data_description="made up"), """\
Working crystal symmetry is not compatible with crystal symmetry from reflection file:
  made up: fake
  Unit cell from file: (10, 10, 12, 90, 90, 120)
    Working unit cell: (11, 11, 13, 90, 90, 120)
  Space group from file: P 6
    Working space group: P 3 1 m""")
  #
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  data = flex.double((1,2,3,4))
  sigmas = flex.double((0.1,0.2,0.3,0.4)).reversed()
  ms = miller.set(xs, mi)
  ma = miller.array(ms, data)
  assert ma.min_f_over_sigma() is None
  ma = miller.array(ms, data, sigmas)
  assert approx_equal(ma.min_f_over_sigma(), 2.5)
  sigmas[1] = 0
  assert approx_equal(ma.min_f_over_sigma(), 0)
  sigmas[2] = -1
  ma = miller.array(ms, data, sigmas)
  try: ma.min_f_over_sigma()
  except AssertionError: pass
  else: raise Exception_expected
  ma = miller.set(xs, indices=flex.miller_index()).array(
    data=flex.double(), sigmas=flex.double())
  assert ma.min_f_over_sigma() is None
  #
  ma = miller.array(ms, data[:2])
  maa = ma.combine(other=ma)
  assert not maa.anomalous_flag()
  ms = miller.set(xs, mi, anomalous_flag=True)
  ma = miller.array(ms, data[:2])
  maa = ma.combine(other=ma)
  assert maa.anomalous_flag()
  #
  maa = ma * ma
  assert approx_equal(maa.data(), [1,4])
  maa = ma * ma.data()
  assert approx_equal(maa.data(), [1,4])
  #
  ma = ms.array(
    data=flex.complex_double([1+2j, 2-3j]), sigmas=flex.double((2,3)))
  maa = ma.as_amplitude_array()
  assert maa.is_xray_amplitude_array()
  assert approx_equal(maa.data(), [5**0.5, 13**0.5])
  assert approx_equal(maa.sigmas(), (2,3))
  mai = ma.as_intensity_array()
  assert mai.is_xray_intensity_array()
  assert approx_equal(mai.data(), [5, 13])
  assert approx_equal(maa.as_intensity_array().data(), [5, 13])
  assert approx_equal(maa.as_intensity_array().sigmas(), mai.sigmas())
  assert approx_equal(mai.as_amplitude_array().data(), [5**0.5, 13**0.5])
  mai.set_sigmas(None)
  assert mai.sigmas() is None
  mai.set_sigmas(maa.as_intensity_array().sigmas())
  assert approx_equal(maa.as_intensity_array().sigmas(), mai.sigmas())
  #
  ms = miller.build_set(
    crystal_symmetry=xs,
    anomalous_flag=False,
    d_min=2)
  sz = ms.indices().size()
  assert sz >= 10
  mt = flex.mersenne_twister(seed=1)
  mc = ms.array(
    data=flex.complex_double(
      mt.random_double(sz)*2-1,
      mt.random_double(sz)*2-1))
  mb = ms.array(data=mt.random_double(sz))
  ma = mc.as_amplitude_array()
  for fc in [mc, ma]:
    sel = ma.f_obs_f_calc_fan_outlier_selection(f_calc=fc)
    assert sel.count(True) == 0
    sel = mb.f_obs_f_calc_fan_outlier_selection(f_calc=fc)
    assert sel.count(True) == 8 # depends on mt seed above
  xs = crystal.symmetry((3,4,5,85,95,105), "P 1")
  mi = flex.miller_index(((1,2,3), (3,0,3), (2,4,1),(0,1,2)))
  ms = miller.set(xs, mi, False)
  ma = miller.array(ms, data=flex.double([1,2,3,4]))
  maa = ma.slice(axis="l", slice_index=3)
  assert approx_equal(maa.data(), [1.0,2.0])
  maa = ma.slice(axis="h", slice_start=1, slice_end=2)
  assert approx_equal(maa.data(), [1.0,3.0])
  #
  ms = miller.set(crystal_symmetry=crystal.symmetry(
    unit_cell=(7.783,8.7364,10.9002,90,102.984,90),
    space_group_symbol='P 1 21 1'),
                  indices=flex.miller_index((
                    (1,0,7),
                    (1,0,8),
                    (1,0,9),
                    (1,0,10),
                    (1,0,11),
                    (1,0,12),
                    (1,1,0),
                    (1,1,1),
                    (1,1,2),
                    (1,1,3),
                    (1,1,4),
                    (1,1,5))))
  fo2 = ms.array(data=flex.double((45.53,-3.10,168.10,-1.51,18.92,-6.10,856.51,2795.86,
                                   6685.00,959.27,1125.99,765.78)),
                 sigmas=flex.double((3.06,3.10,4.50,3.42,4.12,6.72,22.33,662.76,
                                     177.27,16.72,16.36,17.00)))\
      .set_observation_type_xray_intensity()
  fc2 = ms.array(data=flex.double((45.88,0.98,159.37,2.39,22.74,0.03,907.43,
                                   4412.96,7872.73,886.38,1144.82,783.14)))\
      .set_observation_type_xray_intensity()
  s = StringIO()
  r = fo2.show_disagreeable_reflections(fc2, out=s)
  assert r.delta_f_sq_over_sigma.size() == 12
  assert not show_diff(s.getvalue(), """\
  h   k   l       Fo^2      Fc^2   |Fo^2-Fc^2|/sig(F^2)   Fc/max(Fc)  d spacing(A)
  1   1   2    6685.00   7872.73             6.70              1.00          3.60
  1   1   3     959.27    886.38             4.36              0.34          2.81
  1   1   1    2795.86   4412.96             2.44              0.75          4.72
  1   1   0     856.51    907.43             2.28              0.34          5.73
  1   0   9     168.10    159.37             1.94              0.14          1.13
  1   0   8      -3.10      0.98             1.32              0.01          1.26
  1   1   4    1125.99   1144.82             1.15              0.38          2.27
  1   0  10      -1.51      2.39             1.14              0.02          1.02
  1   1   5     765.78    783.14             1.02              0.32          1.89
  1   0  11      18.92     22.74             0.93              0.05          0.93
  1   0  12      -6.10      0.03             0.91              0.00          0.86
  1   0   7      45.53     45.88             0.11              0.08          1.43
""")

def exercise_debye_waller():
  xs = crystal.symmetry((3,4,5,85,95,105), "P 1")
  mi = flex.miller_index(((1,2,3), (-2,1,-3)))
  ms = miller.set(xs, mi, False)
  dw = ms.debye_waller_factors
  assert approx_equal(
    dw(u_iso=0.01).data(),
    [0.8488782, 0.8378216])
  assert approx_equal(
    dw(b_iso=1.01).data(),
    [0.8109245, 0.7974382])
  assert approx_equal(
    dw(u_cart=[0.06,0.04,0.05,0.01,0.02,0.03]).data(),
    [0.2393822, 0.2899416])
  assert approx_equal(
    dw(b_cart=[1.06,1.04,1.05,1.01,1.02,1.03]).data(),
    [0.5460102, 0.6644463])
  assert approx_equal(
    dw(u_cif=[0.04,0.06,0.05,0.01,0.02,0.03]).data(),
    [0.2517064, 0.4104245])
  assert approx_equal(
    dw(u_star=[0.0004,0.0006,0.0005,0.0001,0.0002,0.0003]).data(),
    [0.7813437, 0.8726676])
  #
  ma = ms.array(data=flex.double([2,3]), sigmas=flex.double([4,5]))
  ap = ma.apply_debye_waller_factors
  maa = ap(u_iso=0.01)
  assert approx_equal(maa.data(), [1.697756, 2.513465])
  assert approx_equal(maa.sigmas(), [3.395513, 4.189108])
  maa = ap(u_cart=[0.06,0.04,0.05,0.01,0.02,0.03], apply_to_sigmas=False)
  assert approx_equal(maa.data(), [0.4787644, 0.8698248])
  assert approx_equal(maa.sigmas(), [4, 5])
  #
  ap(b_iso=-250, exp_arg_limit=60)
  ap(b_iso=-250, truncate_exp_arg=True)

def exercise_r1_factor():
  cs = crystal.symmetry((1,2,3), "P21/a")
  mi = flex.miller_index([(1,-1,1), (2,0,-3), (4,1,-2)])
  f_o = miller.array(miller.set(cs, mi),
                     data=flex.double([10, 3, -1]))
  f_c = miller.array(miller.set(cs, mi),
                     data=flex.complex_double([5+3j, 2-1j, 1+2j]))
  assert approx_equal(f_o.r1_factor(f_c, assume_index_matching=True),
                      0.440646)
  f1_o = f_o.select(flex.random_permutation(f_o.size()))
  f1_c = f_c.select(flex.random_permutation(f_c.size()))
  assert approx_equal(f1_o.r1_factor(f1_c), 0.440646)
  assert approx_equal(f1_o.r1_factor(f1_c, scale_factor=Auto), 0.2802183)
  f1_o *= 2
  assert approx_equal(f1_o.r1_factor(f1_c, scale_factor=2), 0.440646)
  assert approx_equal(f1_o.r1_factor(f1_c, scale_factor=Auto), 0.2802183)
  f_o *= 2
  assert approx_equal(
    f_o.r1_factor(f_c, scale_factor=2, assume_index_matching=True), 0.440646)
  f1_c.indices().append((4,5,6))
  f1_c.data().append(0.5)
  try:
    f1_o.r1_factor(f1_c)
    raised = False
  except AssertionError:
    raised = True
  assert raised

def exercise_scale_factor():
  flex.set_random_seed(0)
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,11,12,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=2)
  f_c = miller_set.array(
    data=flex.polar(
      flex.random_double(miller_set.size(), 10)*10-5,
      flex.random_double(miller_set.size(), 10)*10-5))
  f_o = f_c.as_amplitude_array()
  f_sq_o = f_c.as_intensity_array()
  # first with no scale factor
  assert approx_equal(f_o.scale_factor(f_c), 1.)
  assert approx_equal(f_sq_o.scale_factor(f_c), 1.)
  assert approx_equal(
    f_o.scale_factor(f_c, cutoff_factor=0), 1.)
  # now try with scale factor
  scale_factor = flex.random_double()
  f_o = f_c.customized_copy(data=flex.abs(f_c.data())*scale_factor)
  f_o.set_observation_type_xray_amplitude()
  f_sq_o = f_o.as_intensity_array()
  assert approx_equal(f_o.scale_factor(f_c), scale_factor)
  assert approx_equal(f_sq_o.scale_factor(
    f_c, cutoff_factor=0.9), scale_factor*scale_factor)
  # let's make things more random and using weights
  f_o = miller_set.array(data=scale_factor * flex.abs(f_c.data())
                         + (flex.random_double(miller_set.size())*2-1),
                         sigmas=flex.random_double(miller_set.size()))
  f_o.set_observation_type_xray_amplitude()
  f_sq_o = f_o.f_as_f_sq()
  weights = 1./flex.pow2(f_o.sigmas())
  assert approx_equal(f_o.scale_factor(f_c, weights=weights),
                      scale_factor, eps=0.1)
  # using all of the data should give better estimate of true scale factor
  assert abs(f_o.scale_factor(f_c, weights=weights, cutoff_factor=0.5)
             - scale_factor) \
         > abs(f_o.scale_factor(f_c, weights=weights) - scale_factor)
  weights = 1./flex.pow2(f_sq_o.sigmas())
  assert approx_equal(f_sq_o.scale_factor(f_c, weights=weights),
                      scale_factor*scale_factor, eps=0.1)
  assert abs(f_sq_o.scale_factor(f_c, weights=weights, cutoff_factor=0.5)
             - scale_factor*scale_factor) \
         > abs(f_sq_o.scale_factor(f_c, weights=weights)
               - scale_factor*scale_factor)
  f_sq_o.setup_binner(n_bins=10)
  for w in [weights, None]:
    binned_data = f_sq_o.scale_factor(f_c, weights=w, use_binning=True)
    assert isinstance(binned_data, miller.binned_data)
    assert len(binned_data.data) == binned_data.binner.n_bins_all()

def exercise_array_2(space_group_info):
  xs = space_group_info.any_compatible_crystal_symmetry(volume=60)
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
        r = m.redundancies().data().select(p)
        sr = s * flex.sqrt(r.as_double())
        assert approx_equal(sr, sigmas)
      #
      us = ps.select(flex.random_permutation(size=ps.indices().size())) \
        .unique_under_symmetry()
      assert us.indices().size() == m.array().indices().size()
      assert us.map_to_asu().common_set(m.array()).indices().size() \
          == m.array().indices().size()
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
  #
  mi = flex.miller_index([(0,0,0)])
  matches = miller.match_bijvoet_mates(sgtbx.space_group_type(), mi)
  assert matches.pairs().size() == 0
  assert list(matches.singles('+')) == [0]
  assert matches.singles('-').size() == 0

def exercise_map_to_asu(space_group_info):
  crystal_symmetry = space_group_info.any_compatible_crystal_symmetry(
    asu_volume=200)
  for anomalous_flag in [False, True]:
    miller_set = miller.build_set(
      crystal_symmetry=crystal_symmetry,
      anomalous_flag=anomalous_flag,
      d_min=2)
    assert miller_set.indices().size() > 30
    ampl = miller_set.array(
      data=flex.random_double(size=miller_set.indices().size()))
    phases_rad = flex.random_double(
      size=miller_set.indices().size(), factor=2*math.pi)
    cmplx = ampl.phase_transfer(phase_source=phases_rad, deg=False)
    acentric = ~cmplx.centric_flags().data()
    for trig in [flex.cos, flex.sin]:
      assert approx_equal(
        trig(cmplx.phases(deg=False).data().select(acentric)),
        trig(phases_rad.select(acentric)))
    cmplx_exp = cmplx.expand_to_p1().customized_copy(
      crystal_symmetry=cmplx)
    for deg,f in [(False,1), (True,math.pi/180)]:
      phases_exp = cmplx_exp.phases(deg=deg)
      cmplx_exp_asu_phases = cmplx_exp.map_to_asu().phases(deg=deg)
      phases_exp_asu = phases_exp.map_to_asu(deg=deg)
      for trig in [flex.cos, flex.sin]:
        assert approx_equal(
          trig(cmplx_exp_asu_phases.data()*f),
          trig(phases_exp_asu.data()*f))
      phases = cmplx.phases(deg=deg)
      phases_exp = phases.expand_to_p1(phase_deg=deg).customized_copy(
        crystal_symmetry=phases)
      phases_exp_asu = phases_exp.map_to_asu(deg=deg)
      for trig in [flex.cos, flex.sin]:
        assert approx_equal(
          trig(cmplx_exp_asu_phases.data()*f),
          trig(phases_exp_asu.data()*f))

def exercise_complete_array():
  crystal_symmetry = crystal.symmetry((2.1,3,4), "P 2 2 2")
  set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=1)
  ni = set.indices().size()
  for sigmas in [None, flex.random_double(ni)]:
    array = set.array(
      data=flex.random_double(size=set.indices().size()),
      sigmas=sigmas)
    compl = array.complete_array()
    assert compl.indices().all_eq(array.indices())
    sel = flex.random_permutation(size=ni)[:(ni*2)//3]
    selected = array.select(sel)
    compl = selected.complete_array(d_min=1)
    assert compl.indices().size() == array.size()
    new_data_sel = compl.data() >= 0
    assert new_data_sel.count(True) == sel.size()
    if (sigmas is None):
      assert compl.sigmas() is None
    else:
      new_sigmas_sel = compl.sigmas() >= 0
      assert new_sigmas_sel.count(True) == sel.size()
      assert new_sigmas_sel.all_eq(new_data_sel)
    sel = flex.random_permutation(size=ni)[:(ni*2)//3]
    selected = array.select(sel)
    compl_to_4_angstrom = selected.complete_array(d_min=4)
    assert compl_to_4_angstrom.completeness() < 1
    assert approx_equal(
      compl_to_4_angstrom.resolution_filter(d_min=4).completeness(), 1)

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
    assert not fft_map.real_map_unpadded(in_place=False).is_padded()
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
  mpe = f_calc.mean_phase_error(f_calc.phases())
  assert approx_equal(mpe, 0.0)
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

def exercise_lsd_map():
  symmetry = crystal.symmetry(
    unit_cell=(10, 10, 10, 90, 90, 90),
    space_group_symbol="P1")
  d_min = 2
  structure = xray.structure(crystal_symmetry=symmetry)
  # first check with no atoms
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  lsd = f_calc.local_standard_deviation_map(radius=3.5)
  assert flex.max(lsd.real_map_unpadded()) == 0
  assert flex.min(lsd.real_map_unpadded()) == 0
  # now let's add some atoms
  atmrad = flex.double()
  xyzf = flex.vec3_double()
  from cctbx.eltbx import van_der_waals_radii
  a,b,c,_,_,_ = symmetry.unit_cell().parameters()
  for k in xrange(10):
    scatterer = xray.scatterer(
      site=(0, 0.1*math.sin(math.pi*4*k/b), k/b),
      scattering_type="S",
      u=0.02)
    structure.add_scatterer(scatterer)
    atmrad.append(van_der_waals_radii.vdw.table[scatterer.element_symbol()])
    xyzf.append(scatterer.site)
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  lsd = f_calc.local_standard_deviation_map(radius=1.5)
  lsd_map = lsd.real_map_unpadded()
  hist = flex.histogram(data=lsd_map.as_1d(), n_slots=lsd_map.size())
  cutoff = hist.get_cutoff(int(lsd_map.size()*(1-0.85)))
  mask = flex.size_t()
  mask.resize(lsd_map.accessor(), 1)
  mask.set_selected(lsd_map > cutoff, 0)
  # compare to vdw radii-based mask
  from cctbx.masks import around_atoms
  m1 = around_atoms(
    structure.unit_cell(),
    structure.space_group().order_z(),
    structure.sites_frac(),
    atmrad,
    lsd.n_real(),
    solvent_radius=1,
    shrink_truncation_radius=1)
  corr = flex.linear_correlation(
    mask.as_double().as_1d(), m1.data.as_double().as_1d())
  assert corr.coefficient() > 0.90


def exercise_phased_translation_coeff(d_min = 1.0,
                                      algorithm = "direct",
                                      shift = [0.7, 1.2, 1.4],
                                      resolution_factor = 1./3):
  cs = crystal.symmetry((5, 5, 5, 90, 90, 90), "P 1")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer(
    [xray.scatterer("c", (0.2, 0.3, 0.4), u=0.2, occupancy=1)])
  xrs = xray.structure(sp, scatterers)
  f_calc = xrs.structure_factors(d_min = d_min, algorithm = algorithm).f_calc()
  f_obs = abs(f_calc)
  xrs_t = xrs.translate(x = shift[0], y = shift[1], z = shift[2])
  f_calc_t = xrs_t.structure_factors(d_min = d_min,
    algorithm = algorithm).f_calc()
  for fom in [None, f_obs.array(data=flex.double(f_obs.data().size(),1))]:
    result = f_obs.phased_translation_function_coeff(
      phase_source = f_calc,
      f_calc = f_calc_t,
      fom = fom)
    from cctbx import maptbx
    fft_map = result.fft_map(resolution_factor = resolution_factor,
                             symmetry_flags = maptbx.use_space_group_symmetry)
    fft_map.apply_fourier_scaling() # just to exercise this method
    fft_map.apply_sigma_scaling()
    fft_map_data = fft_map.real_map_unpadded()
    crystal_gridding_tags = fft_map.tags()
    cluster_analysis = crystal_gridding_tags.peak_search(
      parameters = maptbx.peak_search_parameters(),
      map = fft_map_data)
    assert cluster_analysis.sites().size() == 1
    expected_shift = [1.-(0.7/5), 1.-(1.2/5), 1.-(1.4/5)]
    assert approx_equal(expected_shift, cluster_analysis.sites()[0], eps=1e-2)

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
      arrays[0] = arrays[0].generate_bijvoet_mates()
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

def generate_random_hl(miller_set, coeff_range=5, max_centric_multiplier=None, set_a=None):
  # added options to skew the phase distribution in one or the other way
  phase_restriction = miller_set.space_group().phase_restriction
  hl = flex.hendrickson_lattman()
  for h in miller_set.indices():
    phase_info = phase_restriction(h)
    if (phase_info.is_centric()):
      min_fom = 0.01
      if set_a is not None:
        min_fom = 1.0-1E-12
      fom = max(min_fom, random.random()*0.95)
      if  max_centric_multiplier is not None:
        fom = fom*max_centric_multiplier
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
      coefs = [one_random_hl(f) for i in xrange(4)]
      if set_a is not None:
        coefs[0] = set_a
      hl.append(coefs)
  return miller.array(miller_set=miller_set, data=hl)

def exercise_average_and_generate_bijvoet_mates_hl(hl):
  hl_ave1 = hl.average_bijvoet_mates()
  hl_gen1 = hl_ave1.generate_bijvoet_mates().adopt_set(hl)
  hl_ave2 = hl_gen1.average_bijvoet_mates().adopt_set(hl_ave1)
  hl_gen2 = hl_ave2.generate_bijvoet_mates().adopt_set(hl_gen1)
  assert approx_equal(hl_ave2.data(), hl_ave1.data())
  assert approx_equal(hl_gen2.data(), hl_gen1.data())

def exercise_phase_integrals(space_group_info):
  crystal_symmetry = space_group_info.any_compatible_crystal_symmetry(
    asu_volume=200)
  is_centric = space_group_info.group().is_centric
  for anomalous_flag in [False, True]:
    miller_set = miller.build_set(
      crystal_symmetry=crystal_symmetry,
      anomalous_flag=anomalous_flag,
      d_min=2)
    assert miller_set.indices().size() > 30
    sg_hl = generate_random_hl(miller_set=miller_set)
    if (anomalous_flag):
      exercise_average_and_generate_bijvoet_mates_hl(sg_hl)
    p1_hl = sg_hl.expand_to_p1()
    sg_phase_integrals = sg_hl.phase_integrals(n_steps=360//5)
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
    for phase_integrator_n_steps in [None, 360//5]:
      with_phases = amplitude_array.phase_transfer(
        phase_source=sg_hl,
        phase_integrator_n_steps=phase_integrator_n_steps)
      assert flex.max(  flex.abs(amplitude_array.data()
                      - abs(with_phases).data())) < 1.e-6
      assert with_phases.mean_weighted_phase_error(
        phase_source=sg_phase_integrals) < 1.e-6


    # test the entropy calculations: maximum uncertainty
    sg_hl = generate_random_hl(miller_set=miller_set, coeff_range=1E-12,max_centric_multiplier=1E-12)
    mean_entropy = sg_hl.phase_entropy(False,False,True)
    assert mean_entropy < 1e-3


    # test the entropy calculations: no uncertainty
    sg_hl = generate_random_hl(miller_set=miller_set, coeff_range=1E-12,set_a=90000.0)
    mean_entropy = sg_hl.phase_entropy(False,False,True)
    if miller_set.space_group() == sgtbx.space_group_info( "P1" ).group():
      assert mean_entropy > 8.2
    if miller_set.space_group() == sgtbx.space_group_info( "P-1" ).group():
      assert mean_entropy > 0.99





def exercise_map_correlation():
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (1,-5,4)))

  data = flex.double((1,2))
  a = miller.array(miller.set(xs, mi), data)
  ph = flex.double((10,20))
  x = a.phase_transfer(ph, deg=True)

  data = flex.double((1,2))
  a = miller.array(miller.set(xs, mi), data)
  ph = flex.double((10,20))
  y = a.phase_transfer(ph, deg=True)
  assert approx_equal(x.map_correlation(y), 1.0)

  data = flex.double((1,2))
  a = miller.array(miller.set(xs, mi), data)
  ph = flex.double((10+180,20+180))
  y = a.phase_transfer(ph, deg=True)
  assert approx_equal(x.map_correlation(y), -1.0)

  data = flex.double((1,2))
  a = miller.array(miller.set(xs, mi), data)
  ph = flex.double((10+90,20+90))
  y = a.phase_transfer(ph, deg=True)
  assert approx_equal(x.map_correlation(y), 0.0)

  data = flex.double((-1,-2))
  a = miller.array(miller.set(xs, mi), data)
  ph = flex.double((10+90,20+90))
  y = a.phase_transfer(ph, deg=True)
  assert approx_equal(x.map_correlation(y), 0.0)

  data = flex.double((-1,-2))
  a = miller.array(miller.set(xs, mi), data)
  ph = flex.double((10+45,20+45))
  y = a.phase_transfer(ph, deg=True)
  assert approx_equal(x.map_correlation(y), -1./math.sqrt(2.0))

def exercise_concatenate():
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  data = flex.double((1,2))
  sigmas = flex.double((0.1,0.2))
  ms = miller.set(xs, mi)
  ma1 = miller.array(ms, data, sigmas)
  #
  mi = flex.miller_index(((2,-4,6), (0,0,-8)))
  data = flex.double((2,4))
  sigmas = flex.double((0.2,0.4))
  ms = miller.set(xs, mi)
  ma2 = miller.array(ms, data, sigmas)
  result = ma1.concatenate(other = ma2)
  assert approx_equal(result.indices(), [(1,-2,3),(0,0,-4),(2,-4,6),(0,0,-8)])
  assert approx_equal(result.data(), [1.0, 2.0, 2.0, 4.0])
  assert approx_equal(result.sigmas(), [0.1, 0.2, 0.2, 0.4])
  mi = flex.miller_index(((2,-4,6), (0,0,-8)))
  data = flex.double((2,4))
  ms = miller.set(xs, mi)
  ma3 = miller.array(ms, data)
  result = ma1.concatenate(other = ma3)
  assert approx_equal(result.indices(), [(1,-2,3),(0,0,-4),(2,-4,6),(0,0,-8)])
  assert approx_equal(result.data(), [1.0, 2.0, 2.0, 4.0])
  assert result.sigmas() is None

def run_call_back(flags, space_group_info):
  exercise_array_2(space_group_info)
  exercise_map_to_asu(space_group_info)
  exercise_squaring_and_patterson_map(space_group_info, verbose=flags.Verbose)
  exercise_array_correlation(space_group_info)
  exercise_as_hendrickson_lattman(space_group_info)
  exercise_phase_integrals(space_group_info)
  exercise_generate_r_free_flag_on_lat_sym(space_group_info)

def exercise_difference_map():
  cs = crystal.symmetry((2,2,2), "P1")
  mi = miller.set(cs, flex.miller_index(((1,-2,3), (0,0,-4))))
  f_o = miller.array(mi, data=flex.double((1,2)))
  f_c = miller.array(mi, data=flex.complex_double((2j, 3j)))
  diff = f_o.f_obs_minus_f_calc(f_obs_factor=1, f_calc=f_c)
  assert approx_equal(tuple(diff.data()), (-1j, -1j))
  cs = crystal.symmetry((2,2,2), "P23")
  mi = miller.set(cs, flex.miller_index(((1,-2,3), (0,0,-4))))
  f_o = miller.array(mi, data=flex.double((1,2)))
  f_c = miller.array(mi, data=flex.complex_double((2j, 3j)))
  diff = f_o.f_obs_minus_f_calc(f_obs_factor=1, f_calc=f_c)
  assert approx_equal(diff.amplitudes().data(), (1,1))
  assert not cs.space_group().phase_restriction(mi.indices()[0]).is_centric()
  assert approx_equal(diff.data()[0], -1j)
  ph_restrict = cs.space_group().phase_restriction(mi.indices()[1])
  assert ph_restrict.is_centric()
  assert ph_restrict.nearest_valid_phase(90, deg=True) == 0
  assert approx_equal(diff.data()[1], -1)
  f_o = miller.array(mi, data=(flex.complex_double((3j, 4j))))
  diff = f_o.f_obs_minus_f_calc(f_obs_factor=1, f_calc=f_c)
  assert approx_equal(diff.data(), ((1j, 1j)))

def exercise_multiscale():
  xrs = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info(number=19),
    elements=["C"]*10)
  f1 = abs(xrs.structure_factors(d_min=3).f_calc())
  f2 = abs(f1.deep_copy())
  assert approx_equal(
    flex.sum(f1.data()*f2.data())/flex.sum(f2.data()*f2.data()), 1)
  f2 = f2.array(data = f2.data()*10)
  assert approx_equal(
    flex.sum(f1.data()*f2.data())/flex.sum(f2.data()*f2.data()), 0.1)
  f1 = f2.multiscale(other = f1)
  assert approx_equal(
    flex.sum(f1.data()*f2.data())/flex.sum(f2.data()*f2.data()), 1)

def exercise_symmetry_agreement_factor():
  from cmath import phase as arg
  uc = uctbx.unit_cell((1, 2, 3, 90, 101, 90))
  cs = crystal.symmetry(uc, 'P 1')
  mi = miller.set(
    cs,
    flex.miller_index(((1,1,1), (-1,1,-1), (-1,1,1), (1,1,-1))),
    anomalous_flag=True)
  cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx('-x,y,-z+1/2'))
  ma = miller.array(
    mi, flex.complex_double((1+1j, -1+1j, 1-1j, -1-1j)))
  cc = [ (1+1j)*(1+1j), (-1+1j)*(-1+1j), (1-1j)*(1-1j), (-1-1j)*(-1-1j) ]
  expected = (3/math.pi**2 * sum([ abs(z)*arg(z)**2 for z in cc ])
                           / sum([ abs(z) for z in cc ]))
  assert approx_equal(ma.symmetry_agreement_factor(cb_op),
                      expected,
                      eps=1e-3)

def exercise_shelxl_extinction_correction():
  cs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  ms = miller.set(cs, mi, anomalous_flag=True)
  data = flex.complex_double((1+2j, 2+3j))
  fc = miller.array(miller_set=ms, data=data)
  exti_val = 0.44
  wavelength = 0.71073
  corr = fc.shelxl_extinction_correction(x=exti_val, wavelength=wavelength)
  assert approx_equal(corr, [0.99965712792728978, 0.99906041805376866])
  ec = xray.shelx_extinction_correction(cs.unit_cell(), value=exti_val,
                                        wavelength=wavelength)
  fc_sq = fc.as_intensity_array().data()
  coef = [math.sqrt(ec.compute(mi[0],fc_sq[0],False)[0]),
          math.sqrt(ec.compute(mi[1],fc_sq[1],False)[0])]
  assert approx_equal(corr, coef)

def exercise_complete_with1():
  xs = crystal.symmetry((3,15,20), "P1")
  mi = flex.miller_index(((0,0,1), ))
  ms = miller.set(xs, mi, anomalous_flag=False)
  cms = ms.complete_set(d_min=2, d_max=20)
  cms_size = cms.indices().size()
  #
  data1 = flex.double(cms_size, 1)
  ma1 = miller.array(cms, data1)
  data2 = flex.double(cms_size, -1)
  ma2 = miller.array(cms, data2)
  #
  sel = flex.random_bool(cms_size, 0.5)
  ma1i = ma1.select(sel)
  ma1del = ma1.select(~sel)
  #
  ma1c = ma1i.complete_with(ma2)
  assert ma1c.data().size() == ma1.data().size()
  x = ma1.sort(by_value="packed_indices")
  y = ma1c.sort(by_value="packed_indices")
  assert x.indices().all_eq(y.indices())
  assert ma1c.common_set(ma1del).data().all_eq(-1)
  assert ma1c.lone_set(ma1del).data().all_eq(1)
  sel = ma1c.data() == -1
  assert ma1c.select(sel).sort(by_value="packed_indices").indices().all_eq(
    ma1del.sort(by_value="packed_indices").indices()  )

def exercise_complete_with2():
  xs = crystal.symmetry((3,15,20), "P1")
  mi1 = flex.miller_index(((0,0,1), (0,0,2)))
  ms1 = miller.set(xs, mi1, anomalous_flag=False)
  d1 = flex.double([1,2])
  ma1 = miller.array(ms1, d1)
  #
  mi2 = flex.miller_index(((0,0,0),(0,0,1), (0,0,3)))
  ms2 = miller.set(xs, mi2, anomalous_flag=False)
  d2 = flex.double([0,-1,-2])
  ma2 = miller.array(ms2, d2)
  #
  ma1c = ma1.complete_with(other=ma2)
  assert list(ma1c.indices()) == [(0, 0, 1), (0, 0, 2), (0, 0, 0), (0, 0, 3)]
  assert list(ma1c.data()) == [1.0, 2.0, 0.0, -2.0]
  #
  ma2c = ma2.complete_with(other=ma1)
  assert list(ma2c.indices()) == [(0, 0, 0), (0, 0, 1), (0, 0, 3), (0, 0, 2)]
  assert list(ma2c.data()) == [0.0, -1.0, -2.0, 2.0]

def exercise_complete_with3():
  xs = crystal.symmetry((3,15,20), "P1")
  mi1 = flex.miller_index(((0,0,1),(0,0,2),(0,0,3), (0,0,4), (0,0,5), (0,0,6),
   (0,0,7),(0,0,8),(0,0,9),(0,0,10),(0,0,11),(0,0,12),(0,0,13),(0,0,14)))
  ms1 = miller.set(xs, mi1, anomalous_flag=False)
  d1 = flex.double(14,1)
  ma1 = miller.array(ms1, d1)
  mi2 = flex.miller_index(((0,0,3), (0,0,4), (0,0,5), (0,0,8), (0,0,9),
    (0,0,10),(0,0,11),(0,0,12)))
  ms2 = miller.set(xs, mi2, anomalous_flag=False)
  d2 = flex.double(8,2)
  ma2 = miller.array(ms2, d2)
  # Complete all
  mi2a = ma2.complete_with(ma1)
  mi2a = mi2a.sort()
  assert mi2a.indices().all_eq(ma1.indices())
  assert list(mi2a.data()) == [1,1,2,2,2,1,1,2,2,2,2,2,1,1]
  # Low resolution only
  mi2a = ma2.complete_with(ma1,d_min=6)
  mi2a = mi2a.sort()
  assert approx_equal(list(mi2a.indices()), ((0,0,1), (0,0,2), (0,0,3), (0,0,4),
    (0,0,5), (0,0,8), (0,0,9), (0,0,10),(0,0,11),(0,0,12)))
  assert list(mi2a.data()) == [1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
  # High resolution only
  mi2a = ma2.complete_with(ma1,d_max=2)
  mi2a = mi2a.sort()
  assert approx_equal(list(mi2a.indices()), ((0,0,3), (0,0,4), (0,0,5),
   (0,0,8), (0,0,9), (0,0,10),(0,0,11),(0,0,12),(0,0,13),(0,0,14)))
  assert list(mi2a.data()) == [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0]
  # Medium resolution only
  mi2a = ma2.complete_with(ma1,d_max=4, d_min=2.5)
  mi2a = mi2a.sort()
  assert approx_equal(list(mi2a.indices()), ((0,0,3), (0,0,4), (0,0,5), (0,0,6),
   (0,0,7), (0,0,8), (0,0,9), (0,0,10),(0,0,11),(0,0,12)))
  assert list(mi2a.data()) == [2.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0]
  # Complete all - reverse
  ma1a = ma1.complete_with(ma2)
  ma1a = ma1a.sort()
  assert ma1a.indices().all_eq(ma1.indices())
  assert approx_equal(ma1a.data(), ma1.data())

def exercise_complete_with_complete_with_bin_average():
  xs = crystal.symmetry((30,15,20), "P1")
  mi = flex.miller_index(((0,0,1),(0,0,2),(0,0,3), (0,0,4), (0,0,5), (0,0,6),
   (0,0,7),(0,0,8),(0,0,9),(0,0,10),(0,0,11),(0,0,12),(0,0,13),(0,0,14)))
  ms = miller.set(xs, mi, anomalous_flag=False)
  ms = ms.complete_set()
  ms = ms.select(selection=flex.random_bool(ms.indices().size(), 0.3))
  d = flex.double(ms.indices().size(),1)
  ma = miller.array(ms, d)
  ma.complete_with_bin_average()

def exercise_merge_equivalents_special_cases():
  xs = crystal.symmetry((10,10,10), "P1")
  mi = flex.miller_index([(1,2,3), (1,2,3)])
  ms = miller.set(xs, mi, anomalous_flag=False)
  d = flex.double([0.0, 165.827774048])
  s = flex.double([8.89633977597e-10, 144.83241272])
  ma = miller.array(ms, data=d, sigmas=s)
  me = ma.merge_equivalents().array()
  assert approx_equal(me.data(), [165.827774048])
  assert approx_equal(me.sigmas(), [144.83241272])

def run(args):
  exercise_complete_with_complete_with_bin_average()
  exercise_complete_with1()
  exercise_complete_with2()
  exercise_complete_with3()
  exercise_lsd_map()
  exercise_shelxl_extinction_correction()
  exercise_symmetry_agreement_factor()
  exercise_multiscale()
  exercise_difference_map()
  exercise_concatenate()
  exercise_phased_translation_coeff()
  exercise_set()
  exercise_union_of_sets()
  exercise_generate_r_free_flags(
    use_lattice_symmetry=False, verbose="--verbose" in args)
  exercise_generate_r_free_flags(
    use_lattice_symmetry=True, verbose="--verbose" in args)
  exercise_enforce_positive_amplitudes()
  exercise_binner()
  exercise_array()
  exercise_debye_waller()
  exercise_r1_factor()
  exercise_scale_factor()
  exercise_complete_array()
  exercise_crystal_gridding()
  exercise_fft_map()
  exercise_map_correlation()
  exercise_merge_equivalents_special_cases()
  debug_utils.parse_options_loop_space_groups(args, run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
