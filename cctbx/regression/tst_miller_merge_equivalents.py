from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO
import random
import sys
from six.moves import range

def exercise(space_group_info, anomalous_flag,
             n_scatterers=8, d_min=2, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers)
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=anomalous_flag).f_calc()
  f = abs(f_calc)
  fs = miller.array(miller_set=f, data=f.data(), sigmas=flex.sqrt(f.data()))
  assert fs.is_unique_set_under_symmetry()
  for a in (f, fs):
    for algorithm in ["gaussian", "shelx"]:
      m = a.merge_equivalents(algorithm=algorithm)
      m.show_summary(out=StringIO())
      j = m.array().adopt_set(a)
      assert flex.linear_correlation(
        j.data(), a.data()).coefficient() > 1-1.e-6
      if (a.sigmas() is not None):
        assert flex.linear_correlation(
          j.sigmas(), a.sigmas()).coefficient() > 1-1.e-6
  redundancies = flex.size_t()
  for i in range(fs.indices().size()):
    redundancies.append(random.randrange(5)+1)
  space_group = space_group_info.group()
  r_indices = flex.miller_index()
  r_data = flex.double()
  r_sigmas = flex.double()
  for i,n in enumerate(redundancies):
    h = fs.indices()[i]
    h_eq = miller.sym_equiv_indices(space_group, h).indices()
    for j in range(n):
      r_indices.append(h_eq[random.randrange(len(h_eq))].h())
      r_data.append(fs.data()[i])
      r_sigmas.append(fs.sigmas()[i])
  r = miller.array(
    miller_set=miller.set(
      crystal_symmetry=fs,
      indices=r_indices,
      anomalous_flag=fs.anomalous_flag()),
    data=r_data,
    sigmas=r_sigmas)
  assert not r.is_unique_set_under_symmetry()
  noise = flex.random_double(size=r.indices().size())
  r = r.sort(by_value=noise)
  for algorithm in ["gaussian", "shelx"]:
    m = r.merge_equivalents(algorithm=algorithm)
    m.show_summary(out=StringIO())
    j = m.array().adopt_set(fs)
    assert j.is_unique_set_under_symmetry()
    assert flex.linear_correlation(
      j.data(), fs.data()).coefficient() > 1-1.e-6
    fssr = fs.sigmas() / flex.sqrt(redundancies.as_double())
    assert flex.linear_correlation(j.sigmas(), fssr).coefficient() > 1-1.e-6
  #
  if (anomalous_flag):
    f_calc_ave = f_calc.average_bijvoet_mates() # uses merge_equivalents
    f_calc_com = f_calc.as_non_anomalous_array().common_set(f_calc_ave)
    assert f_calc_com.indices().all_eq(f_calc_ave.indices())
    for part in [flex.real, flex.imag]:
      assert flex.linear_correlation(
        part(f_calc_com.data()),
        part(f_calc_ave.data())).coefficient() > 1-1.e-6
  # test use_internal_variance=False
  m = r.merge_equivalents(algorithm="gaussian", use_internal_variance=False)
  j = m.array().adopt_set(fs)
  fssr = fs.sigmas() / flex.sqrt(redundancies.as_double())
  assert flex.linear_correlation(j.sigmas(), fssr).coefficient() > 1-1.e-6

def exercise_incompatible_flags_replacement():
  i = flex.miller_index(((1,2,3), (1,2,3), (3,0,3), (3,0,3), (3,0,3), (1,1,2)))
  d = flex.int((1,1,0,1,0,1))
  from cctbx import crystal
  cs = crystal.symmetry(unit_cell=(10,10,10,90,90,90), space_group_symbol="P1")
  ms = miller.set(cs, i)
  ma = miller.array(ms, data=d)
  try: ma.merge_equivalents()
  except Sorry as e: assert "merge_equivalents_exact: incompatible flags" in str(e)
  else: raise Exception_expected
  merging = ma.merge_equivalents(incompatible_flags_replacement=0)
  me = merging.array()
  assert approx_equal(me.data(), (1,1,0))
  merging = ma.merge_equivalents(incompatible_flags_replacement=2)
  me = merging.array()
  assert approx_equal(me.data(), (1,1,2))

def exercise_split_unmerged():
  import random
  random.seed(42)
  flex.set_random_seed(42)

  from cctbx import crystal
  base_set = miller.build_set(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(10,10,10,90,90,90), space_group_symbol="P1"),
    d_min=1.6,
    anomalous_flag=False)
  indices = base_set.indices()
  assert (len(indices) == 510)
  unmerged_hkl = flex.miller_index()
  unmerged_data = flex.double()
  unmerged_sigmas = flex.double()
  redundancies = flex.size_t()
  # XXX grossly overengineered, but I wanted to get a realistic CC to make sure
  # the reflections are being split properly
  for i, hkl in enumerate(indices):
    n_obs = min(8, 1 + i % 12)
    redundancies.append(n_obs)
    intensity_merged = (510 - i) + (510 % 27)
    for j in range(n_obs):
      unmerged_hkl.append(hkl)
      intensity = intensity_merged + 20 * (510 % (7 * (j+1)))
      sigma = max(0.5, i % 10)
      unmerged_data.append(intensity)
      unmerged_sigmas.append(sigma)
  assert (unmerged_hkl.size() == 2877)
  unmerged_array = miller.set(
    crystal_symmetry=base_set,
    indices=unmerged_hkl,
    anomalous_flag=False).array(data=unmerged_data, sigmas=unmerged_sigmas)
  split = miller.split_unmerged(
    unmerged_indices=unmerged_hkl,
    unmerged_data=unmerged_data,
    unmerged_sigmas=unmerged_sigmas)
  assert (split.data_1.size() == split.data_2.size() == 467)
  cc = miller.compute_cc_one_half(unmerged_array)
  assert approx_equal(cc, 0.861, eps=0.001)
  unmerged_array.setup_binner(n_bins=10)
  unmerged_array.set_observation_type_xray_intensity()
  result = unmerged_array.cc_one_half(use_binning=True)
  assert approx_equal(
    result.data[1:-1],
    [0.549, 0.789, 0.843, 0.835, 0.863, 0.860, 0.893, 0.847, 0.875, 0.859],
    eps=0.05)

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(space_group_info, anomalous_flag)

def run():
  exercise_incompatible_flags_replacement()
  exercise_split_unmerged()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print("OK")

if (__name__ == "__main__"):
  run()
