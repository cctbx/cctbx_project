from cctbx import maptbx
from cctbx import miller
from cctbx.development import random_structure
from cctbx.development import structure_factor_utils
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx import fftpack
import sys

def exercise(space_group_info, anomalous_flag, conjugate_flag,
             d_min=3., resolution_factor=0.5, max_prime=5,
             verbose=0):
  structure_factors = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O"),
    random_f_prime_d_min=1,
    random_f_double_prime=anomalous_flag,
    anisotropic_flag=0001,
    random_u_iso=0001,
    random_occupancy=0001
    ).structure_factors_direct(
        anomalous_flag=anomalous_flag,
        d_min=d_min)
  f_calc_array = structure_factors.f_calc_array()
  n_real = f_calc_array.determine_gridding(
    resolution_factor=resolution_factor,
    d_min=d_min,
    max_prime=max_prime)
  if (not anomalous_flag):
    rfft = fftpack.real_to_complex_3d(n_real)
    n_complex = rfft.n_complex()
  else:
    cfft = fftpack.complex_to_complex_3d(n_real)
    n_complex = cfft.n()
  map = maptbx.structure_factors.to_map(
    f_calc_array.space_group(),
    anomalous_flag,
    f_calc_array.indices(),
    f_calc_array.data(),
    n_real,
    flex.grid(n_complex),
    conjugate_flag)
  if (not anomalous_flag):
    real_map = rfft.backward(map.complex_map())
    assert real_map.all() == rfft.m_real()
    complex_map = rfft.forward(real_map)
  else:
    real_map = cfft.backward(map.complex_map())
    assert not real_map.is_padded()
    complex_map = cfft.forward(real_map)
  complex_map /= n_real[0] * n_real[1] * n_real[2]
  assert real_map.focus() == n_real
  assert complex_map.focus() == n_complex
  from_map = maptbx.structure_factors.from_map(
    f_calc_array.unit_cell(),
    f_calc_array.space_group_info().type(),
    anomalous_flag,
    d_min,
    complex_map,
    conjugate_flag)
  match = miller.match_indices(
    f_calc_array.indices(),
    from_map.miller_indices())
  for i in (0,1):
    if (i == 0): m = f_calc_array.indices()
    else: m = from_map.miller_indices()
    for j in match.singles(i):
      assert abs(f_calc_array.unit_cell().d(m[j]) - d_min) < 1.e-5
  structure_factor_utils.check_correlation(
    "from_map-1", f_calc_array.indices(), match,
    f_calc_array.data(), from_map.data(),
    min_corr_ampl=0.9999, max_mean_w_phase_error=.01,
    verbose=verbose)
  from_map = maptbx.structure_factors.from_map(
    anomalous_flag,
    f_calc_array.indices(),
    complex_map,
    conjugate_flag)
  assert from_map.miller_indices().size() == 0
  structure_factor_utils.check_correlation(
    "from_map-2", f_calc_array.indices(), 0,
    f_calc_array.data(), from_map.data(),
    min_corr_ampl=0.9999, max_mean_w_phase_error=.01,
    verbose=verbose)

def exercise_under_sampled(space_group_info, anomalous_flag, conjugate_flag,
                           under_sampling,
                           d_min=2., resolution_factor=0.5, max_prime=5,
                           verbose=0):
  structure_factors = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O"),
    random_f_prime_d_min=1,
    random_f_double_prime=anomalous_flag,
    anisotropic_flag=0001,
    random_u_iso=0001,
    random_occupancy=0001
    ).structure_factors_direct(
        anomalous_flag=anomalous_flag,
        d_min=d_min)
  f_calc_array = structure_factors.f_calc_array()
  n_real = maptbx.determine_gridding(
    unit_cell=f_calc_array.unit_cell(),
    d_min=d_min,
    resolution_factor=resolution_factor,
    max_prime=max_prime,
    mandatory_factors=(under_sampling,)*3)
  if (not anomalous_flag):
    rfft = fftpack.real_to_complex_3d(n_real)
    n_complex = rfft.n_complex()
  else:
    cfft = fftpack.complex_to_complex_3d(n_real)
    n_complex = cfft.n()
  map = maptbx.structure_factors.to_map(
    f_calc_array.space_group(),
    anomalous_flag,
    f_calc_array.indices(),
    f_calc_array.data(),
    n_real,
    flex.grid(n_complex),
    conjugate_flag)
  if (not anomalous_flag):
    real_map = rfft.backward(map.complex_map())
    assert real_map.all() == rfft.m_real()
  else:
    real_map = cfft.backward(map.complex_map())
    assert not real_map.is_padded()
  if (0 or verbose):
    if (not anomalous_flag):
      maptbx.statistics(real_map).show_summary()
      maptbx.statistics(real_map).show_summary()
    else:
      maptbx.statistics(flex.real(real_map)).show_summary()
      maptbx.statistics(flex.imag(real_map)).show_summary()
  n_real_under_sampled = [n/under_sampling for n in n_real]
  if (not anomalous_flag):
    rfft = fftpack.real_to_complex_3d(n_real_under_sampled)
    n_complex_under_sampled = rfft.n_complex()
  else:
    cfft = fftpack.complex_to_complex_3d(n_real_under_sampled)
    n_complex_under_sampled = cfft.n()
  under_sampled_map = maptbx.structure_factors.to_map(
    f_calc_array.space_group(),
    anomalous_flag,
    f_calc_array.indices(),
    f_calc_array.data(),
    n_real_under_sampled,
    flex.grid(n_complex_under_sampled),
    conjugate_flag)
  if (not anomalous_flag):
    under_sampled_map_before_fft = under_sampled_map.complex_map().deep_copy()
    under_sampled_real_map = rfft.backward(under_sampled_map.complex_map())
    assert under_sampled_real_map.all() == rfft.m_real()
  else:
    under_sampled_real_map = cfft.backward(under_sampled_map.complex_map())
    assert not under_sampled_real_map.is_padded()
  if (0 or verbose):
    if (not anomalous_flag):
      maptbx.statistics(under_sampled_real_map).show_summary()
      maptbx.statistics(under_sampled_real_map).show_summary()
    else:
      maptbx.statistics(flex.real(under_sampled_real_map)).show_summary()
      maptbx.statistics(flex.imag(under_sampled_real_map)).show_summary()
  if (0 or verbose):
    print real_map.all(), n_complex
    print under_sampled_real_map.all(), n_complex_under_sampled
  if (not anomalous_flag):
    x_source = real_map
    y_source = under_sampled_real_map
  else:
    x_source = flex.real(real_map)
    y_source = flex.real(under_sampled_real_map)
  x = flex.double()
  n = x_source.focus()
  for i in xrange(0, n[0], under_sampling):
    for j in xrange(0, n[1], under_sampling):
      for k in xrange(0, n[2], under_sampling):
        x.append(x_source[(i,j,k)])
  y = maptbx.copy(y_source, flex.grid(y_source.focus())).as_1d()
  if (0 or verbose):
    print "x:", tuple(x)
    print "y:", tuple(y)
  assert flex.max(flex.abs(x-y)) \
      < (flex.max(flex.abs(x))+flex.max(flex.abs(y)))/2*1.e-6
  if (under_sampling == 1):
    x = maptbx.copy(x_source, flex.grid(x_source.focus())).as_1d()
    c = flex.linear_correlation(x, y)
    assert c.coefficient() >= 0.9999

def run_call_back(flags, space_group_info):
  for anomalous_flag in (00000, 0001)[:]: #SWITCH
    for conjugate_flag in (00000, 0001)[:]: #SWITCH
      exercise(space_group_info, anomalous_flag, conjugate_flag,
               verbose=flags.Verbose)
  for anomalous_flag in (00000, 0001)[:]: #SWITCH
    for conjugate_flag in (00000, 0001)[:]: #SWITCH
      for under_sampling in (1,2,3,4,5):
        exercise_under_sampled(space_group_info,
                               anomalous_flag,
                               conjugate_flag,
                               under_sampling,
                               verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
