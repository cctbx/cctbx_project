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
  n_real = maptbx.determine_grid(
    unit_cell=f_calc_array.unit_cell(),
    d_min=d_min,
    resolution_factor=resolution_factor,
    max_prime=max_prime,
    mandatory_factors=(1,1,1))
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
      assert abs(f_calc_array.unit_cell.d(m[j]) - d_min) < 1.e-5
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

def run_call_back(flags, space_group_info):
  for anomalous_flag in (00000, 0001)[:]: #SWITCH
    for conjugate_flag in (00000, 0001)[:]: #SWITCH
      exercise(space_group_info, anomalous_flag, conjugate_flag,
               verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
