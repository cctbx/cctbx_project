from cctbx import xray
from cctbx import maptbx
from cctbx.development import random_structure
from cctbx.development import structure_factor_utils
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx import fftpack
import sys

def exercise(space_group_info, anomalous_flag, anisotropic_flag,
             d_min=1., resolution_factor=1./3, max_prime=5,
             quality_factor=100, wing_cutoff=1.e-3,
             exp_table_one_over_step_size=-100,
             force_complex=00000,
             verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O", "N", "C", "C", "O"),
    random_f_prime_d_min=1,
    random_f_double_prime=anomalous_flag,
    anisotropic_flag=anisotropic_flag,
    random_u_iso=0001,
    random_occupancy=0001
    )
  f_direct = structure.structure_factors(
    anomalous_flag=anomalous_flag,
    d_min=d_min,
    direct=0001).f_calc()
  n_real = f_direct.crystal_gridding(
    resolution_factor=resolution_factor,
    d_min=d_min,
    symmetry_flags=maptbx.use_space_group_symmetry,
    max_prime=max_prime).n_real()
  rfft = fftpack.real_to_complex_3d(n_real)
  u_extra = xray.calc_u_extra(d_min, resolution_factor, quality_factor)
  electron_density_must_be_positive = 1
  sampled_density = xray.sampled_model_density(
    structure.unit_cell(),
    structure.scatterers(),
    rfft.n_real(),
    rfft.m_real(),
    u_extra,
    wing_cutoff,
    exp_table_one_over_step_size,
    force_complex,
    electron_density_must_be_positive)
  assert sampled_density.anomalous_flag() == (anomalous_flag or force_complex)
  if (0 or verbose):
    print "number of scatterers passed:", \
      sampled_density.n_scatterers_passed()
    print "number of contributing scatterers:", \
      sampled_density.n_contributing_scatterers()
    print "number of anomalous scatterers:", \
      sampled_density.n_anomalous_scatterers()
    print "wing_cutoff:", sampled_density.wing_cutoff()
    print "exp_table_one_over_step_size:", \
      sampled_density.exp_table_one_over_step_size()
    print "exp_table_size:", sampled_density.exp_table_size()
    print "max_shell_radii:", sampled_density.max_shell_radii(),
    print "(%.4f, %.4f, %.4f)" % sampled_density.max_shell_radii_frac()
  tags = maptbx.grid_tags(rfft.n_real())
  symmetry_flags = maptbx.symmetry_flags(use_space_group_symmetry=0001)
  tags.build(structure.space_group_info().type(), symmetry_flags)
  sampled_density.apply_symmetry(tags)
  if (not sampled_density.anomalous_flag()):
    map = sampled_density.real_map()
    assert map.all() == rfft.m_real()
    assert map.focus() == rfft.n_real()
    sf_map = rfft.forward(map)
    assert sf_map.all() == rfft.n_complex()
    assert sf_map.focus() == rfft.n_complex()
    collect_conj = 1
  else:
    cfft = fftpack.complex_to_complex_3d(rfft.n_real())
    map = sampled_density.complex_map()
    assert map.all() == cfft.n()
    assert map.focus() == cfft.n()
    sf_map = cfft.backward(map)
    assert sf_map.all() == cfft.n()
    assert sf_map.focus() == cfft.n()
    collect_conj = 0
  f_fft_data = maptbx.structure_factors.from_map(
    sampled_density.anomalous_flag(),
    f_direct.indices(),
    sf_map,
    collect_conj).data()
  sampled_density.eliminate_u_extra_and_normalize(
    f_direct.indices(),
    f_fft_data)
  structure_factor_utils.check_correlation(
    "direct/fft_regression", f_direct.indices(), 0,
    f_direct.data(), f_fft_data,
    min_corr_ampl=1*0.99, max_mean_w_phase_error=1*3.,
    verbose=verbose)
  f_fft = xray.structure_factors.from_scatterers(
    miller_set=f_direct,
    grid_resolution_factor=resolution_factor,
    quality_factor=quality_factor,
    wing_cutoff=wing_cutoff,
    exp_table_one_over_step_size=exp_table_one_over_step_size,
    max_prime=max_prime)(
      xray_structure=structure,
      miller_set=f_direct,
      fft=0001).f_calc()
  structure_factor_utils.check_correlation(
    "direct/fft_xray", f_direct.indices(), 0,
    f_direct.data(), f_fft.data(),
    min_corr_ampl=1*0.99, max_mean_w_phase_error=1*3.,
    verbose=verbose)

def run_call_back(flags, space_group_info):
  for anomalous_flag in (00000, 0001)[:]: #SWITCH
    for anisotropic_flag in (00000, 0001)[:]: #SWITCH
      exercise(space_group_info, anomalous_flag, anisotropic_flag,
               verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
