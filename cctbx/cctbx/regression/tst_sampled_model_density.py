from cctbx import xray
from cctbx import maptbx
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
from cctbx.development import random_structure
from cctbx.development import structure_factor_utils
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx import fftpack
import random
import sys

def assign_custom_gaussians(structure):
  custom_gaussians = {
    "N": eltbx.xray_scattering.gaussian(
      (11.893779754638672, 3.2774789333343506),
      (0.00015799999528098851, 10.232723236083984), 0),
    "C": eltbx.xray_scattering.gaussian(
      (2.657505989074707, 1.0780789852142334, 1.4909089803695679),
      (14.780757904052734, 0.77677500247955322, 42.086841583251953), 0),
  }
  d = structure.scattering_dict(
    custom_dict=custom_gaussians,
    table="wk1995").dict()
  assert d["N"].gaussian.n_terms() == 2
  assert d["N"].gaussian.c() == 0
  assert d["C"].gaussian.n_terms() == 3
  assert d["C"].gaussian.c() == 0
  assert d["O"].gaussian.n_terms() == 5
  assert d["O"].gaussian.c() != 0
  d = structure.scattering_dict().dict()
  assert d["N"].gaussian.n_terms() == 2

def exercise(space_group_info, const_gaussian,
             anomalous_flag, anisotropic_flag,
             d_min=1., resolution_factor=1./3, max_prime=5,
             quality_factor=100, wing_cutoff=1.e-6,
             exp_table_one_over_step_size=-100,
             force_complex=00000,
             verbose=0):
  if (const_gaussian):
    elements=["const"]*8
  else:
    elements=["N", "C", "C", "O", "N", "C", "C", "O"]
  if (random.random() < 0.5):
    random_f_prime_scale=0.6
  else:
    random_f_prime_scale=0
  structure = random_structure.xray_structure(
    space_group_info,
    elements=elements,
    random_f_prime_d_min=1,
    random_f_prime_scale=random_f_prime_scale,
    random_f_double_prime=anomalous_flag,
    anisotropic_flag=anisotropic_flag,
    random_u_iso=0001,
    random_occupancy=0001)
  if (not const_gaussian and random.random() < 0.5):
    assign_custom_gaussians(structure)
  f_direct = structure.structure_factors(
    anomalous_flag=anomalous_flag,
    d_min=d_min,
    algorithm="direct").f_calc()
  crystal_gridding = f_direct.crystal_gridding(
    resolution_factor=resolution_factor,
    d_min=d_min,
    max_prime=max_prime)
  assert crystal_gridding.symmetry_flags() is None
  rfft = fftpack.real_to_complex_3d(crystal_gridding.n_real())
  u_base = xray.calc_u_base(d_min, resolution_factor, quality_factor)
  electron_density_must_be_positive = 1
  tolerance_positive_definite = 1.e-5
  sampled_density = xray.sampled_model_density(
    structure.unit_cell(),
    structure.scatterers(),
    structure.scattering_dict(),
    rfft.n_real(),
    rfft.m_real(),
    u_base,
    wing_cutoff,
    exp_table_one_over_step_size,
    force_complex,
    electron_density_must_be_positive,
    tolerance_positive_definite)
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
    print "max_sampling_box_edges:", sampled_density.max_sampling_box_edges(),
    print "(%.4f, %.4f, %.4f)" % sampled_density.max_sampling_box_edges_frac()
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
    f_direct.space_group(),
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
      algorithm="fft").f_calc()
  structure_factor_utils.check_correlation(
    "direct/fft_xray", f_direct.indices(), 0,
    f_direct.data(), f_fft.data(),
    min_corr_ampl=1*0.99, max_mean_w_phase_error=1*3.,
    verbose=verbose)

def run_call_back(flags, space_group_info):
  for anomalous_flag in (00000, 0001)[:]: #SWITCH
    for anisotropic_flag in (00000, 0001)[:]: #SWITCH
      exercise(
        space_group_info,
        const_gaussian=random.random()<0.5,
        anomalous_flag=anomalous_flag,
        anisotropic_flag=anisotropic_flag,
        verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
