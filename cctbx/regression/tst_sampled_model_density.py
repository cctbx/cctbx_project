from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import maptbx
from cctbx import crystal
from cctbx import adptbx
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
from cctbx.development import random_structure
from cctbx.development import structure_factor_utils
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx import fftpack
import omptbx
from libtbx.test_utils import approx_equal
import libtbx.utils
import libtbx.introspection
import random
import sys
from six.moves import range

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def assign_custom_gaussians(structure, negative_a=False):
  if (negative_a):
    f = -1
  else:
    f = 1
  custom_gaussians = {
    "N": eltbx.xray_scattering.gaussian(
      (f*11.893779754638672, f*3.2774789333343506),
      (0.00015799999528098851, 10.232723236083984), 0),
    "C": eltbx.xray_scattering.gaussian(
      (f*2.657505989074707, f*1.0780789852142334, f*1.4909089803695679),
      (14.780757904052734, 0.77677500247955322, 42.086841583251953), 0),
  }
  reg = structure.scattering_type_registry(
    custom_dict=custom_gaussians,
    table="wk1995")
  assert reg.gaussian("N").n_terms() == 2
  assert reg.gaussian("N").c() == 0
  assert reg.gaussian("C").n_terms() == 3
  assert reg.gaussian("C").c() == 0
  assert reg.gaussian("O").n_terms() == 5
  assert reg.gaussian("O").c() != 0
  reg = structure.scattering_type_registry()
  assert reg.gaussian("N").n_terms() == 2

def exercise(space_group_info, const_gaussian, negative_gaussian,
             anomalous_flag,
             allow_mix,
             use_u_iso,
             use_u_aniso,
             d_min=1., resolution_factor=1./3, max_prime=5,
             quality_factor=100, wing_cutoff=1.e-6,
             exp_table_one_over_step_size=-100,
             force_complex=False,
             verbose=0):
  if (const_gaussian):
    elements=["const"]*8
  elif (negative_gaussian):
    elements=["H"]*8
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
    use_u_aniso= True,
    use_u_iso= False,
    random_u_cart_scale=0.3,
    random_u_iso=True,
    random_occupancy=True)
  random_structure.random_modify_adp_and_adp_flags_2(
                                 scatterers         = structure.scatterers(),
                                 use_u_iso          = use_u_iso,
                                 use_u_aniso        = use_u_aniso,
                                 allow_mix          = allow_mix,
                                 random_u_iso_scale = 0.3,
                                 random_u_iso_min   = 0.0)
  sampled_density_must_be_positive = True
  if (negative_gaussian):
    reg = structure.scattering_type_registry(
      custom_dict={"H": eltbx.xray_scattering.gaussian(-1)})
    assert reg.gaussian("H").n_terms() == 0
    assert reg.gaussian("H").c() == -1
    sampled_density_must_be_positive = False
  elif (not const_gaussian and random.random() < 0.5):
    if (random.random() < 0.5):
      sampled_density_must_be_positive = False
    assign_custom_gaussians(
      structure,
      negative_a=not sampled_density_must_be_positive)
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
  omptbx.env.num_threads = libtbx.introspection.number_of_processors()
  sampled_density = xray.sampled_model_density(
    unit_cell=structure.unit_cell(),
    scatterers=structure.scatterers(),
    scattering_type_registry=structure.scattering_type_registry(),
    fft_n_real=rfft.n_real(),
    fft_m_real=rfft.m_real(),
    u_base=u_base,
    wing_cutoff=wing_cutoff,
    exp_table_one_over_step_size=exp_table_one_over_step_size,
    force_complex=force_complex,
    sampled_density_must_be_positive=sampled_density_must_be_positive,
    tolerance_positive_definite=1.e-5,
    use_u_base_as_u_extra=False)
  focus = sampled_density.real_map_unpadded().focus()
  all   = sampled_density.real_map_unpadded().all()
  last  = sampled_density.real_map_unpadded().last()
  assert approx_equal(focus, last)
  assert approx_equal(all, last)
  assert sampled_density.anomalous_flag() == (anomalous_flag or force_complex)
  if (0 or verbose):
    print("const_gaussian:", const_gaussian)
    print("negative_gaussian:", negative_gaussian)
    print("number of scatterers passed:", \
      sampled_density.n_scatterers_passed())
    print("number of contributing scatterers:", \
      sampled_density.n_contributing_scatterers())
    print("number of anomalous scatterers:", \
      sampled_density.n_anomalous_scatterers())
    print("wing_cutoff:", sampled_density.wing_cutoff())
    print("exp_table_one_over_step_size:", \
      sampled_density.exp_table_one_over_step_size())
    print("exp_table_size:", sampled_density.exp_table_size())
    print("max_sampling_box_edges:", sampled_density.max_sampling_box_edges(), end=' ')
    print("(%.4f, %.4f, %.4f)" % sampled_density.max_sampling_box_edges_frac())
    if (not sampled_density.anomalous_flag()):
      print("map min:", flex.min(sampled_density.real_map()))
      print("map max:", flex.max(sampled_density.real_map()))
    else:
      print("map min:", flex.min(flex.real(sampled_density.complex_map())), end=' ')
      print(flex.min(flex.imag(sampled_density.complex_map())))
      print("map max:", flex.max(flex.real(sampled_density.complex_map())), end=' ')
      print(flex.max(flex.imag(sampled_density.complex_map())))
  if (not sampled_density.anomalous_flag() and negative_gaussian):
    assert flex.min(sampled_density.real_map()) < 0
    assert flex.max(sampled_density.real_map()) == 0
  if (not sampled_density.anomalous_flag()):
    map = sampled_density.real_map()
    assert map.all() == rfft.m_real()
    assert map.focus() == rfft.n_real()
    sf_map = rfft.forward(map)
    assert sf_map.all() == rfft.n_complex()
    assert sf_map.focus() == rfft.n_complex()
    collect_conj = True
  else:
    cfft = fftpack.complex_to_complex_3d(rfft.n_real())
    map = sampled_density.complex_map()
    assert map.all() == cfft.n()
    assert map.focus() == cfft.n()
    sf_map = cfft.backward(map)
    assert sf_map.all() == cfft.n()
    assert sf_map.focus() == cfft.n()
    collect_conj = False
  f_fft_data = maptbx.structure_factors.from_map(
    space_group=f_direct.space_group(),
    anomalous_flag=sampled_density.anomalous_flag(),
    miller_indices=f_direct.indices(),
    complex_map=sf_map,
    conjugate_flag=collect_conj).data()
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
    sampled_density_must_be_positive=sampled_density_must_be_positive,
    max_prime=max_prime)(
      xray_structure=structure,
      miller_set=f_direct,
      algorithm="fft").f_calc()
  structure_factor_utils.check_correlation(
    "direct/fft_xray", f_direct.indices(), 0,
    f_direct.data(), f_fft.data(),
    min_corr_ampl=1*0.99, max_mean_w_phase_error=1*3.,
    verbose=verbose)

def exercise_negative_parameters(verbose=0):
  structure_default = xray.structure(
    crystal_symmetry = crystal.symmetry(
      unit_cell=((10,13,17,75,80,85)),
      space_group_symbol="P 1"),
    scatterers=flex.xray_scatterer([
      xray.scatterer(label="C", site=(0,0,0), u=0.25)]))
  negative_gaussian = eltbx.xray_scattering.gaussian((1,2), (2,3), -4)
  for i_trial in range(7):
    structure = structure_default.deep_copy_scatterers()
    scatterer = structure.scatterers()[0]
    if (i_trial == 1):
      scatterer.occupancy *= -1
    elif (i_trial == 2):
      structure.scattering_type_registry(custom_dict={"C": negative_gaussian})
    elif (i_trial == 3):
      scatterer.u_iso *= -1
    elif (i_trial == 4):
      u_cart = adptbx.random_u_cart(u_scale=1, u_min=-1.1)
      assert max(adptbx.eigenvalues(u_cart)) < 0
      u_star = adptbx.u_cart_as_u_star(structure.unit_cell(), u_cart)
      scatterer.u_star = u_star
      scatterer.flags.set_use_u_aniso_only()
    elif (i_trial == 5):
      scatterer.fp = -10
    elif (i_trial == 6):
      scatterer.fp = -3
    f_direct = structure.structure_factors(
      d_min=1, algorithm="direct", cos_sin_table=False).f_calc()
    f_fft = structure.structure_factors(
      d_min=1, algorithm="fft",
      quality_factor=1.e8, wing_cutoff=1.e-10).f_calc()
    if (i_trial == 2):
      assert negative_gaussian.at_d_star_sq(f_fft.d_star_sq().data()).all_lt(0)
    if (i_trial in [5,6]):
      f = structure.scattering_type_registry().gaussian_not_optional(
        scattering_type="C").at_d_star_sq(f_fft.d_star_sq().data())
      if (i_trial == 5):
        assert flex.max(f) + scatterer.fp < 0
      else:
        assert flex.max(f) + scatterer.fp > 0
        assert flex.min(f) + scatterer.fp < 0
    cc = flex.linear_correlation(
      abs(f_direct).data(),
      abs(f_fft).data()).coefficient()
    if (cc < 0.999):
      raise AssertionError("i_trial=%d, correlation=%.6g" % (i_trial, cc))
    elif (0 or verbose):
      print("correlation=%.6g" % cc)
    #
    # very simple test of gradient calculations with negative parameters
    structure_factor_gradients = \
      cctbx.xray.structure_factors.gradients(
        miller_set=f_direct,
        cos_sin_table=False)
    target_functor = xray.target_functors.intensity_correlation(
      f_obs=abs(f_direct))
    target_result = target_functor(f_fft, True)
    xray.set_scatterer_grad_flags(scatterers = structure.scatterers(),
                                  site       = True,
                                  u_iso      = True,
                                  u_aniso    = True,
                                  occupancy  = True,
                                  fp         = True,
                                  fdp        = True)
    for algorithm in ["direct", "fft"]:
      grads = structure_factor_gradients(
        xray_structure=structure,
        u_iso_refinable_params=None,
        miller_set=f_direct,
        d_target_d_f_calc=target_result.derivatives(),
        n_parameters = structure.n_parameters(),
        algorithm=algorithm).packed()

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True)[:]: #SWITCH
    for (use_u_iso,use_u_aniso) in [(True,True),
                                    (False,True),
                                    (True,False),
                                    (False,False)]:
      for allow_mix in [False, True]:
        r = random.random()
        exercise(
          space_group_info,
          const_gaussian=r<1/3.,
          negative_gaussian=r>2/3.,
          anomalous_flag=anomalous_flag,
          use_u_iso = use_u_iso,
          use_u_aniso = use_u_aniso,
          verbose=flags.Verbose,
          allow_mix=allow_mix)

def run():
  show_times = libtbx.utils.show_times()
  debug_utils.parse_options_loop_space_groups(
    argv=sys.argv[1:],
    call_back=run_call_back,
    show_cpu_times=False)
  exercise_negative_parameters()
  show_times()

if (__name__ == "__main__"):
  run()
