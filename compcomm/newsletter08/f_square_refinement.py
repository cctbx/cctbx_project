from cctbx import xray
from cctbx import crystal
from cctbx import miller
from cctbx.development import random_structure

indices = miller.build_set(
  crystal_symmetry=crystal.symmetry(unit_cell=(10,11,12, 90,105,90),
                                    space_group_symbol="P21/c"),
  anomalous_flag=False,
  d_min=0.8)
structure = random_structure.xray_structure(indices.space_group_info(),
                                            elements=['C']*6 + ['O']*2 + ['N'],
                                            volume_per_atom=18.6,
                                            random_u_iso=True)
f_ideal = structure.structure_factors(d_min=indices.d_min()).f_calc()

f_obs = f_ideal.amplitudes()
f_obs *= 2
f_obs.set_observation_type_xray_amplitude()
f_obs_square = f_ideal.norm()
f_obs_square *= 3
f_obs_square.set_observation_type_xray_intensity()

ls_against_f = xray.unified_least_squares_residual(f_obs)
ls_against_f_square = xray.unified_least_squares_residual(f_obs_square)

residuals = ls_against_f(f_ideal, compute_derivatives=True)
print "against F: value=%.3f, scale=%.3f" % (residuals.target(),
                                             residuals.scale_factor())
residuals = ls_against_f_square(f_ideal, compute_derivatives=True)
print "against F^2: value=%.3f, scale=%.3f" % (residuals.target(),
                                               residuals.scale_factor())

perturbed_structure = structure.random_shift_sites(max_shift_cart=0.2)
for s in perturbed_structure.scatterers():
  s.flags.set_grad_site(True)

refining_structure = perturbed_structure.deep_copy_scatterers()
optimiser = xray.lbfgs(
    target_functor=ls_against_f_square,
    xray_structure=refining_structure,
    structure_factor_algorithm="direct")
print "Initial L.S. residual:%.3f" % optimiser.first_target_value
structure.show_scatterers()
print "Final L.S. residual:%.3f" % optimiser.final_target_value
refining_structure.show_scatterers()

weighting = xray.weighting_schemes.shelx_weighting()
shelx_weighted_ls_against_f_square = xray.unified_least_squares_residual(
  f_obs_square, weighting=weighting)
refining_structure = perturbed_structure.deep_copy_scatterers()
optimiser = xray.lbfgs(
    target_functor=ls_against_f_square,
    xray_structure=refining_structure,
    structure_factor_algorithm="direct")
print "Initial L.S. residual:%.3f" % optimiser.first_target_value
structure.show_scatterers()
print "Final L.S. residual:%.3f" % optimiser.final_target_value
refining_structure.show_scatterers()
