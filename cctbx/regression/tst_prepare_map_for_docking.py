from __future__ import print_function
from __future__ import division
import math
from iotbx.map_model_manager import map_model_manager
from cctbx.maptbx.prepare_map_for_docking import add_ordered_volume_mask
from cctbx.maptbx.prepare_map_for_docking import assess_cryoem_errors
from cctbx.maptbx.prepare_map_for_docking import get_d_star_sq_step
from cctbx.maptbx.prepare_map_for_docking import write_mtz
from cctbx import adptbx
from scitbx.array_family import flex
from libtbx.utils import format_cpu_times
import random

def get_power_spectrum(mc):

  # Using bins of equal width in d_star_sq matches variation with resolution
  # better, but would have problems with very asymmetric boxes (not used here).
  power_spectrum = flex.double(mc.size(),1.)
  mc_copy = mc.deep_copy()
  d_star_sq_step = get_d_star_sq_step(mc_copy)
  mc_copy.setup_binner_d_star_sq_step(d_star_sq_step=d_star_sq_step)
  for i_bin in mc_copy.binner().range_used():
    sel = mc_copy.binner().selection(i_bin)
    mcsel = mc_copy.select(sel)
    fsq = flex.pow2(flex.abs(mcsel.data()))
    meanfsq = flex.mean_default(fsq, 0)
    power = math.sqrt(meanfsq)
    power_spectrum.set_selected(sel,power)

  return power_spectrum

def exercise():
  """Test prepare_map_for_docking using data with known errors."""

  # Generate two half-maps with same anisotropic signal, independent anisotropic
  # noise. Test to see how well optimal map coefficients are estimated.

  # Start by working out how large the padding will have to be so that
  # starting automatically-generated map will be large enough to contain
  # sphere with room to spare around model.
  # Use reasonably spherical model by choosing residues 69-218 of 3jd6 standard
  n_residues = 150
  start_res = 69
  b_iso = 50
  d_min = 3.0
  from cctbx.development.create_models_or_maps import generate_model
  test_model = generate_model(n_residues=n_residues, start_res=start_res,
      b_iso=b_iso)
  sites_cart = test_model.get_sites_cart()
  cart_min = flex.double(sites_cart.min())
  cart_max = flex.double(sites_cart.max())
  box_centre = (cart_min+cart_max)/2
  dsqrmax = flex.max( (sites_cart - tuple(box_centre)).norms() )**2
  model_radius = math.sqrt(dsqrmax)
  min_model_extent = flex.min(cart_max - cart_min)
  pad_to_allow_cube = model_radius - min_model_extent/2
  # Extra space needed for eventual masking
  boundary_to_smoothing_ratio = 2
  soft_mask_radius = d_min
  padding = soft_mask_radius * boundary_to_smoothing_ratio
  box_cushion = padding + pad_to_allow_cube + d_min # A bit extra

  # Make map in box big enough to cut out cube containing sphere
  mmm = map_model_manager()
  mmm.generate_map(
      n_residues=n_residues, d_min=d_min, k_sol=0.1, b_sol=50.,
      box_cushion=box_cushion)
  # Keep copy of perfect map for tests of success
  mm_start = mmm.map_manager().deep_copy()
  mmm.add_map_manager_by_id(mm_start,'perfect_map')
  model = mmm.model()
  sites_cart = model.get_sites_cart()
  cart_min = flex.double(sites_cart.min())
  cart_max = flex.double(sites_cart.max())

  # Turn starting map into map coeffs for the signal
  ucpars = mmm.map_manager().unit_cell().parameters()
  d_max=max(ucpars[0], ucpars[1], ucpars[2])
  start_map_coeffs = mmm.map_as_fourier_coefficients(
      d_min=d_min, d_max=d_max)

  # Apply anisotropic scaling to map coeffs
  b_target = (100.,200.,300.,-50.,50.,100.)
  u_star_s = adptbx.u_cart_as_u_star(
      start_map_coeffs.unit_cell(), adptbx.b_as_u(b_target))
  b_model = (b_iso,b_iso,b_iso,0.,0.,0.)  # All atoms in model have B=b_iso
  b_expected = list((flex.double(b_target) + flex.double(b_model)))
  scaled_map_coeffs = start_map_coeffs.apply_debye_waller_factors(u_star=u_star_s)

  # Generate map coefficient errors for first half-map from complex normal
  # distribution
  b_target_e = (0.,0.,0.,-50.,-75.,125.) # Anisotropy for error terms
  u_star_e = adptbx.u_cart_as_u_star(
      start_map_coeffs.unit_cell(), adptbx.b_as_u(b_target_e))
  se_target = 100. # Target for SigmaE variance term
  rsigma = math.sqrt(se_target / 2.)
  jj = 0.+1.j  # Define I for generating complex numbers
  random_complexes1 = flex.complex_double()
  ncoeffs=start_map_coeffs.size()
  random.seed(123457) # Make runs reproducible
  for i in range(ncoeffs):
    random_complexes1.append(random.gauss(0.,rsigma) + random.gauss(0.,rsigma)*jj)
  rc1_miller = start_map_coeffs.customized_copy(data=random_complexes1)
  mc1_delta = rc1_miller.apply_debye_waller_factors(u_star=u_star_e)
  map1_coeffs = scaled_map_coeffs.customized_copy(
    data=scaled_map_coeffs.data() + mc1_delta.data())

  # Repeat for second half map with independent errors from same distribution
  random_complexes2 = flex.complex_double()
  for i in range(ncoeffs):
    random_complexes2.append(random.gauss(0.,rsigma) + random.gauss(0.,rsigma)*jj)
  rc2_miller = start_map_coeffs.customized_copy(data=random_complexes2)
  mc2_delta = rc2_miller.apply_debye_waller_factors(u_star=u_star_e)
  map2_coeffs = scaled_map_coeffs.customized_copy(
    data=scaled_map_coeffs.data() + mc2_delta.data())

  # mmm.write_model("fake_map.pdb")
  mmm.add_map_from_fourier_coefficients(
      map1_coeffs, map_id = 'map_manager_1')
  mmm.add_map_from_fourier_coefficients(
      map2_coeffs, map_id = 'map_manager_2')
  # Replace original map_manager with mean of half-maps
  mm_mean_data = (mmm.map_manager_1().map_data() + mmm.map_manager_2().map_data()) / 2
  mmm.map_manager().set_map_data(map_data = mm_mean_data)
  # Add mask map for ordered component of map
  protein_mw = n_residues * 110. # MW from model would be better...
  nucleic_mw = None
  mask_id = 'ordered_volume_mask'
  add_ordered_volume_mask(mmm, d_min,
      protein_mw=protein_mw, nucleic_mw=nucleic_mw,
      map_id_out=mask_id)
  box_centre = tuple(flex.double((ucpars[0],ucpars[1],ucpars[2]))/2)
  # Now refine to assess parameters describing map errors
  results = assess_cryoem_errors(
      mmm, d_min,
      sphere_cent=tuple(box_centre), radius=model_radius+d_min, verbosity=0)

  # resultsdict = results.resultsdict
  # b_refined_a = resultsdict["a_baniso"]
  # print("\nIdeal A tensor as Baniso: ", b_expected)
  # print("Refined A tensor as Baniso", b_refined_a)

  # Note that all maps have been cut out with a spherical mask, so compare using these
  new_mmm = results.new_mmm
  perfect_mapCC = new_mmm.map_model_cc(map_id = 'perfect_map')
  mapCC = new_mmm.map_model_cc(map_id = 'map_manager_wtd') # Achieved map
  start_mapCC = new_mmm.map_model_cc() # Starting map with noise and anisotropy

  mc_perfect = new_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max, map_id='perfect_map')
  mc_achieved = new_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max, map_id='map_manager_wtd')

  # Compare with results using theoretically perfect error parameters to compute
  # ideal map coefficients.
  sigmaS_terms = flex.pow2(get_power_spectrum(mc_perfect)) # Actual signal power before anisotropy
  mc_start = new_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max)
  eE_ideal = mc_start.deep_copy()
  ones_array = flex.double(eE_ideal.size(), 1)
  all_ones = eE_ideal.customized_copy(data = ones_array)
  u_star_s2 = tuple(flex.double(u_star_s)*2.) # Square anisotropy for signal power calc
  sigmaS_terms = sigmaS_terms * all_ones.apply_debye_waller_factors(
    u_star=u_star_s2).data() # Corrected for anisotropy

  u_star_e2 = tuple(flex.double(u_star_e)*2.)
  sigmaE_terms = all_ones.apply_debye_waller_factors(u_star=u_star_e2).data() * se_target

  scale_terms = 1./flex.sqrt(sigmaS_terms + sigmaE_terms/2.)
  dobs_terms = 1./flex.sqrt(1. + sigmaE_terms/(2*sigmaS_terms))
  mc_ideal = eE_ideal.customized_copy(data = eE_ideal.data()*scale_terms*dobs_terms)
  # write_mtz(mc_achieved,"achieved_map.mtz","achieved")
  # write_mtz(mc_ideal,"ideal_map.mtz","ideal")

  mapCC_ideal_achieved = mc_ideal.map_correlation(other=mc_achieved)
  # print("CC between ideal and achieved maps:",mapCC_ideal_achieved)
  new_mmm.add_map_from_fourier_coefficients(
      mc_ideal, map_id = 'ideal_map')
  ideal_mapCC = new_mmm.map_model_cc(map_id = 'ideal_map')
  # print("Perfect, starting, ideal and achieved mapCC: ", perfect_mapCC, start_mapCC, ideal_mapCC, mapCC)
  assert(mapCC_ideal_achieved > 0.9)
  # slightly lower tolerance on older compiler
  import platform, sys
  if sys.platform.startswith('linux') \
    and ('.el6.' in platform.platform()
         or 'centos-6' in platform.platform()
         or ('azure' in platform.platform() and 'glibc2.10' in platform.platform())
         or ('azure' in platform.platform() and 'glibc2.12' in platform.platform())):
    assert(mapCC > 0.948*ideal_mapCC)
  else:
    assert(mapCC > 0.95*ideal_mapCC)


if(__name__ == "__main__"):
  exercise()
  print(format_cpu_times())
  print("OK")
