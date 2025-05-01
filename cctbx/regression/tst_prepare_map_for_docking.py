from __future__ import print_function
from __future__ import division
import math
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager
from cctbx.maptbx.prepare_map_for_docking import add_ordered_volume_mask
from cctbx.maptbx.prepare_map_for_docking import assess_cryoem_errors
#from cctbx.maptbx.prepare_map_for_docking import write_mtz
from cctbx import adptbx
from scitbx.array_family import flex
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
import os
import random

def get_power_spectrum(mc):

  power_spectrum = flex.double(mc.size(),1.)
  mc_copy = mc.deep_copy()
  nref = mc_copy.size()
  num_per_bin = 1000
  max_bins = 100
  min_bins = 6
  n_bins = int(round(max(min(nref / num_per_bin, max_bins), min_bins)))
  mc_copy.setup_binner(n_bins=n_bins)
  for i_bin in mc_copy.binner().range_used():
    sel = mc_copy.binner().selection(i_bin)
    mcsel = mc_copy.select(sel)
    fsq = flex.pow2(flex.abs(mcsel.data()))
    meanfsq = flex.mean_default(fsq, 0)
    power_spectrum.set_selected(sel,meanfsq)

  return power_spectrum

def exercise():
  """Test prepare_map_for_docking using data with known errors."""

  # Generate two half-maps with same anisotropic signal, independent anisotropic
  # noise. Test to see how well optimal map coefficients are estimated.

  # Use reasonably spherical model placed in the centre of a cubic unit cell
  # twice the maximum extent of the model, to simulate the situation with a
  # typical cryo-EM map.

  # Set flag for whether to print out debugging information instead of carrying
  # out regression tests with assert statements.
  # Should be set as False for released version!
  debug = False

  import libtbx.load_env
  iotbx_regression = os.path.join(libtbx.env.find_in_repositories("iotbx"),
      'regression')
  file_name=os.path.join(iotbx_regression,'data', 'big_cube_model.pdb')
  model_radius = 30. # Applies to this chosen model
  protein_mw = 19440.
  nucleic_mw = None

  d_min = 3.
  dm = DataManager()
  start_model = dm.get_model(file_name)

  # Reset B-values from zero to chosen constant
  b_iso = 30.
  b_values=flex.double(start_model.get_sites_cart().size(), b_iso)
  ph = start_model.get_hierarchy()
  ph.atoms().set_b(b_values)

  mmm = map_model_manager()
  mmm.generate_map(
      model=start_model,
      d_min=d_min, k_sol=0.1, b_sol=50.)

  if debug:
    mmm.write_model("fake_map.pdb")

  # Turn starting map into map coeffs for the signal
  start_unit_cell = mmm.map_manager().unit_cell()
  ucpars = start_unit_cell.parameters()
  d_max=max(ucpars[0], ucpars[1], ucpars[2])
  start_vol = ucpars[0] * ucpars[1] * ucpars[2]
  start_map_coeffs = mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max)

  # Keep copy of perfect map (map for B=0 model) for tests of success
  # Add it to mmm so that it has the same masking applied as test map.
  mc_perfect_ori = start_map_coeffs.deep_copy()
  mc_perfect_ori = mc_perfect_ori.apply_debye_waller_factors(b_iso = -b_iso)
  mmm.add_map_from_fourier_coefficients(mc_perfect_ori, map_id='perfect_map')

  # Apply anisotropic scaling to map coeffs
  b_target = (50.,100.,150.,-25.,25.,50.)
  u_star_s = adptbx.u_cart_as_u_star(
      start_unit_cell, adptbx.b_as_u(b_target))
  b_model = (b_iso,b_iso,b_iso,0.,0.,0.)  # All atoms in model have B=b_iso
  b_expected = list((flex.double(b_target) + flex.double(b_model)))
  scaled_map_coeffs = start_map_coeffs.apply_debye_waller_factors(u_star=u_star_s)

  # Generate map coefficient errors for first half-map from complex normal
  # distribution
  b_target_e = (0.,0.,0.,-50.,-50.,50.) # Anisotropy for error terms
  u_star_e = adptbx.u_cart_as_u_star(
      start_unit_cell, adptbx.b_as_u(b_target_e))
  se_target = 10000. # Target for SigmaE variance term
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

  mmm.add_map_from_fourier_coefficients(
      map1_coeffs, map_id = 'map_manager_1')
  mmm.add_map_from_fourier_coefficients(
      map2_coeffs, map_id = 'map_manager_2')
  # Replace original map_manager with mean of half-maps
  mm_mean_data = (mmm.map_manager_1().map_data() + mmm.map_manager_2().map_data()) / 2
  mmm.map_manager().set_map_data(map_data = mm_mean_data)

  # Add mask map for ordered component of map
  ordered_mask_id = 'ordered_volume_mask'
  add_ordered_volume_mask(mmm, d_min,
      protein_mw=protein_mw, nucleic_mw=nucleic_mw,
      ordered_mask_id=ordered_mask_id)
  box_centre = tuple(flex.double((ucpars[0],ucpars[1],ucpars[2]))/2)

  # Now refine to assess parameters describing map errors
  if debug:
    verbosity = 1
  else:
    verbosity = 0
  results = assess_cryoem_errors(
                  mmm=mmm,
                  d_min=d_min,
                  half_maps_provided=True,
                  ordered_mask_id=ordered_mask_id,
                  sphere_cent=box_centre,
                  radius=model_radius+d_min,
                  verbosity=verbosity)

  if debug:
    resultsdict = results.resultsdict
    b_refined_a = resultsdict["a_baniso"]
    print("\nIdeal   A tensor as Baniso: ", b_expected)
    print(  "Refined A tensor as Baniso: ", b_refined_a)

  # Note that all maps have been cut out with a spherical mask, so compare using these
  new_mmm = results.new_mmm
  perfect_mapCC = new_mmm.map_model_cc(map_id = 'perfect_map')
  achieved_mapCC = new_mmm.map_model_cc(map_id = 'map_manager_wtd') # Achieved map
  start_mapCC = new_mmm.map_model_cc() # Starting map with noise and anisotropy
  new_unit_cell = new_mmm.map_manager().crystal_symmetry().unit_cell()
  ucpars = new_unit_cell.parameters()
  d_max = max(ucpars[0], ucpars[1], ucpars[2])
  end_vol = 4./3. * math.pi * model_radius**3 # Random error left only in cut-out sphere
  volume_ratio = end_vol/start_vol
  mc_perfect = new_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max, map_id='perfect_map')
  mc_achieved = new_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max, map_id='map_manager_wtd')

  # Use theoretically perfect error parameters to compute ideal map coefficients.
  # First put back overall b_iso, then model the application of anisotropy and
  # addition of errors to the starting power spectrum
  mc_perfect_unsharp = mc_perfect.apply_debye_waller_factors(b_iso = b_iso)
  sigmaS_terms = get_power_spectrum(mc_perfect_unsharp) # Actual signal power before anisotropy
  ones_array = flex.double(mc_perfect.size(), 1)
  all_ones = mc_perfect.customized_copy(data = ones_array)
  # Recalculate u_star_s and u_star_e appropriate for new cut-out cell
  u_star_s = adptbx.u_cart_as_u_star(new_unit_cell, adptbx.b_as_u(b_target))
  u_star_s2 = tuple(flex.double(u_star_s)*2.) # Anisotropy of signal power
  u_star_e = adptbx.u_cart_as_u_star(new_unit_cell, adptbx.b_as_u(b_target_e))
  u_star_e2 = tuple(flex.double(u_star_e)*2.) # Anisotropy of error power
  sigmaS_terms = sigmaS_terms * all_ones.apply_debye_waller_factors(
      u_star=u_star_s2).data() # Corrected for anisotropy
  sigmaE_terms = all_ones.apply_debye_waller_factors(u_star=u_star_e2).data() * se_target * volume_ratio
  scale_terms = 1./flex.sqrt(sigmaS_terms + sigmaE_terms/2.)
  dobs_terms = flex.sqrt(sigmaS_terms / (sigmaS_terms + sigmaE_terms/2.))

  # Now take starting map, with anisotropy and average errors, and apply ideal
  # correction factors to get E and Dobs
  eE_ideal = new_mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max)
  eE_ideal = eE_ideal.customized_copy(data = eE_ideal.data()*scale_terms)
  mean_Esqr_ideal = flex.mean_default(flex.pow2(flex.abs(eE_ideal.data())),0)
  mc_ideal_wtd = eE_ideal.customized_copy(data = eE_ideal.data()*dobs_terms)

  # Now get results from assessment of signal and errors
  expectE = results.expectE
  mean_Esqr = flex.mean_default(flex.pow2(flex.abs(expectE.data())),0)
  eps_mean_Esqr = 0.15
  if debug:
    print("\n\nRegression tests:")
    print("\nMean value of ideal E**2: ", mean_Esqr_ideal)
    print("Mean value of E**2 for docking map: ", mean_Esqr)
    print("   Target for each is 1 with an allowed deviation of ",eps_mean_Esqr)
  else:
    assert approx_equal(mean_Esqr_ideal, 1.0, eps=eps_mean_Esqr)
    assert approx_equal(mean_Esqr, 1.0, eps=eps_mean_Esqr)

  target_ideal_achieved = 0.95
  mapCC_ideal_achieved = mc_ideal_wtd.map_correlation(other=mc_achieved)
  if debug:
    print("\nCC between ideal and achieved maps:",mapCC_ideal_achieved)
    print("   Target is >= ", target_ideal_achieved)
  new_mmm.add_map_from_fourier_coefficients(
      mc_ideal_wtd, map_id = 'ideal_map')
  ideal_mapCC = new_mmm.map_model_cc(map_id = 'ideal_map')
  target_achieved = 0.95*ideal_mapCC
  if debug:
    print("\nPerfect, starting and ideal mapCC: ", perfect_mapCC, start_mapCC, ideal_mapCC)
    print("Achieved mapCC: ", achieved_mapCC)
    print("   Target is >= ", target_achieved)
  else:
    assert(mapCC_ideal_achieved > target_ideal_achieved)
    assert(achieved_mapCC > target_achieved)

if(__name__ == "__main__"):
  exercise()
  print(format_cpu_times())
  print("OK")
