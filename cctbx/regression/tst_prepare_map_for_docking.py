from __future__ import print_function
from __future__ import division
import math
from iotbx.map_model_manager import map_model_manager
from cctbx.maptbx.prepare_map_for_docking import add_ordered_volume_mask
from cctbx.maptbx.prepare_map_for_docking import run_refine_cryoem_errors
from cctbx import adptbx
from scitbx.array_family import flex
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
import random

def exercise():
  """Test prepare_map_for_docking using data with known errors."""

  # Generate two half-maps with same anisotropic signal,
  # independent anisotropic noise. Test to see if
  # simulation parameters are recovered.

  # Start by working out how large the padding will have to be so that
  # starting automatically-generated map will be large enough to contain
  # sphere with room to spare around model.
  n_residues = 25
  d_min = 2.5
  from cctbx.development.create_models_or_maps import generate_model
  test_model = generate_model(n_residues=n_residues)
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
  b_model = (30.,30.,30.,0.,0.,0.)  # All atoms in model have B=30
  u_star = adptbx.u_cart_as_u_star(
      start_map_coeffs.unit_cell(), adptbx.b_as_u(b_target))
  b_expected = list((flex.double(b_target) + flex.double(b_model)))
  scaled_map_coeffs = start_map_coeffs.apply_debye_waller_factors(u_star=u_star)

  # Generate map coefficient errors for first half-map from complex normal
  # distribution
  b_target_e = (0.,0.,0.,-50.,-50.,100.) # Anisotropy for error terms
  u_star_e = adptbx.u_cart_as_u_star(
      start_map_coeffs.unit_cell(), adptbx.b_as_u(b_target_e))
  se_target = 10. # Target for SigmaE variance term
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
  # mm1 = mmm.map_manager_1()
  # mm2 = mmm.map_manager_2()
  # mm_mean_data = (mm1.map_data() + mm2.map_data()) / 2
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
  results = run_refine_cryoem_errors(
      mmm, d_min,
      sphere_cent=tuple(box_centre), radius=model_radius+d_min, verbosity=0)

  resultsdict = results.resultsdict
  b_refined_a = resultsdict["a_baniso"]
  # print("\nIdeal A tensor as Baniso: ", b_expected)
  # print("Refined A tensor as Baniso", b_refined_a)
  b_refined_e = resultsdict["sigmaE_baniso"] # Note: refers to intensity scale
  # print("\nIdeal SigmaE tensor as Baniso: ",list(flex.double(b_target_e)*2))
  # print("Refined SigmaE tensor as Baniso", b_refined_e)
  for i in range(6):
    assert(approx_equal(b_refined_a[i],b_expected[i], eps=25))
    assert(approx_equal(b_refined_e[i],2*b_target_e[i], eps=25))

if(__name__ == "__main__"):
  exercise()
  print(format_cpu_times())
  print("OK")
