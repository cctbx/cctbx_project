from __future__ import print_function
from __future__ import division
import os
import math
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce_cryoem_resolution
from scitbx.array_family import flex
from libtbx.utils import null_out
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
import numpy as np


def exercise(half_maps=True):
  """Test reduce_cryoem_resolution using data with defined errors."""

  # Generate two half-maps with same signal, independent noise, based on
  # calculated structure factors from a model. Add independent noise to the two
  # half-maps to match FSC curve corresponding to new lower-resolution limit.
  # Check that result is as expected.

  # Use reasonably spherical model placed in the centre of a cubic unit cell
  # twice the maximum extent of the model, to simulate the situation with a
  # typical cryo-EM map.

  # Set flag for whether to print out debugging information instead of carrying
  # out regression tests with assert statements.
  # Should be set as False for released version!
  debug = False
  if debug:
    plot = True
    mute = False
    logger = None
  else:
    plot = False
    mute = True
    logger = null_out()
  import libtbx.load_env
  iotbx_regression = os.path.join(libtbx.env.find_in_repositories("iotbx"),
      'regression')
  file_name=os.path.join(iotbx_regression,'data', 'big_cube_model.pdb')

  d_min = 3.
  target_d_min = 4.
  dm = DataManager()
  start_model = dm.get_model(file_name)

  # Reset B-values from zero to chosen constant
  b_iso = 10.
  b_values=flex.double(start_model.get_sites_cart().size(), b_iso)
  ph = start_model.get_hierarchy()
  ph.atoms().set_b(b_values)

  # Start with map_model_manager containing true map and two half-map copies
  # with no added error
  mmm = map_model_manager()
  mmm.generate_map(
      model=start_model,
      d_min=d_min, k_sol=0.1, b_sol=50.)

  # Turn starting map into map coeffs for the signal
  unit_cell = mmm.map_manager().unit_cell()
  ucpars = unit_cell.parameters()
  d_max = max(ucpars[0], ucpars[1], ucpars[2])
  start_map_coeffs = mmm.map_as_fourier_coefficients(d_min=d_min, d_max=d_max)
  nref_total = start_map_coeffs.size()
  nref_bin1 = 50. # target number of reflections in first bin
  nbins = max(6, round(nref_total / nref_bin1)**(2./3.))
  start_map_coeffs.setup_binner_d_star_sq_bin_size(max_bins=nbins)
  from mmtbx.programs.reduce_cryoem_resolution import apply_spherical_mask
  from mmtbx.programs.reduce_cryoem_resolution import fsc_params_from_d_min
  masked_mmm = apply_spherical_mask(mmm, target_d_min)
  mask_fraction = masked_mmm.mask_info().mean # Use to correct noise for full box
  filenames = []
  if half_maps:
    if debug:
      print("\n\nTesting with half maps\n")
    mc1 = start_map_coeffs.deep_copy()
    mc2 = start_map_coeffs.deep_copy()
    target_params = fsc_params_from_d_min(d_min)
    r_target = target_params.signalratio
    dB_target = target_params.deltaB
    nref_total = mc1.size()
    nref_bin1 = 50. # target number of reflections in first bin
    nbins = max(6, round(nref_total / nref_bin1)**(2./3.))
    mc1.setup_binner_d_star_sq_bin_size(max_bins=nbins)

    for i_bin in mc1.binner().range_used():
      sel = mc1.binner().selection(i_bin)
      mc1sel = mc1.select(sel)
      nrefsel = mc1sel.data().size()
      mean_ssqr = flex.mean(mc1sel.d_star_sq().data())
      sigmaS = flex.mean(flex.pow2(mc1sel.amplitudes().data()))
      sigmaE_factor = math.exp(dB_target * mean_ssqr / 4.) / r_target
      sigmaE_target = sigmaS * sigmaE_factor / mask_fraction
      target_sigma = math.sqrt(sigmaE_target/2.) # For real and imaginary parts
      random_complex = flex.complex_double(
                          np.random.normal(scale=target_sigma,size=nrefsel) +
                       1j*np.random.normal(scale=target_sigma,size=nrefsel))
      mc1.data().set_selected(sel, mc1.data().select(sel) + random_complex )
      random_complex = flex.complex_double(
                          np.random.normal(scale=target_sigma,size=nrefsel) +
                       1j*np.random.normal(scale=target_sigma,size=nrefsel))
      mc2.data().set_selected(sel, mc2.data().select(sel) + random_complex )

    mmm.add_map_from_fourier_coefficients(mc1, map_id='map_manager_1')
    mmm.add_map_from_fourier_coefficients(mc2, map_id='map_manager_2')
    working_mmm = apply_spherical_mask(mmm, target_d_min)
    wmc1 = working_mmm.map_manager_1().map_as_fourier_coefficients(d_min=d_min, d_max=1000.)
    wmc2 = working_mmm.map_manager_2().map_as_fourier_coefficients(d_min=d_min, d_max=1000.)
    from mmtbx.programs.reduce_cryoem_resolution import d_min_from_fsc_analytical
    d_min_results = d_min_from_fsc_analytical(wmc1, wmc2, guess_d_min=d_min)
    d_min_from_fsc = d_min_results.d_min
    if debug:
      print("Requested and actual d_min_from fsc for simulated half-maps: ",
            d_min, d_min_from_fsc)
    else:
      assert approx_equal(d_min, d_min_from_fsc, eps=0.1)

    map_fn1 = "tmp_map1.map"
    map_fn2 = "tmp_map2.map"
    working_mmm.map_manager_1().write_map(map_fn1)
    working_mmm.map_manager_2().write_map(map_fn2)
    args_rcr = [
      map_fn1,
      map_fn2,
      "d_min = " + str(d_min),
      "target_d_min = " + str(target_d_min),
      "file_root = with_halfmaps",
      "plot = " + str(plot),
      "mute = " + str(mute)
    ]
    results = run_program(program_class=reduce_cryoem_resolution.Program,
                          args=args_rcr, logger=logger)
    d_min_from_fsc = results.d_min_from_fsc
    filenames = filenames + results.filenames
    os.remove(map_fn1)
    os.remove(map_fn2)
  else:
    if debug:
      print("\n\nTesting with only full map\n")
    map_fn = "tmp_map.map"
    mmm.map_manager().write_map(map_fn)
    args_rcr = [
      map_fn,
      "target_d_min = " + str(target_d_min),
      "file_root = without_halfmaps",
      "plot = " + str(plot),
      "mute = " + str(mute)
    ]
    results = run_program(program_class=reduce_cryoem_resolution.Program,
                          args=args_rcr, logger=logger)
    d_min_from_fsc = results.d_min_from_fsc
    filenames = filenames + results.filenames
    os.remove(map_fn)

  if debug:
    print("Target and actual d_min_from_fsc for final half-maps: ",
          target_d_min, d_min_from_fsc)
  else:
    assert approx_equal(target_d_min, d_min_from_fsc, eps=0.1)

  if debug:
    print("\nFiles produced during tests:")
    for file in filenames:
      print("    ", file)
  else:
    for file in filenames:
      os.remove(file)


if(__name__ == "__main__"):
  for half_maps in (False, True):
    # print("\n\nTesting with half_maps = ",half_maps)
    exercise(half_maps)
  print(format_cpu_times())
  print("OK")
