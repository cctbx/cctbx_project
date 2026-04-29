from __future__ import division
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
from scitbx.array_family import flex
import sys, time
from libtbx.utils import null_out
from mmtbx.refinement import wrappers
import mmtbx.refinement.refinement_flags
import pydiscamb
from libtbx.test_utils import approx_equal
from cctbx import maptbx, miller
from mmtbx import map_tools

import random
random.seed(0)
flex.set_random_seed(0)

def compute_maps(fmodel, crystal_gridding, map_type):
  map_coefficients = map_tools.electron_density_map(
    fmodel = fmodel).map_coefficients(
      map_type         = map_type,
      isotropize       = True,
      fill_missing     = False)
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = map_coefficients)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded(), map_coefficients

pdb_str = """
CRYST1   20.584   14.037   18.865  90.00  90.00  90.00 P 1
ATOM      1  N   GLN A   5       6.376   6.692   8.279  1.00 20.00           N
ATOM      2  CA  GLN A   5       7.421   7.443   8.965  1.00 20.00           C
ATOM      3  C   GLN A   5       8.806   6.937   8.573  1.00 20.00           C
ATOM      4  O   GLN A   5       9.098   5.747   8.689  1.00 20.00           O
ATOM      5  CB  GLN A   5       7.238   7.355  10.481  1.00 20.00           C
ATOM      6  CG  GLN A   5       5.942   7.968  10.987  1.00 20.00           C
ATOM      7  CD  GLN A   5       5.554   7.459  12.361  1.00 20.00           C
ATOM      8  OE1 GLN A   5       5.133   6.313  12.515  1.00 20.00           O
ATOM      9  NE2 GLN A   5       5.695   8.312  13.369  1.00 20.00           N
ATOM     10  N   ASN A   6       9.655   7.851   8.108  1.00 20.00           N
ATOM     11  CA  ASN A   6      11.017   7.527   7.691  1.00 20.00           C
ATOM     12  C   ASN A   6      11.965   7.928   8.815  1.00 20.00           C
ATOM     13  O   ASN A   6      12.295   9.106   8.975  1.00 20.00           O
ATOM     14  CB  ASN A   6      11.366   8.231   6.384  1.00 20.00           C
ATOM     15  CG  ASN A   6      10.372   7.935   5.278  1.00 20.00           C
ATOM     16  OD1 ASN A   6      10.269   6.802   4.808  1.00 20.00           O
ATOM     17  ND2 ASN A   6       9.634   8.955   4.857  1.00 20.00           N
ATOM     18  N   TYR A   7      12.403   6.944   9.594  1.00 20.00           N
ATOM     19  CA  TYR A   7      13.314   7.194  10.705  1.00 20.00           C
ATOM     20  C   TYR A   7      14.752   7.322  10.215  1.00 20.00           C
ATOM     21  O   TYR A   7      15.093   6.871   9.122  1.00 20.00           O
ATOM     22  CB  TYR A   7      13.209   6.077  11.746  1.00 20.00           C
ATOM     23  CG  TYR A   7      11.853   5.980  12.406  1.00 20.00           C
ATOM     24  CD1 TYR A   7      11.499   6.836  13.441  1.00 20.00           C
ATOM     25  CD2 TYR A   7      10.925   5.032  11.995  1.00 20.00           C
ATOM     26  CE1 TYR A   7      10.260   6.751  14.047  1.00 20.00           C
ATOM     27  CE2 TYR A   7       9.683   4.939  12.595  1.00 20.00           C
ATOM     28  CZ  TYR A   7       9.357   5.801  13.621  1.00 20.00           C
ATOM     29  OH  TYR A   7       8.122   5.713  14.222  1.00 20.00           O
ATOM     30  OXT TYR A   7      15.608   7.881  10.901  1.00 20.00           O
TER
END
"""

def run(d_min            = 1.0,
        scattering_table = "electron",
        algorithm        = "direct",
        macro_cycles     = 10,
        max_iterations   = 50,
        shake_model      = True,
        use_cctbx        = True):
  # Ground truth model
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model.process(make_restraints = True)
  model.setup_scattering_dictionaries(
    scattering_table = scattering_table,
    d_min            = d_min)
  xrs = model.get_xray_structure()
  # Fobs from ground truth model
  f_obs = abs(xrs.structure_factors(d_min=d_min, algorithm=algorithm).f_calc())
  if not use_cctbx:
    wrapper = pydiscamb.DiscambWrapper(xrs, method=pydiscamb.FCalcMethod.TAAM)
    wrapper.set_indices(f_obs.indices())
    data = flex.complex_double(wrapper.f_calc())
    f_obs = f_obs.array(data = flex.abs(data))
  r_free_flags = f_obs.generate_r_free_flags()
  assert scattering_table == xrs.get_scattering_table(), [
    scattering_table, xrs.get_scattering_table()]
  # Generate list of trial sites_cart and shake B factors
  sites_cart_list = []
  xrs_dc = xrs.deep_copy_scatterers()
  if shake_model:
    xrs_dc.shake_sites_in_place(mean_distance=0.1)
    xrs_dc.shake_adp()
  sites_cart_list.append(xrs_dc.sites_cart())
  model.set_xray_structure(xrs_dc)
  model.setup_scattering_dictionaries( # deep_copy_scatterers looses the table!
    scattering_table = scattering_table,
    d_min            = d_min)
  # get fmodel from data manager
  xrs = model.get_xray_structure()
  assert scattering_table == xrs.get_scattering_table(), [
    scattering_table, xrs.get_scattering_table()]
  rf = mmtbx.refinement.refinement_flags.manager(
    individual_sites   = True,
    individual_adp     = True,
    sites_individual   = flex.bool(model.size(), True),
    adp_individual_iso = flex.bool(model.size(), True))
  model.set_refinement_flags(flags = rf)
  #
  sfparams = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfparams.algorithm = algorithm
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    target_name    = "ls_wunit_k1", # ml in reality
    sf_and_grads_accuracy_params = sfparams,
    xray_structure = model.get_xray_structure(),
    discamb_mode = None if use_cctbx else "taam")
  fmodel.update_all_scales()
  print("r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()))
  assert model.get_xray_structure() == fmodel.xray_structure
  #
  # Run refinement
  #
  sc_start = model.get_sites_cart().deep_copy()
  o = wrappers.simple_fsr(
    model           = model,
    fmodel          = fmodel,
    macro_cycles    = macro_cycles,
    max_iterations  = max_iterations,
    sites_cart_list = [model.get_sites_cart(),],
    log             = sys.stdout,
    refine_xyz      = True,
    refine_adp      = True,
    )
  o.run()
  # Sanity checks
  print()
  sc_final = model.get_sites_cart()
  print("shift from start:", flex.mean(flex.sqrt((sc_start - sc_final).dot())))
  assert approx_equal(o.results[-1].r_work, fmodel.r_work())
  fmodel.update_xray_structure(
    xray_structure = model.get_xray_structure(), update_f_calc=True)
  assert approx_equal(o.results[-1].r_work, fmodel.r_work())
  #
  print()
  fmodel.show(show_header=False, show_approx=False)
  #
  with open("refined.pdb", "w") as fo:
    fo.write(model.model_as_pdb())
  # Compute maps
  crystal_gridding = fmodel.f_obs().crystal_gridding(
    d_min             = d_min,
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 1./4)
  fofc_map, fofc = compute_maps(
    fmodel           = fmodel,
    crystal_gridding = crystal_gridding,
    map_type         = "mFo-DFc")
  tfofc_map, tfofc = compute_maps(
    fmodel           = fmodel,
    crystal_gridding = crystal_gridding,
    map_type         = "2mFo-DFc")
  # Output maps as MTZ
  mtz_dataset = fofc.as_mtz_dataset(column_root_label = "mFobs-DFmodel")
  mtz_dataset.add_miller_array(
    miller_array      = tfofc,
    column_root_label = "2mFobs-DFmodel")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "refined.mtz")

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.4f"%(time.time()-t0))
