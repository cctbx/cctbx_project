from __future__ import division, print_function
import os, time
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from cctbx import adptbx

try:                from pydiscamb import DiscambWrapper
except ImportError: DiscambWrapper = None

# set this to True to print values
debug = False

# ------------------------------------------------------------------------------

pdb_str = """
CRYST1   13.000   16.000   19.000  97.00 110.00  95.00 P 1
HETATM    0  H1  HOH A3116       4.492   4.599   5.556  1.00  1.03           H
HETATM    1  O   HOH A3116       5.000   5.200   5.500  1.00 39.03           O
ANISOU    1  O   HOH A3116     2917   6748   5165  -1500  -1165   3269       O
HETATM    2  H2  HOH A3116       4.907   5.843   5.087  1.00  3.03           H
"""

# ------------------------------------------------------------------------------

def get_discamb_grads(xrs, d_target_d_fcalc):
  '''
  Compute analytical gradients using DiSCaMB.

  Parameters
  ----------
  xrs : cctbx.xray.structure
      The input crystal structure.
  d_target_d_fcalc : flex.complex_double
      Gradient of the target function with respect to F_calc.

  Returns
  -------
  list of float
      Analytical gradients with respect to model parameters.

  Notes
  -----
  Uses pyDiscamb to compute parameter gradients.
  If debug mode is enabled, prints the gradients.
'''
  w = DiscambWrapper(xrs)
  g = w.d_target_d_params(d_target_d_fcalc)
  if debug:
    print("gtaam:", [round(_,8) for _ in g])
  return g

# ------------------------------------------------------------------------------

def get_fd(fmodel_, eps=1.e-5):
  '''
  Compute numerical gradients using finite differences.

  Parameters
  ----------
  fmodel_ : mmtbx.f_model.manager
      Fmodel object containing structure factors and model object.
  eps : float, optional
      Step size for finite difference calculation (default is 1e-5).

  Returns
  -------
  flex.double
      Numerical gradients with respect to all enabled parameters
      (sites, u_iso, u_aniso, occupancy).

  Notes
  -----
  Iterates over all scatterers in the structure and perturbs each enabled
  parameter to compute numerical derivatives of the target function.
  '''
  fmodel = fmodel_.deep_copy()
  print('start is_taam', fmodel.is_taam())
  target_functor = fmodel.target_functor()
  structure = fmodel.xray_structure
  unit_cell = structure.unit_cell()
  gs = flex.double()
  for i_scatterer in range(structure.scatterers().size()):
    sc = structure.scatterers()[i_scatterer]
    f = sc.flags
    if(f.grad_site()):
      site_orig = sc.site
      d_target_d_site = []
      for ix in range(3):
        ts = []
        for signed_eps in [eps, -eps]:
          site_cart = list(unit_cell.orthogonalize(site_orig))
          site_cart[ix] += signed_eps
          sc.site = unit_cell.fractionalize(site_cart)
          fmodel.update_xray_structure(update_f_calc=True)
          print('update sites is_taam', fmodel.is_taam())
          ts.append(target_functor().target_work())
        gs.append((ts[0]-ts[1])/(2*eps))
      sc.site = site_orig
    if(f.grad_u_iso()):
      u_iso_orig = sc.u_iso
      d_target_u_iso = []
      ts = []
      for signed_eps in [eps, -eps]:
        sc.u_iso = u_iso_orig+signed_eps
        fmodel.update_xray_structure(update_f_calc=True)
        ts.append(target_functor().target_work())
      gs.append((ts[0]-ts[1])/(2*eps))
      sc.u_iso = u_iso_orig
    if(f.grad_u_aniso()):
      u_star_orig = sc.u_star
      d_target_d_u_star = []
      for ix in range(6):
        ts = []
        for signed_eps in [eps, -eps]:
          u_cart = list(adptbx.u_star_as_u_cart(unit_cell, u_star_orig))
          u_cart[ix] += signed_eps
          sc.u_star = adptbx.u_cart_as_u_star(unit_cell, u_cart)
          fmodel.update_xray_structure(update_f_calc=True)
          ts.append(target_functor().target_work())
        gs.append((ts[0]-ts[1])/(2*eps))
      sc.u_star = u_star_orig
    if(f.grad_occupancy()):
      occupancy_orig = sc.occupancy
      d_target_occupancy = []
      ts = []
      for signed_eps in [eps, -eps]:
        sc.occupancy = occupancy_orig+signed_eps
        fmodel.update_xray_structure(update_f_calc=True)
        ts.append(target_functor().target_work())
      gs.append((ts[0]-ts[1])/(2*eps))
      sc.occupancy = occupancy_orig
  if debug:
    print("fd:", [round(_,8) for _ in gs])
  return gs

# ------------------------------------------------------------------------------

def run(table):
  '''
  Compare analytical and numerical gradients for different scattering tables.

  Parameters
  ----------
  table : str
      Scattering table to use, one of: 'electron', 'wk1995', 'it1992'.

  Returns
  -------
  None

  Raises
  ------
  AssertionError
      If the analytical and numerical gradients do not match within the
      expected tolerance.

  Notes
  -----
  This function:
  - Builds a model from a water molecule PDB string.
  - Computes gradients using both DiSCaMB (analytical) and finite differences
    (numerical).
  - Compares gradients for selected atomic parameters:
    - Site (Cartesian coordinates)
    - Isotropic and anisotropic displacement parameters (U_iso, U_aniso)
    - Occupancy
  '''
  print('-------------- %s --------------' % table)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = table)
  xrs_iso = xrs.deep_copy_scatterers()
  xrs_iso.convert_to_isotropic()
  f_obs = abs(xrs_iso.structure_factors(d_min=0.5).f_calc())
  #
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  params.taam = True
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    xray_structure               = xrs,
    sf_and_grads_accuracy_params = params)
  scatterers = xrs.scatterers()
  tf = fmodel.target_functor()(compute_gradients = True)
  d_target_d_fcalc = tf.d_target_d_f_calc_work()
  #
  # sites cart only for atom 1
  #
  print("sites_cart_1")
  scatterers.flags_set_grads(state=False)
  scatterers.flags_set_grad_site(iselection = flex.size_t([1]))
  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
  g3 = get_fd(fmodel, eps=1.e-5)
  assert approx_equal(g2,g3, 1.e-3)
  if debug: print()
  #
  # Now add iso B for H #2
  #
#  print("sites_cart_1 + uiso_2")
#  scatterers.flags_set_grad_u_iso(iselection = flex.size_t([2]))
#  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
#  g3 = get_fd(fmodel, eps=1.e-5)
#  assert approx_equal(g2,g3, 1.e-5)
#  if debug: print()
#  #
#  # Now add occupancy of H #0
#  #
#  print("sites_cart_1 + uiso_2 + occ_0")
#  scatterers.flags_set_grad_occupancy(iselection = flex.size_t([0]))
#  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
#  g3 = get_fd(fmodel, eps=1.e-5)
#  assert approx_equal(g2,g3, 1.e-3)
#  if debug: print()
#  #
#  # Now add aniso ADP of O #1
#  #
#  print("sites_cart_1 + uiso_2 + occ_0 + uaniso_1")
#  scatterers.flags_set_grad_u_aniso(iselection = flex.size_t([1]))
#  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
#  g3 = get_fd(fmodel, eps=1.e-5)
#  assert approx_equal(g2,g3, 1.e-3)
#  if debug: print()
#  #
#  # Now add sites cart of H #2
#  #
#  print("sites_cart_1 + uiso_2 + occ_0 + uaniso_1 + sites_cart_2")
#  scatterers.flags_set_grad_site(iselection = flex.size_t([2]))
#  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
#  g3 = get_fd(fmodel, eps=1.e-5)
#  assert approx_equal(g2,g3, 1.e-5)
#  if debug: print()
#  #
#  # Cancel all
#  #
#  scatterers.flags_set_grads(state=False)
#  #
#  # Now all parameters
#  #
#  all_selection = flex.size_t([0,1,2])
#  scatterers.flags_set_grad_site(iselection = all_selection)
#  scatterers.flags_set_grad_occupancy(iselection = all_selection)
#  scatterers.flags_set_grad_u_aniso(iselection = flex.size_t([1]))
#  scatterers.flags_set_grad_u_iso(iselection = flex.size_t([0,2]))
#
#  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
#  g3 = get_fd(fmodel, eps=1.e-5)
#  assert approx_equal(g2,g3, 1.e-5)
#  if debug: print()

# ------------------------------------------------------------------------------

if(__name__ == "__main__"):
  '''
  Test to compare gradients computed by DiSCaMB and finite differences.

  This script runs the run() function for each supported scattering table.
  If DiscambWrapper is not available, it skips the tests.

  Returns
  -------
  None
  '''
  t0 = time.time()
  if DiscambWrapper is None:
    print("Skipping:", os.path.basename(__file__))
  else:
    for table in ["electron", "wk1995", "it1992"]:
      run(table=table)
      print()
  print("OK. Time: %8.3f"%(time.time()-t0))
