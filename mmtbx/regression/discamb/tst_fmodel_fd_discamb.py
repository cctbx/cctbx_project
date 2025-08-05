from __future__ import division, print_function
import os, time
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from cctbx import adptbx

try:
  from pydiscamb import DiscambWrapper
  import pydiscamb
except ImportError: DiscambWrapper = None

debug = True

pdb_str = """
CRYST1   13.000   16.000   19.000  97.00 110.00  95.00 P 1
HETATM    0  H1  HOH A3116       4.492   4.599   5.556  1.00  1.03           H
HETATM    1  O   HOH A3116       5.000   5.200   5.500  1.00 39.03           O
ANISOU    1  O   HOH A3116     2917   6748   5165  -1500  -1165   3269       O
HETATM    2  H2  HOH A3116       4.907   5.843   5.087  1.00  3.03           H
"""

# ------------------------------------------------------------------------------

def get_discamb_grads(xrs, fmodel):
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
  w.set_indices(fmodel.f_obs().indices())
  data = flex.complex_double(w.f_calc())
  fc = fmodel.f_obs().customized_copy(data = data)
  fmodel.update(f_calc = fc)
  tf = fmodel.target_functor()(compute_gradients = True)
  d_target_d_fcalc = tf.d_target_d_f_calc_work()
  g = w.d_target_d_params(d_target_d_fcalc)
  if debug:
    print("gdiscamb:", [round(_,8) for _ in g])
  return g

# ------------------------------------------------------------------------------

def get_cctbx_grads(xrs, tf):
  '''
  Compute analytical gradients using CCTBX.

  Parameters
  ----------
  xrs : cctbx.xray.structure
      The crystal structure with enabled gradient flags.
  tf : mmtbx.f_model.target_functor
      Target functor object from fmodel used to compute gradients.

  Returns
  -------
  flex.double
      Packed gradients with respect to atomic parameters.

  Notes
  -----
  Uses CCTBX internal functions to compute gradients for enabled parameters.
  '''
  g = tf.gradients_wrt_atomic_parameters().packed()
  if debug:
    print("n_parameters:", xrs.n_parameters())
    print("cctbx:", [round(_,8) for _ in g])
  return g

# ------------------------------------------------------------------------------

def get_fd(fmodel_, eps=1.e-5, use_discamb=False):
  '''
  Compute numerical gradients via finite differences.

  Parameters
  ----------
  fmodel_ : mmtbx.f_model.manager
      fmodel object containing structure factors and model object.
  eps : float, optional
      Step size for finite difference calculation (default is 1e-5).

  Returns
  -------
  flex.double
      Finite difference approximations of the gradients.

  Notes
  -----
  For each enabled parameter in the structure, this function perturbs the value
  positively and negatively by eps, recomputes the target function, and uses
  central difference to approximate the gradient.
  '''
  fmodel = fmodel_.deep_copy()
  #
  def update_fc():
    if use_discamb:
      pdw = pydiscamb.DiscambWrapper(
        fmodel.xray_structure,
        method = pydiscamb.FCalcMethod.IAM)
      pdw.set_indices(fmodel.f_obs().indices())
      data = flex.complex_double(pdw.f_calc())
      fc = fmodel.f_obs().customized_copy(data = data)
      fmodel.update(f_calc = fc)
    else:
      fmodel.update_xray_structure(update_f_calc=True)
  #
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
          #fmodel.update_xray_structure(update_f_calc=True)
          update_fc()
          ts.append(target_functor().target_work())
        gs.append((ts[0]-ts[1])/(2*eps))
      sc.site = site_orig
    if(f.grad_u_iso()):
      u_iso_orig = sc.u_iso
      d_target_u_iso = []
      ts = []
      for signed_eps in [eps, -eps]:
        sc.u_iso = u_iso_orig+signed_eps
        #fmodel.update_xray_structure(update_f_calc=True)
        update_fc()
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
          #fmodel.update_xray_structure(update_f_calc=True)
          update_fc()
          ts.append(target_functor().target_work())
        gs.append((ts[0]-ts[1])/(2*eps))
      sc.u_star = u_star_orig
    if(f.grad_occupancy()):
      occupancy_orig = sc.occupancy
      d_target_occupancy = []
      ts = []
      for signed_eps in [eps, -eps]:
        sc.occupancy = occupancy_orig+signed_eps
        #fmodel.update_xray_structure(update_f_calc=True)
        update_fc()
        ts.append(target_functor().target_work())
      gs.append((ts[0]-ts[1])/(2*eps))
      sc.occupancy = occupancy_orig
  if debug:
    print("fd:", [round(_,8) for _ in gs])
  return gs

def run(table, fd_delta = 1.e-5):
  """
  Check packed gradinets CCTBX (IAM) vs DiSCaMB (IAM)
  """
  print('-------------- %s --------------' % table)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = table)
  xrs_iso = xrs.deep_copy_scatterers()
  xrs_iso.convert_to_isotropic()
  f_obs = abs(xrs_iso.structure_factors(d_min=1.0).f_calc())
  #
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  params.exp_table_one_over_step_size=0
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    xray_structure               = xrs,
    target_name                  = "ls_wunit_kunit",
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
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, fmodel)
  g3 = get_fd(fmodel, eps=fd_delta)
  g4 = get_fd(fmodel, eps=fd_delta, use_discamb=True)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g1,g3, 1.e-6)
  assert approx_equal(g1,g4, 1.e-3)
  if debug: print()
  #
  # Now add iso B for H #2
  #
  print("sites_cart_1 + uiso_2")
  scatterers.flags_set_grad_u_iso(iselection = flex.size_t([2]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, fmodel)
  g3 = get_fd(fmodel, eps=fd_delta)
  g4 = get_fd(fmodel, eps=fd_delta, use_discamb=True)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g1,g3, 1.e-6)
  assert approx_equal(g1,g4, 1.e-3)
  if debug: print()
  #
  # Now add occupancy of H #0
  #
  print("sites_cart_1 + uiso_2 + occ_0")
  scatterers.flags_set_grad_occupancy(iselection = flex.size_t([0]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, fmodel)
  g3 = get_fd(fmodel, eps=fd_delta)
  g4 = get_fd(fmodel, eps=fd_delta, use_discamb=True)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g1,g3, 1.e-6)
  assert approx_equal(g1,g4, 1.e-3)
  if debug: print()
  #
  # Now add aniso ADP of O #1
  #
  print("sites_cart_1 + uiso_2 + occ_0 + uaniso_1")
  scatterers.flags_set_grad_u_aniso(iselection = flex.size_t([1]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, fmodel)
  g3 = get_fd(fmodel, eps=fd_delta)
  g4 = get_fd(fmodel, eps=fd_delta, use_discamb=True)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g1,g3, 1.e-6)
  assert approx_equal(g1,g4, 1.e-3)
  if debug: print()
  #
  # Now add sites cart of H #2
  #
  print("sites_cart_1 + uiso_2 + occ_0 + uaniso_1 + sites_cart_2")
  scatterers.flags_set_grad_site(iselection = flex.size_t([2]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, fmodel)
  g3 = get_fd(fmodel, eps=fd_delta)
  g4 = get_fd(fmodel, eps=fd_delta, use_discamb=True)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g1,g3, 1.e-6)
  assert approx_equal(g1,g4, 1.e-3)
  if debug: print()
  #
  # Cancell all
  #
  scatterers.flags_set_grads(state=False)
  #
  # Now all parameters
  #
  all_selection = flex.size_t([0,1,2])
  scatterers.flags_set_grad_site(iselection = all_selection)
  scatterers.flags_set_grad_occupancy(iselection = all_selection)
  scatterers.flags_set_grad_u_aniso(iselection = flex.size_t([1]))
  scatterers.flags_set_grad_u_iso(iselection = flex.size_t([0,2]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, fmodel)
  g3 = get_fd(fmodel, eps=fd_delta)
  g4 = get_fd(fmodel, eps=fd_delta, use_discamb=True)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g2,g3, 1.e-6)
  assert approx_equal(g1,g4, 1.e-2)
  if debug: print()

# ------------------------------------------------------------------------------

if(__name__ == "__main__"):
  '''
  Compare gradients from CCTBX, DiSCaMB, and finite differences.

  Parameters
  ----------
  table : str
      Scattering factor table to use. Must be one of: 'electron', 'wk1995', 'it1992'.

  Returns
  -------
  None

  Raises
  -------
  AssertionError
      If gradients from the three methods differ beyond the allowed tolerance.

  Notes
  -----
  For a test water molecule (with anisotropic O and isotropic H), the function:
  - Enables different parameter gradients (site, u_iso, u_aniso, occupancy).
  - Computes gradients using:
      1. CCTBX analytical method,
      2. DiSCaMB analytical method,
      3. Finite difference method.
  - Compares all three methods for consistency.
  '''
  t0 = time.time()
  if DiscambWrapper is None:
    print("Skipping:", os.path.basename(__file__))
  else:
    for table in ["electron", "wk1995", "it1992"]:
      run(table=table)
      print()
  print("OK. Time: %8.3f"%(time.time()-t0))
