from __future__ import division, print_function
import os, time
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from cctbx import adptbx

try:
  from pydiscamb import DiscambWrapper, FCalcMethod
except ImportError: DiscambWrapper = None

pdb_str = """
CRYST1   16.170   14.591   15.187  90.00  90.00  90.00 P 1
SCALE1      0.061843  0.000000  0.000000        0.00000
SCALE2      0.000000  0.068535  0.000000        0.00000
SCALE3      0.000000  0.000000  0.065846        0.00000
ATOM      1  N   TYR A   4       8.357   9.217   8.801  1.00 10.55           N
ATOM      2  CA  TYR A   4       9.150   8.055   9.050  1.00 10.24           C
ATOM      3  C   TYR A   4      10.419   8.399   9.804  1.00  9.86           C
ATOM      4  O   TYR A   4      10.726   9.591  10.050  1.00 11.39           O
ATOM      5  CB  TYR A   4       9.496   7.352   7.737  1.00 30.00           C
ATOM      6  CG  TYR A   4       8.296   6.791   7.006  1.00 30.00           C
ATOM      7  CD1 TYR A   4       7.820   5.517   7.293  1.00 30.00           C
ATOM      8  CD2 TYR A   4       7.642   7.534   6.034  1.00 30.00           C
ATOM      9  CE1 TYR A   4       6.724   5.000   6.628  1.00 30.00           C
ATOM     10  CE2 TYR A   4       6.545   7.024   5.364  1.00 30.00           C
ATOM     11  CZ  TYR A   4       6.091   5.758   5.665  1.00 30.00           C
ATOM     12  OH  TYR A   4       5.000   5.244   5.000  1.00 30.00           O
ATOM     13  OXT TYR A   4      11.170   7.502  10.187  1.00  9.86           O
ATOM     14  H1  TYR A   4       8.254   9.323   7.923  1.00 10.55           H
ATOM     15  H2  TYR A   4       7.560   9.120   9.184  1.00 10.55           H
ATOM     16  H3  TYR A   4       8.764   9.932   9.140  1.00 10.55           H
ATOM     17  HA  TYR A   4       8.637   7.441   9.599  1.00 10.24           H
ATOM     18  HB2 TYR A   4       9.929   7.989   7.148  1.00 30.00           H
ATOM     19  HB3 TYR A   4      10.097   6.615   7.928  1.00 30.00           H
ATOM     20  HD1 TYR A   4       8.245   5.005   7.942  1.00 30.00           H
ATOM     21  HD2 TYR A   4       7.946   8.389   5.830  1.00 30.00           H
ATOM     22  HE1 TYR A   4       6.415   4.146   6.829  1.00 30.00           H
ATOM     23  HE2 TYR A   4       6.116   7.532   4.714  1.00 30.00           H
ATOM     24  HH  TYR A   4       4.712   5.804   4.444  1.00 30.00           H
"""

# ------------------------------------------------------------------------------

def show(g1,g2,g3):
  print('wrapper from fmodel ', [round(_,8) for _ in g1])
  print('new wrapper from xrs', [round(_,8) for _ in g2])
  print('finite diff fmodel  ', [round(_,8) for _ in g3])

# ------------------------------------------------------------------------------

def get_wrapper_grads(xrs, d_target_d_fcalc):
  '''
  '''
  w = DiscambWrapper(xrs,
                     method = FCalcMethod.TAAM,
                     frozen_lcs = True)
  g = w.d_target_d_params(d_target_d_fcalc)
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
  return gs

def run(table, fd_delta = 1.e-5):
  """
  Check packed gradients
  """
  print('-------------- %s --------------' % table)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = table)
  xrs_iso = xrs.deep_copy_scatterers()
  xrs_iso.convert_to_isotropic()
  f_obs = abs(xrs_iso.structure_factors(d_min=0.7).f_calc())
  #
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "taam"
  params.extra.discamb.taam.freeze_local_coordinate_system = True
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
  scatterers.flags_set_grads(state=False)
  #
  # All parameters
  #
  all_selection = flex.size_t([0,1,2])
  scatterers.flags_set_grad_site(iselection = all_selection)
  scatterers.flags_set_grad_occupancy(iselection = all_selection)
  scatterers.flags_set_grad_u_iso(iselection = flex.size_t([0,1,2]))
  g1 = tf.gradients_wrt_atomic_parameters().packed()
  g2 = get_wrapper_grads(xrs, d_target_d_fcalc)
  g3 = get_fd(fmodel, eps=fd_delta)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g1,g3, 3.e-4)
  show(g1,g2,g3)

# ------------------------------------------------------------------------------

if(__name__ == "__main__"):
  '''
  Compare gradients from DiSCaMB and finite differences.

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
      If gradients from the methods differ beyond the allowed tolerance.

  Notes
  -----
  For a test tyr molecule (with anisotropic O and isotropic H), the function
  computes gradients using:
      1. DiSCaMB analytical method via fmodel object,
      2. DiSCaMB analytical method via a new wrapper instance,
      3. Finite difference method based on fmodel.
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
