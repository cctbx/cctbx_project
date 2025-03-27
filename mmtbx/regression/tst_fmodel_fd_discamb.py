from __future__ import division, print_function
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from cctbx import adptbx
import os

try:                from pydiscamb import DiscambWrapper
except ImportError: DiscambWrapper = None

pdb_str = """
CRYST1   13.000   16.000   19.000  97.00 110.00  95.00 P 1
HETATM    0  H1  HOH A3116       4.492   4.599   5.556  1.00  1.03           H
HETATM    1  O   HOH A3116       5.000   5.200   5.500  1.00 39.03           O
ANISOU    1  O   HOH A3116     2917   6748   5165  -1500  -1165   3269       O
HETATM    2  H2  HOH A3116       4.907   5.843   5.087  1.00  3.03           H
"""

def get_discamb_grads(xrs, d_target_d_fcalc):
  w = DiscambWrapper(xrs)
  g = w.d_target_d_params(d_target_d_fcalc)
  print("gdiscamb:", [round(_,8) for _ in g])
  return g

def get_cctbx_grads(xrs, tf):
  print("n_parameters:", xrs.n_parameters())
  g = tf.gradients_wrt_atomic_parameters().packed()
  print("cctbx:", [round(_,8) for _ in g])
  return g

def get_fd(fmodel_, eps=1.e-5):
  fmodel = fmodel_.deep_copy()
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
  print("fd:", [round(_,8) for _ in gs])
  return gs

def run(table):
  """
  Check packed gradinets CCTBX (IAM) vs DiSCaMB (IAM)
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = table)
  xrs_iso = xrs.deep_copy_scatterers()
  xrs_iso.convert_to_isotropic()
  f_obs = abs(xrs_iso.structure_factors(d_min=0.5).f_calc())
  #
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
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
  print("sites_cart")
  scatterers.flags_set_grads(state=False)
  scatterers.flags_set_grad_site(iselection = flex.size_t([1]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
  g3 = get_fd(fmodel, eps=1.e-5)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g1,g3, 1.e-6)
  print()
  #
  # Now add iso B for H #2
  #
  print("sites_cart1+uiso2")
  scatterers.flags_set_grad_u_iso(iselection = flex.size_t([2]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
  g3 = get_fd(fmodel, eps=1.e-5)
  assert approx_equal(g1,g2, 1.e-5)
  assert approx_equal(g1,g3, 1.e-5)
  print()
  #
  # Now add occupancy of H #0
  #
  scatterers.flags_set_grad_occupancy(iselection = flex.size_t([0]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
  g3 = get_fd(fmodel, eps=1.e-5)
  assert approx_equal(g1,g2, 1.e-5)
  assert approx_equal(g1,g3, 1.e-3)
  print()
  #
  # Now add aniso ADP of O #1
  #
  scatterers.flags_set_grad_u_aniso(iselection = flex.size_t([1]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
  g3 = get_fd(fmodel, eps=1.e-5)
  assert approx_equal(g1,g2, 1.e-6)
  assert approx_equal(g1,g3, 1.e-3)
  print()
  #
  # Now add sites cart of H #2
  #
  scatterers.flags_set_grad_site(iselection = flex.size_t([2]))
  g1 = get_cctbx_grads(xrs, tf)
  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
  assert approx_equal(g1,g2, 1.e-6)
  print()
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
  g2 = get_discamb_grads(xrs, d_target_d_fcalc)
  assert approx_equal(g1,g2, 1.e-6)
  print()

if(__name__ == "__main__"):
  if DiscambWrapper is None:
    print("Skipping:", os.path.basename(__file__))
  else:
    for table in ["electron", "wk1995", "it1992"]:
      run(table=table)
