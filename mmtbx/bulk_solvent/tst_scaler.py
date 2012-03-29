from libtbx.test_utils import approx_equal
import mmtbx.f_model
import random, time
from scitbx.array_family import flex
from mmtbx import bulk_solvent
from cctbx import adptbx
from cctbx import sgtbx
from cctbx.development import random_structure
import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")
from mmtbx.bulk_solvent import scaler

if(1):
  random.seed(0)
  flex.set_random_seed(0)

def get_random_and_traceless_b_cart(crystal_symmetry):
  symbol = crystal_symmetry.space_group().type().lookup_symbol()
  point_group = sgtbx.space_group_info(
    symbol=symbol).group().build_derived_point_group()
  adp_constraints = sgtbx.tensor_rank_2_constraints(
    space_group=point_group,
    reciprocal_space=True)
  u_star = adptbx.u_cart_as_u_star(crystal_symmetry.unit_cell(),
    adptbx.random_u_cart(u_scale=1,u_min=0.1))
  u_indep = adp_constraints.independent_params(all_params=u_star)
  u_star = adp_constraints.all_params(independent_params=u_indep)
  b_cart = adptbx.u_as_b(adptbx.u_star_as_u_cart(
    crystal_symmetry.unit_cell(), u_star))
  tr = (b_cart[0]+b_cart[1]+b_cart[2])/3
  b_cart = [b_cart[0]-tr, b_cart[1]-tr, b_cart[2]-tr,
           b_cart[3],b_cart[4],b_cart[5]]
  return b_cart

def run_0(symbol = "C 2"):
  space_group_info = sgtbx.space_group_info(symbol = symbol)
  xrs = random_structure.xray_structure(
    space_group_info  = space_group_info,
    elements          = ["N"]*50,
    volume_per_atom   = 100.0,
    random_u_iso      = True)
  #
  b_cart = get_random_and_traceless_b_cart(crystal_symmetry=xrs.crystal_symmetry())
  u_star = adptbx.u_cart_as_u_star(xrs.unit_cell(), adptbx.b_as_u(b_cart))
  #
  F = xrs.structure_factors(d_min = 1.0).f_calc()
  k_anisotropic = mmtbx.f_model.ext.k_anisotropic(F.indices(), u_star)
  bin_selections = scaler.binning(
    unit_cell      = F.unit_cell(),
    miller_indices = F.indices())
  #
  d_spacings = F.d_spacings().data()
  ss = 1./flex.pow2(d_spacings) / 4.
  k_mask_tmp = mmtbx.f_model.ext.k_mask(ss, 0.35, 80.)
  k_mask = flex.double(F.data().size(), 0)
  k_isotropic = flex.double(F.data().size(), 0)
  for s in bin_selections:
    d = d_spacings.select(s)
    k_mask.set_selected(s, flex.mean(k_mask_tmp.select(s)))
    k_isotropic.set_selected(s, random.randint(1,10))
  #
  fmodel = mmtbx.f_model.manager(
    xray_structure = xrs,
    f_obs          = abs(F),
    k_isotropic    = k_isotropic,
    k_anisotropic  = k_anisotropic,
    k_mask         = k_mask)
  f_calc  = fmodel.f_calc()
  f_masks = fmodel.f_masks()
  f_model = fmodel.f_model()
  f_obs   = abs(f_model)
  r_free_flags = f_obs.generate_r_free_flags()
  #
  assert approx_equal(bulk_solvent.r_factor(f_obs.data(), f_model.data()), 0)
  aso = scaler.run(
    f_obs          = f_obs,
    f_calc         = f_calc,
    f_mask         = f_masks,
    r_free_flags   = r_free_flags,
    bin_selections = bin_selections,
    number_of_cycles = 500,
    auto_convergence_tolerance = 1.e-6,
    ss             = ss,
    try_poly       = True,
    try_expanal    = True,
    try_expmin     = True,
    verbose        = False)
  assert approx_equal(aso.r_final, 0.00030, 0.00001)
  assert approx_equal(aso.r_low,   0.00001, 0.00001)
  assert approx_equal(aso.r_high,  0.00026, 0.00001)
  assert approx_equal(
    bulk_solvent.r_factor(f_obs.data(), abs(aso.core.f_model).data(), 1),
    bulk_solvent.r_factor(f_obs.data(), abs(aso.core.f_model).data()))

if (__name__ == "__main__"):
  t0 = time.time()
  run_0()
  print "Time: %6.4f"%(time.time()-t0)
  print "OK"
