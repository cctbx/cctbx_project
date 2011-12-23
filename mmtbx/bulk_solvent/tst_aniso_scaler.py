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
from mmtbx import bulk_solvent

if(1):
  random.seed(0)
  flex.set_random_seed(0)

def run_00():
  time_aniso_u_scaler = 0
  for symbol in sgtbx.bravais_types.acentric + sgtbx.bravais_types.centric:
    #print symbol, "-"*50
    space_group_info = sgtbx.space_group_info(symbol = symbol)
    xrs = random_structure.xray_structure(
      space_group_info  = space_group_info,
      elements          = ["N"]*100,
      volume_per_atom   = 50.0,
      random_u_iso      = True)
    # XXX ad a method to adptbx to do this
    point_group = sgtbx.space_group_info(
      symbol=symbol).group().build_derived_point_group()
    adp_constraints = sgtbx.tensor_rank_2_constraints(
      space_group=point_group,
      reciprocal_space=True)
    u_star = adptbx.u_cart_as_u_star(xrs.unit_cell(),
      adptbx.random_u_cart(u_scale=1,u_min=0.1))
    u_indep = adp_constraints.independent_params(all_params=u_star)
    u_star = adp_constraints.all_params(independent_params=u_indep)
    b_cart_start=adptbx.u_as_b(adptbx.u_star_as_u_cart(xrs.unit_cell(), u_star))
    #
    tr = (b_cart_start[0]+b_cart_start[1]+b_cart_start[2])/3
    b_cart_start = [b_cart_start[0]-tr,b_cart_start[1]-tr,b_cart_start[2]-tr,
           b_cart_start[3],b_cart_start[4],b_cart_start[5]]
    tr = (b_cart_start[0]+b_cart_start[1]+b_cart_start[2])/3
    #
    #print "Input b_cart :", " ".join(["%8.4f"%i for i in b_cart_start]), "tr:", tr
    F = xrs.structure_factors(d_min = 2.0).f_calc()
    fmodel = mmtbx.f_model.manager(
      f_obs          = abs(F),
      xray_structure = xrs,
      b_cart         = b_cart_start)
    f_obs = abs(fmodel.f_model())
    fmodel = mmtbx.f_model.manager(
        f_obs          = f_obs,
        xray_structure = xrs,
        b_cart         = [0,0,0,0,0,0])
    t0 = time.time()
    obj = bulk_solvent.aniso_u_scaler(
      f_model        = fmodel.f_model().data(),
      f_obs          = fmodel.f_obs().data(),
      miller_indices = fmodel.f_obs().indices(),
      adp_constraint_matrix = adp_constraints.gradient_sum_matrix())
    time_aniso_u_scaler += (time.time()-t0)
    b_cart_final = adptbx.u_as_b(adptbx.u_star_as_u_cart(f_obs.unit_cell(),
      adp_constraints.all_params(tuple(obj.u_star_independent))))
    #print "Output b_cart:", " ".join(["%8.4f"%i for i in b_cart_final])
    assert approx_equal(b_cart_start, b_cart_final, 1.e-4)
  print "Time (aniso_u_scaler only): %6.4f"%time_aniso_u_scaler

def run_01():
  time_aniso_u_scaler = 0
  for symbol in sgtbx.bravais_types.acentric + sgtbx.bravais_types.centric:
    #print symbol, "-"*50
    space_group_info = sgtbx.space_group_info(symbol = symbol)
    xrs = random_structure.xray_structure(
      space_group_info  = space_group_info,
      elements          = ["N"]*100,
      volume_per_atom   = 50.0,
      random_u_iso      = True)
    # XXX ad a method to adptbx to do this
    point_group = sgtbx.space_group_info(
      symbol=symbol).group().build_derived_point_group()
    adp_constraints = sgtbx.tensor_rank_2_constraints(
      space_group=point_group,
      reciprocal_space=True)
    u_star = adptbx.u_cart_as_u_star(xrs.unit_cell(),
      adptbx.random_u_cart(u_scale=1,u_min=0.1))
    u_indep = adp_constraints.independent_params(all_params=u_star)
    u_star = adp_constraints.all_params(independent_params=u_indep)
    b_cart_start=adptbx.u_as_b(adptbx.u_star_as_u_cart(xrs.unit_cell(), u_star))
    #
    tr = (b_cart_start[0]+b_cart_start[1]+b_cart_start[2])/3
    b_cart_start = [b_cart_start[0]-tr,b_cart_start[1]-tr,b_cart_start[2]-tr,
           b_cart_start[3],b_cart_start[4],b_cart_start[5]]
    tr = (b_cart_start[0]+b_cart_start[1]+b_cart_start[2])/3
    #
    #print "Input b_cart :", " ".join(["%8.4f"%i for i in b_cart_start]), "tr:", tr
    F = xrs.structure_factors(d_min = 2.0).f_calc()
    fmodel = mmtbx.f_model.manager(
      f_obs          = abs(F),
      xray_structure = xrs,
      b_cart         = b_cart_start)
    f_obs = abs(fmodel.f_model())
    fmodel = mmtbx.f_model.manager(
        f_obs          = f_obs,
        xray_structure = xrs,
        b_cart         = [0,0,0,0,0,0])
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    obj = bulk_solvent.aniso_u_scaler(
      f_model        = fmodel.f_model().data(),
      f_obs          = fmodel.f_obs().data(),
      miller_indices = fmodel.f_obs().indices(),
      unit_cell      = f_obs.unit_cell())
    a = obj.a
    ####
    #print "Input a :", " ".join(["%7.3f"%i for i in a])
    overall_anisotropic_scale = ext.overall_anisotropic_scale(
      f_obs.indices(), a, f_obs.unit_cell())
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data()*overall_anisotropic_scale)
    fmodel = mmtbx.f_model.manager(
      f_obs          = abs(F),
      xray_structure = xrs)
    f_obs = abs(fmodel.f_model())
    f_obs = f_obs.customized_copy(data = f_obs.data() * overall_anisotropic_scale)
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    fmodel = mmtbx.f_model.manager(
        f_obs          = f_obs,
        xray_structure = xrs)
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    t0 = time.time()
    # XXX try "long doible" to see if this decreases the tolerances
    obj = bulk_solvent.aniso_u_scaler(
      f_model        = fmodel.f_model().data(),
      f_obs          = fmodel.f_obs().data(),
      miller_indices = fmodel.f_obs().indices(),
      unit_cell      = f_obs.unit_cell())
    time_aniso_u_scaler += (time.time()-t0)
    overall_anisotropic_scale = ext.overall_anisotropic_scale(
      f_obs.indices(), obj.a, f_obs.unit_cell())
    assert approx_equal(bulk_solvent.r_factor(f_obs.data(),
      fmodel.f_model().data()*overall_anisotropic_scale), 0.0, 1.e-2) # XXX seems to be low
    #print "Output a:", " ".join(["%7.3f"%i for i in obj.a])
    assert approx_equal(a, obj.a, 1.e-1) # XXX seems to be low
  print "Time (aniso_u_scaler only): %6.4f"%time_aniso_u_scaler

if (__name__ == "__main__"):
  t0 = time.time()
  run_00()
  run_01()
  print "Time: %6.4f"%(time.time()-t0)
  print "OK"
