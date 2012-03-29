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
    u_star = adptbx.u_cart_as_u_star(
      F.unit_cell(), adptbx.b_as_u(b_cart_start))
    fbc = mmtbx.f_model.ext.k_anisotropic(F.indices(), u_star)
    fc = F.structure_factors_from_scatterers(xray_structure=xrs).f_calc()
    f_obs = F.customized_copy(data = flex.abs(fc.data()*fbc))
    t0 = time.time()
    obj = bulk_solvent.aniso_u_scaler(
      f_model        = fc.data(),
      f_obs          = f_obs.data(),
      miller_indices = f_obs.indices(),
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
    F = xrs.structure_factors(d_min = 2.0).f_calc()
    u_star = adptbx.u_cart_as_u_star(
      F.unit_cell(), adptbx.b_as_u(b_cart_start))
    fbc = mmtbx.f_model.ext.k_anisotropic(F.indices(), u_star)
    fc = F.structure_factors_from_scatterers(xray_structure=xrs).f_calc()
    f_obs = F.customized_copy(data = flex.abs(fc.data()*fbc))
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    obj = bulk_solvent.aniso_u_scaler(
      f_model        = fc.data(),
      f_obs          = f_obs.data(),
      miller_indices = f_obs.indices(),
      unit_cell      = f_obs.unit_cell())
    a = obj.a
    ####
    #print "Input a :", " ".join(["%7.3f"%i for i in a])
    overall_anisotropic_scale = mmtbx.f_model.ext.k_anisotropic(
      f_obs.indices(), a, f_obs.unit_cell())
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data()*overall_anisotropic_scale)
    f_obs = abs(fc)
    f_obs = f_obs.customized_copy(data = f_obs.data() * overall_anisotropic_scale)
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    t0 = time.time()
    # XXX try "long double" to see if this decreases the tolerances
    obj = bulk_solvent.aniso_u_scaler(
      f_model        = fc.data(),
      f_obs          = f_obs.data(),
      miller_indices = f_obs.indices(),
      unit_cell      = f_obs.unit_cell())
    time_aniso_u_scaler += (time.time()-t0)
    overall_anisotropic_scale = mmtbx.f_model.ext.k_anisotropic(
      f_obs.indices(), obj.a, f_obs.unit_cell())
    assert approx_equal(bulk_solvent.r_factor(f_obs.data(),
      fc.data()*overall_anisotropic_scale), 0.0, 1.e-2) # XXX seems to be low
    #print "Output a:", " ".join(["%7.3f"%i for i in obj.a])
    assert approx_equal(a, obj.a, 1.e-3) # XXX can it be smaller?
  print "Time (aniso_u_scaler only): %6.4f"%time_aniso_u_scaler

def run_02():
  time_aniso_u_scaler = 0
  for symbol in sgtbx.bravais_types.acentric + sgtbx.bravais_types.centric:
    #print symbol, "-"*50
    space_group_info = sgtbx.space_group_info(symbol = symbol)
    xrs = random_structure.xray_structure(
      space_group_info  = space_group_info,
      elements          = ["N"]*100,
      volume_per_atom   = 50.0,
      random_u_iso      = True)
    xrs.scattering_type_registry(table = "wk1995")
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
    reg = xrs.scattering_type_registry(table="wk1995", d_min=1/12)
    f_000 = reg.sum_of_scattering_factors_at_diffraction_angle_0()
    F = xrs.structure_factors(d_min = 2.0).f_calc()
    i = F.indices()
    i.append([0,0,0])
    d = F.data()
    d.append(f_000)
    F = F.customized_copy(indices = i, data = d)

    u_star = adptbx.u_cart_as_u_star(
      F.unit_cell(), adptbx.b_as_u(b_cart_start))
    fbc = mmtbx.f_model.ext.k_anisotropic(F.indices(), u_star)
    fc = F.structure_factors_from_scatterers(xray_structure=xrs).f_calc()
    f_obs = F.customized_copy(data = flex.abs(fc.data()*fbc))
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    obj = bulk_solvent.aniso_u_scaler(
      f_model        = fc.data(),
      f_obs          = f_obs.data(),
      miller_indices = f_obs.indices(),
      unit_cell      = f_obs.unit_cell())
    a = obj.a
    ####
    #print "Input a :", " ".join(["%7.3f"%i for i in a])
    overall_anisotropic_scale = mmtbx.f_model.ext.k_anisotropic(
      f_obs.indices(), a, f_obs.unit_cell())
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data()*overall_anisotropic_scale)
    f_obs = abs(fc)
    f_obs = f_obs.customized_copy(data = f_obs.data() * overall_anisotropic_scale)
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    #print bulk_solvent.r_factor(f_obs.data(), fmodel.f_model().data())
    t0 = time.time()
    # XXX try "long double" to see if this decreases the tolerances
    obj = bulk_solvent.aniso_u_scaler(
      f_model        = fc.data(),
      f_obs          = f_obs.data(),
      miller_indices = f_obs.indices(),
      unit_cell      = f_obs.unit_cell())
    time_aniso_u_scaler += (time.time()-t0)
    overall_anisotropic_scale = mmtbx.f_model.ext.k_anisotropic(
      f_obs.indices(), obj.a, f_obs.unit_cell())
    assert approx_equal(bulk_solvent.r_factor(f_obs.data(),
      fc.data()*overall_anisotropic_scale), 0.0, 1.e-2) # XXX seems to be low
    #print "Output a:", " ".join(["%7.3f"%i for i in obj.a])
    assert approx_equal(a, obj.a, 1.e-4) # XXX can it be smaller?
    assert overall_anisotropic_scale[len(overall_anisotropic_scale)-1]==1
  print "Time (aniso_u_scaler only): %6.4f"%time_aniso_u_scaler

if (__name__ == "__main__"):
  t0 = time.time()
  run_00()
  run_01()
  run_02() # same as run_01 but with f000 added
  print "Time: %6.4f"%(time.time()-t0)
  print "OK"
