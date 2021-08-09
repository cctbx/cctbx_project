from __future__ import absolute_import, division, print_function
from iotbx import pdb
from cctbx.array_family import flex
import os
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
import random
import mmtbx.f_model
from cctbx import sgtbx
import cctbx.sgtbx.bravais_types
from cctbx.development import random_structure
import libtbx.load_env
from cctbx import adptbx
from six.moves import zip
from six.moves import range

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def get_xray_structure_from_file():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2ERL_noH.pdb", test=os.path.isfile)
  xray_structure = pdb.input(file_name=pdb_file).xray_structure_simple()
  return xray_structure

def get_xray_structure_random(space_group_info):
  xray_structure = random_structure.xray_structure(
      space_group_info  = space_group_info,
      elements          = ["N"]*100,
      volume_per_atom   = 50.0,
      random_u_iso      = True)
  return xray_structure

def get_f_obs_freer(d_min, k_sol, b_sol, b_cart, xray_structure,
    radial_shell_width=None):
  f_dummy = abs(xray_structure.structure_factors(d_min = d_min,
    anomalous_flag = False).f_calc())
  r_free_flags = f_dummy.generate_r_free_flags(fraction = 0.1,
                                               max_free = 99999999)
  mask_params = mmtbx.masks.mask_master_params.extract()
  if( radial_shell_width is not None ):
    mask_params.radial_shell_width = radial_shell_width
  if( type(k_sol) is list):
    mask_params.n_radial_shells = len(k_sol)
  fmodel = mmtbx.f_model.manager(
    mask_params    = mask_params,
    r_free_flags   = r_free_flags,
    f_obs          = f_dummy,
    xray_structure = xray_structure)
  fmodel.update_xray_structure(
    xray_structure = xray_structure,
    update_f_calc  = True,
    update_f_mask  = True)
  fmodel_kbu = fmodel.fmodel_kbu()
  fmodel_kbu.update(
    k_sols = k_sol,
    b_sols = b_sol,
    b_cart = b_cart)
  f_obs = abs(fmodel_kbu.f_model)
  return f_obs, r_free_flags

def get_sf(k_sol, b_sol, b_cart, xrs, miller_set=None, d_min=None, twin_law=None,
           sfg_params=None):
  random.seed(0)
  flex.set_random_seed(0)
  if(miller_set is None):
    assert d_min is not None
    f_dummy = abs(xrs.structure_factors(d_min = d_min,
      anomalous_flag = False).f_calc())
  else:
    f_dummy = miller_set
    assert d_min is None
  r_free_flags = f_dummy.generate_r_free_flags(fraction = 0.1)
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_dummy,
    sf_and_grads_accuracy_params = sfg_params,
    xray_structure = xrs,
    twin_law       = twin_law)
  ss = 1./flex.pow2(r_free_flags.d_spacings().data()) / 4.
  k_mask = mmtbx.f_model.ext.k_mask(ss, k_sol, b_sol)
  u_star = adptbx.u_cart_as_u_star(xrs.unit_cell(), adptbx.b_as_u(b_cart))
  k_anisotropic = mmtbx.f_model.ext.k_anisotropic(r_free_flags.indices(),u_star)
  fmodel.update_xray_structure(
    xray_structure = xrs,
    update_f_calc  = True,
    update_f_mask  = True)
  fmodel.update_core(
    k_mask        = k_mask,
    k_anisotropic = k_anisotropic)
  f_obs = abs(fmodel.f_model())
  return f_obs, r_free_flags

def exercise_00(d_min = 3.0, k_sol = 0.35, b_sol = 50.0,
                b_cart = [2.78, -1.54, 8.55, 0, 0, 0],):
  """
  Recover k_sol, b_sol and b_cart up to eps=1.e-6. Somewhat duplication with
  existing tests.
  """
  xrs = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(symbol="P212121"),
    elements               =(("O","N","C")*100),
    volume_per_atom        = 100,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  f_obs, r_free_flags = get_sf(
    d_min  = d_min,
    k_sol  = k_sol,
    b_sol  = b_sol,
    b_cart = b_cart,
    xrs    = xrs)
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    xray_structure = xrs)
  r_work_start = fmodel.r_work()
  assert r_work_start > 0.3
  params = bss.master_params.extract()
  params.number_of_macro_cycles=4
  r = fmodel.update_all_scales(params=params, fast=False, remove_outliers=False)
  r_work_final = fmodel.r_work()
  assert approx_equal(r_work_final, 0,   1.e-5)
  assert approx_equal(r.k_sol[0], k_sol, 1.e-5)
  assert approx_equal(r.b_sol[0], b_sol, 1.e-1)
  assert approx_equal(r.b_cart, b_cart,  1.e-1)

def exercise_01_general(d_mins = [1.6,],
             solvkb = [(0,0), (0.39,58.0), (0.1,6.),(0.54,87.)],
             b_carts = [(4., 10., -14., 0, 5., 0.),
                        (0., 0., 0., 0., 0., 0.)],
             ):
  xray_structure = get_xray_structure_from_file()
  params = bss.master_params.extract()
  params.number_of_macro_cycles=3
  for fast in [True, False]:
    for d_min in d_mins:
      for kb in solvkb:
        for b_cart in b_carts:
          f_obs, r_free_flags = \
            get_f_obs_freer(d_min  = d_min,
                            k_sol  = [kb[0]],
                            b_sol  = [kb[1]],
                            b_cart = b_cart,
                            xray_structure = xray_structure)
          #
          bin_selections = []
          f_obs.setup_binner(reflections_per_bin=50)
          for i_bin in f_obs.binner().range_used():
            sel = f_obs.binner().selection(i_bin)
            bin_selections.append(sel)
          #
          fmodel = mmtbx.f_model.manager(
            r_free_flags   = r_free_flags,
            f_obs          = f_obs,
            xray_structure = xray_structure,
            bin_selections = bin_selections)
          fmodel.update_all_scales(fast=fast, params=params,
            remove_outliers=False)
          result = bss.bulk_solvent_and_scales(
            fmodel_kbu = fmodel.fmodel_kbu(), params = params)
          if(not fast):
            assert approx_equal(fmodel.r_work(), result.fmodel_kbu.r_factor())
          else:
            assert fmodel.r_work() < 0.036, fmodel.r_work()
          assert approx_equal(result.fmodel_kbu.r_factor(), 0.0, eps = 1.e-4)
          assert approx_equal(result.k_sols()[0],      kb[0],    eps = 1.e-4)
          assert approx_equal(result.b_sols()[0],      kb[1],    eps = 1.e-4)
          assert approx_equal(result.b_cart(),        b_cart,    eps = 1.e-2)

def exercise_02_b_cart_sym_constr(d_min = 2.0, tolerance = 1.e-6):
  for symbol in sgtbx.bravais_types.acentric + sgtbx.bravais_types.centric:
    space_group_info = sgtbx.space_group_info(symbol = symbol)
    xray_structure = get_xray_structure_random(space_group_info)
    sg = xray_structure.space_group()
    uc = xray_structure.unit_cell()
    u_cart_p1 = adptbx.random_u_cart(u_scale=5, u_min=5)
    u_star_p1 = adptbx.u_cart_as_u_star(uc, u_cart_p1)
    b_cart_1 = adptbx.u_star_as_u_cart(uc, u_star_p1)
    b_cart_2 = adptbx.u_star_as_u_cart(uc, sg.average_u_star(u_star = u_star_p1))
    for b_cart in (b_cart_1, b_cart_2):
      f_obs, r_free_flags = \
        get_f_obs_freer(d_min  = d_min,
                        k_sol  = 0,
                        b_sol  = 0,
                        b_cart = b_cart,
                        xray_structure = xray_structure)
      fmodel = mmtbx.f_model.manager(
        r_free_flags   = r_free_flags,
        f_obs          = f_obs,
        xray_structure = xray_structure)
      flag=True
      params = bss.master_params.extract()
      params.number_of_macro_cycles=3
      params.bulk_solvent = False
      params.anisotropic_scaling = True
      params.k_sol_b_sol_grid_search = False
      params.minimization_k_sol_b_sol = False
      params.minimization_b_cart = True
      params.symmetry_constraints_on_b_cart = flag
      params.max_iterations = 50
      params.min_iterations = 50
      result = bss.bulk_solvent_and_scales(
        fmodel_kbu = fmodel.fmodel_kbu(), params = params)
      if(flag == True and approx_equal(b_cart, b_cart_2, out=None)):
        assert approx_equal(result.b_cart(), b_cart, tolerance)
      if(flag == True and approx_equal(b_cart, b_cart_1, out=None)):
        for u2, ufm in zip(b_cart_2, result.b_cart()):
          if(abs(u2) < 1.e-6): assert approx_equal(ufm, 0.0, tolerance)

def exercise_03_do_nothing(d_min = 2.0):
  xray_structure = get_xray_structure_from_file()
  k_sol = 0.98
  b_sol = 127.0
  b_cart = [1,2,3,0,4,0]
  f_obs, r_free_flags = get_f_obs_freer(
    d_min  = d_min,
    k_sol  = k_sol,
    b_sol  = b_sol,
    b_cart = b_cart,
    xray_structure = xray_structure)
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    xray_structure = xray_structure)
  r_work_start = fmodel.r_work()*100.
  params = bss.master_params.extract()
  params.bulk_solvent = False
  params.anisotropic_scaling = False
  fmodel.update_all_scales(params = params, fast=False, remove_outliers=False)
  result = bss.bulk_solvent_and_scales(
    fmodel_kbu = fmodel.fmodel_kbu(), params  = params)
  r_work1 = fmodel.r_work()*100.
  assert r_work_start > 0.0
  assert approx_equal(r_work1, r_work_start, eps = 1.e-6)
  assert approx_equal(result.k_sols()[0], 0, eps = 1.e-6)
  assert approx_equal(result.b_sols()[0], 0, eps = 1.e-6)
  assert approx_equal(result.b_cart(), [0,0,0,0,0,0], eps = 1.e-6)

def exercise_04_fix_k_sol_b_sol_b_cart(d_min = 2.0):
  xray_structure = get_xray_structure_from_file()
  k_sol = 0.98
  b_sol = 127.0
  b_cart = [1,2,3,0,4,0]
  f_obs, r_free_flags = get_f_obs_freer(
    d_min  = d_min,
    k_sol  = k_sol,
    b_sol  = b_sol,
    b_cart = b_cart,
    xray_structure = xray_structure)
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    xray_structure = xray_structure)
  r_work_start = fmodel.r_work()*100.
  params = bss.master_params.extract()
  params.k_sol_b_sol_grid_search = False
  params.minimization_k_sol_b_sol = False
  params.minimization_b_cart = False
  params.fix_k_sol = k_sol
  params.fix_b_sol = b_sol
  params.fix_b_cart.b11 = b_cart[0]
  params.fix_b_cart.b22 = b_cart[1]
  params.fix_b_cart.b33 = b_cart[2]
  params.fix_b_cart.b12 = b_cart[3]
  params.fix_b_cart.b13 = b_cart[4]
  params.fix_b_cart.b23 = b_cart[5]
  result = fmodel.update_all_scales(params = params, fast = False,
    remove_outliers=False)
  r_work = fmodel.r_work()*100.
  assert r_work_start > 0.0
  assert approx_equal(r_work,          0.0, eps = 1.e-6)
  assert approx_equal(result.k_sol[0], k_sol, eps = 1.e-6)
  assert approx_equal(result.b_sol[0], b_sol, eps = 1.e-6)
  assert approx_equal(result.b_cart, b_cart, eps = 1.e-6)

def exercise_05_k_sol_b_sol_only(d_min = 2.0):
  xray_structure = get_xray_structure_from_file()
  k_sol = 0.33
  b_sol = 34.0
  b_cart = [1,2,3,0,4,0]
  f_obs, r_free_flags = get_f_obs_freer(
    d_min  = d_min,
    k_sol  = k_sol,
    b_sol  = b_sol,
    b_cart = b_cart,
    xray_structure = xray_structure)
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    xray_structure = xray_structure)
  params = bss.master_params.extract()
  params.anisotropic_scaling = False
  params.number_of_macro_cycles=5
  u_star = adptbx.u_cart_as_u_star(
    fmodel.f_obs().unit_cell(),adptbx.b_as_u(b_cart))
  fmodel_kbu = mmtbx.f_model.manager_kbu(
    f_obs   = fmodel.f_obs(),
    f_calc  = fmodel.f_calc(),
    f_masks = fmodel.arrays.core.f_masks,
    f_part1 = fmodel.arrays.core.f_part1,
    f_part2 = fmodel.arrays.core.f_part2,
    ss      = fmodel.ss,
    u_star  = u_star)
  r_work_start = fmodel_kbu.r_factor()
  result = bss.bulk_solvent_and_scales(
    fmodel_kbu = fmodel_kbu, params = params)
  r_work = result.fmodel_kbu.r_factor()*100.
  assert r_work_start > 0.05
  #
  assert approx_equal(r_work,              0.0,   eps = 1.e-4)
  assert approx_equal(result.k_sols()[0], k_sol,  eps = 1.e-4)
  assert approx_equal(result.b_sols()[0], b_sol,  eps = 1.e-4)
  assert approx_equal(result.b_cart(),    b_cart, eps = 1.e-4)

def exercise_06_b_cart_only(d_min = 2.0):
  xray_structure = get_xray_structure_from_file()
  k_sol = 0.33
  b_sol = 34.0
  b_cart = [1,2,3,0,4,0]
  f_obs, r_free_flags = get_f_obs_freer(
    d_min  = d_min,
    k_sol  = k_sol,
    b_sol  = b_sol,
    b_cart = b_cart,
    xray_structure = xray_structure)
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    xray_structure = xray_structure)
  fmodel_kbu = mmtbx.f_model.manager_kbu(
    f_obs   = fmodel.f_obs(),
    f_calc  = fmodel.f_calc(),
    f_masks = fmodel.arrays.core.f_masks,
    f_part1 = fmodel.arrays.core.f_part1,
    f_part2 = fmodel.arrays.core.f_part2,
    ss      = fmodel.ss,
    k_sols  = [k_sol],
    b_sols  = [b_sol])
  r_work_start = fmodel_kbu.r_factor()*100.
  params = bss.master_params.extract()
  params.bulk_solvent = False
  result = bss.bulk_solvent_and_scales(
    fmodel_kbu = fmodel_kbu, params  = params)
  r_work = result.fmodel_kbu.r_factor()*100.
  assert r_work_start > 0.0
  assert approx_equal(r_work,               0.0,  eps = 1.e-6)
  assert approx_equal(result.k_sols()[0], k_sol,  eps = 1.e-6)
  assert approx_equal(result.b_sols()[0], b_sol,  eps = 1.e-6)
  assert approx_equal(result.b_cart(),    b_cart, eps = 1.e-6)

def exercise_radial_shells(k_sol=0.33,d_min=1.5,grid_search=False,shell_width=0.6):
  xray_structure = get_xray_structure_from_file()
  b_sol = 34.0
  if( type(k_sol) is list ):
    b_sol = [b_sol,]*len(k_sol)
  b_cart = [1,2,3,0,4,0]
  f_obs, r_free_flags = get_f_obs_freer(
    d_min  = d_min,
    k_sol  = k_sol,
    b_sol  = b_sol,
    b_cart = b_cart,
    xray_structure = xray_structure,
    radial_shell_width=shell_width)
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.radial_shell_width = shell_width
  if( type(k_sol) is list ):
    mask_params.n_radial_shells = len(k_sol)
  else:
    mask_params.n_radial_shells = 2
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    xray_structure = xray_structure,
    mask_params    = mask_params)
  u_star = adptbx.u_cart_as_u_star(
    fmodel.f_obs().unit_cell(),adptbx.b_as_u(b_cart))
  fmodel_kbu = fmodel.fmodel_kbu()
  fmodel_kbu.update(u_star = u_star)
  r_work_start = fmodel_kbu.r_factor()*100.
  msk = fmodel.mask_manager
  print('Solvent content: ', msk.solvent_content_via_mask)
  print('Layer volume fractions: ', msk.layer_volume_fractions)
  if( type(k_sol) is list):
    for i in range(len(k_sol)):
      if( msk.layer_volume_fractions[i] == 0. ):
        k_sol[i] = 0.
  params = bss.master_params.extract()
  params.anisotropic_scaling = False
  params.k_sol_b_sol_grid_search = grid_search
  params.number_of_macro_cycles = 10
  params.k_sol_max = 1.2
  result = bss.bulk_solvent_and_scales(
    fmodel_kbu = fmodel_kbu, params = params)
  r_work = result.fmodel_kbu.r_factor()
  print('R-work: ', r_work)
  print('Solvent radius: ', fmodel.mask_params.solvent_radius)
  assert r_work_start > 0.0
  assert approx_equal(r_work, 0.0, eps = 1.e-3)
  if( type(k_sol) is list ):
    ksols = list(result.fmodel_kbu.k_sols())
    # XXX if layer_volume_fractions=0, then ksol is more or less undefined ?
    # XXX should it be done in bulk_solvent_and_scaling.py ?
    for i in range(len(ksols)):
      if(msk.layer_volume_fractions[i] == 0.):
        ksols[i] = 0.
    assert len(k_sol) == len(ksols)
    for ik in range(len(k_sol)):
      assert approx_equal(ksols[ik], k_sol[ik], eps=0.005),[ksols[ik],k_sol[ik]]
  else:
    for ksol in result.fmodel_kbu.k_sols():
      assert approx_equal(ksol,   k_sol, eps = 1.e-3)
  n=len(result.b_sols())
  if(n>1 and type(b_sol) is float): b_sol = [b_sol,]*n
  assert approx_equal(result.b_sols(), b_sol,  eps = 1.)
  assert approx_equal(result.b_cart(), b_cart, eps = 1.e-6)

def run():
  exercise_00()
  exercise_01_general()
  exercise_02_b_cart_sym_constr()
  exercise_03_do_nothing()
  exercise_04_fix_k_sol_b_sol_b_cart()
  exercise_05_k_sol_b_sol_only()
  exercise_06_b_cart_only()
  exercise_radial_shells()
  exercise_radial_shells(k_sol=[0.33,0.5, 0.9])
  exercise_radial_shells(k_sol=[0.3,0.4,0.5],grid_search=True,shell_width=0.3)
  print("OK: ",format_cpu_times())

if (__name__ == "__main__"):
  run()
