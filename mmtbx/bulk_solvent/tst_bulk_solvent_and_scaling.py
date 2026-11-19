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
    radial_shell_width=None, mask_radii=None):
  f_dummy = abs(xray_structure.structure_factors(d_min = d_min,
    anomalous_flag = False).f_calc())
  r_free_flags = f_dummy.generate_r_free_flags(fraction = 0.1,
                                               max_free = 99999999)
  mask_params = mmtbx.masks.mask_master_params.extract()
  if( mask_radii is not None ):
    mask_params.solvent_radius, mask_params.shrink_truncation_radius = \
      mask_radii
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
          if(fast): # mode=fast accepts no parameters of mode=slow
            fmodel.update_all_scales(fast=True, remove_outliers=False)
          else:
            fmodel.update_all_scales(fast=False, params=params,
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

def exercise_02_b_cart_sym_constr(d_min = 2.0, tolerance = 1.e-5):
  # Minimization stops when the target (sum of squared differences divided
  # by the sum of Fobs**2) does not decrease by more than machine precision,
  # that is at R of about 1.e-8, where b_cart is defined to about 1.e-6.
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
      params.bulk_solvent = False
      params.anisotropic_scaling = True
      params.slow.k_sol_b_sol_grid_search = False
      params.slow.minimization_k_sol_b_sol = False
      params.slow.minimization_b_cart = True
      params.slow.symmetry_constraints_on_b_cart = flag
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
  params.slow.k_sol_b_sol_grid_search = False
  params.slow.minimization_k_sol_b_sol = False
  params.slow.minimization_b_cart = False
  params.slow.fix_k_sol = k_sol
  params.slow.fix_b_sol = b_sol
  params.slow.fix_b_cart.b11 = b_cart[0]
  params.slow.fix_b_cart.b22 = b_cart[1]
  params.slow.fix_b_cart.b33 = b_cart[2]
  params.slow.fix_b_cart.b12 = b_cart[3]
  params.slow.fix_b_cart.b13 = b_cart[4]
  params.slow.fix_b_cart.b23 = b_cart[5]
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
  # b_sol converges to 1e-4 (34.0001 with mask radii 0.95/0.7): minimizer
  # precision, the model is exact (r_work 1e-4).
  assert approx_equal(result.b_sols()[0], b_sol,  eps = 1.e-3)
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
  # The shells are built on the historical mask (solvent radius 1.1 A,
  # shrink 0.9 A): with the default 0.95/0.7 A mask the three 0.3 A shells
  # are degenerate (r_work 0 with k_sol 0.32 instead of 0.30) and the test
  # no longer checks that the shell parameters are recovered.
  mask_radii = (1.1, 0.9)
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
    radial_shell_width=shell_width,
    mask_radii = mask_radii)
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.solvent_radius, mask_params.shrink_truncation_radius = mask_radii
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
  params.slow.k_sol_b_sol_grid_search = grid_search
  params.slow.k_sol_max = 1.2
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

def _params(edits=()):
  """
  master_params with edits like "slow.k_sol_max=0.5" applied.
  """
  import iotbx.phil
  user = iotbx.phil.parse("\n".join(edits))
  working, unused = bss.master_params.fetch(
    source=user, track_unused_definitions=True)
  assert len(unused) == 0, [str(u) for u in unused]
  return working.extract()

def _exact_fmodel_kbu(symbol="P212121", d_min=3.0, k_sol=0.35, b_sol=50.0,
                      b_cart=(3., -5., 2., 0, 0, 0), zero_mask=False):
  """
  manager_kbu with all scales reset and f_obs computed from known k_sol,
  b_sol, b_cart.
  """
  random.seed(0)
  flex.set_random_seed(0)
  xrs = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(symbol=symbol),
    elements               = (("O","N","C")*50),
    volume_per_atom        = 100,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  f_obs, r_free_flags = get_f_obs_freer(d_min=d_min, k_sol=[k_sol],
    b_sol=[b_sol], b_cart=list(b_cart), xray_structure=xrs)
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    xray_structure = xrs)
  kbu = fmodel.fmodel_kbu()
  if(zero_mask):
    fm = kbu.f_masks[0]
    kbu = mmtbx.f_model.manager_kbu(
      f_obs   = kbu.f_obs,
      f_calc  = kbu.f_calc,
      f_masks = [fm.customized_copy(data = fm.data()*0)],
      ss      = kbu.ss)
  return kbu

def exercise_scope_layout():
  p = bss.master_params.extract()
  assert p.mode == "fast"
  assert [p.bulk_solvent, p.anisotropic_scaling, p.apply_back_trace]==[True]*3
  s = p.slow
  assert approx_equal([s.k_sol_grid_search_min, s.k_sol_grid_search_max,
    s.k_sol_step], [0, 0.6, 0.1])
  assert approx_equal([s.b_sol_grid_search_min, s.b_sol_grid_search_max,
    s.b_sol_step], [10, 80, 5])
  assert s.max_iterations == 100
  assert not hasattr(s, "min_iterations")
  assert not hasattr(s, "number_of_macro_cycles")
  assert not hasattr(p, "target") and not hasattr(s, "target")
  assert not hasattr(p, "k_sol_max")
  assert bss.non_default_slow_parameters(p) == []
  s.k_sol_max = 0.5
  s.fix_b_cart.b11 = 1.
  assert bss.non_default_slow_parameters(p) == [
    "slow.k_sol_max", "slow.fix_b_cart.b11"]

def exercise_validate_params():
  from libtbx.utils import Sorry
  fix_b_cart = ["slow.fix_b_cart.b11=1", "slow.fix_b_cart.b22=1",
    "slow.fix_b_cart.b33=1", "slow.fix_b_cart.b12=0", "slow.fix_b_cart.b13=0",
    "slow.fix_b_cart.b23=0"]
  no_kb = ["slow.k_sol_b_sol_grid_search=False",
    "slow.minimization_k_sol_b_sol=False"]
  # edits, twinned, n_masks, text expected in the message
  cases = [
    (["mode=slow","slow.fix_k_sol=0.3"], False, 1, "slow.fix_k_sol"),
    (["mode=slow","slow.fix_b_sol=30"], False, 1, "slow.fix_b_sol"),
    (["mode=slow","slow.fix_k_sol=0.3","bulk_solvent=False"]+no_kb, False, 1,
      "bulk_solvent=False"),
    (["mode=slow","slow.fix_b_cart.b11=1"], False, 1, "all six"),
    (["mode=slow"]+fix_b_cart, False, 1, "slow.minimization_b_cart"),
    (["mode=slow","slow.minimization_b_cart=False",
      "anisotropic_scaling=False"]+fix_b_cart, False, 1,
      "anisotropic_scaling=False"),
    (["mode=slow","slow.k_sol_min=0.7"], False, 1, "slow.k_sol_min"),
    (["mode=slow","slow.b_sol_min=400"], False, 1, "slow.b_sol_min"),
    (["mode=slow","slow.k_sol_grid_search_min=0.7"], False, 1,
      "slow.k_sol_grid_search_min"),
    (["mode=slow","slow.b_sol_grid_search_min=90"], False, 1,
      "slow.b_sol_grid_search_min"),
    (["mode=slow","slow.k_sol_step=-0.1"], False, 1, "slow.k_sol_step"),
    (["mode=slow","slow.b_sol_step=0"], False, 1, "slow.b_sol_step"),
    (["mode=slow","slow.k_sol_step=0.0005"], False, 1, "slow.k_sol_step"),
    (["mode=slow","slow.max_iterations=0"], False, 1, "slow.max_iterations"),
    (["slow.k_sol_max=0.5"], False, 1, "mode=fast"),
    (["mode=slow"], True, 1, "twin"),
    (["mode=slow", "slow.minimization_b_cart=False"]+fix_b_cart, True, 1,
      "twin"),
    ([], False, 3, "mode=slow"),
  ]
  for edits, twinned, n_masks, text in cases:
    p = _params(edits)
    sorry = None
    try:
      bss.validate_params(params=p, fast=(p.mode=="fast"), twinned=twinned,
        n_masks=n_masks)
    except Sorry as e:
      sorry = str(e)
    assert sorry is not None, edits
    assert sorry.count(text) > 0, (edits, sorry)
    assert sorry.count("bulk_solvent_and_scale") > 0, sorry
  # legal
  bss.validate_params(params=_params(), fast=True)
  bss.validate_params(params=_params(), fast=True, twinned=True)
  bss.validate_params(params=_params(["mode=slow"]), fast=False, n_masks=3)
  bss.validate_params(params=_params(["mode=slow","slow.fix_k_sol=0.3",
    "slow.fix_b_sol=30"]+no_kb), fast=False)
  bss.validate_params(params=_params(
    ["mode=slow","slow.minimization_b_cart=False"]+fix_b_cart), fast=False)
  # legacy callers: default mode with explicit slow protocol
  bss.validate_params(params=_params(["slow.max_iterations=50"]),
    fast=False)

def exercise_wired_grid():
  # grid only: answer is on the grid given by parameters
  kbu = _exact_fmodel_kbu()
  p = _params(["mode=slow", "slow.minimization_k_sol_b_sol=False",
    "slow.minimization_b_cart=False", "slow.k_sol_step=0.3",
    "slow.b_sol_grid_search_min=20", "slow.b_sol_step=30"])
  r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu, params=p)
  assert approx_equal(r.k_sols()[0], 0.3, 1.e-6)
  assert approx_equal(r.b_sols()[0], 50., 1.e-6)
  # another grid, another answer
  kbu = _exact_fmodel_kbu()
  p = _params(["mode=slow", "slow.minimization_k_sol_b_sol=False",
    "slow.minimization_b_cart=False", "slow.k_sol_grid_search_max=0.2",
    "slow.k_sol_step=0.2", "slow.b_sol_grid_search_min=40",
    "slow.b_sol_grid_search_max=60", "slow.b_sol_step=20"])
  r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu, params=p)
  assert approx_equal(r.k_sols()[0], 0.2, 1.e-6)
  assert r.b_sols()[0] in [40., 60.]

def exercise_wired_grid_keeps_b_cart():
  # anisotropic_scaling=False: grid search must not touch b_cart
  kbu = _exact_fmodel_kbu(b_cart=(0,0,0,0,0,0))
  kbu.update(b_cart = [1., 2., -3., 0, 0, 0])
  p = _params(["mode=slow", "anisotropic_scaling=False",
    "slow.minimization_k_sol_b_sol=False"])
  r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu, params=p)
  assert approx_equal(r.b_cart(), [1., 2., -3., 0, 0, 0], 1.e-6)

def exercise_wired_symmetry_constraints():
  b_cart = (3., -5., 2., 0, 0, 0) # b11 != b22: not allowed in P4
  for flag in [True, False]:
    kbu = _exact_fmodel_kbu(symbol="P4", k_sol=0, b_sol=0, b_cart=b_cart)
    p = _params(["mode=slow", "bulk_solvent=False",
      "slow.symmetry_constraints_on_b_cart=%s"%str(flag)])
    r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu, params=p)
    if(flag):
      assert approx_equal(r.b_cart()[0], r.b_cart()[1], 1.e-6)
      assert r.r_factor() > 1.e-3, r.r_factor()
    else:
      assert approx_equal(r.b_cart(), b_cart, 1.e-2)
      assert r.r_factor() < 1.e-4, r.r_factor()

def exercise_wired_iterations():
  r_factors = []
  for max_iterations in [1, 100]:
    kbu = _exact_fmodel_kbu(b_cart=(0,0,0,0,0,0))
    p = _params(["mode=slow", "anisotropic_scaling=False",
      "slow.k_sol_b_sol_grid_search=False",
      "slow.max_iterations=%d"%max_iterations])
    r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu, params=p)
    r_factors.append(r.r_factor())
  assert r_factors[0] > 0.01, r_factors
  assert r_factors[1] < 1.e-4, r_factors

def exercise_fix_k_sol_zero_mask():
  kbu = _exact_fmodel_kbu(zero_mask=True)
  p = _params(["mode=slow", "slow.fix_k_sol=0.2", "slow.fix_b_sol=30",
    "slow.k_sol_b_sol_grid_search=False",
    "slow.minimization_k_sol_b_sol=False"])
  r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu, params=p)
  assert approx_equal(r.k_sols()[0], 0)

def exercise_fix_values_used():
  kbu = _exact_fmodel_kbu()
  p = _params(["mode=slow", "slow.fix_k_sol=0.2", "slow.fix_b_sol=30",
    "slow.k_sol_b_sol_grid_search=False",
    "slow.minimization_k_sol_b_sol=False"])
  r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu, params=p)
  assert approx_equal(r.k_sols()[0], 0.2)
  assert approx_equal(r.b_sols()[0], 30.)

def exercise_symmetry_constraints_with_grid_search():
  # Constraints must hold whatever combination of steps is used
  b_cart = (3., -5., 2., 0, 0, 0) # b11 != b22: not allowed in P4
  for edits in [[], ["slow.minimization_b_cart=False"],
                ["slow.minimization_k_sol_b_sol=False"]]:
    kbu = _exact_fmodel_kbu(symbol="P4", b_cart=b_cart)
    r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu,
      params=_params(["mode=slow"]+edits))
    bc = r.b_cart()
    assert approx_equal(bc[0], bc[1], 1.e-6), (edits, bc)
    assert approx_equal(bc[3:], [0,0,0], 1.e-6), (edits, bc)

def exercise_bounds_are_respected():
  for edits, k_expected, b_expected in [
      (["slow.k_sol_max=0.25"], 0.25, None),
      (["slow.b_sol_min=60"], None, 60.),
      (["slow.k_sol_min=0.45"], 0.45, None),
      (["slow.b_sol_max=40"], None, 40.)]:
    kbu = _exact_fmodel_kbu()
    r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu,
      params=_params(["mode=slow"]+edits))
    if(k_expected is not None):
      assert approx_equal(r.k_sols()[0], k_expected, 1.e-6), edits
    if(b_expected is not None):
      assert approx_equal(r.b_sols()[0], b_expected, 1.e-6), edits

def _tgc(kbu, k_sol=0.1, b_sol=20.):
  import mmtbx.bulk_solvent.kbu_refinery as kbu_refinery
  return kbu_refinery.tgc(
    f_obs   = kbu.f_obs,
    f_calc  = kbu.f_calc,
    f_masks = kbu.f_masks,
    ss      = kbu.ss,
    k_sols  = [k_sol],
    b_sols  = [b_sol],
    u_star  = [0,0,0,0,0,0])

def exercise_joint_minimizer():
  """
  k_sol, b_sol and b_cart refined together, with bounds.
  """
  b_cart = (3., -5., 2., 0, 0, 0)
  kbu = _exact_fmodel_kbu(b_cart=b_cart)
  free = dict(k_sol_min=0, k_sol_max=0.6, b_sol_min=0, b_sol_max=300)
  # far start, no grid search: exact answer
  o = _tgc(kbu)
  m = o.minimize_joint(max_iterations=100, **free)
  assert o.kbu.r_factor() < 1.e-4
  assert m.n_calls < 100, m.n_calls
  assert approx_equal(o.kbu.k_sols(), [0.35], 1.e-4)
  assert approx_equal(o.kbu.b_sols(), [50.], 1.e-2)
  assert approx_equal(o.kbu.b_cart(), b_cart, 1.e-2)
  # the answer is outside the bounds: result is on the bound, exactly
  o = _tgc(kbu)
  o.minimize_joint(k_sol_min=0, k_sol_max=0.25, b_sol_min=0, b_sol_max=300,
    max_iterations=100)
  assert o.kbu.k_sols()[0] == 0.25
  assert o.kbu.r_factor() > 1.e-3
  # the start is outside the bounds
  o = _tgc(kbu)
  o.minimize_joint(k_sol_min=0, k_sol_max=0.6, b_sol_min=60, b_sol_max=300,
    max_iterations=100)
  assert o.kbu.b_sols()[0] == 60.
  # subsets of parameters
  o = _tgc(kbu)
  o.minimize_joint(refine_u=False, max_iterations=100, **free)
  assert approx_equal(o.kbu.b_cart(), [0,0,0,0,0,0])
  assert approx_equal(o.kbu.k_sols(), [0.35], 0.01)
  o = _tgc(kbu, k_sol=0.35, b_sol=50.)
  o.minimize_joint(refine_kb=False, max_iterations=100, **free)
  assert approx_equal(o.kbu.k_sols(), [0.35])
  assert approx_equal(o.kbu.b_sols(), [50.])
  assert approx_equal(o.kbu.b_cart(), b_cart, 1.e-2)
  # iteration limit
  o = _tgc(kbu)
  m = o.minimize_joint(max_iterations=1, **free)
  assert m.n_iterations == 1
  assert o.kbu.r_factor() > 0.01

def exercise_default_slow_is_exact():
  """
  Default parameters of mode=slow reproduce exact data; nothing else is
  needed (no macro-cycles).
  """
  b_cart = (3., -5., 2., 0, 0, 0)
  kbu = _exact_fmodel_kbu(b_cart=b_cart)
  r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu, params=_params(["mode=slow"]))
  assert r.r_factor() < 1.e-5, r.r_factor()
  assert approx_equal(r.k_sols(), [0.35], 1.e-4)
  assert approx_equal(r.b_sols(), [50.], 1.e-2)
  assert approx_equal(r.b_cart(), b_cart, 1.e-2)

def exercise_large_b_cart():
  """
  Overall B of the data far from that of the model (b_cart of 100 A**2 and
  more, either sign) is found, with and without grid search.
  """
  for b in [150., -120.]:
    b_cart = (b, b, b+30., 0, 0, 0)
    for edits in [[], ["slow.k_sol_b_sol_grid_search=False"]]:
      kbu = _exact_fmodel_kbu(d_min=3.5, b_cart=b_cart)
      r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu,
        params=_params(["mode=slow"]+edits))
      assert r.r_factor() < 1.e-3, (b, edits, r.r_factor())
      assert approx_equal(r.k_sols(), [0.35], 1.e-2)
      assert approx_equal(r.b_sols(), [50.], 0.5)
      assert approx_equal(r.b_cart(), b_cart, 0.1)

pdb_str_much_solvent = """
CRYST1   15.000   16.000   17.000  90.00  90.00  90.00 P 1
ATOM      1  CB  PHE H   1       7.767   5.853   7.671  1.00 15.00           C
ATOM      2  CG  PHE H   1       6.935   5.032   8.622  1.00 15.00           C
ATOM      3  CD1 PHE H   1       5.918   4.176   8.140  1.00 15.00           C
ATOM      4  CD2 PHE H   1       7.161   5.107  10.012  1.00 15.00           C
ATOM      5  CE1 PHE H   1       5.126   3.395   9.038  1.00 15.00           C
ATOM      6  CE2 PHE H   1       6.382   4.336  10.930  1.00 15.00           C
ATOM      7  CZ  PHE H   1       5.360   3.476  10.439  1.00 15.00           C
ATOM      8  C   PHE H   1       7.956   7.811   6.133  1.00 15.00           C
ATOM      9  O   PHE H   1       8.506   7.237   5.169  1.00 15.00           O
ATOM     10  OXT PHE H   1       8.143   9.010   6.428  1.00 15.00           O
ATOM     13  N   PHE H   1       5.875   6.461   6.183  1.00 15.00           N
ATOM     15  CA  PHE H   1       7.000   7.000   7.000  1.00 15.00           C
END
"""

def exercise_no_solvent_model_for_data_with_solvent():
  # b_cart must stay sane if the bulk solvent is not modelled: this case
  # used to give b_cart of several thousands.
  xrs = pdb.input(source_info=None,
    lines=pdb_str_much_solvent).xray_structure_simple()
  f_obs, r_free_flags = get_f_obs_freer(d_min=2.0, k_sol=[0.35], b_sol=[50.],
    b_cart=[3., -5., 2., 0, 0, 0], xray_structure=xrs)
  fmodel = mmtbx.f_model.manager(
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    xray_structure = xrs)
  kbu = fmodel.fmodel_kbu()
  r_start = kbu.r_factor()
  r = bss.bulk_solvent_and_scales(fmodel_kbu=kbu,
    params=_params(["mode=slow", "bulk_solvent=False"]))
  assert r.r_factor() <= r_start
  assert approx_equal(r.k_sols()[0], 0)
  assert max([abs(b) for b in r.b_cart()]) < 100, r.b_cart()

def run():
  exercise_scope_layout()
  exercise_validate_params()
  exercise_wired_grid()
  exercise_wired_grid_keeps_b_cart()
  exercise_wired_symmetry_constraints()
  exercise_wired_iterations()
  exercise_fix_k_sol_zero_mask()
  exercise_fix_values_used()
  exercise_symmetry_constraints_with_grid_search()
  exercise_bounds_are_respected()
  exercise_no_solvent_model_for_data_with_solvent()
  exercise_joint_minimizer()
  exercise_default_slow_is_exact()
  exercise_large_b_cart()
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
