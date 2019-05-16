from __future__ import absolute_import, division, print_function
from iotbx import pdb
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import mmtbx.f_model
from cctbx import adptbx
import iotbx.pdb
from mmtbx import f_model
from cctbx import maptbx
from mmtbx.bulk_solvent import kbu_refinery
import time
from mmtbx import bulk_solvent
from six.moves import range

pdb_str = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P1
ATOM      1  O   HOH A   1       2.500   2.500   2.500  1.00 10.00           O
ATOM      1  O   HOH A   1       5.500   2.500   2.500  1.00 10.00           O
"""

def get_mask_data(xrs, d_min):
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell          = xrs.unit_cell(),
    d_min              = d_min,
    resolution_factor  = 1./4,
    symmetry_flags     = maptbx.use_space_group_symmetry,
    space_group_info   = xrs.space_group_info())
  mp = mmtbx.masks.mask_master_params.extract()
  return mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xrs,
    p1                       = True,
    solvent_radius           = mp.solvent_radius,
    shrink_truncation_radius = mp.shrink_truncation_radius,
    for_structure_factors    = True,
    n_real                   = crystal_gridding.n_real()).mask_data

def fd_k_sols(kbu, f_obs, eps=1.e-3):
  p = flex.double([0.1, 0.2])
  kbu.update(k_sols = p)
  tg = bulk_solvent.ls_kbp_sol_u_star(
    f_model     = kbu.data,
    f_obs       = f_obs.data(),
    scale       = 1.0,
    kb_sol_grad = True,
    p_sol_grad  = False,
    u_star_grad = True,
    kb_sol_curv = True,
    p_sol_curv  = False)
  g_k_sols_anal = list(tg.grad_k_sols())
  c_k_sols_anal = list(tg.curv_k_sols())
  print("C anal k:", c_k_sols_anal)
  #
  g_k_sols_fd = []
  c_k_sols_fd = []
  t0 = tg.target()
  for shift in [flex.double([eps,0]),
                flex.double([0,eps])]:
    # plus shift
    kbu.update(k_sols = p+shift)
    tg = bulk_solvent.ls_kbp_sol_u_star(
      f_model     = kbu.data,
      f_obs       = f_obs.data(),
      scale       = 1.0,
      kb_sol_grad = True,
      p_sol_grad  = False,
      u_star_grad = True,
      kb_sol_curv = True,
      p_sol_curv  = False)
    t1 = tg.target()
    # minus shift
    kbu.update(k_sols = p-shift)
    tg = bulk_solvent.ls_kbp_sol_u_star(
      f_model     = kbu.data,
      f_obs       = f_obs.data(),
      scale       = 1.0,
      kb_sol_grad = True,
      p_sol_grad  = False,
      u_star_grad = True,
      kb_sol_curv = True,
      p_sol_curv  = False)
    t2 = tg.target()
    g_k_sols_fd.append( (t1-t2)/(2*eps) )
    c_k_sols_fd.append( (t1+t2-2*t0)/eps**2 )
  #
  print("C df:",c_k_sols_fd)
  assert approx_equal(g_k_sols_anal, g_k_sols_fd, 1.e-4)

def fd_b_sols(kbu, f_obs):
  p = flex.double([10, 20])
  kbu.update(b_sols = p)
  tg = bulk_solvent.ls_kbp_sol_u_star(
    f_model     = kbu.data,
    f_obs       = f_obs.data(),
    scale       = 1.0,
    kb_sol_grad = True,
    p_sol_grad  = False,
    u_star_grad = True,
    kb_sol_curv = True,
    p_sol_curv  = False)
  g_b_sols_anal = list(tg.grad_b_sols())
  c_b_sols_anal = list(tg.curv_b_sols())
  print("C anal b:", c_b_sols_anal)
  #
  g_b_sols_fd = []
  c_b_sols_fd = []
  t0 = tg.target()
  for shift in [flex.double([1.e-3,0]),
                flex.double([0,1.e-3])]:
    # plus shift
    kbu.update(b_sols = p+shift)
    tg = bulk_solvent.ls_kbp_sol_u_star(
      f_model     = kbu.data,
      f_obs       = f_obs.data(),
      scale       = 1.0,
      kb_sol_grad = True,
      p_sol_grad  = False,
      u_star_grad = True,
      kb_sol_curv = True,
      p_sol_curv  = False)
    t1 = tg.target()
    # minus shift
    kbu.update(b_sols = p-shift)
    tg = bulk_solvent.ls_kbp_sol_u_star(
      f_model     = kbu.data,
      f_obs       = f_obs.data(),
      scale       = 1.0,
      kb_sol_grad = True,
      p_sol_grad  = False,
      u_star_grad = True,
      kb_sol_curv = True,
      p_sol_curv  = False)
    t2 = tg.target()
    g_b_sols_fd.append( (t1-t2)/(2*1.e-3) )
    c_b_sols_fd.append( (t1+t2-2*t0)/1.e-3**2 *2)
    #
  print("C df:",c_b_sols_fd)
  assert approx_equal(g_b_sols_anal, g_b_sols_fd)

def run(d_min  = 2.0,
        b_cart = [1,2,3,0,4,0],
        k_sols = [0.3, 0.5],
        b_sols = [50., 20.]):
  xrs=iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  #
  f_calc = xrs.structure_factors(d_min=d_min).f_calc()
  #
  mask_data = get_mask_data(xrs=xrs, d_min=d_min)
  n = mask_data.all()
  mask_data1 = flex.double(flex.grid(n), 0)
  mask_data2 = flex.double(flex.grid(n), 0)
  I,J,K = range(n[0]), range(n[1]), range(n[2])
  for i in I:
    for j in J:
      for k in K:
        if(i < n[0]//2 and j < n[1]//2 and k < n[2]//2):
          mask_data1[i,j,k]=mask_data[i,j,k]
        else:
          mask_data2[i,j,k]=mask_data[i,j,k]
  f_mask1 = f_calc.structure_factors_from_map(map=mask_data1,
    use_scale = True, anomalous_flag = False, use_sg = False)
  f_mask2 = f_calc.structure_factors_from_map(map=mask_data2,
    use_scale = True, anomalous_flag = False, use_sg = False)
  #
  u_star = adptbx.u_cart_as_u_star(xrs.unit_cell(), adptbx.b_as_u(b_cart))
  k_total = mmtbx.f_model.ext.k_anisotropic(f_calc.indices(), u_star)
  #
  d_spacings = f_calc.d_spacings().data()
  ss = 1./flex.pow2(d_spacings) / 4.
  k_mask1 = mmtbx.f_model.ext.k_mask(ss, k_sols[0], b_sols[0])
  k_mask2 = mmtbx.f_model.ext.k_mask(ss, k_sols[1], b_sols[1])
  #
  f_obs_data = flex.abs( k_total*(f_calc.data()+k_mask1*f_mask1.data() +
                                                k_mask2*f_mask2.data()) )
  f_obs = f_calc.array(data = f_obs_data)
  #
  kbu = f_model.manager_kbu(
    f_obs   = f_obs,
    f_calc  = f_calc,
    f_masks = [f_mask1, f_mask2],
    ss      = ss,
    k_sols  = k_sols,
    b_sols  = b_sols)
  #
  fd_k_sols(kbu=kbu, f_obs=f_obs)
  fd_b_sols(kbu=kbu, f_obs=f_obs)
  #
  u_star_ini = adptbx.u_cart_as_u_star(xrs.unit_cell(),
    adptbx.b_as_u([0.1,0.2,0.3,0,0.4,0]))
  TGCO = kbu_refinery.tgc(
    f_obs   = f_obs,
    f_calc  = f_calc,
    f_masks = [f_mask1, f_mask2],
    ss      = ss,
    k_sols  = [0.1,0.1],
    b_sols  = [10,10],
    u_star  = u_star_ini)
  TGCO.minimize_kbu(n_cycles=100)
  TGCO.show_kbu()
  assert approx_equal(TGCO.kbu.k_sols(), k_sols, 1.e-3)
  assert approx_equal(TGCO.kbu.b_sols(), b_sols, 1.e-3)
  assert approx_equal(TGCO.kbu.b_cart(), b_cart, 1.e-4)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
