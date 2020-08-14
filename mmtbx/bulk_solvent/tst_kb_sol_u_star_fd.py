from __future__ import absolute_import, division, print_function
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost_adaptbx.boost.python as bp
from six.moves import range
ext = bp.import_ext("cctbx_asymmetric_map_ext")
from cctbx_asymmetric_map_ext import *
from cctbx.array_family import flex
from cctbx import maptbx
import mmtbx.f_model
from libtbx import adopt_init_args
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from libtbx.test_utils import approx_equal
from cctbx import adptbx

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_f_model_ext")

def get_f_masks(xrs, miller_array):
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell          = xrs.unit_cell(),
    d_min              = miller_array.d_min(),
    resolution_factor  = 1./4,
    symmetry_flags     = maptbx.use_space_group_symmetry,
    space_group_info   = xrs.space_group_info())
  mp = mmtbx.masks.mask_master_params.extract()
  mask_data = mmtbx.masks.mask_from_xray_structure(
    xray_structure           = xrs,
    p1                       = True,
    solvent_radius           = mp.solvent_radius,
    shrink_truncation_radius = mp.shrink_truncation_radius,
    for_structure_factors    = True,
    n_real                   = crystal_gridding.n_real()).mask_data
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
  f_mask1 = miller_array.structure_factors_from_map(map=mask_data1,
    use_scale = True, anomalous_flag = False, use_sg = False)
  f_mask2 = miller_array.structure_factors_from_map(map=mask_data2,
    use_scale = True, anomalous_flag = False, use_sg = False)
  return [f_mask1.data(), f_mask2.data()]

class cpp_tg(object):
  def __init__(self, fmodel):
    adopt_init_args(self, locals())
    self.f_masks = get_f_masks(
      xrs=self.fmodel.xray_structure, miller_array=self.fmodel.f_calc())

  def get_tg(self, k_sol, b_sol, u_star):
    fmodel_data = ext.core(
      f_calc        = self.fmodel.f_calc().data(),
      shell_f_masks = self.f_masks,
      k_sols        = k_sol,
      b_sols        = b_sol,
      f_part1       = flex.complex_double(self.fmodel.f_calc().data().size(),0),
      f_part2       = flex.complex_double(self.fmodel.f_calc().data().size(),0),
      u_star        = list(u_star),
      hkl           = self.fmodel.f_calc().indices(),
      ss            = self.fmodel.ss)
    return bss.bulk_solvent.ls_kbp_sol_u_star(
      f_model     = fmodel_data,
      f_obs       = self.fmodel.f_obs().data(),
      scale       = 10.0,
      kb_sol_grad = True,
      p_sol_grad  = False,
      u_star_grad = True,
      kb_sol_curv = True,
      p_sol_curv  = False)

class cpp_tg_u_star_only(object):
  def __init__(self, fmodel):
    adopt_init_args(self, locals())
    self.f_masks = get_f_masks(
      xrs=self.fmodel.xray_structure, miller_array=self.fmodel.f_calc())

  def get_tg(self, k_sol, b_sol, u_star):
    fmodel_data = ext.core(
      f_calc        = self.fmodel.f_calc().data(),
      shell_f_masks = self.f_masks,
      k_sols        = k_sol,
      b_sols        = b_sol,
      f_part1       = flex.complex_double(self.fmodel.f_calc().data().size(),0),
      f_part2       = flex.complex_double(self.fmodel.f_calc().data().size(),0),
      u_star        = list(u_star),
      hkl           = self.fmodel.f_calc().indices(),
      ss            = self.fmodel.ss)
    return bss.bulk_solvent.ls_u_star(
      f_model_abs_no_k_total = flex.abs(fmodel_data.f_model_no_aniso_scale),
      f_obs                  = self.fmodel.f_obs().data(),
      miller_indices         = self.fmodel.f_calc().indices(),
      k_anisotropic          = fmodel_data.k_anisotropic)

def fd(TGO, k_sol, b_sol, u_star, param, e = 1.e-3):
  tg = TGO.get_tg(k_sol=k_sol, b_sol=b_sol, u_star=u_star).target()
  g, c = [], []
  for shift in [flex.double([e,0]),
                flex.double([0,e])]:
    if(param=="k_sol"):
      tg_p = TGO.get_tg(k_sol=k_sol+shift, b_sol=b_sol, u_star=u_star).target()
      tg_m = TGO.get_tg(k_sol=k_sol-shift, b_sol=b_sol, u_star=u_star).target()
    if(param=="b_sol"):
      tg_p = TGO.get_tg(k_sol=k_sol, b_sol=b_sol+shift, u_star=u_star).target()
      tg_m = TGO.get_tg(k_sol=k_sol, b_sol=b_sol-shift, u_star=u_star).target()
    if(param in ["k_sol", "b_sol"]):
      g.append((tg_p-tg_m)/(2*e))
      c.append((tg_p+tg_m-2*tg)/(e**2))
  if(param=="u_star"):
    e=e**2/10. # because u_star is very small
    for shift in [flex.double([e,0,0,0,0,0]),
                  flex.double([0,e,0,0,0,0]),
                  flex.double([0,0,e,0,0,0]),
                  flex.double([0,0,0,e,0,0]),
                  flex.double([0,0,0,0,e,0]),
                  flex.double([0,0,0,0,0,e])]:
      tg_p = TGO.get_tg(k_sol=k_sol, b_sol=b_sol, u_star=u_star+shift).target()
      tg_m = TGO.get_tg(k_sol=k_sol, b_sol=b_sol, u_star=u_star-shift).target()
      g.append((tg_p-tg_m)/(2*e))
  return g, c

def run_group(symbol):
  group = space_group_info(symbol);
  print("\n==")
  elements = ('C', 'N', 'O', 'H')*11
  xrs = random_structure.xray_structure(
    space_group_info = group,
    volume_per_atom = 25.,
    general_positions_only = False,
    elements = elements,
    min_distance = 1.0)
  fo = abs(xrs.structure_factors(d_min=2).f_calc())
  fmodel = mmtbx.f_model.manager(
    f_obs = fo,
    xray_structure = xrs)
  #
  k_sol=flex.double([10.35,5.34])
  b_sol=flex.double([30.0, 24.0])
  b_cart = [10,20,30,40,50,60]
  u_star = flex.double(adptbx.b_as_u(
    adptbx.u_cart_as_u_star(xrs.unit_cell(), b_cart)))
  #
  TGO = cpp_tg(fmodel=fmodel)
  tg = TGO.get_tg(k_sol=k_sol, b_sol=b_sol, u_star=u_star)
  # k_sol
  gk_a=list(tg.grad_k_sols())
  ck_a=list(tg.curv_k_sols())
  gk_fd, ck_fd=fd(TGO=TGO, k_sol=k_sol, b_sol=b_sol, u_star=u_star, param="k_sol")
  # b_sol
  gb_a=list(tg.grad_b_sols())
  cb_a=list(tg.curv_b_sols())
  gb_fd, cb_fd=fd(TGO=TGO, k_sol=k_sol, b_sol=b_sol, u_star=u_star, param="b_sol")
  # u_star
  gu_a=list(tg.grad_u_star())
  gu_fd, junk=fd(TGO=TGO, k_sol=k_sol, b_sol=b_sol, u_star=u_star, param="u_star")
  print("u_star:",gu_a)
  print("u_star:",gu_fd)
  TGO2 = cpp_tg_u_star_only(fmodel=fmodel)
  tg2 = TGO2.get_tg(k_sol=k_sol, b_sol=b_sol, u_star=u_star)
  gu_a2=list(tg2.grad_u_star())
  gu_fd2, junk=fd(TGO=TGO2, k_sol=k_sol, b_sol=b_sol, u_star=u_star, param="u_star")
  print("u_star:",gu_a2)
  print("u_star:",gu_fd2)
  #
  print("k_sol:", gk_a,  ck_a)
  print("k_sol:", gk_fd, ck_fd)
  print("b_sol:", gb_a,  cb_a)
  print("b_sol:", gb_fd, cb_fd)
  #
  assert approx_equal(gk_a, gk_fd, eps=1.e-4)
  assert approx_equal(gb_a, gb_fd, eps=1.e-4)
  assert approx_equal(ck_a, ck_fd, eps=1.e-4)
  assert approx_equal(cb_a, cb_fd, eps=1.e-4)
  assert approx_equal(gu_a, gu_fd, eps=1.e-4)
  assert approx_equal(gu_a2, gu_fd2, eps=1.e-5)

if (__name__ == "__main__"):
  for i in list(range(1,231))[:10]: # XXX Do first 10 (to save time)
    print(i)
    run_group(i)
