from __future__ import absolute_import, division, print_function
from mmtbx import bulk_solvent
import iotbx.phil
from cctbx.array_family import flex
from cctbx import adptbx
from libtbx import group_args

import boost_adaptbx.boost.python as bp
from six.moves import range
ext = bp.import_ext("mmtbx_f_model_ext")

master_params = iotbx.phil.parse("""\
  bulk_solvent = True
    .type = bool
  anisotropic_scaling = True
    .type = bool
  k_sol_b_sol_grid_search = True
    .type = bool
    .expert_level=2
  minimization_k_sol_b_sol = True
    .type = bool
    .expert_level=2
  minimization_b_cart = True
    .type = bool
    .expert_level=2
  target = ls_wunit_k1 *ml
    .type = choice
  symmetry_constraints_on_b_cart = True
    .type = bool
    .expert_level=2
  k_sol_max = 0.6
    .type = float
    .expert_level=2
  k_sol_min = 0.0
    .type = float
    .expert_level=2
  b_sol_max = 300.0
    .type = float
    .expert_level=2
  b_sol_min = 0.0
    .type = float
    .expert_level=2
  k_sol_grid_search_max = 0.6
    .type = float
    .expert_level=2
  k_sol_grid_search_min = 0.0
    .type = float
    .expert_level=2
  b_sol_grid_search_max = 80.0
    .type = float
    .expert_level=2
  b_sol_grid_search_min = 20.0
    .type = float
    .expert_level=2
  k_sol_step = 0.2
    .type = float
    .expert_level=2
  b_sol_step = 20.0
    .type = float
    .expert_level=2
  number_of_macro_cycles = 1
    .type = int
    .expert_level=2
  max_iterations = 25
    .type = int
    .expert_level=3
  min_iterations = 25
    .type = int
    .expert_level=3
  fix_k_sol = None
    .type = float
    .expert_level=2
  fix_b_sol = None
    .type = float
    .expert_level=2
  fix_b_cart
    .expert_level=2
    .style = box
  {
    b11 = None
      .type = float
    b22 = None
      .type = float
    b33 = None
      .type = float
    b12 = None
      .type = float
    b13 = None
      .type = float
    b23 = None
      .type = float
  }
""")

def k_sol_b_sol_b_cart_minimizer(
      fmodel_kbu,
      params        = None,
      refine_k_sol  = False,
      refine_b_sol  = False,
      refine_u_star = False,
      refine_kbu    = False):
  if(params is None): params = master_params.extract()
  fmodel_core_data_work = fmodel_kbu.core_data_work()
  import mmtbx.bulk_solvent.kbu_refinery as kbu_refinery
  obj = kbu_refinery.tgc(
    f_obs   = fmodel_kbu.f_obs,
    f_calc  = fmodel_kbu.f_calc,
    f_masks = fmodel_kbu.f_masks,
    ss      = fmodel_kbu.ss,
    k_sols  = list(fmodel_core_data_work.data.k_sols()),
    b_sols  = list(fmodel_core_data_work.data.b_sols()),
    u_star  = fmodel_core_data_work.data.u_star)
  if(not refine_kbu):
    if(refine_k_sol and refine_b_sol):
      obj.set_refine_kb()
      obj.minimize_kb(use_curvatures_options=[False,True],
        n_cycles=params.number_of_macro_cycles)
    if(refine_u_star):
      obj.set_refine_u()
      obj.minimize_u(n_cycles=params.number_of_macro_cycles)
  else:
    obj.minimize_kbu(n_cycles=5)
  return obj

def _extract_fix_b_cart(fix_b_cart_scope):
  fbs = fix_b_cart_scope
  b_cart = [fbs.b11,fbs.b22,fbs.b33,fbs.b12,fbs.b13,fbs.b23]
  if(b_cart.count(None) > 0): return None
  else: return b_cart

class bulk_solvent_and_scales(object):

  def __init__(self,
               fmodel_kbu = None,
               params = None,
               log    = None):
    self.fmodel_kbu = fmodel_kbu
    start_target = self.fmodel_kbu.r_factor()
    self.params = params
    self.log = log
    if(self.params is None): self.params = master_params.extract()
    if([self.params.bulk_solvent,
        self.params.anisotropic_scaling].count(True) > 0):
      if(not self.params.bulk_solvent):
        self.params.k_sol_b_sol_grid_search = False
        self.params.minimization_k_sol_b_sol = False
      if(not self.params.anisotropic_scaling):
        self.params.minimization_b_cart = False
      params_target = self.params.target
      if(self.fmodel_kbu.check_f_mask_all_zero()):
        self.params.bulk_solvent = False
        self.params.k_sol_b_sol_grid_search = False
        self.params.minimization_k_sol_b_sol = False
      macro_cycles = list(range(1, self.params.number_of_macro_cycles+1))
      mask_ok = not self.fmodel_kbu.check_f_mask_all_zero()
      if(self.params.fix_k_sol is not None and mask_ok):
        assert self.params.bulk_solvent
        assert not self.params.k_sol_b_sol_grid_search
        assert not self.params.minimization_k_sol_b_sol
        self.fmodel_kbu.update(k_sols = self.params.fix_k_sol)
      if(self.params.fix_b_sol is not None and mask_ok):
        assert self.params.bulk_solvent
        assert not self.params.k_sol_b_sol_grid_search
        assert not self.params.minimization_k_sol_b_sol
        self.fmodel_kbu.update(b_sols = self.params.fix_b_sol)
      fix_b_cart = _extract_fix_b_cart(fix_b_cart_scope= self.params.fix_b_cart)
      if(fix_b_cart is not None):
        assert self.params.anisotropic_scaling
        assert not self.params.minimization_b_cart
        self.fmodel_kbu.update(b_cart = fix_b_cart)
      target = self.fmodel_kbu.r_factor()
      for mc in macro_cycles:
        if(self.params.k_sol_b_sol_grid_search and mc == macro_cycles[0] and
           len(self.fmodel_kbu.f_masks)==1):
          ksol,bsol,b_cart,target = self._ksol_bsol_grid_search()
          self.fmodel_kbu.update(k_sols = ksol, b_sols = bsol, b_cart = b_cart)
          assert abs(target-self.fmodel_kbu.r_factor()) < 1.e-6
        if(self.params.minimization_k_sol_b_sol):
          ksol, bsol, target = self._ksol_bsol_cart_minimizer()
          self.fmodel_kbu.update(k_sols = ksol, b_sols = bsol)
          assert abs(target-self.fmodel_kbu.r_factor())<1.e-6
        if(self.params.minimization_b_cart):
          b_cart,target = self._b_cart_minimizer()
          self.fmodel_kbu.update(b_cart = b_cart)
          assert abs(target-self.fmodel_kbu.r_factor())<1.e-6
      ksols = list(self.k_sols()[:])
      do_update = False
      for ik in range(len(ksols)):
        if(abs(ksols[ik]) < 0.01):
          ksols[ik] = 0.
          do_update = True
      if( do_update ):
        if( ksols.count(0) == len(ksols) ):
          bsol = [0.]
        if(len(ksols)>len(bsol) and len(bsol)==1):
          bsol = bsol*len(ksols)
        self.fmodel_kbu.update(k_sols = ksols, b_sols = bsol)

  def k_sols(self):
    return self.fmodel_kbu.k_sols()

  def b_sols(self):
    return self.fmodel_kbu.b_sols()

  def b_cart(self):
    return self.fmodel_kbu.b_cart()

  def format_scale_matrix(self, m=None, log=None):
    sm = m
    if(sm is None): sm = self.b_cart()
    out = log
    if(sm is None):
      print("  k_anisotropic=1", file=log)
      return
    if(len(sm)<=6):
      print("      b_cart(11,22,33,12,13,23):",\
        ",".join([str("%8.4f"%i).strip() for i in sm]), file=out)

  def u_star(self):
    return self.fmodel_kbu.u_star()

  def _ksol_bsol_grid_search(self):
    # XXX HERE
    start_r_factor = self.fmodel_kbu.r_factor()
    final_ksol   = self.fmodel_kbu.data.k_sols()
    final_bsol   = self.fmodel_kbu.data.b_sols()
    final_b_cart = self.fmodel_kbu.b_cart()
    final_r_factor = start_r_factor

    k_sols = kb_range(0.6,
                      0.0,
                      0.1)
    b_sols = kb_range(80.,
                      10,
                      5)
    assert len(self.fmodel_kbu.f_masks)==1
    o = bulk_solvent.k_sol_b_sol_k_anisotropic_scaler_twin(
      f_obs          = self.fmodel_kbu.f_obs.data(),
      f_calc         = self.fmodel_kbu.f_calc.data(),
      f_mask         = self.fmodel_kbu.f_masks[0].data(),
      ss             = self.fmodel_kbu.ss,
      k_sol_range    = flex.double(k_sols),
      b_sol_range    = flex.double(b_sols),
      miller_indices = self.fmodel_kbu.f_obs.indices(),
      r_ref          = start_r_factor)
    if(o.updated()):
      assert o.r() < start_r_factor
      final_r_factor = o.r()
      final_ksol     = o.k_sol()
      final_bsol     = o.b_sol()
      final_b_cart   = adptbx.u_as_b(adptbx.u_star_as_u_cart(
        self.fmodel_kbu.f_obs.unit_cell(), o.u_star()))
      self.fmodel_kbu.update(
        k_sols=[final_ksol],
        b_sols=[final_bsol],
        b_cart=final_b_cart)
      assert self.fmodel_kbu.r_factor() <= start_r_factor
      if(final_ksol < 0.01):
        final_ksol=0
        self.fmodel_kbu.update(
          k_sols=[final_ksol])
        final_r_factor = self.fmodel_kbu.r_factor()
      else:
        o = k_sol_b_sol_b_cart_minimizer(
          fmodel_kbu = self.fmodel_kbu,
          refine_kbu = True)
        if(o.kbu.r_factor() < final_r_factor):
          final_b_cart = adptbx.u_as_b(adptbx.u_star_as_u_cart(
            self.fmodel_kbu.f_obs.unit_cell(), o.kbu.u_star()))
          final_ksol = list(o.kbu.k_sols())
          final_bsol = list(o.kbu.b_sols())
          self.fmodel_kbu.update(
            k_sols=final_ksol,
            b_sols=final_bsol,
            b_cart=final_b_cart)
        else:
          self.fmodel_kbu.update(
            k_sols=[final_ksol],
            b_sols=[final_bsol],
            b_cart=final_b_cart)
        final_r_factor = self.fmodel_kbu.r_factor()
        assert self.fmodel_kbu.r_factor() <= start_r_factor
    self.fmodel_kbu.update(
      k_sols = final_ksol,
      b_sols = final_bsol,
      b_cart = final_b_cart)
    assert abs(self.fmodel_kbu.r_factor()-final_r_factor) < 1.e-6
    assert self.fmodel_kbu.r_factor() <= start_r_factor
    return final_ksol, final_bsol, final_b_cart, final_r_factor

  def _ksol_bsol_cart_minimizer(self):
    start_r_factor = self.fmodel_kbu.r_factor()
    final_ksol     = self.fmodel_kbu.data.k_sols()
    final_bsol     = self.fmodel_kbu.data.b_sols()
    final_r_factor = self.fmodel_kbu.r_factor()
    ksol, bsol     = self._k_sol_b_sol_minimization_helper()
    self.fmodel_kbu.update(k_sols = ksol, b_sols = bsol)
    r_factor = self.fmodel_kbu.r_factor()
    if(r_factor < final_r_factor):
      final_ksol = ksol
      final_bsol = bsol
      final_r_factor = r_factor
    assert final_r_factor <= start_r_factor
    return final_ksol, final_bsol, final_r_factor

  def r_factor(self):
    return self.fmodel_kbu.r_factor()

  def r_all(self):
    return self.r_factor()

  def _b_cart_minimizer(self):
    start_r_factor = self.fmodel_kbu.r_factor()
    final_b_cart   = self.fmodel_kbu.b_cart()
    final_r_factor = self.fmodel_kbu.r_factor()
    b_cart = self._b_cart_minimizer_helper()
    self.fmodel_kbu.update(b_cart = b_cart)
    r_factor = self.fmodel_kbu.r_factor()
    if(r_factor < final_r_factor):
      final_b_cart = b_cart
      final_r_factor = r_factor
    assert final_r_factor <= start_r_factor
    return final_b_cart, final_r_factor

  def _b_cart_minimizer_helper(self, n_macro_cycles = 1):
    r_start = self.fmodel_kbu.r_factor()
    b_start = self.fmodel_kbu.b_cart()
    for u_cycle in range(n_macro_cycles):
      u_min = k_sol_b_sol_b_cart_minimizer(
        fmodel_kbu    = self.fmodel_kbu,
        params        = self.params,
        refine_u_star = True).kbu.u_star()
      b_cart = adptbx.u_as_b(
        adptbx.u_star_as_u_cart(self.fmodel_kbu.f_obs.unit_cell(),u_min))
      self.fmodel_kbu.update(b_cart = b_cart)
    r_final = self.fmodel_kbu.r_factor()
    if(r_final >= r_start):
      self.fmodel_kbu.update(b_cart = b_start)
      return b_start
    else: return b_cart

  def _k_sol_b_sol_minimization_helper(self):
    ksol_orig = self.fmodel_kbu.data.k_sols()
    bsol_orig = self.fmodel_kbu.data.b_sols()
    r_start   = self.fmodel_kbu.r_factor()
    minimizer_obj = k_sol_b_sol_b_cart_minimizer(
      fmodel_kbu   = self.fmodel_kbu,
      params       = self.params,
      refine_k_sol = True,
      refine_b_sol = True)
    ksol = list(minimizer_obj.kbu.k_sols())
    bsol = list(minimizer_obj.kbu.b_sols())
    assert type(ksol) is list
    assert len(ksol) >= 1
    ksol = [max(self.params.k_sol_min, min(self.params.k_sol_max, v))
      for v in ksol]
    for i in range(len(bsol)):
      if(bsol[i] > self.params.b_sol_max): bsol[i] = self.params.b_sol_max
      if(bsol[i] < self.params.b_sol_min): bsol[i] = self.params.b_sol_min
    self.fmodel_kbu.update(k_sols = ksol, b_sols = bsol)
    r_end = self.fmodel_kbu.r_factor()
    if(r_end >= r_start):
      self.fmodel_kbu.update(k_sols = ksol_orig, b_sols = bsol_orig)
      return ksol_orig, bsol_orig
    else: return ksol, bsol

  def k_masks(self):
    return self.fmodel_kbu.k_masks()

  def k_anisotropic(self):
    return self.fmodel_kbu.k_anisotropic()

  def k_isotropic(self):
    return self.fmodel_kbu.k_isotropic()

  def apply_back_trace_of_overall_exp_scale_matrix(self, xray_structure=None):
    if(xray_structure is None): return None
    k_sol, b_sol, b_cart = self.k_sols(), self.b_sols(), self.b_cart()
    assert len(k_sol)==1 # XXX Only one mask!
    k_sol = k_sol[0]
    b_sol = b_sol[0]
    #
    xrs = xray_structure
    if(xrs is None): return
    b_min = min(b_sol, xrs.min_u_cart_eigenvalue()*adptbx.u_as_b(1.))
    if(b_min < 0): xrs.tidy_us()
    b_iso = (b_cart[0]+b_cart[1]+b_cart[2])/3.0
    b_test = b_min+b_iso
    if(b_test < 0.0): b_adj = b_iso + abs(b_test) + 0.001
    else: b_adj = b_iso
    b_cart_new = [b_cart[0]-b_adj,b_cart[1]-b_adj,b_cart[2]-b_adj,
                  b_cart[3],      b_cart[4],      b_cart[5]]
    b_sol_new = b_sol + b_adj
    xrs.shift_us(b_shift = b_adj)
    b_min = min(b_sol_new, xrs.min_u_cart_eigenvalue()*adptbx.u_as_b(1.))
    assert b_min >= 0.0
    xrs.tidy_us()
    #
    assert self.fmodel_kbu
    k_masks = [ext.k_mask(self.fmodel_kbu.ss, k_sol, b_sol_new)]
    u_star=adptbx.u_cart_as_u_star(
      self.fmodel_kbu.f_obs.unit_cell(), adptbx.b_as_u(b_cart_new))
    k_anisotropic = ext.k_anisotropic(self.fmodel_kbu.f_obs.indices(), u_star)
    self.fmodel_kbu = self.fmodel_kbu.update(
      b_sols = [b_sol_new],
      b_cart = b_cart_new)
    return group_args(
      xray_structure = xrs,
      b_adj          = b_adj,
      b_sol          = b_sol_new,
      b_cart         = b_cart_new)

def kb_range(x_max, x_min, step):
  sc = 1000.
  return [i/sc for i in range(int(x_min*sc), int(x_max*sc)+1, int(step*sc))]

class ksol_bsol_result(object):
  def __init__(self, k_sol, b_sol, b_cart, r_factor):
    self.k_sol = k_sol
    self.b_sol = b_sol
    self.b_cart = b_cart
    self.r_factor = r_factor
