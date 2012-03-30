from mmtbx import bulk_solvent
import iotbx.phil
from cctbx.array_family import flex
from scitbx import lbfgs
from libtbx import adopt_init_args
from cctbx import adptbx

import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")

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
  number_of_macro_cycles = 2
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

class kbu_minimizer(object):
  def __init__(self,
               fmodel_core_data,
               f_obs,
               k_initial,
               b_initial,
               u_initial,
               refine_k,
               refine_b,
               refine_u,
               min_iterations,
               max_iterations,
               fmodel_core_data_twin = None,
               twin_fraction = None,
               symmetry_constraints_on_b_cart = True,
               u_min_max = 500.,
               u_min_min =-500.,
               k_sol_max = 10.,
               k_sol_min =-10.,
               b_sol_max = 500.,
               b_sol_min =-500.):
    adopt_init_args(self, locals())
    if(twin_fraction == 0):
      twin_fraction = None
      self.twin_fraction = None
    assert [fmodel_core_data_twin,twin_fraction].count(None) in [0,2]
    self.n_shells = self.fmodel_core_data.data.n_shells()
    if not self.fmodel_core_data_twin is None:
      assert self.fmodel_core_data_twin.data.n_shells() == self.n_shells
    assert self.n_shells > 0  and self.n_shells <= 10
    self.k_min = self.k_initial
    assert len(self.k_min) == self.n_shells
    self.b_min = self.b_initial
    self.u_min = self.u_initial
    self.u_factor = self.fmodel_core_data.uc.volume()**(2/3.)
    if(self.symmetry_constraints_on_b_cart):
      self.adp_constraints = self.f_obs.space_group().adp_constraints()
      u_star = self.f_obs.space_group().average_u_star(u_star = self.u_initial)
      self.dim_u = self.adp_constraints.n_independent_params()
      assert self.dim_u <= 6
      independent_params = self.adp_constraints.independent_params(u_star)
      self.x = self.pack(
        u=independent_params,
        k=self.k_min,
        b=self.b_min,
        u_factor=self.u_factor)
    else:
      self.dim_u = len(self.u_initial)
      assert self.dim_u == 6
      self.x = self.pack(
        u=flex.double(self.u_min),
        k=self.k_min,
        b=self.b_min,
        u_factor=self.u_factor)
    lbfgs_exception_handling_params = lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)
    self.minimizer = lbfgs.run(
      target_evaluator = self,
      core_params = lbfgs.core_parameters(),
      termination_params = lbfgs.termination_parameters(
        min_iterations            = min_iterations,
        max_iterations            = max_iterations),
        exception_handling_params = lbfgs_exception_handling_params)
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, u, k, b, u_factor):
    k = list(k)
    assert type(k) is list
    assert len(k) == self.n_shells
    v = []
    if (self.refine_u): v += [ui*u_factor for ui in u]
    if (self.refine_k): v += k
    if (self.refine_b): v.append(b)
    return flex.double(v)

  def unpack_x(self):
    i = 0
    if(self.refine_u):
      if(self.symmetry_constraints_on_b_cart):
        self.u_min = list(self.adp_constraints.all_params(
          iter(self.x[i:self.dim_u]/self.u_factor)))
      else:
        self.u_min = list(iter(self.x[i:self.dim_u]/self.u_factor))
      i = self.dim_u
    if(self.refine_k):
      self.k_min = list(iter(self.x[i:i+self.n_shells]))
      assert len(self.k_min)==self.n_shells
      i += self.n_shells
    if(self.refine_b):
      self.b_min = self.x[i]

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.u_min = tuple([max(self.u_min_min, min(self.u_min_max, v))
      for v in self.u_min])
    if(self.b_min > self.b_sol_max): self.b_min = self.b_sol_max
    if(self.b_min < self.b_sol_min): self.b_min = self.b_sol_min
    self.k_min = [max(self.k_sol_min, min(self.k_sol_max, v))
      for v in self.k_min]
    if(self.twin_fraction is None):
      self.fmodel_core_data.update(k_sols = self.k_min, b_sol = self.b_min,
        u_star = self.u_min)
      tg = bulk_solvent.bulk_solvent_and_aniso_scale_target_and_grads_ls(
        fm                  = self.fmodel_core_data.data,
        fo                  = self.f_obs.data(),
        compute_k_sol_grad  = self.refine_k,
        compute_b_sol_grad  = self.refine_b,
        compute_u_star_grad = self.refine_u)
    else:
      self.fmodel_core_data.update(k_sols = self.k_min, b_sol = self.b_min,
        u_star = self.u_min)
      self.fmodel_core_data_twin.update(k_sols = self.k_min, b_sol = self.b_min,
        u_star = self.u_min)
      tg = bulk_solvent.bulk_solvent_and_aniso_scale_target_and_grads_ls(
        fm1                 = self.fmodel_core_data.data,
        fm2                 = self.fmodel_core_data_twin.data,
        twin_fraction       = self.twin_fraction,
        fo                  = self.f_obs.data(),
        compute_k_sol_grad  = self.refine_k,
        compute_b_sol_grad  = self.refine_b,
        compute_u_star_grad = self.refine_u)
    self.f = tg.target()
    gk=[0.0]*self.n_shells
    gb=0
    gu=[0,0,0,0,0,0]
    if(self.refine_k or self.refine_b):
      gk = list(tg.grad_k_sols())
      assert len(gk) == self.n_shells
      gb = tg.grad_b_sol()
    if(self.refine_u): gu = list(tg.grad_u_star())
    if(self.symmetry_constraints_on_b_cart and self.refine_u):
      independent_params = flex.double(
        self.adp_constraints.independent_gradients(all_gradients=gu))
      self.g = self.pack(
        u=independent_params,
        k=gk,
        b=gb,
        u_factor=1/self.u_factor)
    else:
      self.g = self.pack(u=gu, k=gk, b=gb, u_factor=1/self.u_factor)
    return self.f, self.g

def k_sol_b_sol_b_cart_minimizer(
      fmodels,
      params        = None,
      refine_k_sol  = False,
      refine_b_sol  = False,
      refine_u_star = False):
  if(params is None): params = master_params.extract()
  fmodel_core_data_work = fmodels.fmodel.core_data_work()
  if(fmodels.fmodel_twin is not None):
    fmodel_core_data_work1 = fmodels.fmodel_twin.core_data_work()
  else: fmodel_core_data_work1=None
  return kbu_minimizer(
    fmodel_core_data = fmodel_core_data_work,
    fmodel_core_data_twin= fmodel_core_data_work1,
    twin_fraction    = fmodels.twin_fraction,
    f_obs            = fmodels.fmodel.f_obs,
    k_initial        = fmodel_core_data_work.data.k_sols(),
    b_initial        = fmodel_core_data_work.data.b_sol,
    u_initial        = fmodel_core_data_work.data.u_star,
    refine_k         = refine_k_sol,
    refine_b         = refine_b_sol,
    refine_u         = refine_u_star,
    max_iterations   = params.max_iterations,
    min_iterations   = params.min_iterations,
    symmetry_constraints_on_b_cart = params.symmetry_constraints_on_b_cart)

def _extract_fix_b_cart(fix_b_cart_scope):
  fbs = fix_b_cart_scope
  b_cart = [fbs.b11,fbs.b22,fbs.b33,fbs.b12,fbs.b13,fbs.b23]
  if(b_cart.count(None) > 0): return None
  else: return b_cart

class fmodels_kbu(object):
  def __init__(self, fmodel, fmodel_twin=None, twin_fraction=None):
    adopt_init_args(self, locals())
    if(self.fmodel_twin is not None):
      assert self.fmodel.f_obs.indices().all_eq(self.fmodel_twin.f_obs.indices())

  def update(self, k_sols=None, b_sol=None, b_cart=None):
    self.fmodel.update(k_sols = k_sols, b_sol = b_sol, b_cart = b_cart)
    if(self.fmodel_twin is not None):
      self.fmodel_twin.update(k_sols = k_sols, b_sol = b_sol, b_cart = b_cart)

  def r_factor(self):
    if(self.fmodel_twin is None):
      return self.fmodel.r_factor()
    else:
      return bulk_solvent.r_factor(
        self.fmodel.f_obs.data(),
        self.fmodel.data.f_model,
        self.fmodel_twin.data.f_model,
        self.twin_fraction)

  def select(self, selection):
    fmodel_twin = None
    if(self.fmodel_twin is not None):
      fmodel_twin = self.fmodel_twin.select(selection = selection)
    return fmodels_kbu(
      fmodel        = self.fmodel.select(selection = selection),
      fmodel_twin   = fmodel_twin,
      twin_fraction = self.twin_fraction)

class bulk_solvent_and_scales(object):

  def __init__(self,
               fmodel_kbu = None,
               fmodel_kbu_twin = None,
               twin_fraction = None,
               params = None,
               log    = None,
               nproc  = None):
    self.fmodels = fmodels_kbu(
      fmodel        = fmodel_kbu,
      fmodel_twin   = fmodel_kbu_twin,
      twin_fraction = twin_fraction)
    start_target = self.fmodels.r_factor()
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
      if(self.params.bulk_solvent):
        assert not self.fmodels.fmodel.check_f_mask_all_zero()
      macro_cycles = range(1, self.params.number_of_macro_cycles+1)
      mask_ok = not self.fmodels.fmodel.check_f_mask_all_zero()
      if(self.params.fix_k_sol is not None and mask_ok):
        assert self.params.bulk_solvent
        assert not self.params.k_sol_b_sol_grid_search
        assert not self.params.minimization_k_sol_b_sol
        self.fmodels.update(k_sols = self.params.fix_k_sol)
      if(self.params.fix_b_sol is not None and mask_ok):
        assert self.params.bulk_solvent
        assert not self.params.k_sol_b_sol_grid_search
        assert not self.params.minimization_k_sol_b_sol
        self.fmodels.update(b_sol = self.params.fix_b_sol)
      fix_b_cart = _extract_fix_b_cart(fix_b_cart_scope = self.params.fix_b_cart)
      if(fix_b_cart is not None):
        assert self.params.anisotropic_scaling
        assert not self.params.minimization_b_cart
        self.fmodels.update(b_cart = fix_b_cart)
      target = self.fmodels.r_factor()
      for mc in macro_cycles:
        if(self.params.k_sol_b_sol_grid_search and mc == macro_cycles[0]):
          ksol,bsol,b_cart,target = self._ksol_bsol_grid_search(nproc=nproc)
          self.fmodels.update(k_sols = ksol, b_sol = bsol, b_cart = b_cart)
          assert abs(target-self.fmodels.r_factor()) < 1.e-6
        if(self.params.minimization_k_sol_b_sol):
          ksol, bsol, target = self._ksol_bsol_cart_minimizer()
          self.fmodels.update(k_sols = ksol, b_sol = bsol)
          assert abs(target-self.fmodels.r_factor())<1.e-6
        if(self.params.minimization_b_cart):
          b_cart,target = self._b_cart_minimizer()
          self.fmodels.update(b_cart = b_cart)
          assert abs(target-self.fmodels.r_factor())<1.e-6
      ksols = list(self.k_sols()[:])
      do_update = False
      for ik in range(len(ksols)):
        if(abs(ksols[ik]) < 0.01):
          ksols[ik] = 0.
          do_update = True
      if( do_update ):
        if( ksols.count(0) == len(ksols) ):
          bsol = 0.
        else:
          bsol = fmodel.b_sol()
        self.fmodels.update(k_sols = ksols, b_sol = bsol)

  def k_sol(self, i):
    return self.fmodels.fmodel.k_sol(i)

  def k_sols(self):
    return self.fmodels.fmodel.k_sols()

  def b_sol(self):
    return self.fmodels.fmodel.b_sol()

  def b_cart(self):
    return self.fmodels.fmodel.b_cart()

  def u_star(self):
    return self.fmodels.fmodel.u_star()

  def _ksol_bsol_grid_search(self, nproc=None):
    start_r_factor = self.fmodels.r_factor()
    final_ksol   = self.fmodels.fmodel.data.k_sols()
    final_bsol   = self.fmodels.fmodel.data.b_sol
    final_b_cart = self.fmodels.fmodel.b_cart()
    final_r_factor = start_r_factor
    self._ksol_len = len(final_ksol)
    k_sols = kb_range(self.params.k_sol_grid_search_max,
                      self.params.k_sol_grid_search_min,
                      self.params.k_sol_step)
    b_sols = kb_range(self.params.b_sol_grid_search_max,
                      self.params.b_sol_grid_search_min,
                      self.params.b_sol_step)
    args = []
    for ksol in k_sols:
      for bsol in b_sols:
        args.append((ksol,bsol))
    results = []
    for args_ in args :
      res = self.try_ksol_bsol(args_)
      results.append(res)
    for result in results :
      if(result.r_factor < final_r_factor):
        final_r_factor = result.r_factor
        final_ksol = result.k_sol
        final_bsol = result.b_sol
        final_b_cart = result.b_cart
        self.fmodels.update(k_sols=final_ksol,b_sol=final_bsol,
          b_cart=final_b_cart)
    self.fmodels.update(
      k_sols = final_ksol,
      b_sol  = final_bsol,
      b_cart = final_b_cart)
    assert abs(self.fmodels.r_factor()-final_r_factor) < 1.e-6
    assert self.fmodels.r_factor() <= start_r_factor
    return final_ksol, final_bsol, final_b_cart, final_r_factor

  def try_ksol_bsol(self, args):
    ksol, bsol = args
    for bc in self.b_cart():
      if(abs(bc) > 300.):
        self.fmodels.update(b_cart = [0,0,0,0,0,0])
        break
    ksol_list = [ksol]*self._ksol_len
    self.fmodels.update(k_sols = ksol_list, b_sol = bsol)
    if(self.params.minimization_k_sol_b_sol):
      ksol_, bsol_, dummy = self._ksol_bsol_cart_minimizer()
      self.fmodels.update(k_sols = ksol_, b_sol = bsol_)
    if(self.params.minimization_b_cart):
      b_cart, dummy = self._b_cart_minimizer()
      self.fmodels.update(b_cart = b_cart)
    return ksol_bsol_result(
      r_factor = self.fmodels.r_factor(),
      k_sol    = self.fmodels.fmodel.data.k_sols(),
      b_sol    = self.fmodels.fmodel.data.b_sol,
      b_cart   = self.fmodels.fmodel.b_cart())

  def _ksol_bsol_cart_minimizer(self):
    start_r_factor = self.fmodels.r_factor()
    final_ksol     = self.fmodels.fmodel.data.k_sols()
    final_bsol     = self.fmodels.fmodel.data.b_sol
    final_r_factor = self.fmodels.r_factor()
    ksol, bsol     = self._k_sol_b_sol_minimization_helper()
    self.fmodels.update(k_sols = ksol, b_sol = bsol)
    r_factor = self.fmodels.r_factor()
    if(r_factor < final_r_factor):
      final_ksol = ksol
      final_bsol = bsol
      final_r_factor = r_factor
    assert final_r_factor <= start_r_factor
    return final_ksol, final_bsol, final_r_factor

  def _b_cart_minimizer(self):
    start_r_factor = self.fmodels.r_factor()
    final_b_cart   = self.fmodels.fmodel.b_cart()
    final_r_factor = self.fmodels.r_factor()
    b_cart = self._b_cart_minimizer_helper()
    self.fmodels.update(b_cart = b_cart)
    r_factor = self.fmodels.r_factor()
    if(r_factor < final_r_factor):
      final_b_cart = b_cart
      final_r_factor = r_factor
    assert final_r_factor <= start_r_factor
    return final_b_cart, final_r_factor

  def _b_cart_minimizer_helper(self, n_macro_cycles = 2):
    r_start = self.fmodels.r_factor()
    b_start = self.fmodels.fmodel.b_cart()
    for u_cycle in xrange(n_macro_cycles):
      u_min = k_sol_b_sol_b_cart_minimizer(
        fmodels       = self.fmodels,
        params        = self.params,
        refine_u_star = True).u_min
      b_cart = adptbx.u_as_b(
        adptbx.u_star_as_u_cart(self.fmodels.fmodel.f_obs.unit_cell(),u_min))
      self.fmodels.update(b_cart = b_cart)
    r_final = self.fmodels.r_factor()
    if(r_final >= r_start):
      self.fmodels.update(b_cart = b_start)
      return b_start
    else: return b_cart

  def _k_sol_b_sol_minimization_helper(self):
    ksol_orig = self.fmodels.fmodel.data.k_sols()
    bsol_orig = self.fmodels.fmodel.data.b_sol
    r_start   = self.fmodels.r_factor()
    minimizer_obj = k_sol_b_sol_b_cart_minimizer(
      fmodels      = self.fmodels,
      params       = self.params,
      refine_k_sol = True,
      refine_b_sol = True)
    ksol, bsol = minimizer_obj.k_min, minimizer_obj.b_min
    assert type(ksol) is list
    assert len(ksol) >= 1
    ksol = [max(self.params.k_sol_min, min(self.params.k_sol_max, v))
      for v in ksol]
    if(bsol < self.params.b_sol_min): bsol = self.params.b_sol_min
    if(bsol > self.params.b_sol_max): bsol = self.params.b_sol_max
    self.fmodels.update(k_sols = ksol, b_sol = bsol)
    r_end = self.fmodels.r_factor()
    if(r_end >= r_start):
      self.fmodels.update(k_sols = ksol_orig, b_sol = bsol_orig)
      return ksol_orig, bsol_orig
    else: return ksol, bsol

def kb_range(x_max, x_min, step):
  sc = 1000.
  return [i/sc for i in range(int(x_min*sc), int(x_max*sc)+1, int(step*sc))]

class ksol_bsol_result(object):
  def __init__ (self, k_sol, b_sol, b_cart, r_factor):
    self.k_sol = k_sol
    self.b_sol = b_sol
    self.b_cart = b_cart
    self.r_factor = r_factor

class u_star_minimizer(object):
  def __init__(self,
               fmodel_core_data,
               f_obs,
               u_initial=[0,0,0,0,0,0],
               refine_u=True,
               min_iterations=500,
               max_iterations=500,
               symmetry_constraints_on_b_cart = True,
               u_min_max = 500.,
               u_min_min =-500.):
    adopt_init_args(self, locals())
    self.u_min = self.u_initial
    self.u_factor = self.fmodel_core_data.uc.volume()**(2/3.)
    if(self.symmetry_constraints_on_b_cart):
      self.adp_constraints = self.f_obs.space_group().adp_constraints()
      u_star = self.f_obs.space_group().average_u_star(u_star = self.u_initial)
      self.dim_u = self.adp_constraints.n_independent_params()
      assert self.dim_u <= 6
      independent_params = self.adp_constraints.independent_params(u_star)
      self.x = self.pack(
        u=independent_params,
        u_factor=self.u_factor)
    else:
      self.dim_u = len(self.u_initial)
      assert self.dim_u == 6
      self.x = self.pack(
        u=flex.double(self.u_min),
        u_factor=self.u_factor)
    lbfgs_exception_handling_params = lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)
    self.minimizer = lbfgs.run(
      target_evaluator = self,
      core_params = lbfgs.core_parameters(),
      termination_params = lbfgs.termination_parameters(
        min_iterations            = min_iterations,
        max_iterations            = max_iterations),
        exception_handling_params = lbfgs_exception_handling_params)
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, u, u_factor):
    v = []
    if (self.refine_u): v += [ui*u_factor for ui in u]
    return flex.double(v)

  def unpack_x(self):
    i = 0
    if(self.refine_u):
      if(self.symmetry_constraints_on_b_cart):
        self.u_min = list(self.adp_constraints.all_params(
          iter(self.x[i:self.dim_u]/self.u_factor)))
      else:
        self.u_min = list(iter(self.x[i:self.dim_u]/self.u_factor))
      i = self.dim_u

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.u_min = tuple([max(self.u_min_min, min(self.u_min_max, v))
      for v in self.u_min])
    self.fmodel_core_data.update(k_sols = [0], b_sol = 0, u_star = self.u_min)
    tg = bulk_solvent.bulk_solvent_and_aniso_scale_target_and_grads_ls(
      fm                  = self.fmodel_core_data.data,
      fo                  = self.f_obs.data(),
      compute_k_sol_grad  = False,
      compute_b_sol_grad  = False,
      compute_u_star_grad = self.refine_u)
    self.f = tg.target()
    gu=[0,0,0,0,0,0]
    if(self.refine_u): gu = list(tg.grad_u_star())
    if(self.symmetry_constraints_on_b_cart and self.refine_u):
      independent_params = flex.double(
        self.adp_constraints.independent_gradients(all_gradients=gu))
      self.g = self.pack(
        u=independent_params,
        u_factor=1/self.u_factor)
    else:
      self.g = self.pack(u=gu, u_factor=1/self.u_factor)
    return self.f, self.g
