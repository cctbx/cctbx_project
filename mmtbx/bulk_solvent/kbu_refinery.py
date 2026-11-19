from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.f_model
from scitbx import lbfgs as scitbx_lbfgs
from libtbx import adopt_init_args
import math
from mmtbx import bulk_solvent
from six.moves import range

def lbfgs_run(target_evaluator,
              min_iterations=0,
              max_iterations=None,
              traditional_convergence_test=1,
              use_curvatures=False):
  ext = scitbx_lbfgs.ext
  minimizer = ext.minimizer(target_evaluator.n)
  minimizer.error = None
  if (traditional_convergence_test):
    is_converged = ext.traditional_convergence_test(target_evaluator.n)
  else:
    raise RuntimeError
    is_converged = ext.drop_convergence_test(min_iterations)
  try:
    icall = 0
    requests_f_and_g = True
    requests_diag = use_curvatures
    while 1:
      if (requests_f_and_g):
        icall += 1
      x, f, g, d = target_evaluator(
        requests_f_and_g=requests_f_and_g,
        requests_diag=requests_diag)
      #if (requests_diag):
      #  print "x,f,d:", tuple(x), f, tuple(d)
      #else:
      #  print "x,f:", tuple(x), f
      if (use_curvatures):
        if (d is None): d = flex.double(x.size())
        have_request = minimizer.run(x, f, g, d)
      else:
        have_request = minimizer.run(x, f, g)
      if (have_request):
        requests_f_and_g = minimizer.requests_f_and_g()
        requests_diag = minimizer.requests_diag()
        continue
      assert not minimizer.requests_f_and_g()
      assert not minimizer.requests_diag()
      if (traditional_convergence_test):
        if (minimizer.iter() >= min_iterations and is_converged(x, g)): break
      else:
        if (is_converged(f)): break
      if (max_iterations is not None and minimizer.iter() >= max_iterations):
        break
      if (use_curvatures):
        have_request = minimizer.run(x, f, g, d)
      else:
        have_request = minimizer.run(x, f, g)
      if (not have_request): break
      requests_f_and_g = minimizer.requests_f_and_g()
      requests_diag = minimizer.requests_diag()
  except RuntimeError as e:
    minimizer.error = str(e)
  minimizer.n_calls = icall
  return minimizer

class refinement_flags(object):
  def __init__(self, refine_k=False, refine_b=False, refine_kb=False,
               refine_u=False):
    adopt_init_args(self, locals())

class minimizer:

  def __init__(self, tgc, min_iterations=0, max_iterations=25):
    adopt_init_args(self, locals())
    self.x = self.tgc.x()
    self.n = self.x.size()

  def run(self, use_curvatures=0):
    self.minimizer = lbfgs_run(
      target_evaluator=self,
      min_iterations=self.min_iterations,
      max_iterations=self.max_iterations,
      use_curvatures=use_curvatures)
    self(requests_f_and_g=True, requests_diag=False)
    return self

  def __call__(self, requests_f_and_g, requests_diag):
    self.tgc.update(x=self.x)
    if (not requests_f_and_g and not requests_diag):
      requests_f_and_g = True
      requests_diag = True
    if (requests_f_and_g):
      self.f = self.tgc.target()
      self.g = self.tgc.gradients()
      self.d = None
    if (requests_diag):
      self.d = self.tgc.curvatures()
      #assert self.d.all_ne(0)
      if(self.d.all_eq(0)): self.d=None
      else:
        self.d = 1 / self.d
    return self.x, self.f, self.g, self.d

class joint_minimizer(object):
  """
  Refine k_sols, b_sols and the overall anisotropic scale together, with
  bounds, by L-BFGS-B. Parameters are scaled so that one unit is 0.1 in
  k_sol, 10 A**2 in b_sol and about 1 A**2 in b_cart.

  The overall scale is the least-squares one at every point, so that the
  gradients with respect to the refined parameters are exact.

  Parameters
  ----------
  tgc : tgc
      Supplies the data and the start; tgc.kbu holds the result.
  k_sol_min, k_sol_max, b_sol_min, b_sol_max : float
      Bounds. The start is moved into the bounds.
  refine_kb, refine_u : bool
      What to refine. The anisotropic scale has no bounds.
  max_iterations : int
      Limit on L-BFGS-B iterations.

  Attributes
  ----------
  n_calls : int
      Number of target and gradient evaluations.
  n_iterations : int
  task : str
      Termination message of L-BFGS-B.
  """
  def __init__(self, tgc, k_sol_min, k_sol_max, b_sol_min, b_sol_max,
               refine_kb=True, refine_u=True, max_iterations=100):
    assert refine_kb or refine_u
    assert k_sol_min <= k_sol_max and b_sol_min <= b_sol_max
    adopt_init_args(self, locals())
    self.kbu = tgc.kbu
    self.f_obs = tgc.f_obs.data()
    self.n_masks = len(self.kbu.f_masks)
    self.constraints = None
    if(tgc.symmetry_constraints): self.constraints = tgc.adp_constraints
    self.scales = flex.double()
    self.lower  = flex.double()
    self.upper  = flex.double()
    self.nbd    = flex.int() # 2: lower and upper bound, 0: no bounds
    if(refine_kb):
      for scale, lower, upper in [(0.1, k_sol_min, k_sol_max),
                                  (10., b_sol_min, b_sol_max)]:
        for i in range(self.n_masks):
          self.scales.append(scale)
          self.lower.append(lower/scale)
          self.upper.append(upper/scale)
          self.nbd.append(2)
    if(refine_u):
      for scale in self._u_star_scales():
        self.scales.append(scale)
        self.lower.append(0)
        self.upper.append(0)
        self.nbd.append(0)
    self.x = self._pack()
    self.n_calls = 0
    self.n_iterations = 0
    self.task = None

  def _u_star_scales(self):
    """
    Change of each refined u_star parameter that changes b_cart by about
    1 A**2.
    """
    a, b, c = self.tgc.f_obs.unit_cell().reciprocal_parameters()[:3]
    f = 1./(8*math.pi**2)
    scales = [f*a*a, f*b*b, f*c*c, f*a*b, f*a*c, f*b*c]
    if(self.constraints is not None):
      scales = list(self.constraints.independent_params(scales))
    return [abs(v) if v != 0 else f*a*a for v in scales]

  def _pack(self):
    x = flex.double()
    if(self.refine_kb):
      x.extend(flex.double(list(self.kbu.k_sols())))
      x.extend(flex.double(list(self.kbu.b_sols())))
    if(self.refine_u):
      u_star = self.kbu.u_star()
      if(self.constraints is not None):
        u_star = self.constraints.independent_params(
          self.tgc.space_group.average_u_star(u_star = u_star))
      x.extend(flex.double(list(u_star)))
    x = x/self.scales
    for i in range(x.size()):
      if(self.nbd[i] == 2):
        x[i] = min(self.upper[i], max(self.lower[i], x[i]))
    return x

  def _unpack(self):
    values = self.x*self.scales
    i = 0
    if(self.refine_kb):
      n = self.n_masks
      self.kbu.update(k_sols = list(values[:n]), b_sols = list(values[n:2*n]))
      i = 2*n
    if(self.refine_u):
      u_star = list(values[i:])
      if(self.constraints is not None):
        u_star = list(self.constraints.all_params(u_star))
      self.kbu.update(u_star = u_star)

  def _target_and_gradients(self):
    self._unpack()
    self.n_calls += 1
    t_g = bulk_solvent.ls_kbp_sol_u_star(
      f_model     = self.kbu.data,
      f_obs       = self.f_obs,
      scale       = bulk_solvent.scale(self.f_obs, self.kbu.data.f_model),
      kb_sol_grad = self.refine_kb,
      p_sol_grad  = False,
      u_star_grad = self.refine_u,
      kb_sol_curv = False,
      p_sol_curv  = False)
    g = flex.double()
    if(self.refine_kb):
      g.extend(t_g.grad_k_sols())
      g.extend(t_g.grad_b_sols())
    if(self.refine_u):
      g_u = t_g.grad_u_star()
      if(self.constraints is not None):
        g_u = self.constraints.independent_gradients(all_gradients = g_u)
      g.extend(flex.double(list(g_u)))
    return t_g.target(), g*self.scales

  def run(self):
    from scitbx import lbfgsb
    m = lbfgsb.minimizer(
      n     = self.x.size(),
      l     = self.lower,
      u     = self.upper,
      nbd   = self.nbd,
      # Stop when the target (order of R**2) decreases by less than
      # machine precision; with data that the model explains exactly this
      # is at R of about 1.e-8.
      factr = 1.,
      pgtol = 1.e-12)
    f, g = self._target_and_gradients()
    # Limit on evaluations guards against an endless line search.
    max_calls = 20*max(1, self.max_iterations)
    while True:
      if(m.process(self.x, f, g)):
        f, g = self._target_and_gradients()
      elif(m.is_terminated()):
        break
      if(m.n_iteration() >= self.max_iterations): break
      if(self.n_calls >= max_calls): break
    self._unpack()
    self.n_iterations = m.n_iteration()
    self.task = m.task()
    return self

class tgc(object):
  def __init__(self,
               f_obs,
               f_calc,
               f_masks,
               ss,
               k_sols=None,
               b_sols=None,
               ps=None,
               u_star=[0,0,0,0,0,0], b_max=300, b_min=0, k_min=0.001, k_max=50,
               symmetry_constraints=True, min_iterations=0,
               max_iterations=25):
    if(ps is not None): assert [k_sols, b_sols].count(None) == 2
    else:               assert [k_sols, b_sols].count(None) == 0
    adopt_init_args(self, locals())
    self.kbu = mmtbx.f_model.manager_kbu(
      f_obs   = self.f_obs,
      f_calc  = self.f_calc,
      f_masks = self.f_masks,
      ss      = self.ss,
      k_sols  = self.k_sols,
      b_sols  = self.b_sols,
      u_star  = self.u_star)
    self.t_g_c = None
    self.use_scale=None
    self.refine_kb=False
    self.refine_k=False
    self.refine_b=False
    self.refine_u=False
    self.refine_p=False
    self.space_group = self.f_obs.space_group()
    self.adp_constraints = self.space_group.adp_constraints()

  def set_refine_kb(self):
    self.refine_kb=True
    self.refine_u=False
    self.refine_k=False
    self.refine_b=False
    self.refine_p=False

  def set_refine_k(self):
    self.refine_k=True
    self.refine_b=False
    self.refine_kb=False
    self.refine_u=False
    self.refine_p=False

  def set_refine_b(self):
    self.refine_k=False
    self.refine_b=True
    self.refine_kb=False
    self.refine_u=False
    self.refine_p=False

  def set_refine_u(self):
    self.refine_k=False
    self.refine_b=False
    self.refine_kb=False
    self.refine_u=True
    self.refine_p=False
    if(self.symmetry_constraints):
      u_star = self.space_group.average_u_star(u_star = self.kbu.u_star())
      self.kbu.update(u_star = u_star)
      assert self.adp_constraints.n_independent_params() <= 6

  def set_refine_p(self):
    self.refine_kb=False
    self.refine_u=False
    self.refine_k=False
    self.refine_b=False
    self.refine_p=True

  def set_use_scale(self, value):
    assert value in [True, False]
    self.use_scale=value

  def normalize(self, parameters, p_min, p_max):
    result = flex.double()
    for p in parameters:
      if(p < p_min): p = p_min
      if(p > p_max): p = p_max
      result.append(p)
    return result

  def x(self):
    if(self.refine_k):
      return self.normalize(self.kbu.k_sols(), self.k_min, self.k_max)
    if(self.refine_b):
      return self.normalize(self.kbu.b_sols(), self.b_min, self.b_max)
    if(self.refine_kb):
      x =      self.normalize(self.kbu.k_sols(), self.k_min, self.k_max)
      x.extend(self.normalize(self.kbu.b_sols(), self.b_min, self.b_max))
      return x
    if(self.refine_u):
      if(not self.symmetry_constraints):
        return flex.double(self.kbu.u_star())
      return flex.double(
        self.adp_constraints.independent_params(self.kbu.u_star()))

  def target(self):
    return self.t_g_c.target()

  def gradients(self):
    if(self.refine_k): return self.t_g_c.grad_k_sols()
    if(self.refine_b): return self.t_g_c.grad_b_sols()
    if(self.refine_kb):
      g=self.t_g_c.grad_k_sols()
      g.extend(self.t_g_c.grad_b_sols())
      return g
    if(self.refine_u):
      if(not self.symmetry_constraints):
        return flex.double(self.t_g_c.grad_u_star())
      return flex.double(
        self.adp_constraints.independent_gradients(all_gradients=self.t_g_c.grad_u_star()))

  def curvatures(self):
    # XXX No curvatures for u_star !
    if(self.refine_k): return self.t_g_c.curv_k_sols()
    if(self.refine_b): return self.t_g_c.curv_b_sols()
    if(self.refine_kb):
      d = self.t_g_c.curv_k_sols()
      d.extend(self.t_g_c.curv_b_sols())
      return d

  def update(self, x):
    if(self.refine_k): self.kbu.update(k_sols=x)
    if(self.refine_b): self.kbu.update(b_sols=x)
    if(self.refine_kb):
      self.kbu.update(
        k_sols=x[:len(x)//2],
        b_sols=x[len(x)//2:])
    if(self.refine_u):
      if(not self.symmetry_constraints):
        u_star = x
      else:
        u_star = self.adp_constraints.all_params(list(x))
      self.kbu.update(u_star = list(u_star))
    if(self.use_scale):
      sc = bulk_solvent.scale(self.f_obs.data(), self.kbu.data.f_model)
    else:
      sc = 1.0
    self.t_g_c = bulk_solvent.ls_kbp_sol_u_star(
      f_model     = self.kbu.data,
      f_obs       = self.f_obs.data(),
      scale       = sc,
      kb_sol_grad = self.refine_k or self.refine_b or self.refine_kb,
      p_sol_grad  = False,
      u_star_grad = self.refine_u,
      kb_sol_curv = self.refine_k or self.refine_b or self.refine_kb,
      p_sol_curv  = False)

  def _minimizer(self):
    return minimizer(tgc = self, min_iterations = self.min_iterations,
      max_iterations = self.max_iterations)

  def minimize_joint(self, k_sol_min, k_sol_max, b_sol_min, b_sol_max,
                     refine_kb=True, refine_u=True, max_iterations=None):
    """
    Refine k_sols, b_sols and u_star together within bounds (see
    joint_minimizer). self.kbu holds the result.

    Returns
    -------
    joint_minimizer
    """
    if(max_iterations is None): max_iterations = self.max_iterations
    return joint_minimizer(
      tgc            = self,
      k_sol_min      = k_sol_min,
      k_sol_max      = k_sol_max,
      b_sol_min      = b_sol_min,
      b_sol_max      = b_sol_max,
      refine_kb      = refine_kb,
      refine_u       = refine_u,
      max_iterations = max_iterations).run()

  def minimize_kb_once(self, use_curvatures):
    self.set_refine_kb()
    return self._minimizer().run(use_curvatures=use_curvatures)

  def minimize_u_once(self):
    self.set_refine_u()
    return self._minimizer().run(use_curvatures=False)

  def minimize_u(self, n_cycles=5):
    #print "minimize_u, r:", self.kbu.r_factor()
    for it in range(n_cycles):
      start_r = self.kbu.r_factor()
      save_b_cart = self.kbu.b_cart()
      self.set_refine_u()
      self.set_use_scale(value = True)
      self._minimizer().run(use_curvatures=False)
      #print "  minimize_u, r:", self.kbu.r_factor()
      r = self.kbu.r_factor()
      bc = list(flex.abs(flex.double(self.kbu.b_cart())))
      if(r>start_r and r>1.e-2 and max(bc)>100):
        self.kbu.update(b_cart = save_b_cart)
        break

  def minimize_kbu(self, n_cycles=10):
    #print "minimize_kbu start r:", self.kbu.r_factor()
    for use_curvatures in [False, True]*n_cycles:
      #print "  minimize_kbu r:", self.kbu.r_factor()
      start_r = self.kbu.r_factor()
      save_k_sols = self.kbu.k_sols()
      save_b_sols = self.kbu.b_sols()
      save_b_cart = self.kbu.b_cart()
      #self.set_use_scale(value = random.choice([True, False]))
      self.set_use_scale(value = True)
      m = self.minimize_kb_once(use_curvatures=use_curvatures)
      r = self.kbu.r_factor()
      if(r>start_r and r>1.e-2 and (flex.min(self.kbu.k_sols())<0 or
         flex.max(self.kbu.k_sols())>1 or flex.min(self.kbu.b_sols())<0 or
         flex.max(self.kbu.k_sols())>100.)):
        self.kbu.update(k_sols = save_k_sols, b_sols = save_b_sols)
#      assert m.minimizer.n_calls == m.minimizer.nfun()
      m = self.minimize_u_once()
 #     assert m.minimizer.n_calls == m.minimizer.nfun()
      r = self.kbu.r_factor()
      bc = list(flex.abs(flex.double(self.kbu.b_cart())))
      if(r>start_r and r>1.e-2 and max(bc)>100):
        self.kbu.update(b_cart = save_b_cart)
        break

  def show_k_sols(self):
    print("k_sols:", [round(k,3) for k in self.kbu.k_sols()], self.kbu.r_factor())

  def show_kbu(self):
    print("k_sols:", [round(k,3) for k in self.kbu.k_sols()])
    print("b_sols:", [round(b,3) for b in self.kbu.b_sols()])
    print("b_cart:", [round(b,3) for b in self.kbu.b_cart()])
