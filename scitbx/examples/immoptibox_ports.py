"""\
http://www2.imm.dtu.dk/~hbn/immoptibox/
"""

from __future__ import absolute_import, division, print_function
import scitbx.math
import scitbx.linalg
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import math
import sys

import scitbx.minpack
import scitbx.lbfgs
import scitbx.lbfgsb
from six.moves import range
from six.moves import zip

try:
  import knitro_adaptbx
except ImportError:
  knitro_adaptbx = None

floating_point_epsilon_double = scitbx.math.floating_point_epsilon_double_get()

def cholesky_decomposition(a, relative_eps=1.e-15):
  assert a.is_square_matrix()
  n = a.focus()[0]
  eps = relative_eps * flex.max(flex.abs(a))
  c = flex.double(a.accessor(), 0)
  for k in range(n):
    sum = 0
    for j in range(k):
      sum += c[(k,j)]**2
    d = a[(k,k)] - sum
    if (d <= eps): return None
    c[(k,k)] = math.sqrt(d)
    for i in range(k+1,n):
      sum = 0
      for j in range(k):
        sum += c[(i,j)] * c[(k,j)]
      c[(i,k)] = (a[(i,k)] - sum) / c[(k,k)]
  return c

def cholesky_solve(c, b):
  assert c.is_square_matrix()
  assert b.is_trivial_1d()
  assert b.size() == c.focus()[0]
  n = c.focus()[0]
  z = flex.double(n, 0)
  for k in range(n):
    sum = 0
    for j in range(k):
      sum += c[(k,j)] * z[j]
    z[k] = (b[k] - sum) / c[(k,k)]
  x = flex.double(n, 0)
  for k in range(n-1,-1,-1):
    sum = 0
    for j in range(k+1,n):
      sum += c[(j,k)] * x[j]
    x[k] = (z[k] - sum) / c[(k,k)]
  return x

class levenberg_marquardt:

  def __init__(self,
        function,
        x0,
        tau=1.e-3,
        eps_1=1.e-16,
        eps_2=1.e-16,
        mu_min=1.e-300,
        k_max=1000):
    self.function = function
    nu = 2
    x = x0
    f_x = function.f(x)
    number_of_function_evaluations = 1
    self.f_x0 = f_x
    j = function.jacobian(x)
    number_of_jacobian_evaluations = 1
    number_of_cholesky_decompositions = 0
    j_t = j.matrix_transpose()
    a = j_t.matrix_multiply(j)
    g = j_t.matrix_multiply(f_x)
    found = flex.max(flex.abs(g)) <= eps_1
    mu = tau * flex.max(a.matrix_diagonal())
    k = 0
    while (not found and k < k_max):
      k += 1
      a_plus_mu = a.deep_copy()
      if (mu > mu_min):
        a_plus_mu.matrix_diagonal_add_in_place(value=mu)
      u = a_plus_mu.matrix_symmetric_as_packed_u()
      gmw = scitbx.linalg.gill_murray_wright_cholesky_decomposition_in_place(u)
      number_of_cholesky_decompositions += 1
      h_lm = gmw.solve(b=-g)
      if (h_lm.norm() <= eps_2 * (x.norm() + eps_2)):
        found = True
      else:
        x_new = x + h_lm
        f_x_new = function.f(x_new)
        number_of_function_evaluations += 1
        rho_denom = 0.5 * h_lm.dot(h_lm * mu - g)
        assert rho_denom != 0
        rho = 0.5*(f_x.norm()**2 - f_x_new.norm()**2) / rho_denom
        if (rho > 0):
          x = x_new
          f_x = f_x_new
          j = function.jacobian(x)
          number_of_jacobian_evaluations += 1
          j_t = j.matrix_transpose()
          a = j_t.matrix_multiply(j)
          g = j_t.matrix_multiply(f_x)
          found = flex.max(flex.abs(g)) <= eps_1
          if (mu > mu_min):
            mu *= max(1/3., 1-(2*rho-1)**3)
          else:
            mu = mu_min
          nu = 2
        else:
          mu *= nu
          nu *= 2
    self.x_star = x
    self.f_x_star = f_x
    self.number_of_iterations = k
    self.number_of_function_evaluations = number_of_function_evaluations
    self.number_of_jacobian_evaluations = number_of_jacobian_evaluations
    self.number_of_cholesky_decompositions = number_of_cholesky_decompositions

  def show_statistics(self):
    print("levenberg_marquardt results:")
    print("  function:", self.function.label())
    print("  x_star:", list(self.x_star))
    print("  0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2)
    print("  0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2)
    print("  number_of_iterations:", self.number_of_iterations)
    print("  number_of_function_evaluations:", \
        self.number_of_function_evaluations)
    print("  number_of_jacobian_evaluations:", \
        self.number_of_jacobian_evaluations)
    print("  number_of_cholesky_decompositions:", \
      self.number_of_cholesky_decompositions)

class damped_newton:

  def __init__(self,
        function,
        x0,
        mu0=None,
        tau=1.e-3,
        eps_1=1.e-16,
        eps_2=1.e-16,
        k_max=1000):
    self.function = function
    x = x0
    f_x = function.f(x=x)
    number_of_function_evaluations = 1
    self.f_x0 = f_x
    fp = function.gradients(x=x, f_x=f_x)
    number_of_gradient_evaluations = 1
    number_of_hessian_evaluations = 0
    number_of_cholesky_decompositions = 0
    mu = mu0
    nu = 2
    k = 0
    while (k < k_max):
      if (flex.max(flex.abs(fp)) <= eps_1):
        break
      fdp = function.hessian(x=x)
      number_of_hessian_evaluations += 1
      if (mu is None):
        mu = tau * flex.max(flex.abs(fdp))
        if (mu == 0):
          mu = fp.norm()/max(x.norm(), floating_point_epsilon_double**0.5)
      while True:
        fdp_plus_mu = fdp.deep_copy()
        fdp_plus_mu.matrix_diagonal_add_in_place(value=mu)
        fdp_plus_mu_cholesky = cholesky_decomposition(a=fdp_plus_mu)
        number_of_cholesky_decompositions += 1
        if (fdp_plus_mu_cholesky is not None):
          break
        mu *= 10
      h_dn = cholesky_solve(c=fdp_plus_mu_cholesky, b=-fp)
      if (h_dn.norm() <= eps_2*(eps_2 + x.norm())):
        break
      x_new = x + h_dn
      f_x_new = function.f(x=x_new)
      number_of_function_evaluations += 1
      fp_new = function.gradients(x=x_new, f_x=f_x_new)
      number_of_gradient_evaluations += 1
      f = flex.sum_sq(f_x)
      fn = flex.sum_sq(f_x_new)
      df = f - fn
      accept = 0
      if (df > 0):
        accept = 1
        dl = 0.5 * h_dn.dot(h_dn * mu - fp)
      elif (fn <= f + abs(f)*(1+100*floating_point_epsilon_double)):
        df = (fp + fp_new).dot(fp - fp_new)
        if (df > 0):
          accept = 2
      if (accept != 0):
        if (accept == 1 and dl > 0):
          mu *= max(1/3, 1 - (2*df/dl - 1)**3)
          nu = 2
        else:
          mu *= nu
          nu *= 2
        x = x_new
        f_x = f_x_new
        fp = fp_new
        fdp = None
        k += 1
      else:
        mu *= nu
        nu *= 2
    self.x_star = x
    self.f_x_star = f_x
    self.number_of_iterations = k
    self.number_of_function_evaluations = number_of_function_evaluations
    self.number_of_gradient_evaluations = number_of_gradient_evaluations
    self.number_of_hessian_evaluations = number_of_hessian_evaluations
    self.number_of_cholesky_decompositions = number_of_cholesky_decompositions

  def show_statistics(self):
    print("damped_newton results:")
    print("  function:", self.function.label())
    print("  x_star:", list(self.x_star))
    print("  0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2)
    print("  0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2)
    print("  number_of_iterations:", self.number_of_iterations)
    print("  number_of_function_evaluations:", \
        self.number_of_function_evaluations)
    print("  number_of_gradient_evaluations:", \
        self.number_of_gradient_evaluations)
    print("  number_of_hessian_evaluations:", \
        self.number_of_hessian_evaluations)
    print("  number_of_cholesky_decompositions:", \
        self.number_of_cholesky_decompositions)

class minpack_levenberg_marquardt_adaptor:

  def __init__(self, function, x0, max_iterations=1000):
    self.function = function
    self.f_x0 = function.f(x=x0)
    x = x0.deep_copy()
    self.number_of_iterations = 0
    self.number_of_function_evaluations = 0
    self.number_of_jacobian_evaluations = 0
    minimizer = scitbx.minpack.levenberg_marquardt(
      m=function.m,
      x=x,
      ftol=1.e-16,
      xtol=1.e-16,
      gtol=0,
      maxfev=0,
      factor=1.0e2,
      call_back_after_iteration=True)
    while (not minimizer.has_terminated()):
      if (minimizer.requests_fvec()):
        self.number_of_function_evaluations += 1
        minimizer.process_fvec(
          fvec=function.f(x=x))
      elif (minimizer.requests_fjac()):
        self.number_of_jacobian_evaluations += 1
        minimizer.process_fjac(
          fjac=function.jacobian(x=x).matrix_transpose().as_1d())
      elif (minimizer.calls_back_after_iteration()):
        self.number_of_iterations += 1
        if (self.number_of_iterations >= max_iterations):
          break
        minimizer.continue_after_call_back_after_iteration()
    self.x_star = x
    self.f_x_star = function.f(x=self.x_star)

  def show_statistics(self):
    print("minpack.levenberg_marquardt results:")
    print("  function:", self.function.label())
    print("  x_star:", list(self.x_star))
    print("  0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2)
    print("  0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2)
    print("  number_of_iterations:", self.number_of_iterations)
    print("  number_of_function_evaluations:", \
        self.number_of_function_evaluations)
    print("  number_of_jacobian_evaluations:", \
        self.number_of_jacobian_evaluations)

class lbfgs_adaptor:

  def __init__(self, function, x0, max_iterations=1000):
    self.function = function
    self.f_x0 = function.f(x=x0)
    self.x = x0.deep_copy()
    self.number_of_function_evaluations = 0
    self.number_of_gradient_evaluations = 0
    minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-16,
        max_iterations=max_iterations))
    self.number_of_iterations = minimizer.iter()
    self.x_star = self.x
    self.f_x_star = function.f(x=self.x_star)
    del self.x

  def compute_functional_and_gradients(self):
    f = 0.5*self.function.f(x=self.x).norm()**2
    self.number_of_function_evaluations += 1
    g = self.function.gradients(x=self.x)
    self.number_of_gradient_evaluations += 1
    return f, g

  def show_statistics(self):
    print("lbfgs results:")
    print("  function:", self.function.label())
    print("  x_star:", list(self.x_star))
    print("  0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2)
    print("  0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2)
    print("  number_of_iterations:", self.number_of_iterations)
    print("  number_of_function_evaluations:", \
        self.number_of_function_evaluations)
    print("  number_of_gradient_evaluations:", \
        self.number_of_gradient_evaluations)

class lbfgsb_adaptor:

  def __init__(self, function, x0, max_iterations=1000):
    self.function = function
    self.f_x0 = function.f(x=x0)
    x = x0.deep_copy()
    self.number_of_function_evaluations = 0
    self.number_of_gradient_evaluations = 0
    minimizer = scitbx.lbfgsb.minimizer(n=x.size(), factr=1.e-16, pgtol=1.e-16)
    f = 0
    g = flex.double(x.size(), 0)
    while True:
      if (minimizer.process(x, f, g, False)):
        f = 0.5*self.function.f(x=x).norm()**2
        self.number_of_function_evaluations += 1
        g = self.function.gradients(x=x)
        self.number_of_gradient_evaluations += 1
      elif (minimizer.is_terminated()):
        break
      elif (minimizer.n_iteration() >= max_iterations):
        break
    self.number_of_iterations = minimizer.n_iteration()
    self.x_star = x
    self.f_x_star = function.f(x=self.x_star)

  def show_statistics(self):
    print("lbfgsb results:")
    print("  function:", self.function.label())
    print("  x_star:", list(self.x_star))
    print("  0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2)
    print("  0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2)
    print("  number_of_iterations:", self.number_of_iterations)
    print("  number_of_function_evaluations:", \
        self.number_of_function_evaluations)
    print("  number_of_gradient_evaluations:", \
        self.number_of_gradient_evaluations)

class MeyerFunctionError(RuntimeError): pass

class test_function:

  def __init__(self, m, n, check_with_finite_differences=True, verbose=1):
    self.m = m
    self.n = n
    self.verbose = verbose
    self.check_with_finite_differences = check_with_finite_differences
    self.check_gradients_tolerance = 1.e-6
    self.check_hessian_tolerance = 1.e-6
    self.initialization()
    assert self.m >= self.n
    if (self.x_star is not None):
      assert approx_equal(
        0.5*self.f(x=self.x_star).norm()**2, self.capital_f_x_star)
    if (1):
      self.exercise_levenberg_marquardt()
    if (1):
      self.exercise_minpack_levenberg_marquardt()
    if (1):
      self.exercise_damped_newton()
    if (1):
      self.exercise_scitbx_minimizers_damped_newton()
    if (1):
      self.exercise_scitbx_minimizers_newton_more_thuente_1994()
    if (1):
      self.exercise_lbfgs()
    if (1):
      try:
        self.exercise_lbfgsb()
      except MeyerFunctionError:
        print("Skipping: exercise_lbfgsb() with", self.__class__.__name__)
    if (1 and knitro_adaptbx is not None):
      self.exercise_knitro_adaptbx()

  def label(self):
    return self.__class__.__name__

  def functional(self, x=None, f_x=None):
    if (f_x is None): f_x = self.f(x=x)
    return 0.5*flex.sum_sq(f_x)

  def jacobian_finite(self, x, relative_eps=1.e-8):
    x0 = x
    result = flex.double()
    for i in range(self.n):
      eps = max(1, abs(x0[i])) * relative_eps
      fs = []
      for signed_eps in [eps, -eps]:
        x = x0.deep_copy()
        x[i] += signed_eps
        fs.append(self.f(x=x))
      result.extend((fs[0]-fs[1])/(2*eps))
    result.resize(flex.grid(self.n, self.m))
    return result.matrix_transpose()

  def jacobian(self, x):
    analytical = self.jacobian_analytical(x=x)
    if (self.check_with_finite_differences):
      finite = self.jacobian_finite(x=x)
      scale = max(1, flex.max(flex.abs(analytical)))
      assert approx_equal(analytical/scale, finite/scale, 1.e-5)
    return analytical

  def gradients_finite(self, x, relative_eps=1.e-7):
    x0 = x
    result = flex.double()
    for i in range(self.n):
      eps = max(1, abs(x0[i])) * relative_eps
      fs = []
      for signed_eps in [eps, -eps]:
        x = x0.deep_copy()
        x[i] += signed_eps
        fs.append(self.functional(x=x))
      result.append((fs[0]-fs[1])/(2*eps))
    return result

  def gradients_analytical(self, x, f_x=None):
    if (f_x is None): f_x = self.f(x=x)
    return self.jacobian_analytical(x=x).matrix_transpose() \
      .matrix_multiply(f_x)

  def gradients(self,x, f_x=None):
    analytical = self.gradients_analytical(x=x, f_x=f_x)
    if (self.check_with_finite_differences):
      finite = self.gradients_finite(x=x)
      scale = max(1, flex.max(flex.abs(analytical)))
      assert approx_equal(analytical/scale, finite/scale,
        self.check_gradients_tolerance)
    return analytical

  def hessian_finite(self, x, relative_eps=1.e-8):
    x0 = x
    result = flex.double()
    for i in range(self.n):
      eps = max(1, abs(x0[i])) * relative_eps
      gs = []
      for signed_eps in [eps, -eps]:
        x = x0.deep_copy()
        x[i] += signed_eps
        gs.append(self.gradients_analytical(x=x))
      result.extend((gs[0]-gs[1])/(2*eps))
    result.resize(flex.grid(self.n, self.n))
    return result

  def hessian(self, x):
    analytical = self.hessian_analytical(x=x)
    if (self.check_with_finite_differences):
      finite = self.hessian_finite(x=x)
      scale = max(1, flex.max(flex.abs(analytical)))
      assert approx_equal(analytical/scale, finite/scale,
        self.check_hessian_tolerance)
    return analytical

  def check_minimized_x_star(self, x_star):
    if (self.x_star is not None):
      assert approx_equal(x_star, self.x_star)

  def check_minimized_capital_f_x_star(self, f_x_star):
    if (self.capital_f_x_star is not None):
      assert approx_equal(0.5*f_x_star.norm()**2, self.capital_f_x_star)

  def check_minimized(self, minimized):
    self.check_minimized_x_star(x_star=minimized.x_star)
    self.check_minimized_capital_f_x_star(f_x_star=minimized.f_x_star)

  def exercise_levenberg_marquardt(self):
    minimized = levenberg_marquardt(function=self, x0=self.x0, tau=self.tau0)
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()

  def exercise_minpack_levenberg_marquardt(self):
    minimized = minpack_levenberg_marquardt_adaptor(function=self, x0=self.x0)
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()

  def exercise_damped_newton(self):
    minimized = damped_newton(function=self, x0=self.x0, tau=self.tau0)
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()

  def exercise_scitbx_minimizers_damped_newton(self):
    import scitbx.minimizers
    minimized = scitbx.minimizers.damped_newton(
      function=self, x0=self.x0, tau=self.tau0)
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()

  def exercise_scitbx_minimizers_newton_more_thuente_1994(self):
    import scitbx.minimizers
    minimized = scitbx.minimizers.newton_more_thuente_1994(
      function=self, x0=self.x0)
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()

  def exercise_lbfgs(self):
    minimized = lbfgs_adaptor(function=self, x0=self.x0)
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()

  def exercise_lbfgsb(self):
    minimized = lbfgsb_adaptor(function=self, x0=self.x0)
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()

  def exercise_knitro_adaptbx(self):
    minimized = knitro_adaptbx.solve(function=self, x0=self.x0)
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()
    minimized = knitro_adaptbx.solve(function=self, x0=self.x0,hessopt="bfgs")
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()
    minimized = knitro_adaptbx.solve(function=self, x0=self.x0,hessopt="lbfgs")
    if (self.verbose): minimized.show_statistics()
    self.check_minimized(minimized=minimized)
    if (self.verbose): print()

class linear_function_full_rank(test_function):

  def initialization(self):
    self.a = flex.double(flex.grid(self.n, self.n), -2./self.m)
    self.a.matrix_diagonal_add_in_place(value=1)
    self.a.resize(flex.grid(self.m,self.n), -2./self.m)
    self.x0 = flex.double(self.n, 1)
    self.tau0 = 1.e-8
    self.delta0 = 10
    self.x_star = flex.double(self.n, -1)
    self.capital_f_x_star = 0.5 * (self.m - self.n)

  def label(self):
    return "%s(m=%d, n=%d)" % (self.__class__.__name__, self.m, self.n)

  def f(self, x):
    return self.a.matrix_multiply(x) - flex.double(self.m, 1)

  def jacobian_analytical(self, x):
    return self.a

  def hessian_analytical(self, x):
    return self.a.matrix_transpose().matrix_multiply(self.a)

class linear_function_rank_1(linear_function_full_rank):

  def initialization(self):
    m = self.m
    n = self.n
    self.a = flex.double(range(1,m+1)).matrix_outer_product(
             flex.double(range(1,n+1)))
    self.x0 = flex.double(n, 1)
    self.tau0 = 1.e-8
    self.delta0 = 10
    self.x_star = None
    self.capital_f_x_star = (m*(m-1))/(4*(2*m+1))

  def check_minimized_x_star(self, x_star):
    assert approx_equal(
      flex.double(range(1,self.n+1)).dot(x_star),
      3/(2*self.m+1))

class linear_function_rank_1_with_zero_columns_and_rows(
        linear_function_full_rank):

  def initialization(self):
    m = self.m
    n = self.n
    self.a = flex.double([0]+list(range(1,m-2+1))+[0]).matrix_outer_product(
             flex.double([0]+list(range(2,n-1+1))+[0]))
    self.x0 = flex.double(n, 1)
    self.tau0 = 1.e-8
    self.delta0 = 10
    self.x_star = None
    self.capital_f_x_star = (m**2+3*m-6)/(4*(2*m-3))

  def check_minimized_x_star(self, x_star):
    assert approx_equal(
      flex.double([0]+list(range(2,self.n-1+1))+[0]).dot(x_star),
      3/(2*self.m-3))

class rosenbrock_function(test_function):

  def initialization(self):
    assert self.m == 2
    assert self.n == 2
    self.x0 = flex.double([-1.2, 1])
    self.tau0 = 1
    self.delta0 = 1
    self.x_star = flex.double([1,1])
    self.capital_f_x_star = 0

  def f(self, x):
    return flex.double([10*(x[1]-x[0]**2), 1-x[0]])

  def jacobian_analytical(self, x):
    return flex.double([[-20*x[0], 10], [-1, 0]])

  def hessian_analytical(self, x):
    f = self.f(x=x)
    j = self.jacobian(x=x)
    return j.matrix_transpose().matrix_multiply(j) \
         + flex.double([[-20*f[0],0],[0,0]])

class helical_valley_function(test_function):

  def initialization(self):
    assert self.m == 3
    assert self.n == 3
    self.x0 = flex.double([-1, 0, 0])
    self.tau0 = 1
    self.delta0 = 1
    self.x_star = flex.double([1,0,0])
    self.capital_f_x_star = 0

  def f(self, x):
    x1,x2,x3 = x
    assert x1 != 0
    t = math.atan(x2/x1)/(2*math.pi)
    if (x1 < 0): t += 0.5
    nx = x[:2].norm()
    return flex.double([10*(x3 - 10*t), 10*(nx - 1), x3])

  def jacobian_analytical(self, x):
    x1,x2,x3 = x
    nx = x[:2].norm()
    nx2 = nx*nx
    k1 = 50/math.pi/nx2
    k2 = 10/nx
    return flex.double([
      [k1*x2, -k1*x1, 10],
      [k2*x1, k2*x2, 0],
      [0, 0, 1]])

  def hessian_analytical(self, x):
    x1,x2,x3 = x
    nx = x[:2].norm()
    nx2 = nx*nx
    k1 = 50/math.pi/nx2
    k2 = 10/nx
    f = self.f(x=x)
    j = self.jacobian_analytical(x=x)
    result = j.matrix_transpose().matrix_multiply(j)
    q1 = x1**2
    q2 = x2**2
    p = x1*x2
    terms = f[0]*k1/nx2*flex.double([[-2*p,q1-q2],[q1-q2,2*p]]) \
          + f[1]*k2/nx2*flex.double([[q2,-p],[-p,q1]])
    for i in [0,1]:
      for j in [0,1]:
        result[i*3+j] += terms[i*2+j]
    return result

class powell_singular_function(test_function):

  def initialization(self):
    assert self.m == 4
    assert self.n == 4
    self.x0 = flex.double([3, -1, 0, 1])
    self.tau0 = 1.e-8
    self.delta0 = 1
    self.x_star = flex.double([0,0,0,0])
    self.capital_f_x_star = 0

  def check_minimized_x_star(self, x_star):
    assert approx_equal(x_star, self.x_star, 1.e-5)

  def f(self, x):
    x1,x2,x3,x4 = x
    return flex.double([
      x1+10*x2,
      5**0.5*(x3-x4),
      (x2-2*x3)**2,
      10**0.5*(x1-x4)**2])

  def jacobian_analytical(self, x):
    x1,x2,x3,x4 = x
    d3 = x2 - 2*x3
    d4 = x1 - x4
    s5 = 5**0.5
    s10 = 10**0.5
    return flex.double([
      [1, 10, 0, 0],
      [0, 0, s5, -s5],
      [0, 2*d3, -4*d3, 0],
      [2*s10*d4, 0, 0, -2*s10*d4]])

  def hessian_analytical(self, x):
    s10 = 10**0.5
    f1,f2,f3,f4 = self.f(x=x)
    j = self.jacobian_analytical(x=x)
    result = j.matrix_transpose().matrix_multiply(j)
    result += flex.double([
     [f4*2*s10,0,0,-f4*2*s10],
     [0,f3*2,-f3*4,0],
     [0,-f3*4,f3*8,0],
     [-f4*2*s10,0,0,f4*2*s10]])
    return result

class freudenstein_and_roth_function(test_function):

  def initialization(self):
    assert self.m == 2
    assert self.n == 2
    self.x0 = flex.double([0.5, -2])
    self.tau0 = 1
    self.delta0 = 1
    self.x_star = (flex.double([53,2]) - 22**0.5*flex.double([4,1]))/3
    self.capital_f_x_star = 24.4921268396

  def f(self, x):
    x1,x2 = x
    return flex.double([
      x1-x2*(2-x2*(5-x2))-13,
      x1-x2*(14-x2*(1+x2))-29])

  def jacobian_analytical(self, x):
    x1,x2 = x
    return flex.double([
      [1, (-2 + x2*(10 - 3*x2))],
      [1, (-14 + x2*(2 + 3*x2))]])

  def hessian_analytical(self, x):
    x1,x2 = x
    f1,f2 = self.f(x=x)
    j = self.jacobian_analytical(x=x)
    result = j.matrix_transpose().matrix_multiply(j)
    result[3] += f1*(10-6*x2) + f2*(2+6*x2)
    return result

class bard_function(test_function):

  ys = [0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58,
        0.73, 0.96, 1.34, 2.10, 4.39]

  def initialization(self):
    assert self.m == 15
    assert self.n == 3
    self.x0 = flex.double([1, 1, 1])
    self.tau0 = 1.e-8
    self.delta0 = 1
    self.x_star = flex.double([0.082411, 1.133036, 2.343695])
    self.capital_f_x_star = 4.10744e-3

  def f(self, x):
    x1,x2,x3 = x
    result = flex.double()
    for i,yi in zip(range(1,15+1),bard_function.ys):
      ui = i
      vi = 16-i
      wi = min(ui, vi)
      denominator = x2*vi + x3*wi
      assert denominator != 0
      result.append(yi - (x1 + ui / denominator))
    return result

  def jacobian_analytical(self, x):
    x1,x2,x3 = x
    result = flex.double()
    for i in range(1,15+1):
      ui = i
      vi = 16-i
      wi = min(ui, vi)
      denominator = (x2*vi + x3*wi)**2
      assert denominator != 0
      result.extend(flex.double([-1, ui*vi/denominator, ui*wi/denominator]))
    result.resize(flex.grid(self.m,self.n))
    return result

  def hessian_analytical(self, x):
    x1,x2,x3 = x
    j = self.jacobian_analytical(x=x)
    result = j.matrix_transpose().matrix_multiply(j)
    for i,fi in zip(range(1,15+1), self.f(x=x)):
      ui = i
      vi = 16-i
      wi = min(ui, vi)
      denominator = (x2*vi + x3*wi)**3
      assert denominator != 0
      term = fi*2*ui/denominator
      result[(1,1)] -= term*vi**2
      result[(1,2)] -= term*vi*wi
      result[(2,2)] -= term*wi**2
    result[(2,1)] = result[(1,2)]
    return result

class kowalik_and_osborne_function(test_function):

  ys = [0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627,
        0.0456, 0.0342, 0.0323, 0.0235, 0.0246]
  us = [4.0000, 2.0000, 1.0000, 0.5000, 0.2500, 0.1670,
        0.1250, 0.1000, 0.0833, 0.0714, 0.0625]

  def initialization(self):
    assert self.m == 11
    assert self.n == 4
    self.x0 = flex.double([0.25, 0.39, 0.415, 0.39])
    self.tau0 = 1
    self.delta0 = 0.1
    self.x_star = flex.double([
      0.1928069346, 0.1912823287, 0.1230565069, 0.1360623307])
    self.capital_f_x_star = 1.53753e-4

  def f(self, x):
    x1,x2,x3,x4 = x
    result = flex.double()
    for yi,ui in zip(kowalik_and_osborne_function.ys,
                     kowalik_and_osborne_function.us):
      denominator = ui*(ui+x3)+x4
      assert denominator != 0
      result.append(yi-x1*ui*(ui+x2)/denominator)
    return result

  def jacobian_analytical(self, x):
    x1,x2,x3,x4 = x
    result = flex.double()
    for ui in kowalik_and_osborne_function.us:
      denominator = ui*(ui+x3)+x4
      assert denominator != 0
      denominator_sq = denominator**2
      assert denominator_sq != 0
      result.extend(flex.double([
        -ui*(ui+x2)/denominator,
        -ui*x1/denominator,
        ui**2*x1*(ui+x2)/denominator_sq,
        ui*x1*(ui+x2)/denominator_sq]))
    result.resize(flex.grid(self.m, self.n))
    return result

  def hessian_analytical(self, x):
    x1,x2,x3,x4 = x
    j = self.jacobian_analytical(x=x)
    result = j.matrix_transpose().matrix_multiply(j)
    for ui,fi in zip(kowalik_and_osborne_function.us, self.f(x=x)):
      denominator = ui*(ui+x3)+x4
      assert denominator != 0
      denominator_sq = denominator**2
      assert denominator_sq != 0
      denominator_cu = denominator**3
      assert denominator_cu != 0
      result[(0,0)] -= 0
      result[(0,1)] -= fi*ui/denominator
      result[(0,2)] -= -fi*(ui**2*(ui+x2))/denominator_sq
      result[(0,3)] -= -fi*(ui*(ui+x2))/denominator_sq
      result[(1,1)] -= 0
      result[(1,2)] -= -fi*ui**2*x1/denominator_sq
      result[(1,3)] -= -fi*ui*x1/denominator_sq
      result[(2,2)] -= fi*2*ui**3*x1*(ui+x2)/denominator_cu
      result[(2,3)] -= fi*2*ui**2*x1*(ui+x2)/denominator_cu
      result[(3,3)] -= fi*2*ui*x1*(ui+x2)/denominator_cu
    for i in range(0,4):
      for j in range(i+1,4):
        result[(j,i)] = result[(i,j)]
    return result

class meyer_function(test_function):

  ys = [34780, 28610, 23650, 19630, 16370, 13720, 11540, 9744,
        8261, 7030, 6005, 5147, 4427, 3820, 3307, 2872]
  ts = [45+5*i for i in range(1,16+1)]

  def initialization(self):
    assert self.m == 16
    assert self.n == 3
    self.x0 = flex.double([0.02, 4000, 250])
    self.tau0 = 1.e-6
    self.delta0 = 100
    self.x_star = flex.double([
      5.60963646990603e-3, 6.181346346e3, 3.452236346e2])
    self.capital_f_x_star = 43.9729275853
    self.check_gradients_tolerance = 1.e-3
    self.check_hessian_tolerance = 1.e-2

  def check_minimized_x_star(self, x_star):
    if (self.verbose):
      print("  expected x_star: %.6g %.6g %.6g" % tuple(self.x_star))
      print("    actual x_star: %.6g %.6g %.6g" % tuple(x_star))

  def check_minimized_capital_f_x_star(self, f_x_star):
    if (self.verbose):
      print("  expected 0.5*f_x_star.norm()**2: %.6g" % self.capital_f_x_star)
      print("    actual 0.5*f_x_star.norm()**2: %.6g" % (0.5*f_x_star.norm()**2))

  def f(self, x):
    x1,x2,x3 = x
    result = flex.double()
    for yi,ti in zip(meyer_function.ys,
                     meyer_function.ts):
      denominator = ti + x3
      assert denominator != 0
      result.append(x1 * math.exp(x2/denominator) - yi)
    if (x[0] > 100): # numerical instability on some platforms
      raise MeyerFunctionError
    return result

  def jacobian_analytical(self, x):
    x1,x2,x3 = x
    result = flex.double()
    for ti in meyer_function.ts:
      denominator = ti + x3
      assert denominator != 0
      denominator_sq = denominator**2
      assert denominator_sq != 0
      term = math.exp(x2/denominator)
      result.extend(flex.double([
        term,
        x1*term/denominator,
        -x1*term*x2/denominator_sq]))
    result.resize(flex.grid(self.m, self.n))
    return result

  def hessian_analytical(self, x):
    x1,x2,x3 = x
    result = flex.double(flex.grid(3,3), 0)
    for ti in meyer_function.ts:
      denominator = ti + x3
      assert denominator != 0
      denominator_sq = denominator**2
      assert denominator_sq != 0
      denominator_cu = denominator**3
      assert denominator_cu != 0
      denominator_qa = denominator**4
      assert denominator_qa != 0
      term = math.exp(x2/denominator)
      result[(0,0)] -= 0
      result[(0,1)] -= term/denominator
      result[(0,2)] -= -x2*term/denominator_sq
      result[(1,1)] -= x1*term/denominator_sq
      result[(1,2)] -= -x1*x2*term/denominator_cu - x1*term/denominator_sq
      result[(2,2)] -= x1*x2**2*term/denominator_qa \
                     + 2*x1*x2*term/denominator_cu
    for i in range(0,3):
      for j in range(i+1,3):
        result[(j,i)] = result[(i,j)]
    j = self.jacobian_analytical(x=x)
    result += j.matrix_transpose().matrix_multiply(j)
    return result

def exercise_cholesky():
  mt = flex.mersenne_twister(seed=0)
  for n in range(1,10):
    a = flex.double(n*n,0)
    a.resize(flex.grid(n, n))
    for i in range(n): a[(i,i)] = 1
    c = cholesky_decomposition(a)
    assert c is not None
    assert approx_equal(c.matrix_multiply(c.matrix_transpose()), a)
    b = mt.random_double(size=n, factor=4)-2
    x = cholesky_solve(c, b)
    assert approx_equal(a.matrix_multiply(x), b)
    d = flex.random_size_t(size=n, modulus=10)
    for i in range(n): a[(i,i)] = d[i]+1
    c = cholesky_decomposition(a)
    assert c is not None
    assert approx_equal(c.matrix_multiply(c.matrix_transpose()), a)
    b = mt.random_double(size=n, factor=4)-2
    x = cholesky_solve(c, b)
    assert approx_equal(a.matrix_multiply(x), b)
  #
  a = flex.double([8, -6, 0, -6, 9, -2, 0, -2, 8])
  a.resize(flex.grid(3,3))
  c = cholesky_decomposition(a)
  assert c is not None
  assert approx_equal(c.matrix_multiply(c.matrix_transpose()), a)
  assert approx_equal(c, [
    2.828427125,          0,              0,
    -2.121320344,    2.121320343,         0,
         0.,        -0.9428090418,   2.666666667])
  #
  a0 = matrix.sym(sym_mat3=[3,5,7,1,2,-1])
  for i_trial in range(100):
    r = scitbx.math.euler_angles_as_matrix(
      mt.random_double(size=3,factor=360), deg=True)
    a = flex.double(r * a0 * r.transpose())
    a.resize(flex.grid(3,3))
    c = cholesky_decomposition(a)
    assert c is not None
    assert approx_equal(c.matrix_multiply(c.matrix_transpose()), a)
    for b in [(0.1,-0.5,2), (-0.3,0.7,-1), (1.3,2.9,4), (-10,-20,17)]:
      b = flex.double(b)
      x = cholesky_solve(c, b)
      assert approx_equal(a.matrix_multiply(x), b)
  #
  for n in range(1,10):
    for i in range(10):
      r = mt.random_double(size=n*n, factor=10)-5
      r.resize(flex.grid(n,n))
      a = r.matrix_multiply(r.matrix_transpose())
      c = cholesky_decomposition(a)
      assert c is not None
      b = mt.random_double(size=n, factor=4)-2
      x = cholesky_solve(c, b)
      assert approx_equal(a.matrix_multiply(x), b)
      a[(i%n,i%n)] *= -1
      c = cholesky_decomposition(a)
      assert c is None

def exercise():
  verbose = "--verbose" in sys.argv[1:]
  exercise_cholesky()
  default_flag = True
  if (0 or default_flag):
    for m in range(1,5+1):
      for n in range(1,m+1):
        linear_function_full_rank(m=m, n=n, verbose=verbose)
  if (0 or default_flag):
    for m in range(1,5+1):
      for n in range(1,m+1):
        linear_function_rank_1(m=m, n=n, verbose=verbose)
  if (0 or default_flag):
    for m in range(3,7+1):
      for n in range(3,m+1):
        linear_function_rank_1_with_zero_columns_and_rows(
          m=m, n=n, verbose=verbose)
  if (0 or default_flag):
    rosenbrock_function(m=2, n=2, verbose=verbose)
  if (0 or default_flag):
    helical_valley_function(m=3, n=3, verbose=verbose)
  if (0 or default_flag):
    powell_singular_function(m=4, n=4, verbose=verbose)
  if (0 or default_flag):
    freudenstein_and_roth_function(m=2, n=2, verbose=verbose)
  if (0 or default_flag):
    bard_function(m=15, n=3, verbose=verbose)
  if (0 or default_flag):
    kowalik_and_osborne_function(m=11, n=4, verbose=verbose)
  if (0 or default_flag):
    meyer_function(m=16, n=3, verbose=verbose)
  print(format_cpu_times())

if (__name__ == "__main__"):
  exercise()
