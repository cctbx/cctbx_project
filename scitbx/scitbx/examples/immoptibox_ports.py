"""\
http://www2.imm.dtu.dk/~hbn/immoptibox/
"""

from __future__ import division
import scitbx.math
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import math

floating_point_epsilon_double = scitbx.math.floating_point_epsilon_double_get()

def cholesky_decomposition(a, relative_eps=1.e-15):
  assert a.nd() == 2
  assert a.focus()[0] == a.focus()[1]
  n = a.focus()[0]
  eps = relative_eps * flex.max(flex.abs(a))
  c = flex.double(a.accessor(), 0)
  for k in xrange(n):
    sum = 0
    for j in xrange(k):
      sum += c[(k,j)]**2
    d = a[(k,k)] - sum
    if (d <= eps): return None
    c[(k,k)] = math.sqrt(d)
    for i in xrange(k+1,n):
      sum = 0
      for j in xrange(k):
        sum += c[(i,j)] * c[(k,j)]
      c[(i,k)] = (a[(i,k)] - sum) / c[(k,k)]
  return c

def cholesky_solve(c, b):
  assert c.nd() == 2
  assert c.focus()[0] == c.focus()[1]
  assert b.nd() == 1
  assert b.size() == c.focus()[0]
  n = c.focus()[0]
  z = flex.double(n, 0)
  for k in xrange(n):
    sum = 0
    for j in xrange(k):
      sum += c[(k,j)] * z[j]
    z[k] = (b[k] - sum) / c[(k,k)]
  x = flex.double(n, 0)
  for k in xrange(n-1,-1,-1):
    sum = 0
    for j in xrange(k+1,n):
      sum += c[(j,k)] * x[j]
    x[k] = (z[k] - sum) / c[(k,k)]
  return x

class levenberg_marquardt:

  def __init__(self,
        function,
        x0,
        tau=1.e-3,
        eps_1=1.e-8,
        eps_2=1.e-12,
        k_max=100):
    nu = 2
    x = x0
    f_x = function.f(x)
    number_of_function_evaluations = 1
    self.f_x0 = f_x
    j = function.jacobian(x)
    number_of_jacobian_evaluations = 1
    number_of_svds = 0
    j_t = j.matrix_transpose()
    a = j_t.matrix_multiply(j)
    g = j_t.matrix_multiply(f_x)
    found = flex.max(flex.abs(g)) <= eps_1
    mu = tau * flex.max(a.matrix_diagonal())
    k = 0
    while (not found and k < k_max):
      k += 1
      a_plus_mu = a.deep_copy()
      a_plus_mu.matrix_diagonal_add_in_place(value=mu)
      from tntbx.generalized_inverse import generalized_inverse
      a_plus_mu_svd = generalized_inverse(square_matrix=a_plus_mu)
      number_of_svds += 1
      h_lm = a_plus_mu_svd.matrix_multiply(-g)
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
          mu *= max(1/3., 1-(2*rho-1)**3)
          nu = 2
        else:
          mu *= nu
          nu *= 2
    self.x_star = x
    self.f_x_star = f_x
    self.number_of_iterations = k
    self.number_of_function_evaluations = number_of_function_evaluations
    self.number_of_jacobian_evaluations = number_of_jacobian_evaluations
    self.number_of_svds = number_of_svds

  def show_statistics(self):
    print "x_star:", list(self.x_star)
    print "0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2
    print "0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2
    print "number_of_iterations:", self.number_of_iterations
    print "number_of_function_evaluations:",self.number_of_function_evaluations
    print "number_of_jacobian_evaluations:",self.number_of_jacobian_evaluations
    print "number_of_svds:", self.number_of_svds

class damped_newton:

  def __init__(self,
        function,
        x0,
        mu0=None,
        tau=1.e-3,
        eps_1=1.e-8,
        eps_2=1.e-12,
        k_max=100):
    x = x0
    f_x = function.f(x=x)
    number_of_function_evaluations = 1
    self.f_x0 = f_x
    fp = function.gradients_analytical(x=x, f_x=f_x)
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
        mu = tau * flex.max(fdp.matrix_diagonal())
        if (mu < floating_point_epsilon_double**0.5):
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
      fp_new = function.gradients_analytical(x=x_new, f_x=f_x_new)
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
    print "x_star:", list(self.x_star)
    print "0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2
    print "0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2
    print "number_of_iterations:", self.number_of_iterations
    print "number_of_function_evaluations:",self.number_of_function_evaluations
    print "number_of_gradient_evaluations:",self.number_of_gradient_evaluations
    print "number_of_hessian_evaluations:",self.number_of_hessian_evaluations
    print "number_of_cholesky_decompositions:", \
      self.number_of_cholesky_decompositions

class test_function:

  def __init__(self, m, n):
    assert m >= n
    self.m = m
    self.n = n
    self.initialization()
    if (self.x_star is not None):
      assert approx_equal(
        0.5*self.f(x=self.x_star).norm()**2, self.capital_f_x_star)
    self.exercise_levenberg_marquardt()
    self.exercise_damped_newton()

  def jacobian_finite(self, x, relative_eps=1.e-5):
    x0 = x
    result = flex.double()
    for i in xrange(self.n):
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
    finite = self.jacobian_finite(x=x)
    assert approx_equal(analytical, finite)
    return analytical

  def gradients_analytical(self, x, f_x=None):
    if (f_x is None): f_x = self.f(x=x)
    return self.jacobian_analytical(x=x).matrix_transpose() \
      .matrix_multiply(f_x)

  def hessian_finite(self, x, relative_eps=1.e-5):
    x0 = x
    result = flex.double()
    for i in xrange(self.n):
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
    finite = self.hessian_finite(x=x)
    assert approx_equal(analytical, finite)
    return analytical

  def check_minimized_x_star(self, x_star):
    assert approx_equal(x_star, self.x_star)

  def exercise_levenberg_marquardt(self):
    minimized = levenberg_marquardt(function=self, x0=self.x0, tau=self.tau0)
    minimized.show_statistics()
    self.check_minimized_x_star(x_star=minimized.x_star)
    assert approx_equal(
      0.5*minimized.f_x_star.norm()**2, self.capital_f_x_star)
    print

  def exercise_damped_newton(self):
    minimized = damped_newton(function=self, x0=self.x0, tau=self.tau0)
    minimized.show_statistics()
    self.check_minimized_x_star(x_star=minimized.x_star)
    assert approx_equal(
      0.5*minimized.f_x_star.norm()**2, self.capital_f_x_star)
    print

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
    self.a = flex.double(xrange(1,m+1)).matrix_outer_product(
             flex.double(xrange(1,n+1)))
    self.x0 = flex.double(n, 1)
    self.tau0 = 1.e-8
    self.delta0 = 10
    self.x_star = None
    self.capital_f_x_star = (m*(m-1))/(4*(2*m+1))

  def check_minimized_x_star(self, x_star):
    assert approx_equal(
      flex.double(xrange(1,self.n+1)).dot(x_star),
      3/(2*self.m+1))

class linear_function_rank_1_with_zero_columns_and_rows(
        linear_function_full_rank):

  def initialization(self):
    m = self.m
    n = self.n
    self.a = flex.double([0]+range(1,m-2+1)+[0]).matrix_outer_product(
             flex.double([0]+range(2,n-1+1)+[0]))
    self.x0 = flex.double(n, 1)
    self.tau0 = 1.e-8
    self.delta0 = 10
    self.x_star = None
    self.capital_f_x_star = (m**2+3*m-6)/(4*(2*m-3))

  def check_minimized_x_star(self, x_star):
    assert approx_equal(
      flex.double([0]+range(2,self.n-1+1)+[0]).dot(x_star),
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
    x1, x2, x3 = x
    assert x1 != 0
    t = math.atan(x2/x1)/(2*math.pi)
    if (x1 < 0): t += 0.5
    nx = x[:2].norm()
    return flex.double([10*(x3 - 10*t), 10*(nx - 1), x3])

  def jacobian_analytical(self, x):
    x1, x2, x3 = x
    nx = x[:2].norm()
    nx2 = nx*nx
    k1 = 50/math.pi/nx2
    k2 = 10/nx
    return flex.double([
      [k1*x2, -k1*x1, 10],
      [k2*x1, k2*x2, 0],
      [0, 0, 1]])

  def hessian_analytical(self, x):
    x1, x2, x3 = x
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

def exercise_cholesky():
  mt = flex.mersenne_twister(seed=0)
  for n in xrange(1,10):
    a = flex.double(n*n,0)
    a.resize(flex.grid(n, n))
    for i in xrange(n): a[(i,i)] = 1
    c = cholesky_decomposition(a)
    assert c is not None
    assert approx_equal(c.matrix_multiply(c.matrix_transpose()), a)
    b = mt.random_double(size=n, factor=4)-2
    x = cholesky_solve(c, b)
    assert approx_equal(a.matrix_multiply(x), b)
    d = flex.random_size_t(size=n, modulus=10)
    for i in xrange(n): a[(i,i)] = d[i]+1
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
  a0 = matrix.sym([3,5,7,1,2,-1])
  for i_trial in xrange(100):
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
  for n in xrange(1,10):
    for i in xrange(10):
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
  exercise_cholesky()
  try: import tntbx
  except ImportError: pass
  else:
    if (1):
      for m in xrange(1,5+1):
        for n in xrange(1,m+1):
          linear_function_full_rank(m=m, n=n)
    if (1):
      for m in xrange(1,5+1):
        for n in xrange(1,m+1):
          linear_function_rank_1(m=m, n=n)
    if (1):
      for m in xrange(3,7+1):
        for n in xrange(3,m+1):
          linear_function_rank_1_with_zero_columns_and_rows(m=m, n=n)
    if (1):
      rosenbrock_function(m=2, n=2)
    if (1):
      helical_valley_function(m=3, n=3)
  print "OK"

if (__name__ == "__main__"):
  exercise()
