from __future__ import division
import scitbx.math
from scitbx.array_family import flex

floating_point_epsilon_double = scitbx.math.floating_point_epsilon_double_get()

class damped_newton:

  def __init__(self,
        function,
        x0,
        mu0=None,
        tau=1.e-3,
        eps_1=1.e-16,
        eps_2=1.e-16,
        mu_min=1.e-300,
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
      fdp_plus_mu = fdp.deep_copy()
      if (mu > mu_min):
        fdp_plus_mu.matrix_diagonal_add_in_place(value=mu)
      u = fdp_plus_mu.matrix_symmetric_as_packed_u()
      gmw = u.matrix_cholesky_gill_murray_wright_decomposition_in_place()
      number_of_cholesky_decompositions += 1
      h_dn = gmw.solve(b=-fp)
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
          if (mu > mu_min):
            mu *= max(1/3, 1 - (2*df/dl - 1)**3)
          else:
            mu = mu_min
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
    print "scitbx.minimizers.damped_newton results:"
    print "  function:", self.function.label()
    print "  x_star:", list(self.x_star)
    print "  0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2
    print "  0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2
    print "  number_of_iterations:", self.number_of_iterations
    print "  number_of_function_evaluations:", \
        self.number_of_function_evaluations
    print "  number_of_gradient_evaluations:", \
        self.number_of_gradient_evaluations
    print "  number_of_hessian_evaluations:", \
        self.number_of_hessian_evaluations
    print "  number_of_cholesky_decompositions:", \
        self.number_of_cholesky_decompositions

class newton_more_thuente_1994:

  def __init__(self,
        function,
        x0,
        xtol=None,
        gtol=None,
        ftol=None,
        stpmin=None,
        stpmax=None,
        eps_1=1.e-16,
        eps_2=1.e-16,
        k_max=1000):
    self.function = function
    x = x0.deep_copy()
    f_x = function.f(x=x)
    number_of_function_evaluations = 1
    self.f_x0 = f_x
    fp = function.gradients(x=x, f_x=f_x)
    number_of_gradient_evaluations = 1
    number_of_hessian_evaluations = 0
    number_of_cholesky_decompositions = 0
    line_search = scitbx.math.line_search_more_thuente_1994()
    if (xtol is not None): line_search.xtol = xtol
    if (ftol is not None): line_search.ftol = ftol
    if (gtol is not None): line_search.gtol = gtol
    if (stpmin is not None): line_search.stpmin = stpmin
    if (stpmax is not None): line_search.stpmax = stpmax
    k = 0
    while (k < k_max):
      if (flex.max(flex.abs(fp)) <= eps_1):
        break
      fdp = function.hessian(x=x)
      number_of_hessian_evaluations += 1
      u = fdp.matrix_symmetric_as_packed_u()
      gmw = u.matrix_cholesky_gill_murray_wright_decomposition_in_place()
      number_of_cholesky_decompositions += 1
      h_dn = gmw.solve(b=-fp)
      line_search.start(
        x=x,
        functional=function.functional(f_x=f_x),
        gradients=fp,
        search_direction=h_dn,
        initial_estimate_of_satisfactory_step_length=1)
      while (line_search.info_code == -1):
        f_x = function.f(x=x)
        number_of_function_evaluations += 1
        fp = function.gradients(x=x, f_x=f_x)
        number_of_gradient_evaluations += 1
        line_search.next(
          x=x,
          functional=function.functional(f_x=f_x),
          gradients=fp)
      h_dn *= line_search.stp
      k += 1
      if (h_dn.norm() <= eps_2*(eps_2 + x.norm())):
        break
    self.x_star = x
    self.f_x_star = f_x
    self.number_of_iterations = k
    self.number_of_function_evaluations = number_of_function_evaluations
    self.number_of_gradient_evaluations = number_of_gradient_evaluations
    self.number_of_hessian_evaluations = number_of_hessian_evaluations
    self.number_of_cholesky_decompositions = number_of_cholesky_decompositions
    self.line_search_info = line_search.info_meaning

  def show_statistics(self):
    print "scitbx.minimizers.newton_more_thuente_1994 results:"
    print "  function:", self.function.label()
    print "  x_star:", list(self.x_star)
    print "  0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2
    print "  0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2
    print "  number_of_iterations:", self.number_of_iterations
    print "  number_of_function_evaluations:", \
        self.number_of_function_evaluations
    print "  number_of_gradient_evaluations:", \
        self.number_of_gradient_evaluations
    print "  number_of_hessian_evaluations:", \
        self.number_of_hessian_evaluations
    print "  number_of_cholesky_decompositions:", \
        self.number_of_cholesky_decompositions
    print "  line_search_info:", \
        self.line_search_info
