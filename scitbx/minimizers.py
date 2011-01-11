from __future__ import division
import scitbx.math
from scitbx.array_family import flex

floating_point_epsilon_double = scitbx.math.floating_point_epsilon_double_get()

class FunctionalException(RuntimeError): pass

class damped_newton(object):

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
      gmw = scitbx.linalg.gill_murray_wright_cholesky_decomposition_in_place(u)
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

class newton_more_thuente_1994(object):

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
        matrix_symmetric_relative_epsilon=1.e-12,
        k_max=1000,
        constant_hessian=False):
    self.function = function
    self.constant_hessian = constant_hessian
    fdp = None
    gmw = None
    u = None
    x = x0.deep_copy()
    function_f, callback_after_step = [getattr(function, attr, None)
      for attr in ["f", "callback_after_step"]]
    if (function_f is not None):
      f_x = function_f(x=x)
      functional = function.functional(f_x=f_x)
      fp = function.gradients(x=x, f_x=f_x)
    else:
      f_x = None
      functional = function.functional(x=x)
      fp = function.gradients(x=x)
    self.f_x0 = f_x
    self.functional_x0 = functional
    number_of_function_evaluations = 1
    number_of_functional_exceptions = 0
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
      if (fdp is None) or ( not self.constant_hessian ):
        fdp = function.hessian(x=x)
        number_of_hessian_evaluations += 1
        u = fdp.matrix_symmetric_as_packed_u(
          relative_epsilon=matrix_symmetric_relative_epsilon)
        gmw = scitbx.linalg \
          .gill_murray_wright_cholesky_decomposition_in_place(u)
        number_of_cholesky_decompositions += 1
      h_dn = gmw.solve(b=-fp)
      initial_step_length = 1
      backup_x = x.deep_copy()
      if (function_f is not None):
        backup_f_x = f_x.deep_copy()
      backup_functional = functional
      backup_fp = fp.deep_copy()
      while True:
        try:
          line_search.start(
            x=x,
            functional=functional,
            gradients=fp,
            search_direction=h_dn,
            initial_estimate_of_satisfactory_step_length=initial_step_length)
          while (line_search.info_code == -1):
            if (function_f is not None):
              f_x = function_f(x=x)
              functional = function.functional(f_x=f_x)
              fp = function.gradients(x=x, f_x=f_x)
            else:
              functional = function.functional(x=x)
              fp = function.gradients(x=x)
            number_of_function_evaluations += 1
            number_of_gradient_evaluations += 1
            line_search.next(x=x, functional=functional, gradients=fp)
        except FunctionalException:
          number_of_functional_exceptions += 1
          x = backup_x.deep_copy()
          if (function_f is not None):
            f_x = backup_f_x.deep_copy()
          functional = backup_functional
          fp = backup_fp.deep_copy()
          initial_step_length *= 0.5
        else:
          break
      h_dn *= line_search.stp
      k += 1
      if (callback_after_step): callback_after_step(x=x)
      if (h_dn.norm() <= eps_2*(eps_2 + x.norm())):
        break
    self.x_star = x
    self.f_x_star = f_x
    self.functional_x_star = functional
    self.number_of_iterations = k
    self.number_of_function_evaluations = number_of_function_evaluations
    self.number_of_functional_exceptions = number_of_functional_exceptions
    self.number_of_gradient_evaluations = number_of_gradient_evaluations
    self.number_of_hessian_evaluations = number_of_hessian_evaluations
    self.number_of_cholesky_decompositions = number_of_cholesky_decompositions
    self.line_search_info = line_search.info_meaning

  def show_statistics(self):
    print "scitbx.minimizers.newton_more_thuente_1994 results:"
    get_label = getattr(self.function, "label", None)
    if (get_label is not None):
      print "  function:", get_label()
    print "  x_star:", list(self.x_star)
    if (self.f_x0 is not None):
      print "  0.5*f_x0.norm()**2:", 0.5*self.f_x0.norm()**2
      print "  0.5*f_x_star.norm()**2:", 0.5*self.f_x_star.norm()**2
    print "  functional_x0:", self.functional_x0
    print "  functional_x_star:", self.functional_x_star
    print "  number_of_iterations:", self.number_of_iterations
    print "  number_of_function_evaluations:", \
        self.number_of_function_evaluations
    print "  number_of_functional_exceptions:", \
        self.number_of_functional_exceptions
    print "  number_of_gradient_evaluations:", \
        self.number_of_gradient_evaluations
    print "  number_of_hessian_evaluations:", \
        self.number_of_hessian_evaluations
    print "  number_of_cholesky_decompositions:", \
        self.number_of_cholesky_decompositions
    print "  line_search_info:", \
        self.line_search_info
