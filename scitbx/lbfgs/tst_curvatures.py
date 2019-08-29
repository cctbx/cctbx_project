"""Test case dramatically improves LBFGS multiparameter refinement by
using curvatures to provide relative weightings for the contributions of
each parameter to the target functional.
This implementation is based on code from Ralf Grosse-Kunstleve
in the module scitbx/lbfgs/dev/twisted_gaussian.py.
Here, the re-usable part of the code is abstracted to a mix-in class
that can be used by any other application wishing to use curvatures."""
from __future__ import absolute_import, division, print_function
from six.moves import range

def lbfgs_run(target_evaluator,
              min_iterations=0,
              max_iterations=None,
              traditional_convergence_test=1,
              traditional_convergence_test_eps=None,
              use_curvatures=False,verbose=False):
  from scitbx import lbfgs as scitbx_lbfgs
  ext = scitbx_lbfgs.ext
  minimizer = ext.minimizer(target_evaluator.n)
  minimizer.error = None
  if (traditional_convergence_test and traditional_convergence_test_eps is not None):
    is_converged = ext.traditional_convergence_test(n = target_evaluator.n,
      eps = traditional_convergence_test_eps)
  elif (traditional_convergence_test):
    is_converged = ext.traditional_convergence_test(target_evaluator.n)
  else:
    raise RuntimeError
    is_converged = ext.drop_convergence_test(min_iterations)
  callback_after_step = getattr(target_evaluator, "callback_after_step", None)
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
      if verbose:
        if (requests_diag):
          print("x,f,d:", tuple(x), f, tuple(d))
        else:
          print("x,f:", tuple(x), f)
      if (use_curvatures):
        from scitbx.array_family import flex
        if (d is None): d = flex.double(x.size())
        have_request = minimizer.run(x, f, g, d)
      else:
        have_request = minimizer.run(x, f, g)
      if (callback_after_step is not None):
        callback_after_step(minimizer)
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

class lbfgs_with_curvatures_mix_in(object):
  """Mix in this class with your application-specific minimizer class.
     The app-specific class is required to implement:
       1. In the constructor
            a) define the free parameters as self.x = flex.double()
            b) last line: call the lbfgs_with_curvatures_mix_in constructor
       2. compute_functional_and_gradients() returns (double functional, flex.double gradients)
       3. curvatures() returns flex.double curvatures == diagonal elements of Hessian matrix

     Exposed variable traditional_convergence_test_eps gives some flexibility
     to the application program to specify the epsilon value for hitting
     convergence.  Choosing the right value avoids unnecessary minimization
     iterations that optimize the target functional to the last decimal place.
  """

  def __init__(self, min_iterations=0, max_iterations=1000,
        traditional_convergence_test_eps=None,
        use_curvatures=True):
    self.n = len(self.x)
    self.minimizer = lbfgs_run(
        target_evaluator=self,
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        traditional_convergence_test_eps=traditional_convergence_test_eps,
        use_curvatures=use_curvatures)

  def __call__(self, requests_f_and_g, requests_diag):
    if (not requests_f_and_g and not requests_diag):
      requests_f_and_g = True
      requests_diag = True
    if (requests_f_and_g):
      self.f, self.g = self.compute_functional_and_gradients()
      self.d = None
    if (requests_diag):
      self.d = self.curvatures()
      #assert self.d.all_gt(0) # conceptually curvatures must be positive to be within convergence well
      sel = (self.g != 0)
      self.d.set_selected(~sel,1000) # however, if we've decided to not refine a certain parameter, we
                                     # can indicate this to LBFGS by setting the gradient to zero.
                                     # Then we can set the curvature to an arbitrary positive value that
                                     # will be tested for positivity but otherwise ignored.
      assert self.d.select(sel).all_gt(0)
      self.d = 1 / self.d
    return self.x, self.f, self.g, self.d


class fit_xy_translation(lbfgs_with_curvatures_mix_in):
  """Elementary test case.  Six groups of xy observations generated by
     random Gaussian with mean = 0 and sigma = 1.  Fit the population means with
     six xy model means.  However, since the sample sizes are
     extremely different, convergence is slow if using a global least-squares
     target function.  Convergence speeds up considerably if curvatures
     are added.
  """

  def __init__(self,use_curvatures,verbose):
    self.verbose=verbose
    self.create_test_data()
    from scitbx.array_family import flex

    self.count_iterations = 0
    self.x = flex.double([0.]*12) # x & y displacements for each of 6 data groups
    lbfgs_with_curvatures_mix_in.__init__(self,
      min_iterations=0,
      max_iterations=1000,
      use_curvatures=use_curvatures)
    if self.verbose:
      print(["%8.5f"%a for a in self.x[0::2]])
      print(["%8.5f"%a for a in self.x[1::2]])

  def curvatures(self):
    from scitbx.array_family import flex
    curvs = flex.double([0.]*12)
    for x in range(6):
      selection = (self.master_groups==x)
      curvs[2*x] = 2. * selection.count(True)
      curvs[2*x+1]=2. * selection.count(True)
    return curvs

  def compute_functional_and_gradients(self):
    from scitbx.array_family import flex
    import math

    self.model_mean_x = flex.double(len(self.observed_x))
    self.model_mean_y = flex.double(len(self.observed_x))

    for x in range(6):
      selection = (self.master_groups==x)
      self.model_mean_x.set_selected(selection, self.x[2*x])
      self.model_mean_y.set_selected(selection, self.x[2*x+1])

    delx = self.observed_x - self.model_mean_x
    dely = self.observed_y - self.model_mean_y
    delrsq = delx*delx + dely*dely

    f = flex.sum(delrsq)

    gradients = flex.double([0.]*12)
    for x in range(6):
      selection = (self.master_groups==x)
      gradients[2*x] = -2. * flex.sum( delx.select(selection) )
      gradients[2*x+1]=-2. * flex.sum( dely.select(selection) )
    if self.verbose:
      print("Functional ",math.sqrt(flex.mean(delrsq)))
    self.count_iterations += 1
    return f,gradients

  def create_test_data(self):
    from scitbx.array_family import flex
    self.observed_x = flex.double()
    self.observed_y = flex.double()
    self.master_groups = flex.int()
    igroup = 0
    import random
    random.seed(0.0)

    for group_N in [30000, 25000, 20000, 10000, 2000, 5]:  # six data groups of different sizes
      for x in range(group_N):
        self.observed_x.append( random.gauss(0.,1.) )
        self.observed_y.append( random.gauss(0.,1.) )
        self.master_groups.append(igroup)
      igroup+=1

    self.known_mean_x = flex.double()
    self.known_mean_y = flex.double()
    for x in range(6):
      selection = (self.master_groups==x)
      self.known_mean_x.append( flex.mean( self.observed_x.select(selection) ) )
      self.known_mean_y.append( flex.mean( self.observed_y.select(selection) ) )
    if self.verbose:
      print(["%8.5f"%a for a in self.known_mean_x])
      print(["%8.5f"%a for a in self.known_mean_y])

def run(verbose=False):

  C = fit_xy_translation(use_curvatures=True,verbose=verbose)
  assert C.count_iterations < 10 # should be 7 on Linux
  C = fit_xy_translation(use_curvatures=False,verbose=verbose)
  assert C.count_iterations > 40 # should be 50 on Linux

  return None

if (__name__ == "__main__"):

  result = run(verbose=False)
  print("OK")
