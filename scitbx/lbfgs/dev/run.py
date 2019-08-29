from __future__ import absolute_import, division, print_function
from scitbx import lbfgs as scitbx_lbfgs
from scitbx.array_family import flex
from libtbx import adopt_init_args

# Rosenbrock's function, gradients and curvatures

def target(x,y):
  return (1-x)**2+100*(y-x**2)**2

def grad_x(x,y):
  return -2*(1-x) + 400*(-x*y+x**3)

def grad_y(x,y):
  return 2*100*(y-x**2)

def curv_xx(x,y):
  return 2 + 400*(-y+3*x**2)

def curv_yy(x,y):
  return 200

#

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
      if (requests_diag):
        print("x,f,d:", tuple(x), f, tuple(d))
      else:
        print("x,f:", tuple(x), f)
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

class minimizer:

  def __init__(self, xx=-3, yy=-4, min_iterations=0, max_iterations=10000):
    adopt_init_args(self, locals())
    self.x = flex.double([xx, yy])
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
    self.xx, self.yy = self.x
    if (not requests_f_and_g and not requests_diag):
      requests_f_and_g = True
      requests_diag = True
    if (requests_f_and_g):
      self.f = target(self.xx,self.yy)
      self.g = flex.double(
        (grad_x(self.xx, self.yy),
         grad_y(self.xx, self.yy)))
      self.d = None
    if (requests_diag):
      self.d = flex.double(
        (curv_xx(self.xx, self.yy),
         curv_yy(self.xx, self.yy)))
      assert self.d.all_ne(0)
      self.d = 1 / self.d
    return self.x, self.f, self.g, self.d

def run():
  for use_curvatures in (False, True):
    print("use_curvatures:", use_curvatures)
    m = minimizer().run(use_curvatures=use_curvatures)
    print(tuple(m.x), "final")
    if (abs(m.x[0]) > 1.e-4 or abs(m.x[1]) > 1.e-4):
      print(tuple(m.x), "failure, use_curvatures="+str(use_curvatures))
    print("iter,exception:", m.minimizer.iter(), m.minimizer.error)
    print("n_calls:", m.minimizer.n_calls)
    assert m.minimizer.n_calls == m.minimizer.nfun()
    print()

if (__name__ == "__main__"):
  run()
