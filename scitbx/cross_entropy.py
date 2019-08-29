from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.python_utils import random_transform
import sys
from six.moves import range
from six.moves import zip


class cross_entropy_optimizer(object):
  def __init__(self, function, mean, sigma, alpha, beta, q=7, elite_size=10, sample_size=50, eps=1e-7, inject_eps=1e-3, n_max=50000, monitor_cycle=50):
    self.function = function
    self.mean = mean
    self.sigma = sigma

    self.alpha = alpha
    self.beta = beta
    self.q = q
    self.n = function.n
    self.n_in = elite_size
    self.sample_size = sample_size

    self.eps = eps
    self.inject_eps = inject_eps

    self.monitor_cycle = monitor_cycle
    self.n_max = n_max

    self.count = 0.0
    self.monitor_cycle=monitor_cycle
    self.last_mean = self.mean.deep_copy()
    self.last_sigma = self.sigma.deep_copy()

    self.best_sol = self.mean.deep_copy()
    self.best_score = 1e90
    self.best_score = self.compute_target( self.best_sol )
    self.best_sol_found = 0

    converged = False
    while not converged:
      self.count += 1.0
      self.generate_new_means_and_sigmas()
      self.inject()
      converged = self.convergence_test()

  def inject(self,c=3.0):
     mdelta = flex.max(self.mean - self.last_mean)
     sdelta = flex.min(self.sigma)
     if sdelta < self.inject_eps:
       self.sigma = self.sigma + max(mdelta*c,c*self.inject_eps)
     self.last_mean = self.mean.deep_copy()
     self.last_sigma = self.sigma.deep_copy()


  def compute_target(self,x):
    t = self.function.target(x)
    if t < self.best_score:
      self.best_score = t
      self.best_sol = x.deep_copy()
      self.best_sol_found = self.count
    return t

  def generate_and_score_samples(self):
    sample_list = []
    target_list = flex.double()
    for ii in range(self.sample_size):
      x = random_transform.t_variate(a=max(2,self.n-1),N=self.n)
      x = x*self.sigma + self.mean
      t = self.compute_target(x )
      sample_list.append( x )
      target_list.append( t )

    order = flex.sort_permutation( flex.double(target_list) )
    return sample_list, t, order

  def generate_new_means_and_sigmas(self):
    s,t,o = self.generate_and_score_samples()
    nm = self.mean*0.0
    nv = self.mean*0.0

    for ii in range(self.n_in):
      nm = nm + s[o[ii]]
      nv = nv + s[o[ii]]*s[o[ii]]
    nm = nm/self.n_in
    nv = nv/self.n_in - nm*nm
    nv =  nv
    self.mean  = self.mean*(1.0-self.alpha)  + self.alpha*nm
    beta = self.beta -self.beta*((1.0-1.0/self.count)**self.q)
    self.sigma = flex.sqrt( self.sigma*self.sigma*(1.0-beta) + beta*nv)
    self.compute_target( self.mean )

  def print_status(self,out=None):
    if out is None:
      out = sys.stdout
    print(" Cycle:  %i"%self.count, file=out)
    for mm in self.best_sol:
      print("%5.3e "%mm, file=out)
    print("Target : %5.3e"%self.best_score, file=out)
    print(file=out)

  def convergence_test(self):
    max_var = flex.max( self.sigma )
    if max_var < self.eps:
      return True
    if self.count - self.best_sol_found > self.monitor_cycle:
      return True
    if self.count ==  self.n_max:
      return True

    return False





class test_rosenbrock_function(object):
  def __init__(self, dim=4):
    self.n = dim*2
    self.dim = dim
    self.means = flex.double( self.n, 2.0 )
    self.sigmas = flex.double( self.n, 5.0 )


    self.target_count=0
    self.optimizer =  cross_entropy_optimizer(self,
                                              mean=self.means,
                                              sigma=self.sigmas,
                                              alpha=0.75,
                                              beta=0.75,
                                              q=8.5,
                                              elite_size=10,
                                              sample_size=50, inject_eps=1e-4,
                                              monitor_cycle=500)
    self.sol = self.optimizer.best_sol
    for ii in self.sol:
      assert abs( ii-1.0 ) < 1e-2

  def target(self, vector):
    self.target_count += 1
    x_vec = vector[0:self.dim]
    y_vec = vector[self.dim:]
    result=0
    for x,y in zip(x_vec,y_vec):
      result+=100.0*((y-x*x)**2.0) + (1-x)**2.0
    return result




def run():
  flex.set_random_seed(0)
  test_rosenbrock_function(1)

if __name__ == "__main__":
  run()
  print("OK")
