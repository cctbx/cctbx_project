from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import random
import scitbx.direct_search_simulated_annealing
from six.moves import range
from six.moves import zip

class test_rosenbrock_function(object):
  def __init__(self, dim=4):
    self.n = dim*2
    self.dim = dim
    self.x = flex.double( self.n, 2.0 )

    self.starting_simplex=[]
    for ii in range(self.n+1):
      self.starting_simplex.append(flex.random_double(self.n) + self.x)

    self.optimizer = scitbx.direct_search_simulated_annealing.dssa(
                          dimension=self.n,
                          matrix = self.starting_simplex,
                          evaluator = self,
                          tolerance=1e-8,
                          further_opt=True,
                          coolfactor=0.6)
    self.x = self.optimizer.get_solution()
  #  m = self.optimizer.get_candi()
  #  for mm in m:
  #    print list(mm)
    for x in self.x:
      assert abs(x-1.0)<1e-2


  def target(self, vector):
    x_vec = vector[0:self.dim]
    y_vec = vector[self.dim:]
    result=0
    for x,y in zip(x_vec,y_vec):
      result+=100.0*((y-x*x)**2.0) + (1-x)**2.0
    return result

def run():
  random.seed(0)
  flex.set_random_seed(0)
  test_rosenbrock_function(3)
  random.seed(0)
  flex.set_random_seed(0)
  test_rosenbrock_function(4)

if __name__ == "__main__":
  run()
