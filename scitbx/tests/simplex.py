from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import scitbx.simplex
from six.moves import range
from six.moves import zip

class test_function(object):
  def __init__(self,n):
    self.n = n
    self.starting_simplex=[]
    for ii in range(self.n+1):
      self.starting_simplex.append(flex.random_double(self.n))
    self.optimizer = scitbx.simplex.simplex_opt( dimension=self.n,
                                  matrix  = self.starting_simplex,
                                  evaluator = self,
                                  tolerance=1e-10)
    self.x = self.optimizer.get_solution()
    for ii in range(self.n):
      assert approx_equal(self.x[ii],ii+1,1e-5)

  def target(self, vector):
    result = 0.0
    for ii in range(self.n):
      result += (vector[ii]-ii-1)**2.0
    result +=1.0
    return result

class test_rosenbrock_function(object):
  def __init__(self, dim=2):
    self.dim = dim
    self.n = dim * 2
    self.starting_simplex=[]
    self.x = flex.double(self.n, 2)
    for ii in range(self.n+1):
      self.starting_simplex.append(flex.random_double(self.n)+ self.x)
    self.optimizer = scitbx.simplex.simplex_opt( dimension=self.n,
                                  matrix  = self.starting_simplex,
                                  evaluator = self,
                                  tolerance=1e-10)
    self.x = self.optimizer.get_solution()
    assert (abs(self.x[0]-1.0) < 1e-6)
    assert (abs(self.x[1]-1.0) < 1e-6)

  def target(self, vector):
    result = 0
    for x, y in zip( vector[0:self.dim], vector[self.dim:]):
      result = (1-x)**2.0 + 100.0*(y-x*x)**2.0
    #print(x, y, result)
    return result

def run():
  flex.set_random_seed(0)
  for ii in range(10):
    test_rosenbrock_function(1)
    test_function(1)
    test_function(2)
    test_function(3)
    test_function(4)

if __name__ == "__main__":
  run()
