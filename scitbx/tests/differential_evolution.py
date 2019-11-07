from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import random
import scitbx.differential_evolution
from six.moves import zip

class test_rosenbrock_function(object):
  def __init__(self, dim=5):
    self.x = None
    self.n = 2*dim
    self.dim = dim
    self.domain = [ (1,3) ]*self.n
    self.optimizer = scitbx.differential_evolution.differential_evolution_optimizer(self,population_size=min(self.n*10,40),n_cross=self.n,cr=0.9, eps=1e-8, show_progress=True)
    print(list(self.x))
    for x in self.x:
      assert abs(x-1.0)<1e-2


  def target(self, vector):
    tmp = vector.deep_copy()
    x_vec = vector[0:self.dim]
    y_vec = vector[self.dim:]
    result=0
    for x,y in zip(x_vec,y_vec):
      result+=100.0*((y-x*x)**2.0) + (1-x)**2.0
    #print list(x_vec), list(y_vec), result
    return result

  def print_status(self, mins,means,vector,txt):
    print(txt,mins, means, list(vector))


def run():
  random.seed(0)
  flex.set_random_seed(0)
  test_rosenbrock_function(1)
  print("OK")


if __name__ == "__main__":
  run()
