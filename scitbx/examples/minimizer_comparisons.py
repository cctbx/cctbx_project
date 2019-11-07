from __future__ import absolute_import, division, print_function
from scitbx.examples.functions import Function

from scitbx import simplex
from scitbx import direct_search_simulated_annealing
from scitbx import differential_evolution
from scitbx import cross_entropy

from scitbx import lbfgs

from scitbx.array_family import flex

import libtbx.load_env
from six.moves import range

dim=2
start = flex.double( [4,4] )

funct_count=0

def derivative(vector,f,h=1e-8):
  ori = vector.deep_copy()
  tmp = vector.deep_copy()
  f0 = f(ori)
  f_vector=tmp*0.0
  for ii,var in enumerate(tmp):
    tmp[ii] = var + h
    f_vector[ii]=f(tmp)
    tmp[ii]=var
  f_vector = (f_vector-f0)/h
  return f_vector

class test_lbfgs(object):
  def __init__(self, name):
    self.x = start.deep_copy()
    self.n = 2
    self.name = name
    self.fcount=0
    self.minimizer = lbfgs.run(target_evaluator=self)
    print("LBFGS ITERATIONS", self.minimizer.iter(), self.fcount," SOLUTION", list(self.x), self.function(self.x))

  def function(self,vector):
    self.fcount+=1
    result = Function(self.name)(dim).eval(vector)
    return result

  def compute_functional_and_gradients(self):
    f = self.function(self.x)
    g = derivative(self.x, self.function,h=1e-5)
    return f, g

class test_simplex(object):
  def __init__(self, name):
    self.n = 2
    self.x = start.deep_copy()
    self.name=name
    self.starting_simplex=[]
    self.fcount=0
    for ii in range(self.n+1):
        self.starting_simplex.append((flex.random_double(self.n)/2-1)*0.0003+ self.x)
    self.optimizer = simplex.simplex_opt(
                                  dimension=self.n,
                                  matrix  = self.starting_simplex,
                                  evaluator = self,
                                  tolerance=1e-10)
    self.x = self.optimizer.get_solution()
    print("SIMPLEX ITRATIONS", self.optimizer.count,self.fcount, "SOLUTION", list(self.x), self.target(self.x))

  def target(self,vector):
    self.fcount += 1
    return  Function(self.name)(dim).eval(vector)


class test_dssa(object):
  def __init__(self,name):
    self.name=name
    self.fcount=0
    self.n = 2
    self.x = start.deep_copy()
    self.starting_simplex=[]
    for ii in range(self.n+1):
      self.starting_simplex.append(3.1*(flex.random_double(self.n)/2-1) + self.x)
    self.optimizer = direct_search_simulated_annealing.dssa(
                          dimension=self.n,
                          matrix = self.starting_simplex,
                          evaluator = self,
                          further_opt=True,
                          coolfactor=0.5, simplex_scale=1
                                          )
    self.x = self.optimizer.get_solution()
    print("DSSA ITERATIONS", self.optimizer.count, self.fcount, "SOLUTION", list(self.x), self.target(self.x))

  def function(self,vector):
    self.fcount+=1
    return Function(self.name)(dim).eval(vector)

  def target(self,vector):
    return self.function(vector)

  def compute_functional_and_gradients(self):
    f = self.function(self.x)
    g = derivative(self.x, self.function,h=1e-5)
    return f, g



class test_cross_entropy(object):
  def __init__(self,name):
    self.name=name
    self.fcount=0
    self.n = 2
    self.x = start.deep_copy()
    self.means = flex.double( self.n, 4.0 )
    self.sigmas = flex.double( self.n, 2.0 )
    self.optimizer =  cross_entropy.cross_entropy_optimizer(self,
                                              mean=self.means,
                                              sigma=self.sigmas,
                                              alpha=0.95,
                                              beta=0.75,
                                              q=8.5,
                                              elite_size=10,
                                              sample_size=50, inject_eps=1e-4,
                                              monitor_cycle=150,eps=1e-8)
    print("CROSS ENTROPY ITERATIONS", self.optimizer.count, self.fcount, "SOLUTION", list(self.optimizer.best_sol), self.target(self.optimizer.best_sol))

  def function(self,vector):
    self.fcount += 1
    result = Function(self.name)(dim).eval(vector)
    return result

  def target(self, vector):
    return self.function(vector)

  def compute_functional_and_gradients(self):
    f = self.function(self.x)
    g = derivative(self.x, self.function,h=1e-5)
    return f, g





class test_differential_evolution(object):
  def __init__(self,name,npop=20):
    self.name=name
    self.fcount=0
    self.n = 2
    self.x = None #flex.double(self.n, 2)
    self.domain = [(start[0]-1,start[0]+1),(start[1]-1, start[1]+1)]
    self.optimizer =  differential_evolution.differential_evolution_optimizer(self,population_size=npop,cr=0.9,n_cross=2,eps=1e-12,show_progress=False)
    print("DIFFERENTIAL EVOLUTION ITERATIONS", self.optimizer.count, self.fcount, "SOLUTION", list(self.x), self.target(self.x))

  def function(self,vector):
    self.fcount+=1
    result = Function(self.name)(dim).eval(vector)
    return result

  def target(self,vector):
    return self.function(vector)

  def compute_functional_and_gradients(self):
    f = self.function(self.x)
    g = derivative(self.x, self.function,h=1e-5)
    return f, g

  def print_status(self, mins,means,vector,txt):
    print(txt,mins, means, list(vector))



class test_cma_es(object):
  def __init__(self,name,l=0):
    self.m = start.deep_copy()
    self.s = flex.double( [2,2])
    self.l = l
    self.name = name
    self.fcount = 0
    from cma_es import cma_es_interface
    self.minimizer = cma_es_interface.cma_es_driver( 2, self.m, self.s, self.my_function, self.l )
    print("CMA-ES ITERATIONS", self.minimizer.count, self.fcount,"SOLUTION",  list(self.minimizer.x_final), self.my_function( self.minimizer.x_final ))

    self.x = self.minimizer.x_final.deep_copy()


  def compute_functional_and_gradients(self):
    f = self.my_function(self.x)
    g = derivative(self.x, self.my_function,h=1e-5)
    return f, g

  def my_function(self,vector):
    self.fcount+=1
    tmp = Function(self.name)(dim)
    return tmp.eval(list(vector))


def run(args):
  assert len(args) == 0

  have_cma_es = libtbx.env.has_module("cma_es")
  if (not have_cma_es):
    print("Skipping some tests: cma_es module not available or not configured.")
    print()

  names = ['easom','rosenbrock','ackley','rastrigin']
  for name in names:
    print("****", name, "****")
    if name == 'easom':
      start = flex.double( [0.0,0.0] )
    else:
      start = flex.double( [4,4] )
    for ii in range(1):
      test_lbfgs(name)
      if (have_cma_es):
        test_cma_es(name)
      test_differential_evolution(name)
      test_cross_entropy(name)
      test_simplex(name)
      test_dssa(name)
      print()
    print()
  from libtbx.utils import format_cpu_times
  print(format_cpu_times())

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
