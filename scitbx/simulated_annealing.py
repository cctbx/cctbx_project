"""
a basic simulated annealing engine.
as always, the user supplies the function 'target' in the main driver,
the rest is taken care of here

This implementation is based on

http://nl.wikipedia.org/wiki/Simulated_annealing

"""

from scitbx.array_family import flex
import os,sys,math


class sa_optimizer(object):
  def __init__(self,parent,start_solution, scales,start_t=100, end_t=1, burn_in=100, burn_out=10000, steps=100):
    self.parent = parent
    self.scales = scales
    self.n = scales.size()
    self.sol = start_solution
    self.current_score = self.parent.target(  self.sol )


    self.temp_scheme = []
    self.start_t = start_t
    self.end_t = end_t
    self.burn_in = burn_in
    self.burn_out = burn_out
    self.steps = steps
    self.make_anneal_profile()

    self.best_sol = self.sol.deep_copy()
    self.best_score = float( self.current_score )

    self.run()

  def get_solution(self):
    return self.best_sol, self.best_score

  def run(self):
    #simple now really
    n = len(self.temp_scheme)
    for ii in xrange(n):
      self.temp=self.temp_scheme[ii]
      tmp = self.move()
      if tmp:
        if self.current_score < self.best_score:
          self.best_score = self.current_score
          self.best_sol = self.sol

  def make_anneal_profile(self):
    for ii in xrange(self.burn_in):
      self.temp_scheme.append( self.start_t )
    dt = self.start_t-self.end_t
    dt = dt / (1+self.steps)
    for ii in xrange(self.steps):
      self.temp_scheme.append( self.start_t-ii*dt )
    for ii in xrange(self.burn_out):
      self.temp_scheme.append( self.end_t )


  def perturb(self, vector):
    delta = flex.random_double(self.n)*2-1
    delta = delta*self.scales
    new_vector = vector + delta
    return new_vector

  def move(self):
    new_vector = self.perturb( self.sol )
    new_score = self.parent.target( new_vector )
    if new_score <= self.current_score:
      self.sol = new_vector
      self.current_score = new_score
      return True
    else:
      delta = self.current_score - new_score
      eps = math.exp( delta/self.temp )
      r = flex.random_double(1)[0]
      if eps > r:
        self.sol = new_vector
        self.current_score = new_score
        return True
      return False






class tst_function(object):
  def __init__(self):
    self.a = 1.2
    self.b = 4.5
    self.scales = flex.double([0.1,0.1])
    self.x = flex.double( [12,12] )
    self.optimizer = sa_optimizer(self,self.x, self.scales)

    self.sol, self.score = self.optimizer.get_solution()
    assert (self.score<1e-3)

  def target(self, vector):
    result = (vector[0]-self.a)**2.0
    result = result + abs((vector[1]-self.b)**3.0)
    return result*1000



def tst():
  tst_function()
  print "OK"


if __name__ == """__main__""":
  tst()
