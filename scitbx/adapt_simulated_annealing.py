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
  def __init__(self,parent,start_solution, scales,start_t=1000, end_t=1, burn_in=200, burn_out=100, steps=100, coolfactor=0.96):
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
    self.coolfactor = coolfactor
    self.maxStep_in = 20
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
      for jj in xrange(self.maxStep_in):
        tmp = self.move()
        if tmp:
          if self.current_score < self.best_score:
            self.best_score = self.current_score
            self.best_sol = self.sol

  def make_anneal_profile(self):
    ENE = flex.double()
    for ii in range(25):
      new_vector = self.perturb( self.sol )
      ENE.append( self.parent.target( new_vector ) )
    dE = ( ENE - flex.mean(ENE) ).norm()/ 5.0

    self.start_t = dE
    for ii in xrange(self.burn_in):
      self.temp_scheme.append( self.start_t )
    tmp_t = self.start_t
    while(tmp_t > self.end_t):
      tmp_t = tmp_t * self.coolfactor
      self.temp_scheme.append( tmp_t )
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
    print list(new_vector), new_score, self.current_score, self.temp
    if new_score <= self.current_score:
      self.sol = new_vector
      self.current_score = new_score
      return True
    else:
      delta = self.current_score - new_score
      eps = flex.exp(flex.double( [delta/self.temp] ))[0]
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
    self.optimizer = sa_optimizer(self,self.x, self.scales, coolfactor=0.9)

    self.sol, self.score = self.optimizer.get_solution()
    #assert (self.score<1e-3)
    if(self.score <1e-3):
      print "good"

  def target(self, vector):
    result = (vector[0]-self.a)**2.0
    result = result + abs((vector[1]-self.b)**3.0)
    return result*1000



def tst():
  for ii in range(10):
    tst_function()
  print "OK"


if __name__ == """__main__""":
  tst()
