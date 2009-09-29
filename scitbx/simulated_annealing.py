"""
a basic simulated annealing engine.
as always, the user supplies the function 'target' in the main driver,
the rest is taken care of here

This implementation is based on

http://nl.wikipedia.org/wiki/Simulated_annealing

and some ideas adapated from Adaptive Simulated Annealing

http://www.ingber.com/asa96_lessons.pdf


The main drive behind this routine is to have a simulated annealing routine that is relatively self guiding.
The input parameters are:
               parent:  The parent object containing the target function def target(self,vector)
               start_solution: the best starting solution available at this time
               scales: initial scales defining the domain in which the solution will be found
               burn_in: constant temperature moves before decreasing the temperature
               burn_out: constant temperature moves after annealing is done
               steps: number of steps
               coolfactor=0.96: cooling rate : T(k+1)=T*coolfactor
               max_step_in: perform n steps at a given temperature
               spring_constant=0.95: a constant introducing bias towards the best solution
               reanneal: maxiumum number of re-anealing steps
               n_screen: size of initial random population
               eps=1e-5: epsilon for convergence test (maximum change in varaibles)

In total, reanneal annealing runs are run, each time starting with the last best solution.
The bias towards the best solution encountered is done as a factor of temperature.
The lower T, the more likely one is to get a sample close to the best solution.
Furthermore, we only sample in a known range, the range of samping is biased by the last number
of know good solutions. Without this bias, it is very easy to walk away from the right solution.

although some tweaking is needed, typically set the spring_constant to about 0.8 to 0.95.
The cooling rate should be aggresive, especially when the number of steps is small.


The algorithm works well on the Rosenbrock functions:

http://en.wikipedia.org/wiki/Rosenbrock_function


I suggest playing with number of steps and cooling rate to get reasonable behavoir for a given problem.

"""

from scitbx.array_family import flex
import os,sys,math


class trial_queue(object):
  def __init__(self, n=5, min_span=1e-1):
    self.queue = []
    self.targets = flex.double()
    self.order = flex.double()
    self.n = n
    self.worst = None
    self.best = None
    self.span = None
    self.min_span=min_span

  def set_worst_best_oldest(self):
    self.worst = flex.min_index( self.targets )
    self.best = flex.max_index( self.targets )
    self.oldest = flex.min_index( flex.double(self.order) )

  def add_solution(self, sol, score):
    m = len(self.queue)
    if m < self.n:
      self.queue.append( sol.deep_copy() )
      self.targets.append( score )
      self.order.append(self.order.size())
      self.set_worst_best_oldest()
    else:
      # replace the oldest solution
      self.order = self.order+1
      self.queue[ self.oldest ] = sol.deep_copy()
      self.targets[ self.oldest ] = score
      self.order[ self.oldest  ] = flex.max(self.order+1)
      self.order = self.order-1
      self.set_worst_best_oldest()
      self.get_span()

  def get_span(self):
    s_result = []
    m = self.queue[0].size()
    for ii in range(m):
      s_result.append( [] )

    for index in range(m):
      d = flex.double()
      for sol in self.queue:
        d.append( sol[index] )
      low  = flex.min(d)
      high = flex.max(d)
      if (high-low)<self.min_span:
        low  = low-0.5*self.min_span
        high = high+0.5*self.min_span
      s_result[ index ] = [low,high]
    self.span = s_result

  def span_as_scales(self, min_scale=0.1):
    result = flex.double()
    for lh in self.span:
      result.append( max(lh[1]-lh[0], min_scale) )
    return result








class adaptive_sa_optimizer(object):
  def __init__(self,
               parent,
               start_solution,
               scales,
               burn_in=20,
               burn_out=20,
               steps=150,
               coolfactor=0.96,
               max_step_in=30,
               spring_constant=0.95,
               reanneal=100,
               n_screen=100,
               eps=1e-3):
    self.eps = eps
    self.parent = parent
    self.scales = scales
    self.n = scales.size()

    self.sol = start_solution
    self.current_score = self.parent.target(  self.sol )
    self.best_sol = self.sol.deep_copy()
    self.best_score = float( self.current_score )
    self.last_best_sol = self.sol.deep_copy()

    self.spring_constant = spring_constant
    self.reanneal=reanneal
    self.n_screen = n_screen


    self.memory = None

    self.temp_scheme = []
    self.temp = 1.0e5
    self.start_t = None
    self.end_t = None
    self.burn_in = burn_in
    self.burn_out = burn_out
    self.steps = steps
    self.coolfactor = coolfactor
    self.max_step_in = max_step_in

    self.count = 0
    self.round = 0

    converged=False
    while not converged:
      converged = self.meta_run()
      self.round+=1
      if self.round == self.reanneal:
        converged=True

  def convergence_test(self):
    delta = flex.abs(self.last_best_sol - self.best_sol)
    delta = flex.max(delta)
    if delta> self.eps:
      return False
    return True

  def get_solution(self):
    return self.best_sol, self.best_score


  def meta_run(self, blow_up_factor=1.5):
    self.temp_scheme = []
    self.memory = trial_queue( n=self.n )
    self.make_anneal_profile(self.n_screen)
    self.run()
    self.scales = self.memory.span_as_scales()*blow_up_factor
    conv_test = self.convergence_test()
    self.last_best_sol =  self.best_sol.deep_copy()
    return conv_test

  def run(self):
    self.memory.add_solution( self.best_sol, self.best_score )
    #simple now really
    n = len(self.temp_scheme)
    for ii in xrange(n):
      self.temp=self.temp_scheme[ii]
      for jj in xrange(self.max_step_in):
        tmp = self.move()


  def make_anneal_profile(self,n_sample=25):
    ENE = flex.double()
    for ii in range(n_sample):
      new_vector = self.perturb( self.best_sol, False )
      new_score  = self.parent.target(new_vector)
      self.memory.add_solution(new_vector,new_score)
      if new_score < self.best_score:
        self.best_score = new_score
        self.best_sol = new_vector.deep_copy()

      ENE.append( self.parent.target( new_vector ) )
    dE = ( ENE - flex.mean(ENE) ).norm()/ math.sqrt(n_sample)
    self.memory.get_span()
    self.start_t = dE
    for ii in xrange(self.burn_in):
      self.temp_scheme.append( self.start_t )
    tmp_t = self.start_t

    count = 0
    while(count < self.steps):
      tmp_t = tmp_t * self.coolfactor
      self.temp_scheme.append( tmp_t )
      count += 1
    self.end_t = tmp_t
    for ii in xrange(self.burn_out):
      self.temp_scheme.append( self.end_t )


  def build_delta_vector(self, use_memory=True):
    delta_uniform = flex.random_double(self.n)
    if not use_memory:
      delta_uniform=delta_uniform*2.0-1.0
      return delta_uniform*self.scales

    span_vector = self.memory.span_as_scales()
    y = flex.double()
    for u in delta_uniform:
      this_sign = (u-0.5)
      if this_sign > 0:
        this_sign = 1.0
      else:
        this_sign = -1.0
      this_y = this_sign*self.temp*((1+1.0/self.temp)**(abs(u*2.0-1))-1)
      y.append( this_y )
    delta =  y * span_vector
    return delta

  def get_weight(self):
    x = self.temp
    w = ( 1-x*((1+1/x)**0.5-1.0)*2.0)*self.spring_constant
    return w

  def perturb(self, vector, use_memory=True ):
    delta = self.build_delta_vector( use_memory)
    new_vector = vector.deep_copy()
    w = self.get_weight()
    new_vector = w*self.best_sol + (1.0-w)*new_vector
    new_vector = new_vector + delta
    return new_vector

  def move(self):
    new_vector = self.perturb( self.sol, True)
    new_score = self.parent.target( new_vector )
    self.count+=1
    #print self.count, self.temp, self.best_score, list(self.best_sol), list(self.sol)
    if new_score <= self.current_score:
      self.sol = new_vector
      self.current_score = new_score
      if new_score < self.best_score:
        self.best_score = new_score
        self.best_sol = self.sol
        self.memory.add_solution( self.best_sol, self.best_score )
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
    self.scales = flex.double([5.5,5.5])
    self.x = flex.double( [3,-3] )
    self.optimizer = adaptive_sa_optimizer(self,
                                  self.x,
                                  self.scales,
                                  coolfactor=0.95,
                                  steps=250,
                                  max_step_in=10,
                                  burn_in=10,
                                  spring_constant=0.95,
                                  reanneal=5,
                                  n_screen=100)

    self.x, self.score = self.optimizer.get_solution()
    assert (self.score<1e-3)

  def target(self, vector):
    x1 = vector[0]
    x2 = vector[1]
    result = 100.0*(x2-x1*x1)**2.0 + (1-x1)**2.0
    return result

class test_rosenbrock_function(object):
  def __init__(self, dim=5):
    self.n = 2*dim
    self.x = flex.double( self.n, 2.0 )
    self.dim = dim
    self.scales = flex.double( self.n, 2.0 )

    self.optimizer = adaptive_sa_optimizer(self,
                                  self.x,
                                  self.scales,
                                  coolfactor=0.75,
                                  steps=100,
                                  max_step_in=30,
                                  burn_in=10,
                                  spring_constant=0.95,
                                  reanneal=1000,
                                  n_screen=max(self.n*50,100),
                                  eps=1e-5)
    self.x, self.score = self.optimizer.get_solution()
    for x in self.x:
      assert abs(x-1.0)<1e-2


  def target(self, vector):
    tmp = vector.deep_copy()
    x_vec = vector[0:self.dim]
    y_vec = vector[self.dim:]
    result=0
    for x,y in zip(x_vec,y_vec):
      result+=100.0*((y-x*x)**2.0) + (1-x)**2.0
    return result


def tst():
  flex.set_random_seed(11)
  tst_function()
  test_rosenbrock_function(2)
  print "OK"

if __name__ == """__main__""":
  tst()
