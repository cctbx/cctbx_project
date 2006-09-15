from scitbx.array_family import flex
class differential_evolution_optimizer(object):
  """
This is a python implementation of differential evolution
It assumes an evaluator class is passed in that has the following
functionality
data members:
 n              :: The number of parameters
 domain         :: a  list [(low,high)]*n
                   with approximate upper and lower limits for each parameter
 x              :: a place holder for a final solution

 also a function called 'target' is needed.
 This function should take a parameter vector as input and return a the function to be minimized.

 The code below was implemented on the basis of the following sources of information:
 1. http://www.icsi.berkeley.edu/~storn/code.html
 2. http://www.daimi.au.dk/~krink/fec05/articles/JV_ComparativeStudy_CEC04.pdf
 3. http://ocw.mit.edu/NR/rdonlyres/Sloan-School-of-Management/15-099Fall2003/A40397B9-E8FB-4B45-A41B-D1F69218901F/0/ses2_storn_price.pdf


 The developers of the differential evolution method have this advice:
 (taken from ref. 1)

 If you are going to optimize your own objective function with DE,
 try the following settings for the input file first: Choose method
 e.g. DE/rand/1/exp, set the number of parents NP to 10 times the
 number of parameters, select weighting factor F=0.8 and crossover
 constant CR=0.9. Make sure that you initialize your parameter vectors
 by exploiting their full numerical range, i.e. if a parameter is allowed
 to exhibit values in the range [-100, 100] it's a good idea to pick the
 initial values from this range instead of unnecessarily restricting diversity.
 If you experience misconvergence you usually have to increase the value for NP,
 but often you only have to adjust F to be a little lower or higher than 0.8.
 If you increase NP and simultaneously lower F a little, convergence is more likely
 to occur but generally takes longer, i.e. DE is getting more robust (there is always
 a convergence speed/robustness tradeoff).

 Note: NP is called population size in the routine below.)
  """

  def __init__(self,
               evaluator,
               population_size=50,
               f=0.8,
               cr=0.9,
               eps=1e-2,
               n_cross=1,
               max_iter=10000,
               monitor_cycle=200,
               out=None,
               show_progress=False,
               show_progress_nth_cycle=25):
    self.show_progress=show_progress
    self.show_progress_nth_cycle=show_progress_nth_cycle
    self.evaluator = evaluator
    self.population_size = population_size
    self.f = f
    self.cr = cr
    self.n_cross = n_cross
    self.max_iter = max_iter
    self.monitor_cycle = monitor_cycle
    self.vector_length = evaluator.n
    self.eps = eps
    self.population = []
    for ii in xrange(self.population_size):
      self.population.append( flex.double(self.vector_length,0) )
    self.scores = flex.double(self.population_size,1000)
    self.optimize()
    self.best_score = flex.min( self.scores )
    self.best_vector = self.population[ flex.min_index( self.scores ) ]
    self.evaluator.x = self.best_vector
    if self.show_progress:
      self.evaluator.print_status(
            flex.min(self.scores),
            flex.mean(self.scores),
            self.population[ flex.min_index( self.scores ) ],
            'Final')


  def optimize(self):
    # initialise the population please
    self.make_random_population()
    # score the population please
    self.score_population()
    converged = False
    monitor_score = flex.min( self.scores )
    count = 0
    while not converged:
      self.evolve()
      location = flex.min_index( self.scores )
      if self.show_progress:
        if count%self.show_progress_nth_cycle==0:
          # make here a call to a custom print_status function in the evaluator function
          # the function signature should be (min_target, mean_target, best vector)
          self.evaluator.print_status(
            flex.min(self.scores),
            flex.mean(self.scores),
            self.population[ flex.min_index( self.scores ) ],
            count)

      count += 1
      if count%self.monitor_cycle==0:
        if (monitor_score-flex.min(self.scores) ) < self.eps:
          converged = True
        else:
         monitor_score = flex.min(self.scores)
      rd = (flex.mean(self.scores) - flex.min(self.scores) )
      rd = rd*rd/(flex.min(self.scores)*flex.min(self.scores))
      if ( rd < self.eps ):
        converged = True


      if count>=self.max_iter:
        converged =True

  def make_random_population(self):
    for ii in xrange(self.vector_length):
      delta  = self.evaluator.domain[ii][1]-self.evaluator.domain[ii][0]
      offset = self.evaluator.domain[ii][0]
      random_values = flex.random_double(self.population_size)
      random_values = random_values*delta+offset
      # now please place these values ni the proper places in the
      # vectors of the population we generated
      for vector, item in zip(self.population,random_values):
        vector[ii] = item

  def score_population(self):
    for vector,ii in zip(self.population,xrange(self.population_size)):
      tmp_score = self.evaluator.target(vector)
      self.scores[ii]=tmp_score

  def evolve(self):
    for ii in xrange(self.population_size):
      rnd = flex.random_double(self.population_size-1)
      permut = flex.sort_permutation(rnd)
      # make parent indices
      i1=permut[0]
      if (i1>=ii):
        i1+=1
      i2=permut[1]
      if (i2>=ii):
        i2+=1
      i3=permut[2]
      if (i3>=ii):
        i3+=1
      #
      x1 = self.population[ i1 ]
      x2 = self.population[ i2 ]
      x3 = self.population[ i3 ]
      vi = x1 + self.f*(x2-x3)
      # prepare the offspring vector please
      rnd = flex.random_double(self.vector_length)
      permut = flex.sort_permutation(rnd)
      test_vector = self.population[ii].deep_copy()
      # first the parameters that sure cross over
      for jj in xrange( self.vector_length  ):
        if (jj<self.n_cross):
          test_vector[ permut[jj] ] = vi[ permut[jj] ]
        else:
          if (rnd[jj]>self.cr):
            test_vector[ permut[jj] ] = vi[ permut[jj] ]
      # get the score please
      test_score = self.evaluator.target( test_vector )
      # check if the score if lower
      if test_score < self.scores[ii] :
        self.scores[ii] = test_score
        self.population[ii] = test_vector


class test_function(object):
  def __init__(self):
    self.x = None
    self.n = 9
    self.domain = [ (-100,100) ]*self.n
    self.optimizer =  differential_evolution_optimizer(self,population_size=100,n_cross=5)
    assert flex.sum(self.x*self.x)<1e-9

  def target(self, vector):
    tmp = vector.deep_copy()
    result = (flex.sum(flex.cos(tmp*10))+self.n+1)*flex.sum( (tmp)*(tmp) )
    return result



def run():
  test_function()
  print "OK"


if __name__ == "__main__":
  run()
