from scitbx.array_family import flex
from cma_es import cma_es

class cma_es_driver(object):
  """
  This object provides one with a easy interface to cma_es optimisation.
  For now, no options can be set, this will be added in the future.
  """
  def __init__(self, N, mean, sigma, evaluator):
    self.N = N
    self.x = mean
    self.sigma = sigma
    self.evaluator = evaluator

    self.optimizer = cma_es(self.N, self.x, self.sigma)
    self.count = 0

    while (not self.optimizer.converged() ):
      # get sample population
      p = self.optimizer.sample_population()
      pop_size = p.accessor().all()[0]
      # update objective function
      v = flex.double(pop_size)
      for i in xrange(pop_size):
        vector = p[(i*N):(i*N + N)]
        v[i] = self.evaluator( vector )
      self.optimizer.update_distribution(v)
    self.x_final = self.optimizer.get_result()
    self.score_final = self.evaluator( self.x_final )


def tst_it():
  def function(vector):
    x = vector[0]
    y = vector[1]
    result =100.0*((y-x*x)**2.0) + (1-x)**2.0
    return result

  m = flex.double( [5,5] )
  s = flex.double( [3,3] )
  obj = cma_es_driver( 2, m, s, function )
  assert abs(obj.x_final[0]-1)<1e-8
  assert abs(obj.x_final[1]-1)<1e-8




# =============================================================================
if (__name__ == '__main__'):
  tst_it()
  print 'Ok'
