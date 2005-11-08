import scitbx.lbfgs
import scitbx.math
from cctbx.array_family import flex
import math

"""  A simple polynomial fit """

class polynomial_fit:
  def __init__(self, x_obs, y_obs, w_obs,n):
    assert x_obs.size() == y_obs.size()
    assert n < x_obs.size()
    self.x_obs = x_obs
    self.y_obs = y_obs
    self.w_obs = w_obs*w_obs
    self.n = n
    self.x = flex.double(self.n,0)
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.a = self.x
    del self.x

  def compute_functional_and_gradients(self):
    self.a = self.x
    y_calc = flex.double(self.x_obs.size(),0)
    for i in range(self.n):
      y_calc = y_calc + (self.a[i])*flex.pow(self.x_obs,i)
    y_diff = self.y_obs - y_calc
    f = flex.sum(y_diff*y_diff/self.w_obs)
    g = flex.double(self.n,0)
    for i in range(self.n):
      g[i] = -flex.sum( 2.0*(y_diff/self.w_obs)*flex.pow(self.x_obs,i) )
    return f, g


class fake_data(object):
  def __init__(self,
               x_data,
               y_data):

    self.x_data = x_data
    self.y_data = y_data

    ## Make a permuation reference, this allows one to
    ## do non parametric resampling of multidimensional data
    self.permutation_reference = flex.double( range( x_data.size() ) )

    self.non_para_bootstrap = scitbx.math.non_parametric_bootstrap(
      self.permutation_reference, 0 )

  def fake_it(self, size):
    selection_set = self.non_para_bootstrap.draw( size )
    isel = flex.int()

    for element in selection_set:
      isel.append( int(element) )

    new_x = flex.double( flex.select(self.x_data, isel ) )
    new_y = flex.double( flex.select(self.y_data, isel ) )
    return new_x, new_y


def example():
  x_obs = flex.double( range(20) )
  a = flex.double([1,2,3])
  w_obs = flex.double(20,100.0)
  y_ideal = a[0] + a[1]*x_obs + a[2]*x_obs*x_obs
  y_obs = y_ideal + flex.random_double(size=x_obs.size())*1.5

  faker = fake_data(  x_obs,  y_obs)


  fit = polynomial_fit(x_obs=x_obs,y_obs=y_obs,w_obs=w_obs,n=3)
  print "------------------------------------------- "
  print "       True and fitted coeffcients"
  print "------------------------------------------- "
  for i in range(a.size()):
    print i, a[i], fit.a[i]
  print "------------------------------------------- "
  print " Bootstrapped mean and standard deviations"
  print "------------------------------------------- "
  mean=[0,0,0]
  std=[0,0,0]

  for trial in range(100):
    x_new, y_new =  faker.fake_it(20)
    fit = polynomial_fit(x_obs=x_new,y_obs=y_new,w_obs=w_obs,n=3)
    for i in range(a.size()):
      mean[i]+=fit.a[i]
      std[i]+=fit.a[i]*fit.a[i]

  for i in range(3):
    mean[i]/=100.0
    std[i]/=100.0
    std[i] -= mean[i]*mean[i]
    std[i] = math.sqrt( std[i] )
    print i, mean[i], std[i]


if (__name__ == "__main__"):
  example()
