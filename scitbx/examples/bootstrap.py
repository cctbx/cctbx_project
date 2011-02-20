import scitbx.lbfgs
import scitbx.math
from scitbx.array_family import flex
import math

"""
This example shows and easy way of obtaining reasonable estimates of
standard deviations of refinable parameters.

A sample run is shown below:

generated data:
x_obs  y_obs
0.0 1.89158676083
1.0 7.26573570715
2.0 17.4139845949
3.0 34.7247015887
4.0 57.2279161206
5.0 86.7943470166
6.0 121.658728013
7.0 163.13278756
8.0 209.591863552
9.0 262.234445997
10.0 321.838939765
11.0 386.669697396
12.0 457.737267767
13.0 534.242112317
14.0 618.481785819
15.0 706.341196433
16.0 801.835917893
17.0 903.488114864
18.0 1010.44961062
19.0 1122.7567835

final residual of fit on 'true' data
2.54570348084

Resulting fit and ESDs.

-------------------------------------------
       True and fitted coeffcients
-------------------------------------------
a 1.0 1.90288360993
b 2.0 1.92697354706
c 3.0 3.00469428576
-------------------------------------------
 Bootstrapped mean and standard deviations
-------------------------------------------
a 1.91734755932 0.225087363373
b 1.92587646901 0.0510724881016
c 3.00476495957 0.00291754050371


Cross-check with GNUPLOT fitting shows a good correspondence with bootstrap results:

=================================================================
After 5 iterations the fit converged.
final sum of squares of residuals : 2.5457
rel. change during last iteration : -4.58051e-07

degrees of freedom (ndf) : 17
rms of residuals      (stdfit) = sqrt(WSSR/ndf)      : 0.386972
variance of residuals (reduced chisquare) = WSSR/ndf : 0.149747

Final set of parameters            Asymptotic Standard Error
=======================            ==========================

a               = 1.90288          +/- 0.2356       (12.38%)
b               = 1.92697          +/- 0.05748      (2.983%)
c               = 3.00469          +/- 0.002921     (0.0972%)


correlation matrix of the fit parameters:

               a      b      c
a               1.000
b              -0.840  1.000
c               0.706 -0.965  1.000
=================================================================
"""




class polynomial_fit:
  def __init__(self, x_obs, y_obs, w_obs,n):
    assert x_obs.size() == y_obs.size()
    assert n < x_obs.size()
    self.x_obs = x_obs
    self.y_obs = y_obs
    self.w_obs = w_obs*w_obs*0+1.0
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
    print f
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

  for ii in range(20):
    print x_obs[ii], y_obs[ii]

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
