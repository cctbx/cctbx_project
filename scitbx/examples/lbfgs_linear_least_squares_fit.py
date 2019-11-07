from __future__ import absolute_import, division, print_function
import scitbx.lbfgs
from scitbx.array_family import flex
import math
from six.moves import zip

class linear_least_squares_fit(object):

  def __init__(self, x_obs, y_obs):
    self.x_obs = x_obs
    self.y_obs = y_obs
    self.x = flex.double([1, 0]) # start with slope=1, y_intercept=0
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.slope = self.x[0]
    self.y_intercept = self.x[1]
    del self.x

  def compute_functional_and_gradients(self):
    slope = self.x[0]
    y_intercept = self.x[1]
    y_calc = slope * self.x_obs + y_intercept
    y_diff = self.y_obs - y_calc
    f = flex.sum(flex.pow2(y_diff))
    g = flex.double([
      flex.sum(-2 * y_diff * self.x_obs),
      flex.sum(-2 * y_diff)])
    return f, g

def example():
  x_obs = flex.double([1,2,3,4,5,6,7,8,9,10])
  slope = -math.pi
  y_intercept = math.sqrt(3)
  y_ideal = slope * x_obs + y_intercept
  y_obs = y_ideal + flex.random_double(size=x_obs.size())*0.1
  fit = linear_least_squares_fit(x_obs=x_obs, y_obs=y_obs)
  print("fit.slope:", fit.slope)
  print("fit.y_intercept:", fit.y_intercept)
  y_calc = fit.slope * x_obs + fit.y_intercept
  print(" x_obs  y_obs y_calc  diff")
  for xo,yo,yc in zip(x_obs, y_obs, y_calc):
    print("%6.2f %6.2f %6.2f %6.2f" % (xo,yo,yc,yo-yc))

if (__name__ == "__main__"):
  example()
