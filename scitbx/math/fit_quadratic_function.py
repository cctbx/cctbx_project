from scitbx import lbfgs
from scitbx.array_family import flex


class fit_quadratic_function_2d_data(object):
  """ Fitting the function y=a(x1-x1m)^2 +b(x2-x2m)^2 + 2c(x1-x1m)(x2-x2m) """
  def __init__(self,
               x1_obs,
               x2_obs,
               y_obs):

    self.x1_obs = x1_obs
    self.x2_obs = x2_obs
    self.y_obs = y_obs
    self.x = flex.double([0,0,0,0,0])
    self.minimizer = lbfgs.run(target_evaluator=self)
    self.x1m=self.x[0]
    self.x2m=self.x[1]
    self.a=self.x[2]
    self.b=self.x[3]
    self.c=self.x[4]
    del self.x

  def calc_predict_function_values(self, observation_k):
    x = self.x1_obs[observation_k]-self.x[0]
    y = self.x2_obs[observation_k]-self.x[1]
    a = self.x[2]
    b = self.x[3]
    c = self.x[4]
    result = a*x*x+b*y*y+2.0*c*x*y
    return result

  def calc_f(self):
    result = 0
    for ii in range(len(self.x1_obs)):
      tmp = self.calc_predict_function_values(ii)
      result += (tmp-self.y_obs[ii])**2.0
    return(result)

  def calc_g(self):
    g = [0,0,0,0,0]
    xm=self.x[0]
    ym=self.x[1]
    a=self.x[2]
    b=self.x[3]
    c=self.x[4]
    for jj in range(len(self.x1_obs)): ## loop over all data points
      x = self.x1_obs[jj]
      y = self.x2_obs[jj]
      yobs=self.y_obs[jj]

      g[0] += 2.0*(-2.0*a*(x-xm)-2.0*c*(y-ym))*(a*(x-xm)**2.0 + 2.0*c*(x-xm)*(y-ym) \
                                                   +b*(y-ym)**2.0-yobs)

      g[1] += 2.0*(-2.0*c*(x-xm)-2.0*b*(y-ym))*(a*(x-xm)**2.0 + 2.0*c*(x-xm)*(y-ym) \
                                                   +b*(y-ym)**2.0-yobs)

      g[2] += 2.0*(x-xm)**2.0 * (a*(x-xm)**2.0 + 2.0*c*(x-xm)*(y-ym) + b*(y-ym)**2.0 -yobs)

      g[3] += 2.0*(y-ym)**2.0 * (a*(x-xm)**2.0 + 2.0*c*(x-xm)*(y-ym) + b*(y-ym)**2.0 -yobs)

      g[4] += 4.0*(y-ym)*(x-xm) * (a*(x-xm)**2.0 + 2.0*c*(x-xm)*(y-ym) + b*(y-ym)**2.0 -yobs)

    return(flex.double(g))

  def compute_functional_and_gradients(self):
    f = self.calc_f()
    g = self.calc_g()
    return f, g

## I hope to do something more general here at the appropriate juncture.
## The above implementation is valid, but not as general as desired, but is
## sufficient for the moment.
## Peter Zwart, Aug 2nd, 2005.
