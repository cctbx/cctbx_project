from __future__ import absolute_import, division, print_function
from six.moves import range
import math
import scitbx.math
from scitbx import lbfgs
from scitbx.array_family import flex

class fit_cdf(object):
  """
  =============================================================================
  Class fits a distribution according to its cumulative distribution function

  Arguments:
    x_data - measured property (list)
    y_data - fraction of points with measured property (cdf) (list)
    distribution - type of distribution to be fit

  Useful accessible attributes:
    self.x - the final parameters (list)
    self.distribution - the distribution being modeled (see notes)

  Notes:
    The available distributions are,
      gaussian
      rayleigh
      rice
    The distribution argument should not have quotes since it is used directly
    to access the class.
    The lbfgs method is used to minimize, f = (predicted - observed)^2
  -----------------------------------------------------------------------------
  """
  def __init__(self,x_data=None,y_data=None,distribution=None,
               max_iterations=50):
    # setup data
    assert(len(x_data) == len(y_data))
    self.n = len(x_data)
    self.x_data = flex.double(x_data)
    self.y_data = flex.double(y_data)
    if (flex.max(self.y_data) > 1.0):
      print("The cumulative distribution function (y_data) should only have values be between 0 and 1.")
      exit()

    # intialize distribution with guess
    self.distribution = distribution()
    self.distribution.estimate_parameters_from_cdf(x_data=self.x_data,
                                                   y_data=self.y_data)
    self.x = self.distribution.get_parameters()

    # optimize parameters
    self.minimizer = lbfgs.run(target_evaluator=self)

    # set optimized parameters
    self.distribution.set_parameters(p=self.x)

  def compute_functional_and_gradients(self):
    # caculate difference between predicted and observed values
    self.distribution.set_parameters(p=self.x)
    is_cpp_ = getattr(self.distribution,"interface","Python")=="C++"
    if is_cpp_:
      predicted = self.distribution.cdf(x=self.x_data)
    else:
      predicted = flex.double(self.n)
      for i in range(self.n):
        predicted[i] = self.distribution.cdf(x=self.x_data[i])
    difference = predicted - self.y_data

    # target function for minimization is sum of rmsd
    f = flex.sum(flex.sqrt(difference*difference))
    if is_cpp_:
      gradients = self.distribution.gradients(x=self.x_data, nparams=len(self.x), difference=difference)
      return f,gradients
    gradients = flex.double(len(self.x))
    for i in range(self.n):
      g_i = self.distribution.cdf_gradients(x=self.x_data[i])
      for j in range(len(self.x)):
        gradients[j] = gradients[j] + difference[i]*g_i[j]
    gradients = 2.0*gradients
    return f,gradients

# =============================================================================
class rayleigh(object):
  """
  =============================================================================
  Class models a 1-d Rayleigh distribution using one parameter, sigma.

              x                x^2
    pdf = --------- exp(- ------------)
           sigma^2         2 sigma^2

                        x^2
    cdf = 1 - exp(- -----------)
                     2 sigma^2

  The derivative of the cdf with respect to sigma is,

      d(cdf)          x^2               x^2             x
    ---------- = - --------- exp( - -----------) = - ------- pdf
     d(sigma)       sigma^3          2 sigma^2        sigma

  Methods:
    set_parameters
    get_parameters
    estimate_parameters_from_cdf
    pdf
    cdf
    d_cdf_d_sigma
    d_cdf_d_sigma_finite
    cdf_gradients
  -----------------------------------------------------------------------------
  """
  def __init__(self,mean=None,sigma=None):

    # parameter for Rayleigh distribution
    if (sigma is None):
      sigma = 1.0
    self.sigma = sigma

  def set_parameters(self,p=None):
    """
    Function sets the all the parameters
    """
    assert(len(p) == 1)
    self.sigma = float(p[0])

  def get_parameters(self):
    """
    Function returns all the parameters
    """
    return flex.double([self.sigma])

  def estimate_parameters_from_cdf(self,x_data=None,y_data=None):
    """
    Function estimates the parameter values based on the data (cdf)
    """
    # sigma is the mode of the distribution
    # approximate with the median (cdf = 0.5)
    midpoint = None
    for i in range(len(x_data)):
      if (y_data[i] > 0.5):
        midpoint = i
        break
    if (midpoint is None):
      midpoint = len(x_data) - 1
    self.sigma = x_data[midpoint]

  def pdf(self,x=None):
    """
    Function returns the probability density function at x
    """
    x_sigma = x/self.sigma
    f = (x_sigma/self.sigma)*math.exp(-0.5*x_sigma*x_sigma)
    return f

  def cdf(self,x=None):
    """
    Function returns the cumulative distribution function at x
    """
    x_sigma = x/self.sigma
    f = 1.0 - math.exp(-0.5*x_sigma*x_sigma)
    return f

  def inv_cdf(self,cdf=None):
    """
    Function returns the inverse cumulative distribution function at cdf
    """
    return math.sqrt(-2.*self.sigma*self.sigma*math.log(1.-cdf))

  def d_cdf_d_sigma(self,x=None):
    """
    Function returns the derivative of the cdf at x with respect to the
    standard deviation
    """
    df = -(x/self.sigma)*self.pdf(x=x)
    return df

  def d_cdf_d_sigma_finite(self,x=None,delta=0.00001):
    """
    Function returns the derivative of the cdf at x with respect to the
    standard deviation
    """
    sigma0 = self.sigma
    self.sigma = sigma0 + delta
    f1 = self.cdf(x=x)
    self.sigma = sigma0 - delta
    f2 = self.cdf(x=x)
    self.sigma = sigma0
    return (f1-f2)/(2.0*delta)

  def cdf_gradients(self,x=None):
    """
    Function returns a flex.double containing all derivatives
    """
    return flex.double([self.d_cdf_d_sigma(x)])
