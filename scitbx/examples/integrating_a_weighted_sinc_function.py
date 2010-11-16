from scitbx.array_family import flex
from scitbx.math import chebyshev_polynome
import scitbx.math
import math


class function(object):
  def __init__(self, n, m=100, k=2.5, d_max=45.0):
    self.n = n
    self.m = m
    self.k = k
    self.d_max=d_max
    self.x = 1.0-2.0*(flex.double(range(m+1))/m)
    self.r = 0.5*(1+self.x)*self.d_max
    self.r[0] = 1e-8

    self.coefs = (flex.random_double(self.n)-0.5)*0.0
    self.load_coefs()
    self.polynome = chebyshev_polynome(self.n, -1.0, +1.0, self.coefs)

  def show(self):
    result = get_p_of_r(self.x)
    for r,y in zip(self.r, result):
      print r, y

  def load_coefs(self, coefs=None):
    if coefs is None:
      self.coefs = (flex.random_double(self.n)-0.5)*2.0
    else:
      assert len(coefs)==self.n
      self.coefs = coefs
    # no means to refresh the coefficients yet in an elegant manner
    self.polynome = chebyshev_polynome(self.n, -1.0, +1.0, self.coefs)


  def get_p_of_r(self,x):
    base = flex.pow((1.0-x*x),self.k)
    exp_pol = flex.exp( self.polynome.f( x ) )
    result = exp_pol*base
    return result

  def get_sinc(self, q, x ):
    r = 0.5*(x+1)*self.d_max
    sinc = flex.sin( r*q )/(r*q)
    return sinc

  def integrate(self, q, ni):
    gle = scitbx.math.gauss_legendre_engine(ni)
    x_int = gle.x()
    w_int = gle.w()
    p_of_r = self.get_p_of_r(x_int)
    sinc = self.get_sinc( q, x_int )
    tbi = p_of_r*sinc
    wtbi = tbi*w_int
    result = flex.sum(wtbi)
    return result

def example():
  f = function(5,100)
  q_trials = [0.001, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

  ref_integrals = []
  for q in q_trials:
    ref_integrals.append( f.integrate(q,90) )

  for q, jj in zip(q_trials, range(len(q_trials))):
    print q,jj,
    for ii in range(2,90):
      print 100.0*abs(f.integrate(q,ii)-ref_integrals[jj])/abs(ref_integrals[jj]+1e-13),
    print





if (__name__ == "__main__"):
  example()
