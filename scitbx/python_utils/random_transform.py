## This file conatins some usefull transformations
## for the generation of specific random variables
import random
import math
from scitbx.array_family import flex


## This function generates approximately
## 1 million Gaussian variates in about 0.35 seconds
## on longnose.
##
## Theorectical (central) moments are (with mu=0 and sigma=1)
## given below, as well as the result from sampling, in paranthesis.
## mu1 = 0 (0.00041)
## mu2 = 1 (0.99954)
## mu3 = 0 (6.08714e-05)
## mu4 = 3 (2.99848)
## It behaves as expected.
## Similar results are obtained when changing mu and sigma.
def normal_variate(mu=0.0,sigma=1.0,N=100):
  "Normal variate via Box-Muller transform"
  U1 = flex.random_double(size=N)
  U2 = flex.random_double(size=N)
  return flex.sqrt(-2.0*flex.log(U1))*flex.cos(2.0*math.pi*U2)*sigma+mu




## This function generates approximately
## 1 million t variates in about 0.6 seconds
## on longnose.
##
## for something_very_big < a,
## the t-vbaraites will go to Gaussian variates
## Note that the variance is undefined if a<=2
## for a=6, mu=0, sigma=1, the results are
## mu1 = 0    ( 0.00055754)
## mu2 = 1.5  ( 1.49899419)
## mu3 = 0    (-0.01446477)
## mu3 = 13.5 (13.29293869)
def t_variate(a=1.0,mu=0.0,sigma=1.0,N=100):
  "T-variate via Baley's one-liner"
  U1 = flex.random_double(size=N)
  U2 = flex.random_double(size=N)
  return ( flex.sqrt(a*(flex.pow(U1,-2.0/a)-1.0))
           *flex.cos(2.0*math.pi*U2)*sigma+mu )




## This function generates approximately
## 1 million normalised-Wilson (amplitudes) variates
## per 0.78 seconds.
##
## Expected and simulated raw moments:
## mu1 : 0.886  (0.88592)
## mu2 : 1      (0.99946)
## mu4 : 2      (1.99826)

def wilson_amplitude_variate(N=100):
  "Wilson amplitude variate; The Rayleigh distribution"
  ## Get wilson variate via two gaussians with half a variance
  A = normal_variate(mu=0,sigma=math.sqrt(0.5),N=N)
  B = normal_variate(mu=0,sigma=math.sqrt(0.5),N=N)
  return flex.sqrt(A*A+B*B)

# As before, this takes about 0.81 seconds for 1 million r.v.'s
def wilson_intensity_variate(N=100):
  "Wilson intensity variate; The exponetial distribution"
  return flex.pow2(wilson_amplitude_variate(N=N))

# This random variate is distributed as the (error free) |I+-I-|
# normalised by 4*( (sum f") * (sum  f_nought) )^(1/2)
def pseudo_normalized_abs_delta_i(N=100):
  x = flex.random_double(size=N)
  x = -0.5*flex.log( 1.0-x )
  return(x)
