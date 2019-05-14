from __future__ import absolute_import, division, print_function
import math,copy
from scitbx.array_family import flex
class SALight:
  """Lightweight implementation of a simulated annealing protocol"""

  def random_increment(self):
    random = 2.0*flex.random_double(len(self.x)) - 1.0 # values from -1 to 1
    sum_sq = flex.sum(random*random)
    normalized = random / math.sqrt(sum_sq)
    return normalized * self.L

  def cycles(self,k=600):
    self.ktotal = k
    self.counter = copy.copy(self.ktotal)
    while (self.counter > 0):
      self.counter -=1
      decreasing_increment = (self.counter/self.ktotal)*self.random_increment()
      #print self.format%tuple(decreasing_increment)
      self.x += decreasing_increment
      print(self.format%tuple(self.x))

if __name__=="__main__":
  SA = SALight()
  # Half mosaicity in degrees
  # Mid-wavelength adjustment factor
  # Bandpass fractional full width
  # adjustment angle in degrees
  # adjustment angle in degrees
  # adjustment angle in degrees

  # starting values; likely expected values
  SA.x = flex.double([0.1,1.00,0.006,0.0,0.0,0.0])
  SA.initial = SA.x.deep_copy()

  # reasonable length scale (expected interval, half width)
  SA.L = flex.double([0.04,0.001,0.002,0.1,0.1,0.1])

  SA.format = "Mosaicity %6.3f Wave mean %7.4f bandpass %7.4f Angles %8.5f %8.5f %8.5f"
  #print list(SA.random_increment())
  print(SA.format%tuple(SA.initial))
  SA.cycles()
  print()
  print(SA.format%tuple(SA.initial))
