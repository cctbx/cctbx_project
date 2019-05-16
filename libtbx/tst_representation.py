from __future__ import absolute_import, division, print_function

from six.moves import range
if (__name__ == "__main__"):
  for x in range(30):
    J = 7205759403792790 + x
    print(J * 10**30 // 2**56)

  T = 0.1 # Target
  # T ~= J / (2**N)
  # where J has exactly R bits
  # J ~= T * (2**N)
  print("BEst",T*(2**52))
  print("    ",T*(2**53))

  import math


  print(52 - math.log(T,2))
  print(53 - math.log(T,2))

  def precision_approx_equal(self,other,precision=24):
    # Use concepts from IEEE-754 to determine if the difference between
    # two numbers is within floating point error.  Not within scope to
    # do this for double precision literals; only interested in the case
    # where the data are from a ~single precision digital-analog converter.
    if self==other:
      return True
    if (self > 0.) != (other > 0.):
      return False
    #compute the exponent
    import math
    T = abs(self)
    Np = math.floor(precision-math.log(T,2))
    significand = int(T * 2**Np)
    val1 = significand/(2**Np) # nearest floating point representation of self
    val2 = (significand+1)/(2**Np) # next-nearest
    return abs(T-abs(other)) <= abs(val1-val2)
  print(precision_approx_equal(0.799999,0.800004,precision=17))
