# James Farey was an English surveyor/geologist/writer (ca. 1816)
# he was fascinated by fractions and their decimal equivalents
# distances were measured with rods and chains and the poor man
# had to put up with fact a chain measured 4 rods or 66 feet
# http://en.wikipedia.org/wiki/Farey_series

def farey(v, lim):
    """
    determine closest fraction for floating point value v
    such that lim is the maximum denominator
    the function returns a (numerator, denominator) tuple
    """
    if v < 0:
        n, d = farey(-v, lim)
        return -n, d
    # gives a zero z of the same type as the lim argument
    z = lim - lim
    lower, upper = (z, z+1), (z+1, z)
    while True:
        mediant = (lower[0] + upper[0]), (lower[1] + upper[1])
        if v * mediant[1] > mediant[0]:
            if lim < mediant[1]:
                return upper
            lower = mediant
        elif v * mediant[1] == mediant[0]:
            if lim >= mediant[1]:
                return mediant
            if lower[1] < upper[1]:
                return lower
            return upper
        else:
            if lim < mediant[1]:
                return lower
            upper = mediant

def tst():
  import math as smath
  # farey() returns (numerator, denominator)
  # crude approximation of pi = 22/7.0
  n, d = farey(smath.pi, 100)
  assert n==22
  assert d==7

  # gives 'Chinese' approximation 355/113.0
  # 113355 is easy to remember
  n, d = farey(smath.pi, 1000)
  assert n == 355
  assert d == 113

  n, d = farey(smath.pi, 100000)
  assert n==312689
  assert d==99532


if __name__ == "__main__":
  tst()
