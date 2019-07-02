from __future__ import absolute_import, division, print_function
import math
from six.moves import range

def prime_sieve( n,
           include_unity=False
         ):
  numbers  = []
  marks = []
  primes = []
  for ii in range(n):
    numbers.append( ii+1 )
    marks.append( False )
  marks[0]=True
  k = 1
  if include_unity:
    primes.append( 1 )

  while k <=  math.sqrt(float(n)):
    m=None
    for ii in range(n):
      if m == None:
        if numbers[ii]>k:
          if not marks[ii]:
            m=numbers[ii]
            primes.append( m )

      if m != None:
        if numbers[ii]%m == 0:
          marks[ii]=True
    k = m

  for ii in range(n):
    if not marks[ii]:
      primes.append( numbers[ii] )

  return primes


def tst_prime_sieve():
  known = [2,3,5,7]
  test = prime_sieve( 10 )
  assert ( known ==  test )

if (__name__ == "__main__"):
  tst_prime_sieve()
  print("OK")
