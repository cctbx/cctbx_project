from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from six.moves import range

def demo():
  # toy array
  a = flex.double([10,11,12])
  # one way of looping over the array
  for i in range(a.size()):
    print(a[i])
  print()
  # a better way of looping over the array
  for ai in a:
    print(ai)
  print()
  # another good way of looping over the array
  for i,ai in enumerate(a):
    print(i, ai)
  print()
  # modify the elements one-by-one
  for i in range(a.size()):
    a[i] *= 10
  for ai in a:
    print(ai)
  print()
  # a better way of modifying all elements
  a += 100 # this works at C++ speed
  for ai in a:
    print(ai)
  print()
  print("OK")

if (__name__ == "__main__"):
  demo()
