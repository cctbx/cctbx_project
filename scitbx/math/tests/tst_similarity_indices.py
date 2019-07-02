from __future__ import absolute_import, division, print_function
import scitbx.math

def exercise():
  x = [0.232, 0.002, 0.45, 1.2, 0.233, 1.2, 0.5, 0.231, 0,0,0,0,0]
  assert scitbx.math.similarity_indices(x=x, eps=0.01) == \
    [0, 1, 2, 3, 0, 3, 4, 0, 1, 1, 1, 1, 1]
  x = [0,2,0,2]
  assert scitbx.math.similarity_indices(x=x, eps=0.01) == [0, 1, 0, 1]
  x = [0,0,0,0]
  assert scitbx.math.similarity_indices(x=x, eps=0.01) == [0, 0, 0, 0]
  x = [1,2,3,4]
  assert scitbx.math.similarity_indices(x=x, eps=0.01) == [0, 1, 2, 3]

if (__name__ == "__main__"):
  exercise()
  print("OK")
