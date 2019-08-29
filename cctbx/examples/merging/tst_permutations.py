"""A completely standalone example that explains what a permutation
matrix is and how it is applied in practice.
"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
def stuff_about_permutations():
    trial = sqr((1,2,3,4,5,6,7,8,9))
    for x in range(3):
      for y in range(3):
        print(trial(x,y), end=' ')
      print()

    ordering = col((3,1,2))
    matcode = flex.int(9)
    for islow in range(3):
      matcode[3*islow + ordering[islow]-1] = 1
    matcode = sqr(matcode)
    print("matcode")
    for x in range(3):
      for y in range(3):
        print(matcode(x,y), end=' ')
      print()
    prod = matcode * trial *(matcode.inverse())
    print("product",prod.elems)
    for x in range(3):
      for y in range(3):
        print(prod(x,y), end=' ')
      print()

stuff_about_permutations()
print("OK")
