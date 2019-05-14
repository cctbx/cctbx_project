from __future__ import absolute_import, division, print_function
from six.moves import range
import math,random
from libtbx.test_utils import approx_equal
from scitbx import matrix

def generate_random_rotation(count):
  for x in range(count):

    u1=random.random()
    u2=random.random()
    u3=random.random()
    two_pi=2.0*math.pi

    #random unit quaternion, Steven M. LaValle, Planning Algorithms
    # Cambridge University Press, 2006, p. 198.
    # http://planning.cs.uiuc.edu
    q = [math.sqrt(1.-u1)*math.sin(two_pi*u2),
         math.sqrt(1.-u1)*math.cos(two_pi*u2),
         math.sqrt(u1)*math.sin(two_pi*u3),
         math.sqrt(u1)*math.cos(two_pi*u3),
        ]

    #use quaternion formula to find angle of rotation
    angle =2.* math.acos(q[0])

    #find vector of rotation
    denom = math.sin(angle/2.)
    vector = (q[1]/denom, q[2]/denom, q[3]/denom)
    assert approx_equal(matrix.col(q).length() , 1.0)
    assert approx_equal(matrix.col(vector).length(),1.0)
    yield vector, angle


if __name__=="__main__":
  for vector, angle in generate_random_rotation(50):
    V = matrix.col(vector)
    assert approx_equal(V.length(),1.0)
  print("OK")
