from math import pi, cos, asin, sqrt
import pickle, StringIO
from cctbx.array_family import flex
from cctbx import uctbx
from cctbx.crystal_orientation import crystal_orientation
from scitbx import matrix
from libtbx.test_utils import approx_equal, not_approx_equal
import random
import math
import sys

def exercise_functions():
  orthorhombic = (1,0,0,0,0.5,0.,0.,0.,0.25)
  R = crystal_orientation(orthorhombic,True)
  A = R.rotate_thru((1,0,0,),math.pi/2.)
  assert approx_equal( A.direct_matrix(), (1,0, 0, 0,0,2,0,-4,0) )
  B = R.rotate_thru((0,1,0,),math.pi/2.)
  assert approx_equal( B.direct_matrix(), (0,0,-1, 0,2,0,4, 0,0) )
  C = R.rotate_thru((0,0,1,),math.pi/2.)
  assert approx_equal( C.direct_matrix(), (0,1, 0,-2,0,0,0, 0,4) )

def exercise_basic():
  identity = (1,0,0,0,1,0,0,0,1)
  I = crystal_orientation(identity,False) #direct space
  orthorhombic = (1,0,0,0,0.5,0.,0.,0.,0.25)
  R = crystal_orientation(orthorhombic,True) #reciprocal space
  assert R.direct_matrix() == (1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 4.0)
  assert R.reciprocal_matrix() == orthorhombic
  inverse = (-1,0,0,0,-1,0,0,0,-1)
  negative = crystal_orientation(inverse,False)
  negative.make_positive()
  assert I == negative
  assert R.unit_cell().parameters() == (1.0,2.0,4.0,90.,90.,90.)
  assert approx_equal( R.unit_cell_inverse().parameters(), (1.0,0.5,0.25,90.,90.,90.) )

O1 = crystal_orientation((-0.015553395334476732, -0.0028287158782335244, 0.018868416534039902,
                            -0.0016512962184570643, -0.020998220575299865, 0.0012056332661160732,
                            0.015789188025134133, -0.011166135863736777, 0.013045365404272641),True)
def exercise_change_basis():
  assert approx_equal(O1.unit_cell().parameters(),
                      (47.659, 47.6861, 49.6444, 62.9615, 73.8222, 73.5269),1E-3)
  reindex = (0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0) # swap a & b and take inverse
  O1.change_basis(reindex)
  assert approx_equal(O1.unit_cell().parameters(),
                      (47.6861, 47.659, 49.6444, 73.8222, 62.9615, 73.5269),1E-3)

def exercise_compare():
  pass

def exercise_pickle():
  p = pickle.dumps(O1)
  v = pickle.loads(p)
  assert O1 == v

def exercise_exceptions():
  pass

def run():
  exercise_functions()
  exercise_basic()
  exercise_change_basis()
  exercise_compare()
  exercise_pickle()
  exercise_exceptions()
  print "OK"

if (__name__ == "__main__"):
  run()
