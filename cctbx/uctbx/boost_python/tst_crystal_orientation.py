from math import pi
import pickle
from cctbx import uctbx, sgtbx
from cctbx.crystal_orientation import crystal_orientation, basis_type
from libtbx.test_utils import approx_equal
import math

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
  assert I == negative.make_positive()
  assert R.unit_cell().parameters() == (1.0,2.0,4.0,90.,90.,90.)
  assert approx_equal( R.unit_cell_inverse().parameters(), (1.0,0.5,0.25,90.,90.,90.) )

O1 = crystal_orientation((-0.015553395334476732, -0.0028287158782335244, 0.018868416534039902,
                            -0.0016512962184570643, -0.020998220575299865, 0.0012056332661160732,
                            0.015789188025134133, -0.011166135863736777, 0.013045365404272641),True)
def exercise_change_basis():
  assert approx_equal(O1.unit_cell().parameters(),
                      (47.659, 47.6861, 49.6444, 62.9615, 73.8222, 73.5269),1E-3)
  reindex = (0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0) # swap a & b and take inverse
  O2 = O1.change_basis(reindex)
  assert approx_equal(O2.unit_cell().parameters(),
                      (47.6861, 47.659, 49.6444, 73.8222, 62.9615, 73.5269),1E-3)

  rhombohedral_test = crystal_orientation((
    0.002737747939994224, -0.0049133768326561278, 0.0023634556852316566,
    0.0062204242383498082, 0.006107332242442573, 0.0047036576949967112,
    -0.0057640198753891566, -0.0025891042237382953, 0.0023674924674260264),basis_type.reciprocal)
  rhombohedral_reference = crystal_orientation((
    -0.0076361080646872997, 0.0049061665572297979, 0.0023688116121433865,
    -0.00011109895272056645, -0.0061110173438898583, 0.0047062769738302939,
     0.0031790372319626674, 0.0025876279220667518, 0.0023669727051432361),basis_type.reciprocal)
  # Find a similarity transform that maps the two cells onto each other
  c_inv_r_best = rhombohedral_test.best_similarity_transformation(
      other = rhombohedral_reference, fractional_length_tolerance = 1.00,
      unimodular_generator_range=1)
  c_inv_r_int = tuple([int(round(ij,0)) for ij in c_inv_r_best])
  assert c_inv_r_int == (-1, 0, 0, 1, -1, 0, 0, 0, 1)
  c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r_int))
  cb_op = sgtbx.change_of_basis_op(c_inv)
  rhombohedral_reindex = rhombohedral_test.change_basis(cb_op)
  assert rhombohedral_reindex.difference_Z_score(rhombohedral_reference) < 0.40
  assert rhombohedral_reindex.direct_mean_square_difference(rhombohedral_reference) < 0.1

  #an alternative test from ana that should fail (gives high msd~0.22; cell axes don't match):
  ana_reference = crystal_orientation((
    0.0023650364919947241, 0.012819317075171401, 0.003042762222847376,
    0.0081242553464681254, 0.0050052660206998077, -0.01472465697193685,
   -0.01373896574061278, 0.0083781530252581681, -0.0035301340829149005),basis_type.reciprocal)
  ana_current = crystal_orientation((
   -0.014470153848927263, 0.0095185368147633793, 0.00087746490483763798,
   -0.0049989006573928716, -0.0079714727432991222, 0.014778692772530192,
    0.0051268914129933571, 0.010264066188909109, 0.0044244589492769002),basis_type.reciprocal)
  c_inv_r_best = ana_current.best_similarity_transformation(
      other = ana_reference,
      fractional_length_tolerance = 200.0,
      unimodular_generator_range=1)
  c_inv_r_int = tuple([int(round(ij,0)) for ij in c_inv_r_best])
  c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r_int))
  cb_op = sgtbx.change_of_basis_op(c_inv)
  ana_reindex = ana_reference.change_basis(cb_op.inverse())
  assert 200.0 > ana_reindex.difference_Z_score(ana_current) > 20.

  u = uctbx.unit_cell((10., 10., 10., 90., 90., 90.))
  CO = crystal_orientation(u.fractionalization_matrix(),True)
  assert approx_equal(
    CO.unit_cell().parameters(),
    CO.change_basis((1,0,0, 0,1,0, 0,0,1)).unit_cell().parameters())
  u = uctbx.unit_cell((2,3,5))
  CO = crystal_orientation(u.fractionalization_matrix(),True)
  assert approx_equal(
    CO.change_basis((0,1,0, 0,0,1, 1,0,0)).unit_cell().parameters(),
    (5,2,3,90,90,90))
  cb_op = sgtbx.change_of_basis_op("y,z,x")
  assert approx_equal(
    CO.change_basis(cb_op).unit_cell().parameters(),
    (5,2,3,90,90,90))

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
