from __future__ import division
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.utils.rotations as rt
from scitbx import matrix
import math
import sys


class test_rotation_angles_conversion(object):
  """
  Consider R = Rx(alpha)Ry(beta)Rz(gamma)

  Test the conversion of rotation angles to rotation matrix
  and the rotation matrix to angle
  """
  def __init__(self):
    # set_test_matrix
    self.rot1 = flex.vec3_double([
      (-0.317946, -0.173437, 0.932111),
      ( 0.760735, -0.633422, 0.141629),
      ( 0.565855,  0.754120, 0.333333)])
    self.rot2 = flex.vec3_double([
      (0       ,  0       , 1),
      (0.784042, -0.620708, 0),
      (0.620708,  0.784042, 0)])
    self.rot3 = flex.vec3_double([
      ( 0       ,  0       , -1),
      ( 0.097445, -0.995241,  0),
      (-0.995241, -0.097445,  0)])
    # Angles for rot, in radians
    self.rot_angles1 = flex.double(
      (-0.4017753, 1.2001985, 2.6422171))
    self.rot_angles2 = flex.double(
      (-0.4017753, math.pi/2, 2.6422171))
    self.rot_angles3 = flex.double(
      (-0.4017753, -math.pi/2, 2.6422171))

  def test_angles_to_matrix(self):
    """
    Verify derivation of rotation matrix R = Rx Ry Rz
    from alpha:
    rotation around x, beta: rotation around y, gamma: around z
    """
    print 'Running ',sys._getframe().f_code.co_name
    angles = self.rot_angles1
    expected = self.rot1.as_double()
    result = rt.angles_to_rotation(angles_xyz=angles,deg=False)
    assert approx_equal(expected,result,1e-4)
    # convert to Degrees
    angles = angles/math.pi*180
    result = rt.angles_to_rotation(angles_xyz=angles,deg=True)
    assert approx_equal(expected,result,1e-4)
    # test the rotations with sin(beta)==0
    angles = self.rot_angles2
    expected = self.rot2.as_double()
    result = rt.angles_to_rotation(angles_xyz=angles,deg=False)
    assert approx_equal(expected,result,1e-4)
    angles = self.rot_angles3
    expected = self.rot3.as_double()
    result = rt.angles_to_rotation(angles_xyz=angles,deg=False)
    assert approx_equal(expected,result,1e-4)


  def test_matrix_to_angles(self):
    """
    Note that there are two possible sets of angles for a rotation
    matrix.
    Also note that for the cases where cos(beta)=0, there is no unique
    answer
    """
    print 'Running ',sys._getframe().f_code.co_name
    r = self.rot1.as_double()
    expected_angles = self.rot_angles1
    angles = rt.rotation_to_angles(rotation=r, deg=False)
    assert approx_equal(expected_angles,angles,1e-3)
    # Test cos(beta)=0
    # sin(beta) = 1
    r = self.rot2.as_double()
    # when sin(beta) = 1 the (alpha + gamma) is the solution
    expected_angles_sum = self.rot_angles2[0] + self.rot_angles2[2]
    angles = rt.rotation_to_angles(rotation=r, deg=False)
    angles_sum = angles[0] + angles[2]
    assert approx_equal(expected_angles_sum,angles_sum,1e-3)
    # sin(beta) =  -1
    # when sin(beta) = -1 the (alpha - gamma) is the solution
    expected_angles_sum = self.rot_angles2[0] - self.rot_angles2[2]
    r = self.rot3.as_double()
    expected_angles = self.rot_angles3
    angles = rt.rotation_to_angles(rotation=r, deg=False)
    angles_sum = angles[0] - angles[2]
    assert approx_equal(expected_angles_sum,angles_sum,1e-3)

  def test_rotations_are_good(self):
    """
    Make sure that our rotation matrices are good
    """
    print 'Running ',sys._getframe().f_code.co_name
    for rm in [self.rot1,self.rot2,self.rot3]:
      r = matrix.sqr(rm.as_double())
      assert r.is_r3_rotation_matrix(rms_tolerance=1e-3)

  def test_working_with_tuples(self):
    """
    When working with scitbx matrix.rec or matrix.sqr
    (the form rotation matrices are in)
    the elements of those matrices are available as tuple.

    Verify that we process tuple well
    """
    print 'Running ',sys._getframe().f_code.co_name
    r = tuple(self.rot1.as_double())
    expected_angles = self.rot_angles1
    angles = rt.rotation_to_angles(rotation=r, deg=False)
    assert approx_equal(expected_angles,angles,1e-3)

if __name__=='__main__':
  t = test_rotation_angles_conversion()
  t.test_rotations_are_good()
  t.test_angles_to_matrix()
  t.test_matrix_to_angles()
  t.test_working_with_tuples()



