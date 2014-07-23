from __future__ import division
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.utils.ncs_utils as nu
import iotbx.ncs
from scitbx import matrix
import scitbx.rigid_body
from iotbx import pdb
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
    self.rotation1 = matrix.sqr(self.rot1.as_double())
    self.rotation2 = matrix.sqr(self.rot2.as_double())
    self.translation1 = matrix.rec((0.5,-0.5,0),(3,1))
    self.translation2 = matrix.rec((0,0,0),(3,1))


  def test_concatenate_rot_tran(self):
    """ Verify correct concatenation of rotation and translations """
    print 'Running ',sys._getframe().f_code.co_name
    pdb_obj = pdb.hierarchy.input(pdb_string=test_pdb_str)
    transforms_obj = iotbx.ncs.input(
      pdb_hierarchy_inp=pdb_obj,
      rotations=[self.rotation1,self.rotation2],
      translations=[self.translation1,self.translation2])
    results = nu.concatenate_rot_tran(transforms_obj)
    expected = flex.double([
      -0.40177529, 1.20019851, 2.64221706, 0.5, -0.5, 0.0,
      2.24044161,  1.57079633, 0.0,        0.0,  0.0, 0.0])
    assert approx_equal(results,expected,1.0e-4)
    s = 2.0
    results = nu.concatenate_rot_tran(transforms_obj,s=s)
    expected = flex.double([
      -0.40177529, 1.20019851, 2.64221706, 0.5/s, -0.5/s, 0.0,
      2.24044161,  1.57079633, 0.0,        0.0,  0.0, 0.0])
    assert approx_equal(results,expected,1.0e-4)

  def test_separate_rot_tran(self):
    """
    Verify correct conversion from angles and translation
    to rotation matrices and translations """
    print 'Running ',sys._getframe().f_code.co_name
    pdb_obj = pdb.hierarchy.input(pdb_string=test_pdb_str)
    transforms_obj = iotbx.ncs.input(
      pdb_hierarchy_inp=pdb_obj,
      rotations=[self.rotation1,self.rotation2],
      translations=[self.translation1,self.translation2])
    s = 2.0
    x = flex.double([
      -0.40177529, 1.20019851, 2.64221706, 0.5/s, -0.5/s, 0.0,
      2.24044161,  1.57079633, 0.0,        0.0,  0.0, 0.0])
    transforms_obj = nu.separate_rot_tran(
      x=x,transforms_obj=transforms_obj,s=s)
    rot_results, tran_results = nu.get_rotation_translation_as_list(
      transforms_obj=transforms_obj)
    rot_expected = [self.rotation1, self.rotation2]
    tran_expected = [self.translation1,self.translation2]
    assert approx_equal(tran_results,tran_expected,1.0e-4)
    assert approx_equal(rot_results,rot_expected,1.0e-4)


  def test_angles_to_matrix(self):
    """
    Verify derivation of rotation matrix R = Rx Ry Rz
    from alpha:
    rotation around x, beta: rotation around y, gamma: around z
    """
    print 'Running ',sys._getframe().f_code.co_name
    angles = self.rot_angles1
    expected = self.rot1.as_double()
    result = nu.angles_to_rotation(angles_xyz=angles,deg=False)
    assert approx_equal(expected,result,1e-4)
    # convert to Degrees
    angles = angles/math.pi*180
    result = nu.angles_to_rotation(angles_xyz=angles,deg=True)
    assert approx_equal(expected,result,1e-4)
    # test the rotations with sin(beta)==0
    angles = self.rot_angles2
    expected = self.rot2.as_double()
    result = nu.angles_to_rotation(angles_xyz=angles,deg=False)
    assert approx_equal(expected,result,1e-4)
    angles = self.rot_angles3
    expected = self.rot3.as_double()
    result = nu.angles_to_rotation(angles_xyz=angles,deg=False)
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
    angles = nu.rotation_to_angles(rotation=r, deg=False)
    assert approx_equal(expected_angles,angles,1e-3)
    # Test cos(beta)=0
    # sin(beta) = 1
    r = self.rot2.as_double()
    # when sin(beta) = 1 the (alpha + gamma) is the solution
    expected_angles_sum = self.rot_angles2[0] + self.rot_angles2[2]
    angles = nu.rotation_to_angles(rotation=r, deg=False)
    angles_sum = angles[0] + angles[2]
    assert approx_equal(expected_angles_sum,angles_sum,1e-3)
    # sin(beta) =  -1
    # when sin(beta) = -1 the (alpha - gamma) is the solution
    expected_angles_sum = self.rot_angles2[0] - self.rot_angles2[2]
    r = self.rot3.as_double()
    expected_angles = self.rot_angles3
    angles = nu.rotation_to_angles(rotation=r, deg=False)
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
    angles = nu.rotation_to_angles(rotation=r, deg=False)
    assert approx_equal(expected_angles,angles,1e-3)


  def test_update_x(self):
    """    Verify that transforms are getting updated    """
    print 'Running ',sys._getframe().f_code.co_name
    pdb_obj = pdb.hierarchy.input(pdb_string=test_pdb_str)
    transforms_obj = iotbx.ncs.input(
      pdb_hierarchy_inp=pdb_obj,
      rotations=[self.rotation1,self.rotation2],
      translations=[self.translation1,self.translation2])

    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
    x1 = nu.concatenate_rot_tran(transforms_obj)
    x2 = nu.shake_transformations(
      x = x1,
      shake_angles_sigma=0.035,
      shake_translation_sigma=0.5)
    # update with shaken parameters
    transforms_obj = nu.separate_rot_tran(
      x=x2, transforms_obj=transforms_obj)
    ncs_restraints_group_list = nu.separate_rot_tran(
      x=x2, ncs_restraints_group_list=ncs_restraints_group_list)

    x3 = nu.concatenate_rot_tran(
      ncs_restraints_group_list=ncs_restraints_group_list)
    x4 = nu.concatenate_rot_tran(transforms_obj=transforms_obj)
    assert abs(sum(list(x3-x4))) < 1.0e-3

    the,psi,phi =x2[6:9]
    rot = scitbx.rigid_body.rb_mat_xyz(
      the=the, psi=psi, phi=phi, deg=False)
    a3 = rot.rot_mat()

    the,psi,phi =x4[6:9]
    rot = scitbx.rigid_body.rb_mat_xyz(
      the=the, psi=psi, phi=phi, deg=False)
    a4 = rot.rot_mat()
    assert abs(sum(list(a3-a4))) < 1.0e-3

  def test_ncs_selection(self):
    """
    verify that extended_ncs_selection, which include the master ncs copy and
    the portion of the protein we want to refine.
    """
    print 'Running ',sys._getframe().f_code.co_name
    transforms_obj = iotbx.ncs.input(
      pdb_string=test_pdb_str,
      rotations=[self.rotation1,self.rotation2],
      translations=[self.translation1,self.translation2])
    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
    refine_selection = flex.size_t(range(30))
    result = nu.get_extended_ncs_selection(
      ncs_restraints_group_list=ncs_restraints_group_list,
      refine_selection=refine_selection)
    expected = [0, 1, 2, 3, 4, 5, 6, 21, 22, 23, 24, 25, 26, 27, 28, 29]
    assert list(result) == expected

  def test_ncs_related_selction(self):
    print 'Running ',sys._getframe().f_code.co_name
    transforms_obj = iotbx.ncs.input(
      pdb_string=test_pdb_str,
      rotations=[self.rotation1,self.rotation2],
      translations=[self.translation1,self.translation2])
    ncs_restraints_group_list = transforms_obj.get_ncs_restraints_group_list()
    result = nu.get_ncs_related_selection(
      ncs_restraints_group_list=ncs_restraints_group_list,
      asu_size=25)
    # Note that ASU length is set to 25 (not 21)
    expected = [True, True, True, True, True, True, True,
                True, True, True, True, True, True, True,
                True, True, True, True, True, True, True,
                False, False, False, False]
    assert list(result) == expected


test_pdb_str = '''\
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   1       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   1       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   1      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   1      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   1      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   1      10.523   7.209   9.055  1.00 20.00
'''

if __name__=='__main__':
  t = test_rotation_angles_conversion()
  # TODO: Consider adding more tests
  # run tests
  t.test_rotations_are_good()
  t.test_angles_to_matrix()
  t.test_matrix_to_angles()
  t.test_working_with_tuples()
  t.test_concatenate_rot_tran()
  t.test_separate_rot_tran()
  t.test_update_x()
  t.test_ncs_selection()
  t.test_ncs_related_selction()

