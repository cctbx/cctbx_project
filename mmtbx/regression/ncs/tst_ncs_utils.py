from __future__ import division
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.ncs.ncs_utils as nu
import mmtbx.maps.correlation
from scitbx import matrix
import scitbx.rigid_body
from iotbx import pdb
import iotbx.ncs as ncs
import unittest
import math

__author__ = 'Youval'

test_pdb_str = '''\
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   1       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   1       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   1      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   1      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   1      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   1      10.523   7.209   9.055  1.00 20.00
'''

pdb_answer_0 = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1      11.782  12.419   4.645  1.00 10.00           N
ATOM      2  CA  THR A   1      11.671  11.061   4.125  1.00 10.00           C
ATOM      3  C   THR A   1      11.746  10.033   5.249  1.00 10.00           C
ATOM      4  O   THR A   1      12.561  10.157   6.163  1.00 10.00           O
ATOM      5  CB  THR A   1      12.772  10.760   3.092  1.00 10.00           C
ATOM      6  OG1 THR A   1      12.672  11.682   2.000  1.00 10.00           O
ATOM      7  CG2 THR A   1      12.635   9.340   2.565  1.00 10.00           C
TER
ATOM      1  N   THR B   1       6.768   9.093   9.237  1.00 10.00           N
ATOM      2  CA  THR B   1       7.284   8.654   7.945  1.00 10.00           C
ATOM      3  C   THR B   1       8.638   7.968   8.097  1.00 10.00           C
ATOM      4  O   THR B   1       9.495   8.426   8.852  1.00 10.00           O
ATOM      5  CB  THR B   1       7.423   9.832   6.963  1.00 10.00           C
ATOM      6  OG1 THR B   1       6.144  10.446   6.765  1.00 10.00           O
ATOM      7  CG2 THR B   1       7.962   9.350   5.625  1.00 10.00           C
TER
ATOM      1  N   THR C   1       9.093   2.000  10.493  1.00 10.00           N
ATOM      2  CA  THR C   1       8.879   2.702   9.233  1.00 10.00           C
ATOM      3  C   THR C   1      10.081   3.570   8.875  1.00 10.00           C
ATOM      4  O   THR C   1      10.652   4.241   9.734  1.00 10.00           O
ATOM      5  CB  THR C   1       7.618   3.584   9.284  1.00 10.00           C
ATOM      6  OG1 THR C   1       6.472   2.770   9.559  1.00 10.00           O
ATOM      7  CG2 THR C   1       7.417   4.305   7.960  1.00 10.00           C
END
"""

pdb_poor_0 = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1      11.120  12.388   4.399  1.00 10.00           N
ATOM      2  CA  THR A   1      11.623  11.024   4.280  1.00 10.00           C
ATOM      3  C   THR A   1      12.335  10.587   5.556  1.00 10.00           C
ATOM      4  O   THR A   1      13.094  11.353   6.149  1.00 10.00           O
ATOM      5  CB  THR A   1      12.588  10.881   3.089  1.00 10.00           C
ATOM      6  OG1 THR A   1      11.911  11.230   1.876  1.00 10.00           O
ATOM      7  CG2 THR A   1      13.099   9.452   2.986  1.00 10.00           C
TER
ATOM      8  N   THR B   1       7.221   8.778   9.110  1.00 10.00           N
ATOM      9  CA  THR B   1       7.661   8.575   7.734  1.00 10.00           C
ATOM     10  C   THR B   1       9.058   7.965   7.687  1.00 10.00           C
ATOM     11  O   THR B   1       9.944   8.359   8.445  1.00 10.00           O
ATOM     12  CB  THR B   1       7.660   9.895   6.941  1.00 10.00           C
ATOM     13  OG1 THR B   1       6.338  10.446   6.927  1.00 10.00           O
ATOM     14  CG2 THR B   1       8.122   9.659   5.511  1.00 10.00           C
TER
ATOM     15  N   THR C   1       8.639   1.605   9.684  1.00 10.00           N
ATOM     16  CA  THR C   1       8.751   2.735   8.769  1.00 10.00           C
ATOM     17  C   THR C   1      10.144   3.354   8.827  1.00 10.00           C
ATOM     18  O   THR C   1      10.716   3.520   9.904  1.00 10.00           O
ATOM     19  CB  THR C   1       7.704   3.820   9.080  1.00 10.00           C
ATOM     20  OG1 THR C   1       6.388   3.265   8.969  1.00 10.00           O
ATOM     21  CG2 THR C   1       7.842   4.986   8.114  1.00 10.00           C
END
"""

test_pdb_str_2 = '''\
ATOM     45  N   PHEAa   6     219.693 144.930 112.416  1.00 50.00           N
ATOM     46  CA  PHEAa   6     218.871 146.020 112.886  1.00 50.00           C
ATOM     47  C   PHEAa   6     217.413 145.628 112.926  1.00 50.00           C
ATOM     48  O   PHEAa   6     216.730 145.905 113.908  1.00 50.00           O
TER
ATOM   1244  N   ARGAb   6     303.367 160.705 103.239  1.00 50.00           N
ATOM   1245  CA  ARGAb   6     302.396 160.991 104.331  1.00 50.00           C
ATOM   1246  C   ARGAb   6     302.285 162.473 104.586  1.00 50.00           C
ATOM   1247  O   ARGAb   6     302.837 163.292 103.851  1.00 50.00           O
TER
ATOM   2754  N   PHEAc   6     242.472 151.067 115.352  1.00 50.00           N
ATOM   2755  CA  PHEAc   6     241.314 151.789 115.823  1.00 50.00           C
ATOM   2756  C   PHEAc   6     240.094 150.900 115.864  1.00 50.00           C
ATOM   2757  O   PHEAc   6     239.358 150.912 116.847  1.00 50.00           O
TER
ATOM   3953  N   ARGAd   6     314.882 195.854 106.123  1.00 50.00           N
ATOM   3954  CA  ARGAd   6     313.875 195.773 107.215  1.00 50.00           C
ATOM   3955  C   ARGAd   6     313.239 197.116 107.471  1.00 50.00           C
ATOM   3956  O   ARGAd   6     313.460 198.078 106.736  1.00 50.00           O
TER
ATOM   5463  N   PHEAe   6     261.525 165.024 118.275  1.00 50.00           N
ATOM   5464  CA  PHEAe   6     260.185 165.283 118.746  1.00 50.00           C
ATOM   5465  C   PHEAe   6     259.365 164.014 118.787  1.00 50.00           C
ATOM   5466  O   PHEAe   6     258.673 163.762 119.769  1.00 50.00           O
TER
ATOM   6662  N   ARGAf   6     313.035 232.818 109.051  1.00 50.00           N
ATOM   6663  CA  ARGAf   6     312.124 232.379 110.143  1.00 50.00           C
ATOM   6664  C   ARGAf   6     311.048 233.405 110.399  1.00 50.00           C
ATOM   6665  O   ARGAf   6     310.909 234.383 109.665  1.00 50.00           O
END
'''

phil_str = '''\
ncs_group {
  master_selection = chain 'Aa'
  copy_selection = chain 'Ac'
  copy_selection = chain 'Ae'
}

ncs_group {
  master_selection = chain 'Ab'
  copy_selection = chain 'Ad'
  copy_selection = chain 'Af'
}
'''

class Test_ncs_utils(unittest.TestCase):
  """
  Consider R = Rx(alpha)Ry(beta)Rz(gamma)

  Test the conversion of rotation angles to rotation matrix
  and the rotation matrix to angle
  """
  def setUp(self):
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
    self.rot_angles1_deg = flex.double(
      (-0.4017753, 1.2001985, 2.6422171)) * 180/math.pi
    self.rotation1 = matrix.sqr(self.rot1.as_double())
    self.rotation2 = matrix.sqr(self.rot2.as_double())
    self.rotation3 = matrix.sqr(self.rot3.as_double())
    self.translation1 = matrix.rec((0.5,-0.5,0),(3,1))
    self.translation2 = matrix.rec((0,0,0),(3,1))
    self.translation3 = matrix.rec((0,1,2),(3,1))
    self.r_t = [[self.rotation1, self.translation1],
                [self.rotation2, self.translation2],
                [self.rotation3, self.translation3]]

    self.pdb_obj = pdb.hierarchy.input(pdb_string=test_pdb_str)
    self.tr_obj1 = ncs.input(
      pdb_hierarchy_inp=self.pdb_obj,
      rotations=[self.rotation1,self.rotation2],
      translations=[self.translation1,self.translation2])
    self.tr_obj2 = ncs.input(
      pdb_hierarchy_inp=self.pdb_obj,
      rotations=[self.rotation1,self.rotation2,self.rotation3],
      translations=[self.translation1,self.translation2,self.translation3])
    self.ncs_restraints_group_list = \
      self.tr_obj1.get_ncs_restraints_group_list()

  def test_concatenate_rot_tran(self):
    """ Verify correct concatenation of rotation and translations """
    # print sys._getframe().f_code.co_name
    results = nu.concatenate_rot_tran(self.tr_obj1)
    expected = flex.double([
      -0.40177529, 1.20019851, 2.64221706, 0.5, -0.5, 0.0,
      2.24044161,  1.57079633, 0.0,        0.0,  0.0, 0.0])
    assert approx_equal(results,expected,1.0e-4)

  def test_update_rot_tran(self):
    """
    Verify correct conversion from angles and translation
    to rotation matrices and translations """
    # print sys._getframe().f_code.co_name
    x = flex.double([
      -0.40177529, 1.20019851, 2.64221706, 0.5, -0.5, 0.0,
      2.24044161,  1.57079633, 0.0,        0.0,  0.0, 0.0])
    self.tr_obj1 = nu.update_rot_tran(
      x=x,transforms_obj=self.tr_obj1)
    rot_results, tran_results = nu.get_rotation_translation_as_list(
      transforms_obj=self.tr_obj1)
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
    # print sys._getframe().f_code.co_name
    angles = self.rot_angles1
    expected = self.rot1.as_double()
    result = nu.angles_to_rotation(angles_xyz=angles,deg=False)
    assert approx_equal(expected,result,1e-4)
    # convert to Degrees
    angles = self.rot_angles1/math.pi*180
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
    # print sys._getframe().f_code.co_name
    r = self.rot1.as_double()
    expected_angles = self.rot_angles1
    angles = nu.rotation_to_angles(rotation=r, deg=False)
    assert approx_equal(expected_angles,angles,1e-3)
    expected_angles = self.rot_angles1_deg
    angles = nu.rotation_to_angles(rotation=r, deg=True)
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
    angles = nu.rotation_to_angles(rotation=r, deg=False)
    angles_sum = angles[0] - angles[2]
    assert approx_equal(expected_angles_sum,angles_sum,1e-3)

  def test_rotations_are_good(self):
    """
    Make sure that our rotation matrices are good
    """
    # print sys._getframe().f_code.co_name
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
    # print sys._getframe().f_code.co_name
    r = tuple(self.rot1.as_double())
    expected_angles = self.rot_angles1
    angles = nu.rotation_to_angles(rotation=r, deg=False)
    assert approx_equal(expected_angles,angles,1e-3)


  def test_update_x(self):
    """    Verify that transforms are getting updated    """
    # print sys._getframe().f_code.co_name
    x1 = nu.concatenate_rot_tran(self.tr_obj1)
    x2 = nu.shake_transformations(
      x = x1,
      shake_angles_sigma=0.035,
      shake_translation_sigma=0.5)
    # update with shaken parameters
    self.tr_obj1 = nu.update_rot_tran(
      x=x2, transforms_obj=self.tr_obj1)
    self.ncs_restraints_group_list = nu.update_rot_tran(
      x=x2, ncs_restraints_group_list=self.ncs_restraints_group_list)

    x3 = nu.concatenate_rot_tran(
      ncs_restraints_group_list=self.ncs_restraints_group_list)
    x4 = nu.concatenate_rot_tran(
      transforms_obj=self.tr_obj1)
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

    # test that update_rot_tran dose not unintentionally change x
    round_val = 3
    r_elems = []
    for rec in self.ncs_restraints_group_list:
      for c in rec.copies:
        r = c.r.round(round_val)
        r_elems.append(r.elems)

  def test_ncs_selection(self):
    """
    verify that extended_ncs_selection, which include the master ncs copy and
    the portion of the protein we want to refine.
    """
    # print sys._getframe().f_code.co_name
    refine_selection = flex.size_t(range(30))
    result = nu.get_extended_ncs_selection(
      ncs_restraints_group_list=self.ncs_restraints_group_list,
      refine_selection=refine_selection)
    expected = [0, 1, 2, 3, 4, 5, 6, 21, 22, 23, 24, 25, 26, 27, 28, 29]
    assert list(result) == expected

  def test_ncs_related_selection(self):
    # print sys._getframe().f_code.co_name
    result = nu.get_ncs_related_selection(
      ncs_restraints_group_list=self.ncs_restraints_group_list,
      asu_size=25)
    # ASU length is set to 25 (not 21)
    expected = [True, True, True, True, True, True, True,
                True, True, True, True, True, True, True,
                True, True, True, True, True, True, True,
                False, False, False, False]
    assert list(result) == expected

  def test_center_of_coordinates_shift(self):
    """
    test shifting translation to and from the center of coordinates of the
    master ncs copy
    """
    # print sys._getframe().f_code.co_name

    xrs = self.pdb_obj.xray_structure_simple()
    nrg = self.ncs_restraints_group_list

    shifts = nu.get_ncs_gorups_centers(
      xray_structure = xrs,
      ncs_restraints_group_list=nrg)

    xyz = self.pdb_obj.hierarchy.atoms().extract_xyz()
    center_of_coor = (flex.vec3_double([xyz.sum()]) * (1/xyz.size())).round(8)
    # test shifts
    t1 = shifts[0].round(8)
    t2 = shifts[1].round(8)
    d1 = flex.sqrt((center_of_coor-t1).dot()).min_max_mean().as_tuple()
    d2 = flex.sqrt((center_of_coor-t2).dot()).min_max_mean().as_tuple()
    assert (d1 == d2) and (d1 == (0,0,0))

    # test shift to center
    new_nrg = nu.shift_translation_to_center(
      shifts = shifts,
      ncs_restraints_group_list=nrg)
    expected = (-4.62169, -5.42257, 5.288)
    assert (new_nrg[0].copies[0].t.round(5)).elems == expected
    # back to original coordinates system
    old_nrg = nu.shift_translation_back_to_place(
      shifts=shifts,
      ncs_restraints_group_list=new_nrg)
    expected = (old_nrg[0].copies[0].t.round(5)).elems
    result = (nrg[0].copies[0].t.round(5)).elems
    assert result == expected

  def test_nrg_selection(self):
    """
    test that a atom selection propagates correctly to ncs_restraints_group_list
    """
    # print sys._getframe().f_code.co_name

    nrg = self.ncs_restraints_group_list
    m1 = list(nrg[0].master_iselection)
    c1 = list(nrg[0].copies[0].iselection)
    c2 = list(nrg[0].copies[1].iselection)

    assert len(m1) == len(c1)
    assert m1 == [0,   1,  2,  5,  6,  3,  4]
    assert c1 == [7,   8,  9, 12, 13, 10, 11]
    assert c2 == [14, 15, 16, 19, 20, 17, 18]

    selection1 = flex.size_t([0,1,5,3,100,101])
    selection2 = flex.size_t([0,1,5,3,7,8,9,12,100,101])
    selection3 = flex.size_t([0,1,5,3,7,8,9,12,14,15,19,17,100,101])

    new_nrg = nu.ncs_groups_selection(
      ncs_restraints_group_list=nrg,
      selection=selection1)
    # only atoms in master are selected
    mt = list(new_nrg[0].master_iselection)
    c1t = list(new_nrg[0].copies[0].iselection)

    assert mt == []
    assert c1t == []

    # atoms selected in both master and copies
    new_nrg = nu.ncs_groups_selection(
      ncs_restraints_group_list=nrg,
      selection=selection2)
    # only atoms in master are selected
    mt = list(new_nrg[0].master_iselection)
    c1t = list(new_nrg[0].copies[0].iselection)

    assert mt == []
    assert c1t == []

    new_nrg = nu.ncs_groups_selection(
      ncs_restraints_group_list=nrg,
      selection=selection3)
    # only atoms in master are selected
    mt = list(new_nrg[0].master_iselection)
    c1t = list(new_nrg[0].copies[0].iselection)
    c2t = list(new_nrg[0].copies[1].iselection)

    assert mt == [0, 1, 5]
    assert c1t == [7, 8, 12]
    assert c2t == [14, 15, 19]

  def test_selected_positions(self):
    # print sys._getframe().f_code.co_name

    a=flex.size_t([1,2,5,6,4])
    pos={0,3,4}
    s, d = nu.selected_positions(a,pos)
    assert list(s) == [1,6,4]
    assert list(d) == [2,5]

  def test_remove_items_from_selection(self):
    # print sys._getframe().f_code.co_name
    a=flex.size_t([1,2,5,6,4])
    r = flex.size_t([2,5])
    s = nu.remove_items_from_selection(a,r)
    assert list(s) == [1,6,4]

  def test_ncs_restraints_group_list_switching(self):
    """
    Testing switching of master ncs copy and a copy (turning the copy to master)
    and then switching to the original configuration
    """
    # print sys._getframe().f_code.co_name
    nrg = self.tr_obj2.get_ncs_restraints_group_list()
    master = nrg[0].master_iselection
    copy_1 = nrg[0].copies[0].iselection
    # switch master with the first copy
    nrg_new = nu.change_ncs_groups_master(
      ncs_restraints_group_list=nrg,
      new_masters=[1])
    # test that selection have switched
    assert list(copy_1) == list(nrg_new[0].master_iselection)
    assert list(master) == list(nrg_new[0].copies[0].iselection)
    # test that rotation matrices and translation vectors are properly converted
    R1 = self.rotation1.transpose()
    R2 = self.rotation2 * R1
    R3 = self.rotation3 * R1

    T1 = - R1 * self.translation1
    T2 = self.rotation2 * T1 + self.translation2
    T3 = self.rotation3 * T1 + self.translation3

    assert R1.round(8).elems == nrg_new[0].copies[0].r.round(8).elems
    assert R2.round(8).elems == nrg_new[0].copies[1].r.round(8).elems
    assert R3.round(8).elems == nrg_new[0].copies[2].r.round(8).elems

    assert T1.round(8).elems == nrg_new[0].copies[0].t.round(8).elems
    assert T2.round(6).elems == nrg_new[0].copies[1].t.round(6).elems
    assert T3.round(6).elems == nrg_new[0].copies[2].t.round(6).elems

    # flip back to original formation
    nrg_original = nu.change_ncs_groups_master(
      ncs_restraints_group_list=nrg_new,
      new_masters=[1])
    # test that selection have switched
    r1 = self.rotation1
    r2 = self.rotation2
    r3 = self.rotation3

    t1 = self.translation1
    t2 = self.translation2
    t3 = self.translation3

    assert list(master) == list(nrg_original[0].master_iselection)
    assert list(copy_1) == list(nrg_original[0].copies[0].iselection)
    # test that rotation matrices and translation vectors are properly converted
    assert r1.round(4).elems == nrg_original[0].copies[0].r.round(4).elems
    assert r2.round(4).elems == nrg_original[0].copies[1].r.round(4).elems
    assert r3.round(4).elems == nrg_original[0].copies[2].r.round(4).elems

    assert t1.round(4).elems == nrg_original[0].copies[0].t.round(4).elems
    assert t2.round(4).elems == nrg_original[0].copies[1].t.round(4).elems
    assert t3.round(4).elems == nrg_original[0].copies[2].t.round(4).elems


  def test_get_list_of_best_ncs_copy_map_correlation(self):
    """
    Verifying that we get a list of chain index for the chain with the best
    map correlation
    """
    # print sys._getframe().f_code.co_name

    d_min = 1.0
    ncs_inp = ncs.input(pdb_string=pdb_poor_0)
    ncs_restraints_group_list = ncs_inp.get_ncs_restraints_group_list()

    pdb_inp_poor = pdb.input(lines=pdb_poor_0,source_info=None)
    ph_poor = pdb_inp_poor.construct_hierarchy()
    ph_poor.atoms().reset_i_seq()
    xrs_poor = pdb_inp_poor.xray_structure_simple()

    pdb_inp_answer = pdb.input(lines=pdb_answer_0,source_info=None)
    ph_answer = pdb_inp_answer.construct_hierarchy()
    ph_answer.atoms().reset_i_seq()
    xrs_answer = pdb_inp_answer.xray_structure_simple()

    fc = xrs_answer.structure_factors(d_min=d_min, algorithm="direct").f_calc()
    fft_map = fc.fft_map(resolution_factor = 0.25)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()

    selections = [flex.size_t([0,1,2,3,4,5,6]),flex.size_t([7,8,9,10,11,12,13]),
                flex.size_t([14,15,16,17,18,19,20])]

    mp = mmtbx.maps.correlation.from_map_and_xray_structure_or_fmodel(
      xray_structure = xrs_poor,
      map_data       = map_data,
      d_min          = d_min)

    cc = mp.cc(selections=selections)

    best_list = nu.get_list_of_best_ncs_copy_map_correlation(
      ncs_restraints_group_list,
      xray_structure = xrs_poor,
      map_data       = map_data,
      d_min          = d_min)

    assert [round(x,1) for x in cc] == [0.5,0.8,0.5]
    assert best_list == [1]

  def test_iselection_ncs_to_asu(self):
    # print sys._getframe().f_code.co_name
    pdb_inp = pdb.input(lines=pdb_answer_0,source_info=None)
    ph = pdb_inp.construct_hierarchy()
    isel_asu = flex.size_t([8,9,13])
    isel_ncs = isel_asu - 7
    results = nu.iselection_ncs_to_asu(isel_ncs,'B',ph)
    self.assertEqual(list(isel_asu),list(results))

  def test_iselection_asu_to_ncs(self):
    # print sys._getframe().f_code.co_name
    pdb_inp = pdb.input(lines=pdb_answer_0,source_info=None)
    ph = pdb_inp.construct_hierarchy()
    isel_asu = flex.size_t([8,9,13])
    isel_ncs = isel_asu - 7
    results = nu.iselection_asu_to_ncs(isel_asu,'B',ph)
    self.assertEqual(list(isel_ncs),list(results))

  def test_change_ncs_groups_master(self):
    """ test change_ncs_groups_master when we have more than one group """
    ncs_obj_phil = ncs.input(
      pdb_string=test_pdb_str_2,
      ncs_phil_string=phil_str)
    nrgl = ncs_obj_phil.get_ncs_restraints_group_list()
    self.assertAlmostEqual(len(nrgl),2)
    # check that the masters coordinates are the same after several flips
    # select masters and some copies as references
    master_1 = nrgl[0].master_iselection
    master_2 = nrgl[1].master_iselection
    copy_1_1 = nrgl[0].copies[0].iselection
    copy_2_1 = nrgl[1].copies[0].iselection
    copy_2_2 = nrgl[1].copies[1].iselection
    # Test rotation and translations
    pdb_inp = pdb.input(lines=test_pdb_str_2,source_info=None)
    ph = pdb_inp.construct_hierarchy()
    verify_transforms(nrgl,ph)
    # change only the second master
    nrgl_new = nu.change_ncs_groups_master(
      ncs_restraints_group_list=nrgl,
      new_masters=[0,2])
    self.assertEqual(list(master_1),list(nrgl_new[0].master_iselection))
    self.assertEqual(list(copy_1_1),list(nrgl_new[0].copies[0].iselection))
    self.assertEqual(list(master_2),list(nrgl_new[1].copies[1].iselection))
    self.assertEqual(list(copy_2_1),list(nrgl_new[1].copies[0].iselection))
    self.assertEqual(list(copy_2_2),list(nrgl_new[1].master_iselection))
    verify_transforms(nrgl_new,ph)

def verify_transforms(nrgl,ph):
  """
  Check that all copies relate correctly to master via the transforms
  Will raise an error if necessary

  Args:
    nrgl: ncs_restraints_group_list
    ph: Hierarchy object
  """
  for gr in nrgl:
    m_xyz = ph.select(gr.master_iselection).atoms().extract_xyz()
    for cp in gr.copies:
      c_xyz = ph.select(cp.iselection).atoms().extract_xyz()
      xyz = cp.r.elems * m_xyz + cp.t
      assert approx_equal(c_xyz,xyz,1)

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_change_ncs_groups_master']
  suite = unittest.TestSuite(map(Test_ncs_utils,tests))
  return suite

if __name__=='__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main()
