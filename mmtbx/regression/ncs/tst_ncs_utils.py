from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.ncs.ncs_utils as nu
import mmtbx.maps.correlation
from scitbx import matrix
import iotbx.ncs as ncs
import iotbx.pdb
import unittest
import math
from six.moves import map

__author__ = 'Youval'

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

  def test_get_list_of_best_ncs_copy_map_correlation(self):
    """
    Verifying that we get a list of chain index for the chain with the best
    map correlation
    """
    # print sys._getframe().f_code.co_name

    d_min = 1.0
    pdb_inp = iotbx.pdb.input(lines=pdb_poor_0,source_info=None)
    ncs_inp = ncs.input(hierarchy=pdb_inp.construct_hierarchy())
    ncs_restraints_group_list = ncs_inp.get_ncs_restraints_group_list()

    pdb_inp_poor = iotbx.pdb.input(lines=pdb_poor_0,source_info=None)
    ph_poor = pdb_inp_poor.construct_hierarchy(sort_atoms=False)
    ph_poor.atoms().reset_i_seq()
    xrs_poor = pdb_inp_poor.xray_structure_simple()

    pdb_inp_answer = iotbx.pdb.input(lines=pdb_answer_0,source_info=None)
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

    nu.get_list_of_best_ncs_copy_map_correlation(
      ncs_groups     = ncs_restraints_group_list,
      xray_structure = xrs_poor,
      map_data       = map_data,
      d_min          = d_min)

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_transform_update']
  suite = unittest.TestSuite(list(map(Test_ncs_utils,tests)))
  return suite

if __name__=='__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
