from __future__ import division
import math
from cStringIO import StringIO
from libtbx.test_utils import approx_equal, show_diff
from scitbx import matrix
from cctbx import crystal, sgtbx, uctbx
from dxtbx.model.crystal import crystal_model, \
     crystal_model_from_mosflm_matrix

def random_rotation():
  import random
  from scitbx.math import euler_angles_as_matrix
  return euler_angles_as_matrix([random.uniform(0,360) for i in xrange(3)])

def exercise_crystal_model_from_mosflm_matrix():
  mosflm_matrix = map(float, ''' -0.00495480 -0.01491776  0.00238445
  0.01505572 -0.00661190 -0.00149401
  0.00585043  0.00438127  0.00586415
       0.000       0.000       0.000
  -0.2932645  -0.8829514   0.3665960
   0.8911171  -0.3913446  -0.2296951
   0.3462750   0.2593185   0.9015806
     57.7822     57.7822    150.0931     90.0000     90.0000     90.0000
       0.000       0.000       0.000'''.split())
  A = mosflm_matrix[:9]
  unit_cell = uctbx.unit_cell(mosflm_matrix[21:27])
  cm = crystal_model_from_mosflm_matrix(A, unit_cell=unit_cell)
  assert approx_equal(cm.get_unit_cell().parameters(),
                      unit_cell.parameters(), eps=1.0e-2)

def exercise_crystal_model():

  real_space_a = matrix.col((10,0,0))
  real_space_b = matrix.col((0,11,0))
  real_space_c = matrix.col((0,0,12))
  model = crystal_model(real_space_a=(10,0,0),
                        real_space_b=(0,11,0),
                        real_space_c=(0,0,12),
                        space_group_symbol="P 1")
  assert isinstance(model.get_unit_cell(), uctbx.unit_cell)
  assert model.get_unit_cell().parameters() == (
    10.0, 11.0, 12.0, 90.0, 90.0, 90.0)
  assert approx_equal(model.get_A(), (1/10, 0, 0, 0, 1/11, 0, 0, 0, 1/12))
  assert approx_equal(model.get_A().inverse(), (10, 0, 0, 0, 11, 0, 0, 0, 12))
  assert approx_equal(model.get_B(), model.get_A())
  assert approx_equal(model.get_U(), (1, 0, 0, 0, 1, 0, 0, 0, 1))
  assert approx_equal(model.get_real_space_vectors(),
                      (real_space_a, real_space_b, real_space_c))
  assert approx_equal(model.get_mosaicity(), 0)

  model2 = crystal_model(real_space_a=(10,0,0),
                         real_space_b=(0,11,0),
                         real_space_c=(0,0,12),
                         space_group_symbol="P 1")
  assert model == model2
  model2.set_mosaicity(0.01)
  assert model != model2
  # rotate 45 degrees about x-axis
  R1 = matrix.sqr((1, 0, 0,
                   0, math.cos(math.pi/4), -math.sin(math.pi/4),
                   0, math.sin(math.pi/4), math.cos(math.pi/4)))
  # rotate 30 degrees about y-axis
  R2 = matrix.sqr((math.cos(math.pi/6), 0, math.sin(math.pi/6),
                   0, 1, 0,
                   -math.sin(math.pi/6), 0, math.cos(math.pi/6)))
  # rotate 60 degrees about z-axis
  R3 = matrix.sqr((math.cos(math.pi/3), -math.sin(math.pi/3), 0,
                   math.sin(math.pi/3), math.cos(math.pi/3), 0,
                   0, 0, 1))
  R = R1 * R2 * R3
  model.set_U(R)
  # B is unchanged
  assert approx_equal(model.get_B(), (1/10, 0, 0, 0, 1/11, 0, 0, 0, 1/12))
  assert approx_equal(model.get_U(), R)
  assert approx_equal(model.get_A(), model.get_U() * model.get_B())
  a_, b_, c_ = model.get_real_space_vectors()
  assert approx_equal(a_, R * real_space_a)
  assert approx_equal(b_, R * real_space_b)
  assert approx_equal(c_, R * real_space_c)
  s = StringIO()
  print >> s, model
  assert not show_diff(s.getvalue().replace("-0.0000", " 0.0000"), """\
Crystal:
    Unit cell: (10.000, 11.000, 12.000, 90.000, 90.000, 90.000)
    Space group: P 1
    U matrix:  {{ 0.4330, -0.7500,  0.5000},
                { 0.7891,  0.0474, -0.6124},
                { 0.4356,  0.6597,  0.6124}}
    B matrix:  {{ 0.1000,  0.0000,  0.0000},
                { 0.0000,  0.0909,  0.0000},
                { 0.0000,  0.0000,  0.0833}}
    A = UB:    {{ 0.0433, -0.0682,  0.0417},
                { 0.0789,  0.0043, -0.0510},
                { 0.0436,  0.0600,  0.0510}}

""")
  model.set_B((1/12, 0, 0, 0, 1/12, 0, 0, 0, 1/12))
  assert approx_equal(model.get_unit_cell().parameters(),
                      (12, 12, 12, 90, 90, 90))

  model3 = crystal_model(
    real_space_a=(10,0,0),
    real_space_b=(0,11,0),
    real_space_c=(0,0,12),
    space_group=sgtbx.space_group_info("P 222").group(),
    mosaicity=0.01)
  assert model3.get_space_group().type().hall_symbol() == " P 2 2"
  assert model != model3
  #
  sgi_ref = sgtbx.space_group_info(number=230)
  model_ref = crystal_model(
    real_space_a=(44,0,0),
    real_space_b=(0,44,0),
    real_space_c=(0,0,44),
    space_group=sgi_ref.group())
  assert approx_equal(model_ref.get_U(), (1,0,0,0,1,0,0,0,1))
  assert approx_equal(model_ref.get_B(), (1/44,0,0,0,1/44,0,0,0,1/44))
  assert approx_equal(model_ref.get_A(), model_ref.get_B())
  assert approx_equal(model_ref.get_unit_cell().parameters(),
                      (44,44,44,90,90,90))
  a_ref, b_ref, c_ref = model_ref.get_real_space_vectors()
  cb_op_to_primitive = sgi_ref.change_of_basis_op_to_primitive_setting()
  model_primitive = model_ref.change_basis(cb_op_to_primitive)
  cb_op_to_reference = model_primitive.get_space_group().info()\
    .change_of_basis_op_to_reference_setting()
  a_prim, b_prim, c_prim = model_primitive.get_real_space_vectors()
  #print cb_op_to_primitive.as_abc()
  ##'-1/2*a+1/2*b+1/2*c,1/2*a-1/2*b+1/2*c,1/2*a+1/2*b-1/2*c'
  assert approx_equal(a_prim, -1/2 * a_ref + 1/2 * b_ref + 1/2 * c_ref)
  assert approx_equal(b_prim, 1/2 * a_ref - 1/2 * b_ref + 1/2 * c_ref)
  assert approx_equal(c_prim, 1/2 * a_ref + 1/2 * b_ref - 1/2 * c_ref)
  #print cb_op_to_reference.as_abc()
  ##b+c,a+c,a+b
  assert approx_equal(a_ref, b_prim + c_prim)
  assert approx_equal(b_ref, a_prim + c_prim)
  assert approx_equal(c_ref, a_prim + b_prim)
  assert approx_equal(
    model_primitive.get_U(),
  [-0.5773502691896258, 0.40824829046386285, 0.7071067811865476,
   0.5773502691896257, -0.4082482904638631, 0.7071067811865476,
   0.5773502691896257, 0.8164965809277259, 0.0])
  assert approx_equal(
    model_primitive.get_B(),
    [0.0262431940540739, 0.0, 0.0, 0.00927837023781507, 0.02783511071344521,
     0.0, 0.01607060866333063, 0.01607060866333063, 0.03214121732666125])
  assert approx_equal(
    model_primitive.get_A(), (0, 1/44, 1/44, 1/44, 0, 1/44, 1/44, 1/44, 0))
  assert approx_equal(
    model_primitive.get_unit_cell().parameters(),
    [38.1051177665153, 38.1051177665153, 38.1051177665153,
     109.47122063449069, 109.47122063449069, 109.47122063449069])
  assert model_ref != model_primitive
  model_ref_recycled = model_primitive.change_basis(cb_op_to_reference)
  assert approx_equal(model_ref.get_U(), model_ref_recycled.get_U())
  assert approx_equal(model_ref.get_B(), model_ref_recycled.get_B())
  assert approx_equal(model_ref.get_A(), model_ref_recycled.get_A())
  assert approx_equal(model_ref.get_unit_cell().parameters(),
                      model_ref_recycled.get_unit_cell().parameters())
  assert model_ref == model_ref_recycled
  #
  uc = uctbx.unit_cell((58.2567, 58.1264, 39.7093, 46.9077, 46.8612, 62.1055))
  sg = sgtbx.space_group_info(symbol="P1").group()
  cs = crystal.symmetry(unit_cell=uc, space_group=sg)
  cb_op_to_minimum = cs.change_of_basis_op_to_minimum_cell()
  # the reciprocal matrix
  B = matrix.sqr(uc.fractionalization_matrix()).transpose()
  U = random_rotation()
  direct_matrix = (U * B).inverse()
  model = crystal_model(direct_matrix[:3],
                        direct_matrix[3:6],
                        direct_matrix[6:9],
                        space_group=sg)
  assert uc.is_similar_to(model.get_unit_cell())
  uc_minimum = uc.change_basis(cb_op_to_minimum)
  model_minimum = model.change_basis(cb_op_to_minimum)
  assert uc_minimum.is_similar_to(model_minimum.get_unit_cell())
  assert model_minimum != model
  model_minimum.update(model)
  assert model_minimum == model
  #
  from scitbx.math import euler_angles
  A_static = model.get_A()
  A_as_scan_points = [A_static]
  num_scan_points = 11
  for i in range(num_scan_points-1):
    A_as_scan_points.append(
      A_as_scan_points[-1] * matrix.sqr(euler_angles.xyz_matrix(0.1,0.2,0.3)))
  model.set_A_at_scan_points(A_as_scan_points)
  model_minimum = model.change_basis(cb_op_to_minimum)
  assert model.num_scan_points == model_minimum.num_scan_points == num_scan_points
  M = matrix.sqr(cb_op_to_minimum.c_inv().r().transpose().as_double())
  M_inv = M.inverse()
  for i in range(num_scan_points):
    A_orig = model.get_A_at_scan_point(i)
    A_min = model_minimum.get_A_at_scan_point(i)
    assert A_min == A_orig * M_inv

def exercise_similarity():

  model_1 = crystal_model(real_space_a=(10,0,0),
                          real_space_b=(0,11,0),
                          real_space_c=(0,0,12),
                          space_group_symbol="P 1",
                          mosaicity=0.5)
  model_2 = crystal_model(real_space_a=(10,0,0),
                          real_space_b=(0,11,0),
                          real_space_c=(0,0,12),
                          space_group_symbol="P 1",
                          mosaicity=0.5)
  assert model_1.is_similar_to(model_2)

  # mosaicity tests
  model_1.set_mosaicity(-1)
  model_2.set_mosaicity(-0.5)
  assert model_1.is_similar_to(model_2) # test ignores negative mosaicity
  model_1.set_mosaicity(0.5)
  model_2.set_mosaicity(0.63) # outside tolerance
  assert not model_1.is_similar_to(model_2)
  model_2.set_mosaicity(0.625) #just inside tolerance
  assert model_1.is_similar_to(model_2)

  # orientation tests
  R = model_2.get_U()
  dr1 = matrix.col((1, 0, 0)).axis_and_angle_as_r3_rotation_matrix(0.0101, deg=True)
  dr2 = matrix.col((1, 0, 0)).axis_and_angle_as_r3_rotation_matrix(0.0099, deg=True)
  model_2.set_U(dr1 * R)
  assert not model_1.is_similar_to(model_2) # outside tolerance
  model_2.set_U(dr2 * R)
  assert model_1.is_similar_to(model_2) # inside tolerance

  # unit_cell.is_similar_to is tested elsewhere
  return

def run():
  exercise_crystal_model()
  exercise_crystal_model_from_mosflm_matrix()
  exercise_similarity()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()
