from __future__ import absolute_import, division, print_function

import math
import random
import pytest

from cctbx import crystal, sgtbx, uctbx
from dxtbx.model import (
    Crystal,
    MosaicCrystalKabsch2010,
    MosaicCrystalSauter2014,
    CrystalFactory,
)
from libtbx.test_utils import approx_equal
from scitbx import matrix


def random_rotation():
    from scitbx.math import euler_angles_as_matrix

    return euler_angles_as_matrix([random.uniform(0, 360) for i in range(3)])


def test_crystal_model_from_mosflm_matrix():
    mosflm_matrix = map(
        float,
        """ -0.00495480 -0.01491776  0.00238445
  0.01505572 -0.00661190 -0.00149401
  0.00585043  0.00438127  0.00586415
       0.000       0.000       0.000
  -0.2932645  -0.8829514   0.3665960
   0.8911171  -0.3913446  -0.2296951
   0.3462750   0.2593185   0.9015806
     57.7822     57.7822    150.0931     90.0000     90.0000     90.0000
       0.000       0.000       0.000""".split(),
    )
    A = mosflm_matrix[:9]
    unit_cell = uctbx.unit_cell(mosflm_matrix[21:27])
    cm = CrystalFactory.from_mosflm_matrix(A, unit_cell=unit_cell)
    assert approx_equal(
        cm.get_unit_cell().parameters(), unit_cell.parameters(), eps=1.0e-2
    )


def test_crystal_model():
    real_space_a = matrix.col((10, 0, 0))
    real_space_b = matrix.col((0, 11, 0))
    real_space_c = matrix.col((0, 0, 12))
    model = Crystal(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    # This doesn't work as python class uctbx.unit_cell(uctbx_ext.unit_cell)
    # so C++ and python classes are different types
    # assert isinstance(model.get_unit_cell(), uctbx.unit_cell)
    assert model.get_unit_cell().parameters() == (10.0, 11.0, 12.0, 90.0, 90.0, 90.0)
    assert approx_equal(model.get_A(), (1 / 10, 0, 0, 0, 1 / 11, 0, 0, 0, 1 / 12))
    assert approx_equal(
        matrix.sqr(model.get_A()).inverse(), (10, 0, 0, 0, 11, 0, 0, 0, 12)
    )
    assert approx_equal(model.get_B(), model.get_A())
    assert approx_equal(model.get_U(), (1, 0, 0, 0, 1, 0, 0, 0, 1))
    assert approx_equal(
        model.get_real_space_vectors(), (real_space_a, real_space_b, real_space_c)
    )
    assert (
        model.get_crystal_symmetry().unit_cell().parameters()
        == model.get_unit_cell().parameters()
    )
    assert model.get_crystal_symmetry().space_group() == model.get_space_group()

    model2 = Crystal(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    assert model == model2 and not (model != model2)

    model2a = Crystal(model.get_A(), model.get_space_group())
    assert model == model2a and not (model != model2a)

    model2b = Crystal(
        matrix.sqr(model.get_A()).inverse().elems,
        model.get_space_group().type().lookup_symbol(),
        reciprocal=False,
    )
    assert model == model2b and not (model != model2b)

    # rotate 45 degrees about x-axis
    R1 = matrix.sqr(
        (
            1,
            0,
            0,
            0,
            math.cos(math.pi / 4),
            -math.sin(math.pi / 4),
            0,
            math.sin(math.pi / 4),
            math.cos(math.pi / 4),
        )
    )
    # rotate 30 degrees about y-axis
    R2 = matrix.sqr(
        (
            math.cos(math.pi / 6),
            0,
            math.sin(math.pi / 6),
            0,
            1,
            0,
            -math.sin(math.pi / 6),
            0,
            math.cos(math.pi / 6),
        )
    )
    # rotate 60 degrees about z-axis
    R3 = matrix.sqr(
        (
            math.cos(math.pi / 3),
            -math.sin(math.pi / 3),
            0,
            math.sin(math.pi / 3),
            math.cos(math.pi / 3),
            0,
            0,
            0,
            1,
        )
    )
    R = R1 * R2 * R3
    model.set_U(R)
    # B is unchanged
    assert approx_equal(model.get_B(), (1 / 10, 0, 0, 0, 1 / 11, 0, 0, 0, 1 / 12))
    assert approx_equal(model.get_U(), R)
    assert approx_equal(
        model.get_A(), matrix.sqr(model.get_U()) * matrix.sqr(model.get_B())
    )
    a_, b_, c_ = model.get_real_space_vectors()
    assert approx_equal(a_, R * real_space_a)
    assert approx_equal(b_, R * real_space_b)
    assert approx_equal(c_, R * real_space_c)
    assert (
        str(model).replace("-0.0000", " 0.0000")
        == """\
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
"""
    )
    model.set_B((1 / 12, 0, 0, 0, 1 / 12, 0, 0, 0, 1 / 12))
    assert approx_equal(model.get_unit_cell().parameters(), (12, 12, 12, 90, 90, 90))

    U = matrix.sqr(
        (0.3455, -0.2589, -0.9020, 0.8914, 0.3909, 0.2293, 0.2933, -0.8833, 0.3658)
    )
    B = matrix.sqr((1 / 13, 0, 0, 0, 1 / 13, 0, 0, 0, 1 / 13))
    model.set_A(U * B)
    assert approx_equal(model.get_A(), U * B)
    assert approx_equal(model.get_U(), U, 1e-4)
    assert approx_equal(model.get_B(), B, 1e-5)

    model3 = Crystal(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group=sgtbx.space_group_info("P 222").group(),
    )
    assert model3.get_space_group().type().hall_symbol() == " P 2 2"
    assert model != model3
    #
    sgi_ref = sgtbx.space_group_info(number=230)
    model_ref = Crystal(
        real_space_a=(44, 0, 0),
        real_space_b=(0, 44, 0),
        real_space_c=(0, 0, 44),
        space_group=sgi_ref.group(),
    )
    assert approx_equal(model_ref.get_U(), (1, 0, 0, 0, 1, 0, 0, 0, 1))
    assert approx_equal(model_ref.get_B(), (1 / 44, 0, 0, 0, 1 / 44, 0, 0, 0, 1 / 44))
    assert approx_equal(model_ref.get_A(), model_ref.get_B())
    assert approx_equal(
        model_ref.get_unit_cell().parameters(), (44, 44, 44, 90, 90, 90)
    )
    a_ref, b_ref, c_ref = map(matrix.col, model_ref.get_real_space_vectors())
    cb_op_to_primitive = sgi_ref.change_of_basis_op_to_primitive_setting()
    model_primitive = model_ref.change_basis(cb_op_to_primitive)
    cb_op_to_reference = (
        model_primitive.get_space_group()
        .info()
        .change_of_basis_op_to_reference_setting()
    )
    a_prim, b_prim, c_prim = map(matrix.col, model_primitive.get_real_space_vectors())
    # print cb_op_to_primitive.as_abc()
    ##'-1/2*a+1/2*b+1/2*c,1/2*a-1/2*b+1/2*c,1/2*a+1/2*b-1/2*c'
    assert approx_equal(a_prim, -1 / 2 * a_ref + 1 / 2 * b_ref + 1 / 2 * c_ref)
    assert approx_equal(b_prim, 1 / 2 * a_ref - 1 / 2 * b_ref + 1 / 2 * c_ref)
    assert approx_equal(c_prim, 1 / 2 * a_ref + 1 / 2 * b_ref - 1 / 2 * c_ref)
    # print cb_op_to_reference.as_abc()
    ##b+c,a+c,a+b
    assert approx_equal(a_ref, b_prim + c_prim)
    assert approx_equal(b_ref, a_prim + c_prim)
    assert approx_equal(c_ref, a_prim + b_prim)
    assert approx_equal(
        model_primitive.get_U(),
        [
            -0.5773502691896258,
            0.40824829046386285,
            0.7071067811865476,
            0.5773502691896257,
            -0.4082482904638631,
            0.7071067811865476,
            0.5773502691896257,
            0.8164965809277259,
            0.0,
        ],
    )
    assert approx_equal(
        model_primitive.get_B(),
        [
            0.0262431940540739,
            0.0,
            0.0,
            0.00927837023781507,
            0.02783511071344521,
            0.0,
            0.01607060866333063,
            0.01607060866333063,
            0.03214121732666125,
        ],
    )
    assert approx_equal(
        model_primitive.get_A(),
        (0, 1 / 44, 1 / 44, 1 / 44, 0, 1 / 44, 1 / 44, 1 / 44, 0),
    )
    assert approx_equal(
        model_primitive.get_unit_cell().parameters(),
        [
            38.1051177665153,
            38.1051177665153,
            38.1051177665153,
            109.47122063449069,
            109.47122063449069,
            109.47122063449069,
        ],
    )
    assert model_ref != model_primitive
    model_ref_recycled = model_primitive.change_basis(cb_op_to_reference)
    assert approx_equal(model_ref.get_U(), model_ref_recycled.get_U())
    assert approx_equal(model_ref.get_B(), model_ref_recycled.get_B())
    assert approx_equal(model_ref.get_A(), model_ref_recycled.get_A())
    assert approx_equal(
        model_ref.get_unit_cell().parameters(),
        model_ref_recycled.get_unit_cell().parameters(),
    )
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
    model = Crystal(
        direct_matrix[:3], direct_matrix[3:6], direct_matrix[6:9], space_group=sg
    )
    assert uc.is_similar_to(model.get_unit_cell())
    uc_minimum = uc.change_basis(cb_op_to_minimum)
    model_minimum = model.change_basis(cb_op_to_minimum)
    assert uc_minimum.is_similar_to(model_minimum.get_unit_cell())
    assert model_minimum != model
    model_minimum.update(model)
    assert model_minimum == model
    #
    from scitbx.math import euler_angles

    A_static = matrix.sqr(model.get_A())
    A_as_scan_points = [A_static]
    num_scan_points = 11
    for i in range(num_scan_points - 1):
        A_as_scan_points.append(
            A_as_scan_points[-1] * matrix.sqr(euler_angles.xyz_matrix(0.1, 0.2, 0.3))
        )
    model.set_A_at_scan_points(A_as_scan_points)
    model_minimum = model.change_basis(cb_op_to_minimum)
    assert model.num_scan_points == model_minimum.num_scan_points == num_scan_points
    M = matrix.sqr(cb_op_to_minimum.c_inv().r().transpose().as_double())
    M_inv = M.inverse()
    for i in range(num_scan_points):
        A_orig = matrix.sqr(model.get_A_at_scan_point(i))
        A_min = matrix.sqr(model_minimum.get_A_at_scan_point(i))
        assert approx_equal(A_min, A_orig * M_inv)

    mosaic_model = MosaicCrystalKabsch2010(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    mosaic_model2 = MosaicCrystalKabsch2010(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    assert approx_equal(mosaic_model.get_mosaicity(), 0)
    assert mosaic_model == mosaic_model2
    mosaic_model2.set_mosaicity(0.01)
    assert mosaic_model != mosaic_model2
    # FIXME Crystal == MosaicCrystal gives unexpected result, depending on
    # parameter order
    # model4 = Crystal(real_space_a=(10,0,0),
    #                  real_space_b=(0,11,0),
    #                  real_space_c=(0,0,12),
    #                  space_group_symbol="P 1")
    # print "model4 == mosaic_model", model4 == mosaic_model # True
    # print "mosaic_model == model4", mosaic_model == model4 # False
    # print "mosaic_model2 == model4", mosaic_model2 == model4 # False
    # print "model4 == mosaic_model2", model4 == mosaic_model2 # True

    mosaic_model = MosaicCrystalSauter2014(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    mosaic_model2 = MosaicCrystalSauter2014(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    assert approx_equal(mosaic_model.get_half_mosaicity_deg(), 0)
    assert approx_equal(mosaic_model.get_domain_size_ang(), 0)
    assert mosaic_model == mosaic_model2
    mosaic_model2.set_half_mosaicity_deg(0.01)
    assert mosaic_model != mosaic_model2
    mosaic_model2.set_half_mosaicity_deg(0)
    assert mosaic_model == mosaic_model2
    mosaic_model2.set_domain_size_ang(1000)
    assert mosaic_model != mosaic_model2


def test_similarity():
    model_1 = MosaicCrystalKabsch2010(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    model_1.set_mosaicity(0.5)
    model_2 = MosaicCrystalKabsch2010(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    model_2.set_mosaicity(0.5)
    assert model_1.is_similar_to(model_2)
    model_1.set_mosaicity(-1)
    model_2.set_mosaicity(-0.5)
    assert model_1.is_similar_to(model_2)  # test ignores negative mosaicity
    model_1.set_mosaicity(0.5)
    model_2.set_mosaicity(0.63)  # outside tolerance
    assert not model_1.is_similar_to(model_2)
    model_2.set_mosaicity(0.62)  # just inside tolerance

    # orientation tests
    R = matrix.sqr(model_2.get_U())
    dr1 = matrix.col((1, 0, 0)).axis_and_angle_as_r3_rotation_matrix(0.0101, deg=True)
    dr2 = matrix.col((1, 0, 0)).axis_and_angle_as_r3_rotation_matrix(0.0099, deg=True)
    model_2.set_U(dr1 * R)
    assert not model_1.is_similar_to(model_2)  # outside tolerance
    model_2.set_U(dr2 * R)
    assert model_1.is_similar_to(model_2)  # inside tolerance

    # unit_cell.is_similar_to is tested elsewhere


@pytest.mark.xfail
def test_change_basis_mosaic_crystal():
    from cctbx.sgtbx import change_of_basis_op

    mosaic_model = MosaicCrystalSauter2014(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )
    mosaic_model.set_half_mosaicity_deg(0.01)
    mosaic_model.set_domain_size_ang(1000)
    assert mosaic_model.get_half_mosaicity_deg() == 0.01
    assert mosaic_model.get_domain_size_ang() == 1000

    # Currently change_basis downgrades the MosaicCrystalSauter2014 to
    # a Crystal object so this will fail
    mosaic_model2 = mosaic_model.change_basis(change_of_basis_op())
    assert mosaic_model2.get_half_mosaicity_deg() == 0.01
    assert mosaic_model2.get_domain_size_ang() == 1000
    assert mosaic_model.is_similar_to(mosaic_model2)


def test_check_old_vs_new():
    from dxtbx.tests.model.crystal_model_old import crystal_model_old

    model_1 = Crystal(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )

    model_2 = crystal_model_old(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )

    cov_B = matrix.sqr([1] * (9 * 9))

    model_1.set_B_covariance(cov_B)
    model_2.set_B_covariance(cov_B)

    A_list = [model_1.get_A() for i in range(20)]

    model_1.set_A_at_scan_points(A_list)
    model_2.set_A_at_scan_points(A_list)

    cell_sd_1 = model_1.get_cell_parameter_sd()
    cell_sd_2 = model_2.get_cell_parameter_sd()
    cell_volume_sd_1 = model_1.get_cell_volume_sd()
    cell_volume_sd_2 = model_2.get_cell_volume_sd()
    covB1 = model_1.get_B_covariance()
    covB2 = model_1.get_B_covariance()

    A1 = model_1.get_A()
    A2 = model_2.get_A()
    U1 = model_1.get_U()
    U2 = model_2.get_U()
    B1 = model_1.get_B()
    B2 = model_2.get_B()
    UC1 = model_1.get_unit_cell()
    UC2 = model_2.get_unit_cell()
    RSV1 = model_1.get_real_space_vectors()
    RSV2 = model_2.get_real_space_vectors()
    SG1 = model_1.get_space_group()
    SG2 = model_2.get_space_group()

    assert model_1.num_scan_points == model_2.num_scan_points

    A_list_1 = [
        model_1.get_A_at_scan_point(i) for i in range(model_1.get_num_scan_points())
    ]
    A_list_2 = [
        model_2.get_A_at_scan_point(i) for i in range(model_1.get_num_scan_points())
    ]
    B_list_1 = [
        model_1.get_B_at_scan_point(i) for i in range(model_1.get_num_scan_points())
    ]
    B_list_2 = [
        model_2.get_B_at_scan_point(i) for i in range(model_1.get_num_scan_points())
    ]
    U_list_1 = [
        model_1.get_U_at_scan_point(i) for i in range(model_1.get_num_scan_points())
    ]
    U_list_2 = [
        model_2.get_U_at_scan_point(i) for i in range(model_1.get_num_scan_points())
    ]

    cell_sd_1 = model_1.get_cell_parameter_sd()
    cell_sd_2 = model_2.get_cell_parameter_sd()
    cell_volume_sd_1 = model_1.get_cell_volume_sd()
    cell_volume_sd_2 = model_2.get_cell_volume_sd()
    covB1 = model_1.get_B_covariance()
    covB2 = model_1.get_B_covariance()

    assert approx_equal(A1, A2)
    assert approx_equal(B1, B2)
    assert approx_equal(U1, U2)
    assert approx_equal(UC1.parameters(), UC2.parameters())
    assert approx_equal(RSV1[0], RSV2[0])
    assert approx_equal(RSV1[1], RSV2[1])
    assert approx_equal(RSV1[2], RSV2[2])
    assert str(SG1.info()) == str(SG2.info())

    for i in range(model_1.get_num_scan_points()):
        assert approx_equal(A_list_1[i], A_list_2[i])
        assert approx_equal(B_list_1[i], B_list_2[i])
        assert approx_equal(U_list_1[i], U_list_2[i])

    assert approx_equal(covB1, covB2)
    assert approx_equal(cell_volume_sd_1, cell_volume_sd_2)
    assert approx_equal(cell_sd_1, cell_sd_2)


@pytest.mark.parametrize(
    "crystal_class", [Crystal, MosaicCrystalKabsch2010, MosaicCrystalSauter2014]
)
def test_set_scan_varying_B_covariance(crystal_class):

    xl = crystal_class(
        real_space_a=(10, 0, 0),
        real_space_b=(0, 11, 0),
        real_space_c=(0, 0, 12),
        space_group_symbol="P 1",
    )

    from scitbx.array_family import flex

    cov_B = flex.double(
        (
            [8e-14, -1e-29, 3e-30, 3e-14, 8e-14, 3e-30, 2e-15, -7e-15, 2e-14],
            [-1e-29, 2e-45, -4e-46, -4e-30, -1e-29, -4e-46, 1e-30, 4e-30, -3e-30],
            [3e-30, -4e-46, 1e-46, 9e-31, 3e-30, 1e-46, 2e-31, -2e-31, 9e-31],
            [3e-14, -4e-30, 9e-31, 4e-14, 3e-14, 9e-31, -2e-15, -1e-15, 7e-15],
            [8e-14, -1e-29, 3e-30, 3e-14, 1e-13, 3e-30, 6e-15, -1e-15, 3e-14],
            [3e-30, -4e-46, 1e-46, 9e-31, 3e-30, 1e-46, 2e-31, -2e-31, 9e-31],
            [2e-15, 1e-30, 2e-31, -2e-15, 6e-15, 2e-31, 2e-14, 1e-14, 2e-15],
            [-7e-15, 4e-30, -2e-31, -1e-15, -1e-15, -2e-31, 1e-14, 3e-14, -2e-15],
            [2e-14, -3e-30, 9e-31, 7e-15, 3e-14, 9e-31, 2e-15, -2e-15, 8e-15],
        )
    )
    xl.set_B_covariance(cov_B)

    cov_B_array = flex.double(flex.grid(20, 9, 9))
    cov_B_2d = cov_B.as_scitbx_matrix()
    cov_B.reshape(flex.grid(1, 9, 9))
    for i in range(20):
        cov_B_array[i : (i + 1), :, :] = cov_B

    # This should fail: set_A_at_scan_points should be called first
    with pytest.raises(
        RuntimeError,
        message="Setting B covariance before A at " "scan points is not allowed",
    ):
        xl.set_B_covariance_at_scan_points(cov_B_array)

    A_list = [xl.get_A() for i in range(20)]
    xl.set_A_at_scan_points(A_list)

    # Now the setter should work
    xl.set_B_covariance_at_scan_points(cov_B_array)

    # Check getters for B covariance and cell sd
    cell_sd = xl.get_cell_parameter_sd()
    for i in range(20):
        cov_B_at_scan_point = xl.get_B_covariance_at_scan_point(i)
        assert cov_B_at_scan_point == cov_B_2d
        cell_sd_at_scan_point = xl.get_cell_parameter_sd_at_scan_point(i)
        assert cell_sd_at_scan_point == pytest.approx(cell_sd)
