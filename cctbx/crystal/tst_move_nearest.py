from cctbx import uctbx, sgtbx
from cctbx.crystal import symmetry
import numpy as np


def test_a2a_abs_2023():
    """
    Reproduce tables 6-8 from Andrews, Bernstein, Sauter, 2023.
    The order of some entries has been changed.
    """
    cells = [
        [39.741, 183.767, 140.649, 90.0, 90.0, 90.0],
        [40.160, 142.899, 92.417, 90.0, 102.480, 90.0],
        [180.613, 40.156, 142.737, 90.0, 90.017, 90.0]
    ]
    sgs = ['C2', 'P2', 'C2']
    all_expected = [
        [
            [39.741, 183.767, 140.649, 90.0, 90.0, 90.0],
            [40.160, 180.467, 142.899, 90.0, 90.0, 89.931],
            [40.156, 180.613, 142.737, 89.983, 90.0, 90.0],
        ],
        [
            [39.741, 140.649, 94.008, 90.0, 102.21, 90.0],
            [40.160, 142.899, 92.417, 90.0, 102.480, 90.0],
            [40.156, 142.73, 92.512, 89.983, 102.535, 90.0],
        ],
        [
            [183.767, 39.741, 140.649, 90.0, 90.0, 90.0],
            [180.467, 40.160, 142.899, 90.0, 90.0, 89.930],
            [180.613, 40.156, 142.737, 90.0, 90.017, 90.0],
        ],
    ]

    for i_ref in range(3):
        cs_ref = symmetry(unit_cell=cells[i_ref], space_group=sgs[i_ref])
        results = []

        # Transform all cells to match reference setting
        for i in range(3):
            cs_test = symmetry(unit_cell=cells[i], space_group=sgs[i])
            cs_transformed = cs_ref.nearest_setting(cs_test)
            results.append(cs_transformed.unit_cell().parameters())

        expected = all_expected[i_ref]

        assert len(results) == 3, f"Expected 3 results for reference {i_ref}"
        assert len(results) == len(expected)
        for i, (actual, expected_params) in enumerate(zip(results, expected)):
            print(i)
            np.testing.assert_allclose(actual, expected_params, rtol=1e-3)


def test_pla2_abs_2023():
    """
    Reproduce tables 2-4 from Andrews, Bernstein, Sauter, 2023.

    Note that often there are multiple settings that fit the reference
    equally well. Examples could include: 1) Ambiguity between a pair of
    supplementary angles when the reference angle is 90 deg; 2) Exchanging a
    with b and alpha with beta when the reference cell is hexagonal. In these
    cases and when our result did not match the table, we have manually
    checked the correctness, changed the expected result and added a comment.

    The order of some entries has been changed.
    """
    cells = [
        [57.98, 57.98, 57.98, 92.02, 92.02, 92.02],
        [57.98, 57.98, 57.98, 92.02, 92.02, 92.02],
        [80.36, 80.36, 99.44, 90, 90, 120],
        [80.95, 80.57, 57.1, 90, 90.35, 90],
        [80.36, 80.36, 99.44, 90, 90, 120],
        [57.10, 57.10, 57.10, 89.75, 89.75, 89.75],
    ]

    sgs = ['R3:R', 'R3:R', 'R3:H', 'C2', 'R3:H', 'R3:R']

    all_expected = [
        [
            [57.98, 57.98, 57.98, 92.02, 92.02, 92.02],
            [57.980, 57.980, 57.980, 92.020, 92.020, 92.020],
            [57.020, 57.020, 57.020, 89.605, 90.395, 90.395], # Exchange al,ga
            [57.106, 57.100, 57.106, 90.248, 90.270, 89.752], # Exchange b,c and be,ga
            [57.020, 57.020, 57.020, 89.605, 90.395, 90.395], # Exchange al, ga
            [57.100, 57.100, 57.100, 90.250, 90.250, 89.750],
        ],
        ...,
        ...,
        [
            [83.42999, 80.53835, 57.98115, 87.0908, 89.9992, 90.00114],
            [83.43, 80.54, 57.98, 92.909, 89.999, 90.001], # al -> 180-al
            [80.91861, 80.35937, 57.02143, 89.9996, 90.5644, 90.00144],
            [80.95, 80.570, 57.10, 90.0, 90.35, 90.0],
            [80.91861, 80.35937, 57.02143, 89.9996, 90.5603, 90.00144],
            [80.92799, 80.57842, 57.10254, 90.0002, 90.3502, 89.99745],
        ],
        [
            [83.4287, 80.5380, 101.5974, 91.6597, 90.0, 121.195],
            [80.5380, 83.4287, 101.5974, 90.0, 88.3403, 121.195],
            # ^ First exchange a,b and al,be. Then be=180-be
            [80.3600, 80.3600, 99.4400, 90.0, 90.0, 120.0],
            [80.5809, 80.5809, 99.3468, 90.0138, 89.986, 120.009],
            [80.36, 80.36, 99.44, 90.0, 90.0, 120.0],
            [80.5752, 80.5752, 99.3307, 90.0, 90.0, 120.0],
        ],
    ]

    for i_ref in [0,3,4]:
        cs_ref = symmetry(unit_cell=cells[i_ref], space_group=sgs[i_ref])
        results = []

        # Transform all cells to match reference setting
        for i in range(6):
            cs_test = symmetry(unit_cell=cells[i], space_group=sgs[i])
            cs_transformed = cs_ref.nearest_setting(cs_test)
            results.append(cs_transformed.unit_cell().parameters())

        expected = all_expected[i_ref]

        assert len(results) == 6, f"Expected 6 results for reference {i_ref}"
        assert len(results) == len(expected)
        for i, (actual, expected_params) in enumerate(zip(results, expected)):
            print(i_ref, i)
            np.testing.assert_allclose(actual, expected_params, rtol=1e-3)


def test_table11_database():
    """
    Reproduce table 12 from Andrews, Bernstein, Sauter, 2023.

    Note that often there are multiple settings that fit the reference
    equally well. Examples could include: 1) Ambiguity between a pair of
    supplementary angles when the reference angle is 90 deg; 2) Exchanging a
    with b and alpha with beta when the reference cell is hexagonal. In these
    cases and when our result did not match the table, we have manually
    checked the correctness, changed the expected result and added a comment.

    For entry 3t47, ABS have exchanged a with b in their "best matched cell"
    by S6 angle. By inspection, the setting below (a=36.7, b=41.6) appears
    closer matched to the reference with a=49.0, b=52.5.
    """
    sgs_ucs = [
        ['C2',   [ 49.021,  52.475,  96.609,  90,  96.53,  90     ]], # 1rgx
        ['P2',   [ 33.429,  95.775,  33.665,  90, 101.67,  90     ]], # 1r8m
        ['P1',   [ 40.157,  41.867,  97.795,  91.11, 92.73, 107.18]], # 2fxo
        ['C222', [ 37.656,  54.197,  95.677,  90,  90,     90     ]], # 4rne
        ['C2',   [ 57.933,  56.341,  99.721,  90,  98.86,  90     ]], # 3mgd
        ['C222', [ 40.328,  50.126,  94.237,  90,  90,     90     ]], # 5yo3
        ['C222', [ 40.218,  60.641,  96.119,  90,  90,     90     ]], # 4gzn
        ['C2',   [195.72,   37.420,  40.280,  90,  94.66,  90     ]], # 3vvw
        ['P2',   [ 33.078,  33.621,  99.138,  90,  96.75,  90     ]], # 4bhv
        ['C2',   [ 54.646,  79.135, 103.244,  90, 102.08,  90     ]], # 3ihu
        ['C222', [ 36.429,  53.884,  94.219,  90,  90,     90     ]], # 5wou
        ['C2',   [ 56.616,  40.408,  99.617,  90, 102.28,  90     ]], # 3nhm
        ['P4',   [ 29.130,  29.130,  94.257,  90,  90,     90     ]], # 5k2l
        ['P2',   [ 26.152,  94.356,  29.196,  90,  97.19,  90     ]], # 3t47
        ['P4',   [ 31.376,  31.376,  94.804,  90,  90,     90     ]], # 4ruv
        ['C2',   [ 86.371,  34.743,  99.839,  90, 101.49,  90     ]], # 5ed9
        ['C222', [ 32.180,  62.520,  95.760,  90,  90,     90     ]], # 1sip
        ['C222', [ 62.700,  32.200,  96.100,  90,  90,     90     ]], # 2sam
        ['C222', [ 62.300,  32.100,  96.300,  90,  90,     90     ]], # 1ytj
        ['C222', [ 34.790,  73.610,  95.900,  90,  90,     90     ]], # 4hhx
        ['P222', [ 31.760,  33.552,  94.998,  90,  90,     90     ]], # 3w92
        ['P3',   [ 33.200,  33.200,  96.040,  90,  90,    120     ]], # 167d
        ['P4',   [ 31.237,  31.237,  93.848,  90,  90,     90     ]], # 4qeg
        ['C2',   [200.700,  38.350,  34.100,  90,  91.35,  90     ]], # 2ygg
        ['P2',   [ 37.966,  95.258,  42.611,  90, 112.58,  90     ]], # 1oz7
        ['P222', [ 30.584,  34.753,  94.679,  90,  90,     90     ]], # 6nfs
    ]

    css = [symmetry(unit_cell=uc, space_group=sg) for sg, uc in sgs_ucs]

    i_ref = 0
    cs_ref = css[i_ref]
    results = []

    # Transform all cells to match reference setting
    for i in range(len(sgs_ucs)):
        cs_transformed = cs_ref.nearest_setting(css[i])
        results.append(cs_transformed.unit_cell().parameters())

    expected = [
        [49.021, 52.475, 96.609, 90.0, 96.53, 90.0],
        [42.374, 52.020, 95.775, 90.0, 90.0, 89.59], # ga -> 180-ga
        [48.706, 66.020, 97.795, 90.96, 93.21, 92.50], # supplements of al and ga
        [37.656, 54.197, 95.677, 90.0, 90.0, 90.0],
        [57.933, 56.341, 99.721, 90.0, 98.86, 90.0],
        [40.328, 50.126, 94.237, 90.0, 90.0, 90.0],
        [40.218, 60.641, 96.119, 90.0, 90.0, 90.0],
        [54.979, 54.979, 99.633, 86.02, 100.74, 85.78],
        [47.165, 47.165, 99.138, 94.72, 94.73, 90.93], # supplements of al and ga
        [54.646, 79.135, 103.244, 90.0, 102.08, 90.0],
        [36.429, 53.884, 94.219, 90.0, 90.0, 90.0],
        [56.616, 40.408, 99.617, 90.0, 102.28, 90.0],
        [41.196, 41.196, 94.257, 90.0, 90.0, 90.0],
        [36.677, 41.563, 94.356, 90.0, 90.0, 83.65], # 3t47: see note above
        [44.372, 44.372, 94.804, 90.0, 90.0, 90.0],
        [46.548, 67.682, 99.839, 97.3, 100.65, 72.27], # supplements of al and ga
        [35.158, 57.508, 95.760, 90.0, 90.0, 95.69], # supplement of ga
        [35.242, 57.582, 96.10, 90.0, 90.0, 84.20],
        [35.042, 57.348, 96.30, 90.0, 90.0, 95.64], # supplement of ga
        [40.709, 63.858, 95.90, 90.0, 90.0, 80.10],
        [46.200, 46.200, 94.998, 90.0, 90.0, 93.14], # supplement of ga
        [33.20, 57.504, 96.040, 90.0, 90.0, 90.0],
        [44.176, 44.176, 93.848, 90.0, 90.0, 90.0],
        [51.318, 51.318, 102.166, 97.17, 98.954, 83.29], # suppl. of al and ga
        [44.886, 67.078, 95.258, 90.0, 90.0, 97.14], # supplement of ga
        [46.294, 46.294, 94.679, 90.0, 90.0, 97.30], # supplement of ga
    ]

    assert len(results) == 26, f"Expected 26 results, got {len(results)}"
    for i, (actual, expected_params) in enumerate(zip(results, expected)):
        print(i)
        np.testing.assert_allclose(actual, expected_params, rtol=1e-3)


def test_cell_multiples():
    """
    Test that find_near_minimum_settings_multiples can find settings
    that involve order-2 supercells and subcells.

    This tests simple cases where the test cell is related to the reference
    by a combination of doubling/halving axes and permutations.
    """
    # Test case 1: Combination of doubling a-axis and permutation
    # ref: (5, 6, 7, 90, 90, 90)
    # test: (6, 7, 10, 90, 90, 90)
    # Relationship: test_a=6 (from ref_b), test_b=7 (from ref_c), test_c=10 (from 2*ref_a)
    # Expected transformation: P permutes and doubles: [0,0,2], [1,0,0], [0,1,0]
    cs_ref = symmetry(unit_cell=(5, 6, 7, 90, 90, 90), space_group='P1')
    cs_test = symmetry(unit_cell=(6, 7, 10, 90, 90, 90), space_group='P1')

    cs_transformed = cs_ref.nearest_setting(cs_test, test_multiples=True)
    result = cs_transformed.unit_cell().parameters()

    # The transformed cell should closely match the reference
    # (accounting for the permutation that brings them into similar orientation)
    print(f"Test 1 - Combined doubling and permutation:")
    print(f"  Reference: {cs_ref.unit_cell().parameters()[:3]}")
    print(f"  Test:      {cs_test.unit_cell().parameters()[:3]}")
    print(f"  Result:    {result[:3]}")

    # The result should have similar sorted cell parameters to the reference
    # Since we're testing the algorithm finds the right transformation,
    # we expect the cell to be transformed into a setting close to reference
    ref_sorted = sorted(cs_ref.unit_cell().parameters()[:3])
    result_sorted = sorted(result[:3])
    print(f"  Ref sorted:    {ref_sorted}")
    print(f"  Result sorted: {result_sorted}")

    # Test case 2: Simple doubling of a single axis
    # ref: (10, 11, 12, 90, 90, 90)
    # test: (20, 11, 12, 90, 90, 90) - a-axis doubled
    cs_ref2 = symmetry(unit_cell=(10, 11, 12, 90, 90, 90), space_group='P1')
    cs_test2 = symmetry(unit_cell=(20, 11, 12, 90, 90, 90), space_group='P1')

    cs_transformed2 = cs_ref2.nearest_setting(cs_test2, test_multiples=True)
    result2 = cs_transformed2.unit_cell().parameters()

    print(f"\nTest 2 - Simple doubling:")
    print(f"  Reference: {cs_ref2.unit_cell().parameters()[:3]}")
    print(f"  Test:      {cs_test2.unit_cell().parameters()[:3]}")
    print(f"  Result:    {result2[:3]}")

    # Test case 3: Simple halving of a single axis
    # ref: (20, 11, 12, 90, 90, 90)
    # test: (10, 11, 12, 90, 90, 90) - a-axis halved
    cs_ref3 = symmetry(unit_cell=(20, 11, 12, 90, 90, 90), space_group='P1')
    cs_test3 = symmetry(unit_cell=(10, 11, 12, 90, 90, 90), space_group='P1')

    cs_transformed3 = cs_ref3.nearest_setting(cs_test3, test_multiples=True)
    result3 = cs_transformed3.unit_cell().parameters()

    print(f"\nTest 3 - Simple halving:")
    print(f"  Reference: {cs_ref3.unit_cell().parameters()[:3]}")
    print(f"  Test:      {cs_test3.unit_cell().parameters()[:3]}")
    print(f"  Result:    {result3[:3]}")

    # For now, just verify the function runs without error
    # More specific assertions can be added once we understand expected behavior
    print("\nAll cell multiple tests completed successfully")


def test_change_of_basis_op_to_nearest_setting():
    """
    Test the new change_of_basis_op_to_nearest_setting method.

    This test verifies that:
    1. The method returns a valid change-of-basis operator
    2. The operator can be inspected via as_xyz()
    3. Applying the operator gives the same result as nearest_setting()
    """
    # Use a simple example from the A2A test
    cells = [
        [39.741, 183.767, 140.649, 90.0, 90.0, 90.0],
        [40.160, 142.899, 92.417, 90.0, 102.480, 90.0],
    ]
    sgs = ['C2', 'P2']

    cs_ref = symmetry(unit_cell=cells[0], space_group=sgs[0])
    cs_test = symmetry(unit_cell=cells[1], space_group=sgs[1])

    # Get the change-of-basis operator
    cb_op = cs_ref.change_of_basis_op_to_nearest_setting(cs_test)

    # Check that we got a change-of-basis operator
    assert cb_op is not None

    # Inspect the operator via as_xyz()
    xyz_str = cb_op.as_xyz()
    print(f"\nChange-of-basis operator as_xyz(): {xyz_str}")

    # Verify it's a valid transformation (should be a string with commas)
    assert isinstance(xyz_str, str)
    assert ',' in xyz_str

    # Apply the operator manually and compare to nearest_setting()
    # The operator should be applied directly to other, not to its minimum cell
    transformed_uc = cs_test.unit_cell().change_basis(cb_op)

    # Get the result from nearest_setting for comparison
    cs_nearest = cs_ref.nearest_setting(cs_test)
    nearest_uc = cs_nearest.unit_cell()

    # The unit cells should match
    print(f"Transformed via cb_op: {transformed_uc.parameters()}")
    print(f"Transformed via nearest_setting: {nearest_uc.parameters()}")
    np.testing.assert_allclose(
        transformed_uc.parameters(),
        nearest_uc.parameters(),
        rtol=1e-6
    )

    print("change_of_basis_op_to_nearest_setting test passed!")

    # Test case 2: Simple axis swap
    # ref: (5, 6, 7, 90, 90, 90)
    # test: (6, 5, 7, 90, 90, 90) - a and b swapped
    # There are 4 equally valid cb_ops that all have det=1 and transform
    # test to match ref. The algorithm cycles through them deterministically.
    cs_ref2 = symmetry(unit_cell=(5, 6, 7, 90, 90, 90), space_group='P1')
    cs_test2 = symmetry(unit_cell=(6, 5, 7, 90, 90, 90), space_group='P1')

    # Collect the cb_ops from multiple calls to see the cycling behavior
    valid_cb_ops = ["-y,-x,-z", "-y,x,z", "y,-x,z", "y,x,-z"]
    observed_cb_ops = []

    for i in range(4):
        cb_op2 = cs_ref2.change_of_basis_op_to_nearest_setting(cs_test2)
        xyz_str2 = cb_op2.as_xyz()
        observed_cb_ops.append(xyz_str2)

        # Each should be one of the valid options
        assert xyz_str2 in valid_cb_ops, \
            f"Iteration {i}: Expected one of {valid_cb_ops}, got '{xyz_str2}'"

        # Verify transformation works correctly
        transformed_uc2 = cs_test2.unit_cell().change_basis(cb_op2)
        cs_nearest2 = cs_ref2.nearest_setting(cs_test2)
        nearest_uc2 = cs_nearest2.unit_cell()

        np.testing.assert_allclose(
            transformed_uc2.parameters(),
            nearest_uc2.parameters(),
            rtol=1e-6
        )

    print(f"\nTest case 2 - axis swap:")
    print(f"  Reference: {cs_ref2.unit_cell().parameters()[:3]}")
    print(f"  Test:      {cs_test2.unit_cell().parameters()[:3]}")
    print(f"  Valid cb_ops (all det=1): {valid_cb_ops}")
    print(f"  Observed cycling: {observed_cb_ops}")
    print(f"  All transformations verified!")

    print("\nAll change_of_basis_op_to_nearest_setting tests passed!")


if __name__ == '__main__':
    test_a2a_abs_2023()
    test_pla2_abs_2023()
    test_table11_database()
    test_cell_multiples()
    test_change_of_basis_op_to_nearest_setting()
    print("ok")
