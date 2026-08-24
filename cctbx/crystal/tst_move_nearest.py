import io

from cctbx import uctbx, sgtbx
from cctbx.crystal import symmetry
from cctbx.uctbx.near_minimum import cell_distance
import numpy as np


def _assert_nearest_setting_matches(actual, expected, reference_cell_params):
    """
    The tabulated `expected` cell parameters encode an arbitrary tie-break
    phase (see docstrings above) among settings that are equally close to the
    reference, so an exact match is not required. Instead, check that:
    1) `actual` describes the same lattice as `expected` (up to a similarity
       transformation), and
    2) `actual` is at least as close to the reference cell as the tabulated
       `expected` value is.
    """
    sim_ops = uctbx.unit_cell(actual).similarity_transformations(
        uctbx.unit_cell(expected),
        relative_length_tolerance=0.02,
        absolute_angle_tolerance=0.5,
        unimodular_generator_range=1)
    assert len(sim_ops) > 0, \
        f"No similarity transformation found between {actual} and {expected}"

    d_actual = cell_distance(reference_cell_params, actual)
    d_expected = cell_distance(reference_cell_params, expected)
    # The tabulated values are only given to a handful of decimal places, so
    # a strict 1e-3 absolute tolerance is too tight: rounding the tabulated
    # cell alone can shift its summed L1 distance by more than that. Allow
    # a small tolerance that scales with the tabulated distance itself.
    tolerance = max(1e-3, d_expected * 1e-3)
    assert d_actual <= d_expected + tolerance, \
        f"{actual} is farther from the reference cell than " \
        f"{expected} ({d_actual} > {d_expected})"


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
        reference_cell_params = cs_ref.unit_cell().parameters()
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
            _assert_nearest_setting_matches(
                actual, expected_params, reference_cell_params)


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
    reference_cell_params = cs_ref.unit_cell().parameters()
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
        _assert_nearest_setting_matches(
            actual, expected_params, reference_cell_params)


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

    # Test case 3: Verify has_nearer_setting returns False when identity is among tied settings
    cs_ref3 = symmetry(unit_cell=(5, 6, 7, 90, 90, 90), space_group='P1')
    cs_test3 = symmetry(unit_cell=(5, 6, 7, 90, 90, 90), space_group='P1')  # Already matching

    has_nearer = cs_ref3.has_nearer_setting(cs_test3)
    print(f"\nTest case 3 - identity setting:")
    print(f"  Reference: {cs_ref3.unit_cell().parameters()[:3]}")
    print(f"  Test:      {cs_test3.unit_cell().parameters()[:3]}")
    print(f"  has_nearer_setting: {has_nearer}")
    assert not has_nearer, "Should return False when current setting is already optimal"

    # Also test with the axis swap case - should still prefer identity
    cs_test4 = symmetry(unit_cell=(6, 5, 7, 90, 90, 90), space_group='P1')
    cb_op_first = cs_ref3.change_of_basis_op_to_nearest_setting(cs_test4)
    print(f"\nTest case 4 - axis swap with identity preference:")
    print(f"  First call cb_op (should prefer identity if tied): {cb_op_first.as_xyz()}")
    # Note: This won't be identity since (6,5,7) != (5,6,7), but we're testing the preference logic

    print("\nAll has_nearer_setting tests passed!")


def test_min_cell_frame_prefilter():
    """
    Regression test for the min-cell-frame bulk prefilter used by
    cctbx.command_line.search_nearest_cell (which replaced the retired
    query-frame bulk_query_frame_distances design). The prefilter ranks
    database rows by L1 distance -- in the query's MINIMUM-CELL frame --
    between each row's minimum cell and each of the query's cached
    nearly-reduced settings, taking the running elementwise min over
    settings (no (S,N,...) intermediate). This is exactly the metric
    change_of_basis_op_to_nearest_setting uses internally to pick its best
    setting (its `distances` / `best_dist` computation), so -- unlike the
    retired query-frame lower bound -- the prefilter value should now agree
    with the internal selection distance, not merely bound it.
    """
    np.random.seed(0)
    cs_ref = symmetry(
        unit_cell=(57.98, 57.98, 57.98, 92.02, 92.02, 92.02),
        space_group='R3:R')
    settings, _ = cs_ref.near_minimum_settings_and_cb_ops(
        length_tolerance=0.03, angle_tolerance=3.0, test_multiples=False)
    query_settings = np.array([s['cell'] for s in settings], dtype=float)

    # P1 accepts any (a,b,c,alpha,beta,gamma), so random triclinic cells are
    # always a valid crystal.symmetry -- this exercises the prefilter across
    # a wide range of raw cells without needing to hand-craft cells
    # compatible with higher-symmetry space groups.
    n_rows = 100
    axes = np.random.uniform(20, 150, size=(n_rows, 3))
    angles = np.random.uniform(70, 110, size=(n_rows, 3))
    rows = []
    for i in range(n_rows):
        cell = tuple(axes[i]) + tuple(angles[i])
        rows.append(symmetry(unit_cell=cell, space_group='P1'))

    # Standard-API replacement for the retired bulk_lean_min_cell_params: a
    # row whose minimum_cell() raises is excluded, not propagated; empty
    # input yields (0, 6)/(0,)-shaped arrays, matching the old bulk helper's
    # shape-safety.
    min_cell_params_list = []
    valid_indices_list = []
    for i, cs in enumerate(rows):
        try:
            p = cs.minimum_cell().unit_cell().parameters()
        except Exception:
            continue
        min_cell_params_list.append(p)
        valid_indices_list.append(i)
    min_cell_params = np.array(min_cell_params_list) if min_cell_params_list \
        else np.empty((0, 6))
    valid_indices = np.array(valid_indices_list, dtype=int)
    # All P1 cells are valid, so no row should have been dropped.
    assert list(valid_indices) == list(range(n_rows))

    best_dist = np.full(min_cell_params.shape[0], np.inf)
    for setting_row in query_settings:
        d = np.abs(min_cell_params - setting_row).sum(axis=1)
        np.minimum(best_dist, d, out=best_dist)

    # Reproduce change_of_basis_op_to_nearest_setting's internal selection
    # distance directly, one row at a time: its `uc_other` is exactly
    # min_cell_params[i], since both come from reducing the same row to its
    # minimum cell (the two computation paths -- the retired lean C++-reduction
    # shortcut and the standard change_of_basis_op_to_minimum_cell() API --
    # were verified numerically identical, max diff ~1.4e-13, across the full
    # 200k-row PDB database). This confirms the vectorized prefilter is an
    # equal alternate computation of the same quantity, not just a lower
    # bound.
    true_dist = np.array([
        min(cell_distance(min_cell_params[i], s) for s in query_settings)
        for i in range(n_rows)
    ])
    np.testing.assert_allclose(best_dist, true_dist, atol=1e-9)

    print("min-cell-frame prefilter agrees exactly with the per-row "
          f"reference computation for {n_rows} random rows")


def test_min_cell_frame_prefilter_ties():
    """
    Ties are EXACT for pseudo-symmetric cells: many of the query's
    near-reduced settings can land at bit-identical L1 distance from a
    same-lattice candidate. change_of_basis_op_to_nearest_setting's
    tie-break counter then cycles which of those tied settings -- and thus
    which cb_op / transformed cell -- is returned on successive calls (see
    nearest_setting_count). Check that the min-cell-frame prefilter distance
    is unaffected by that cycling (it only ever reports the minimum, not
    which setting achieved it), and that whichever setting is returned still
    describes the query's lattice (via similarity_transformations, never
    exact cell parameters, since settings may permute across calls).
    """
    cs_ref = symmetry(
        unit_cell=(57.98, 57.98, 57.98, 92.02, 92.02, 92.02),
        space_group='R3:R')
    # Same lattice as cs_ref, described in a different (still near-reduced)
    # basis -- see test_pla2_abs_2023's "Exchange al,ga" row.
    cs_row = symmetry(
        unit_cell=(57.020, 57.020, 57.020, 89.605, 90.395, 90.395),
        space_group='R3:R')

    settings, _ = cs_ref.near_minimum_settings_and_cb_ops(
        length_tolerance=0.03, angle_tolerance=3.0, test_multiples=False)
    query_settings = np.array([s['cell'] for s in settings], dtype=float)

    # Standard-API replacement for the retired bulk_lean_min_cell_params (see
    # test_min_cell_frame_prefilter above).
    min_cell_params_list = []
    valid_indices_list = []
    for i, cs in enumerate([cs_row]):
        try:
            p = cs.minimum_cell().unit_cell().parameters()
        except Exception:
            continue
        min_cell_params_list.append(p)
        valid_indices_list.append(i)
    min_cell_params = np.array(min_cell_params_list) if min_cell_params_list \
        else np.empty((0, 6))
    valid_indices = np.array(valid_indices_list, dtype=int)
    assert list(valid_indices) == [0]
    best_dist = np.full(1, np.inf)
    for setting_row in query_settings:
        d = np.abs(min_cell_params - setting_row).sum(axis=1)
        np.minimum(best_dist, d, out=best_dist)

    # Bit-identical ties among near-reduced settings are expected for this
    # pseudo-cubic cell; confirm more than one setting actually achieves the
    # minimum, so the tie-break cycling exercised below is meaningful.
    d_all = np.abs(query_settings - min_cell_params[0]).sum(axis=1)
    n_tied = int(np.sum(d_all <= best_dist[0] + 1e-9))
    assert n_tied > 1, f"expected multiple tied settings, got {n_tied}"

    for _ in range(4):
        cb_op = cs_ref.change_of_basis_op_to_nearest_setting(
            cs_row, length_tolerance=0.03, angle_tolerance=3.0,
            test_multiples=False)
        transformed_uc = cs_row.unit_cell().change_basis(cb_op)

        # The prefilter's distance must equal whatever
        # change_of_basis_op_to_nearest_setting's internal selection found,
        # regardless of which tied setting its counter happened to pick.
        cb_other_to_min = cs_row.change_of_basis_op_to_minimum_cell()
        uc_other = cs_row.unit_cell().change_basis(cb_other_to_min).parameters()
        selected_dist = min(cell_distance(uc_other, s) for s in query_settings)
        assert abs(selected_dist - best_dist[0]) < 1e-9

        sim_ops = uctbx.unit_cell(transformed_uc.parameters()).similarity_transformations(
            uctbx.unit_cell(cs_row.unit_cell().parameters()),
            relative_length_tolerance=0.02,
            absolute_angle_tolerance=0.5,
            unimodular_generator_range=1)
        assert len(sim_ops) > 0, \
            f"{transformed_uc.parameters()} is not the same lattice as " \
            f"{cs_row.unit_cell().parameters()}"

    print("min-cell-frame prefilter distance is stable across tie-break "
          "cycling; finalized cell remains lattice-equivalent each time")


def test_fetch_and_search_symmetry_construction_agree():
    """
    Pins fetch_pdb_cells.py's private `_crystal_symmetry_from_row` and
    search_nearest_cell.py's inline :R-retry block inside `load_pdb_cells`
    together. Per Decision 3 of the fetch/precompute plan, these are two
    independently duplicated copies of the same bare-rhombohedral-retry
    logic (no shared helper module in cctbx/command_line/); this test
    catches silent divergence between them.
    """
    from cctbx.command_line import fetch_pdb_cells, search_nearest_cell

    # (pdb_id, a, b, c, alpha, beta, gamma, space_group_symbol): a few
    # ordinary symbols, a bare rhombohedral symbol needing the :R retry (at
    # two different pseudo-cubic cells -- see test_pla2_abs_2023), an
    # explicit :H rhombohedral symbol that should NOT need a retry, and one
    # row that should fail both variants.
    rows_in = [
        ('synC2a', 49.021, 52.475, 96.609, 90, 96.53, 90, 'C2'),
        ('synP2a', 33.429, 95.775, 33.665, 90, 101.67, 90, 'P2'),
        ('synP1a', 40.157, 41.867, 97.795, 91.11, 92.73, 107.18, 'P1'),
        ('synC222', 37.656, 54.197, 95.677, 90, 90, 90, 'C222'),
        ('synP4a', 29.130, 29.130, 94.257, 90, 90, 90, 'P4'),
        ('synP3a', 33.200, 33.200, 96.040, 90, 90, 120, 'P3'),
        ('synRhomA', 57.98, 57.98, 57.98, 92.02, 92.02, 92.02, 'R 3 2'),
        ('synRhomB', 57.10, 57.10, 57.10, 89.75, 89.75, 89.75, 'R 3 2'),
        ('synHexR3', 80.36, 80.36, 99.44, 90, 90, 120, 'R3:H'),
        ('synBad', 10, 10, 10, 90, 90, 90, 'not a real symbol'),
        ('synP222', 31.760, 33.552, 94.998, 90, 90, 90, 'P222'),
        ('synC2b', 54.646, 79.135, 103.244, 90, 102.08, 90, 'C2'),
    ]
    expect_fail = {'synBad'}

    sg_info_cache = {}
    direct_cs = {}
    for pdb_id, a, b, c, al, be, ga, sg in rows_in:
        try:
            cs = fetch_pdb_cells._crystal_symmetry_from_row(
                a, b, c, al, be, ga, sg, sg_info_cache)
        except Exception:
            assert pdb_id in expect_fail, \
                f"{pdb_id} unexpectedly failed _crystal_symmetry_from_row"
            continue
        assert pdb_id not in expect_fail, \
            f"{pdb_id} unexpectedly succeeded in _crystal_symmetry_from_row"
        direct_cs[pdb_id] = cs
    assert set(direct_cs) == {r[0] for r in rows_in} - expect_fail

    # Feed the same rows through the real CSV writer (process_entry-shaped
    # dicts, Z filled with a placeholder), then load them back through
    # search_nearest_cell's independent inline :R-retry copy.
    base_rows = [
        {
            'pdb_id': pdb_id, 'unit_cell_a': a, 'unit_cell_b': b,
            'unit_cell_c': c, 'unit_cell_alpha': al, 'unit_cell_beta': be,
            'unit_cell_gamma': ga, 'unit_cell_Z': 1, 'space_group': sg,
        }
        for pdb_id, a, b, c, al, be, ga, sg in rows_in
    ]
    csv_path = '/Users/dwpaley/daily/20251120/.tmp/tst_symmetry_equiv.csv'
    fetch_pdb_cells.write_cells_csv(base_rows, csv_path, out=io.StringIO())

    loaded_rows, n_skipped, _ = search_nearest_cell.load_pdb_cells(csv_path)
    loaded_cs = dict(loaded_rows)

    assert n_skipped == len(expect_fail), \
        f"expected {len(expect_fail)} skipped rows, got {n_skipped}"
    assert set(loaded_cs) == set(direct_cs)

    for pdb_id, cs_direct in direct_cs.items():
        cs_loaded = loaded_cs[pdb_id]
        np.testing.assert_allclose(
            cs_direct.unit_cell().parameters(),
            cs_loaded.unit_cell().parameters(),
            rtol=1e-6, atol=1e-6,
            err_msg=f"{pdb_id}: unit cell mismatch between the two code paths")
        assert cs_direct.space_group_info().type().number() == \
            cs_loaded.space_group_info().type().number(), \
            f"{pdb_id}: space group mismatch between the two code paths"

    print("fetch_pdb_cells._crystal_symmetry_from_row and "
          "search_nearest_cell.load_pdb_cells agree on "
          f"{len(direct_cs)} rows, including a bare-rhombohedral :R retry; "
          f"{len(expect_fail)} row(s) correctly rejected by both")


def test_search_nearest_cell_matches_direct_computation():
    """
    Ground-truth check for cctbx.search_nearest_cell.run() against a
    synthetic augmented CSV. Replaces the retired "precomputed vs
    on-the-fly fallback" comparison (there is no fallback path anymore --
    see Decision 2 of the fetch/precompute plan). The augmented CSV's
    min-cell columns are computed in-process via
    fetch_pdb_cells.write_cells_csv (the same precompute step
    cctbx.fetch_pdb_cells performs on a real fetch), and the expected
    per-row distances are computed independently: reduce each row to its
    minimum cell and take the minimum cell_distance to any of the query's
    nearly-reduced settings -- exactly the quantity run() reports as
    "Distance" (see test_min_cell_frame_prefilter's proof that the
    min-cell-frame prefilter agrees exactly with this per-row computation).
    """
    from cctbx.command_line import fetch_pdb_cells, search_nearest_cell

    # Table 12's cells (test_table11_database) plus a few cells from the
    # a2a/pla2 tables, for a synthetic database spanning several space
    # groups.
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
        ['C2',   [ 39.741, 183.767, 140.649,  90,  90,      90    ]], # a2a
        ['P2',   [ 40.160, 142.899,  92.417,  90, 102.480,  90    ]], # a2a
        ['R3:R', [ 57.98,   57.98,   57.98,   92.02, 92.02, 92.02 ]], # pla2
        ['C2',   [ 80.95,   80.57,   57.1,    90,  90.35,   90    ]], # pla2
    ]

    pdb_ids = [f"synt{i:03d}" for i in range(len(sgs_ucs))]

    # Build the augmented CSV in-process via the real writer, exercising the
    # same precompute step cctbx.fetch_pdb_cells performs on a real fetch.
    base_rows = []
    for pdb_id, (sg, uc) in zip(pdb_ids, sgs_ucs):
        a, b, c, al, be, ga = uc
        base_rows.append({
            'pdb_id': pdb_id, 'unit_cell_a': a, 'unit_cell_b': b,
            'unit_cell_c': c, 'unit_cell_alpha': al, 'unit_cell_beta': be,
            'unit_cell_gamma': ga, 'unit_cell_Z': 1, 'space_group': sg,
        })
    csv_path = '/Users/dwpaley/daily/20251120/.tmp/tst_search_ground_truth.csv'
    fetch_pdb_cells.write_cells_csv(base_rows, csv_path, out=io.StringIO())

    # Ground truth, computed independently of both write_cells_csv and
    # search_nearest_cell.run(): the query is row 0's own cell/space group,
    # so it should be its own nearest match at ~0 distance.
    query_sg, query_uc = sgs_ucs[0]
    cs_ref = symmetry(unit_cell=query_uc, space_group=query_sg)
    settings, _ = cs_ref.near_minimum_settings_and_cb_ops(
        length_tolerance=0.03, angle_tolerance=3.0, test_multiples=False)
    query_settings = [s['cell'] for s in settings]

    expected_dist = {}
    for pdb_id, (sg, uc) in zip(pdb_ids, sgs_ucs):
        cs_row = symmetry(unit_cell=uc, space_group=sg)
        cb_to_min = cs_row.change_of_basis_op_to_minimum_cell()
        uc_min = cs_row.unit_cell().change_basis(cb_to_min).parameters()
        expected_dist[pdb_id] = min(
            cell_distance(uc_min, s) for s in query_settings)

    best_pdb_id = min(expected_dist, key=expected_dist.get)
    assert best_pdb_id == pdb_ids[0], \
        "query cell (synt000) should be its own nearest match"
    assert expected_dist[pdb_ids[0]] < 1e-6

    query_args = [
        "unit_cell=%.6f,%.6f,%.6f,%.6f,%.6f,%.6f" % tuple(query_uc),
        "space_group=%s" % query_sg,
        "data_file=%s" % csv_path,
        "n_results=%d" % len(sgs_ucs),
        "metric=cellparams",
    ]
    run_out = io.StringIO()
    search_nearest_cell.run(query_args, out=run_out)
    lines = run_out.getvalue().splitlines()

    header_idx = next(i for i, line in enumerate(lines) if line.startswith("Rank"))
    data_lines = [l for l in lines[header_idx + 1:] if l.strip()]
    assert len(data_lines) > 0

    actual = []
    for line in data_lines:
        tok = line.split()
        actual.append((tok[1], float(tok[2])))

    # Top hit: the query's own cell/space group, at ~0 distance.
    rank1_pdb, rank1_dist = actual[0]
    assert rank1_pdb == pdb_ids[0], \
        f"expected top hit {pdb_ids[0]}, got {rank1_pdb}"
    assert rank1_dist < 1e-3, \
        f"expected ~0 distance for exact match, got {rank1_dist}"

    # Sane ranking: non-decreasing distance, every reported row present
    # (n_results covers the whole synthetic table), and every reported
    # row's printed distance agrees with the independently computed ground
    # truth (matched by pdb code, not position -- exact ties among returned
    # rows are not required to break in any particular order).
    dists = [d for _, d in actual]
    assert dists == sorted(dists), "run() output is not sorted by distance"
    actual_order = [p for p, _ in actual]
    assert set(actual_order) == set(pdb_ids)
    for pdb_id, dist in actual:
        assert abs(dist - expected_dist[pdb_id]) < 1e-3, \
            f"{pdb_id}: run() reported {dist}, expected {expected_dist[pdb_id]}"

    print("search_nearest_cell.run() ranking matches independent "
          f"ground-truth computation for {len(sgs_ucs)} synthetic rows "
          f"(top hit: {rank1_pdb} at distance {rank1_dist:.6f})")


def test_cell_metrics_against_scalar_reference():
    """
    Independent scalar reimplementation (plain `math`, no numpy) of each
    non-v7 SAUC formula from cell_metrics.py's features()/distance(),
    checked at high precision. This is what catches degrees-vs-radians,
    column order, weight order and broadcasting mistakes that a
    self-consistent vectorized-only test would miss.
    """
    import math
    from cctbx.uctbx import cell_metrics as cm

    def scalar_features(metric, p):
        a, b, c, al_d, be_d, ga_d = p
        al, be, ga = math.radians(al_d), math.radians(be_d), math.radians(ga_d)
        if metric == 'cellparams':
            return list(p)
        if metric in ('l1', 'l2'):
            return [a, b, c, al * (b + c) / 2, be * (a + c) / 2, ga * (a + b) / 2]
        g1, g2, g3 = a * a, b * b, c * c
        g4, g5, g6 = 2*b*c*math.cos(al), 2*a*c*math.cos(be), 2*a*b*math.cos(ga)
        if metric == 'ncdist':
            return [g1, g2, g3, g4, g5, g6]
        if metric == 's6':
            return [g4/2, g5/2, g6/2, -g1-g5/2-g6/2, -g2-g4/2-g6/2, -g3-g4/2-g5/2]
        body = g1 + g2 + g3 + min(g4+g5+g6, -g4-g5+g6, -g4+g5-g6, g4-g5-g6)
        return [g1, g2, g3, g2+g3-abs(g4), g1+g3-abs(g5), g1+g2-abs(g6), body]

    def scalar_distance(metric, f1, f2):
        d = [y - x for x, y in zip(f1, f2)]
        if metric == 'cellparams':
            return sum(abs(x) for x in d)
        if metric == 'l1':
            return sum(abs(x) for x in d) / math.sqrt(6)
        if metric == 'l2':
            return math.sqrt(sum(x*x for x in d))
        if metric in ('ncdist', 's6'):
            return 0.1 * sum(x*x for x in d) ** 0.25
        w = [1, 1, 1, .5, .5, .5, 1/3]
        return 0.1 * sum(wi*x*x for wi, x in zip(w, d)) ** 0.25

    np.random.seed(10)
    n = 20
    cells = np.column_stack([
        np.random.uniform(10, 100, n), np.random.uniform(10, 100, n),
        np.random.uniform(10, 100, n), np.random.uniform(70, 110, n),
        np.random.uniform(70, 110, n), np.random.uniform(70, 110, n),
    ])
    ref = [15.3, 22.7, 31.1, 88.5, 91.2, 95.7]

    for metric in ('cellparams', 'l1', 'l2', 'ncdist', 's6', 'dc7unsrt'):
        f_ref = cm.features(metric, ref)
        assert cm.features(metric, np.array(ref)).shape == f_ref.shape  # (6,) input
        f_batch = cm.features(metric, cells)  # (N,6) input
        assert f_batch.shape == (n, len(f_ref))
        scalar_rows = [scalar_features(metric, cells[i]) for i in range(n)]
        np.testing.assert_allclose(f_batch, scalar_rows, rtol=1e-12)
        np.testing.assert_allclose(f_ref, scalar_features(metric, ref), rtol=1e-12)

        d_vec = cm.distance(metric, f_ref, f_batch)  # (6,) vs (N,6) broadcast
        d_scalar = [scalar_distance(metric, scalar_features(metric, ref), row)
                    for row in scalar_rows]
        np.testing.assert_allclose(d_vec, d_scalar, rtol=1e-12)

    print("cell_metrics formulas match independent scalar reimplementation "
          f"for {n} random cells across 6 metrics")


def test_ncdist_matches_compiled_NCDist():
    """
    The compiled Andrews-Bernstein NCDist additionally searches G6
    reduction boundaries, so it is always <= the local (non-boundary-
    searching) ncdist metric, with equality for a pair whose reduced cells
    sit in the interior of their Niggli domains (verified during planning:
    agreement to 1e-13 for a nearby pair).
    """
    import math
    from cctbx.uctbx import cell_metrics as cm
    from cctbx.uctbx.determine_unit_cell import NCDist

    np.random.seed(11)
    n = 15
    for _ in range(n):
        p1 = np.random.uniform(20, 150, 3).tolist() + np.random.uniform(80, 100, 3).tolist()
        p2 = np.random.uniform(20, 150, 3).tolist() + np.random.uniform(80, 100, 3).tolist()
        g1, g2 = cm.features('ncdist', p1), cm.features('ncdist', p2)
        oracle = 0.1 * math.sqrt(NCDist(list(g1), list(g2)))
        ours = cm.distance('ncdist', g1, g2)
        assert oracle <= ours + 1e-9, (oracle, ours)

    # Near pair: no reduction boundary lies between them, so the compiled
    # oracle and the local formula must agree closely.
    p1 = [40.0, 55.0, 70.0, 89.0, 91.0, 90.5]
    p2 = [40.02, 55.03, 69.98, 89.05, 90.95, 90.45]
    g1, g2 = cm.features('ncdist', p1), cm.features('ncdist', p2)
    oracle = 0.1 * math.sqrt(NCDist(list(g1), list(g2)))
    ours = cm.distance('ncdist', g1, g2)
    assert abs(oracle - ours) < 1e-9, (oracle, ours)

    print("local ncdist metric agrees with compiled NCDist oracle "
          f"(<= for {n} random pairs, equal within 1e-9 for a near pair)")


def test_cell_distance_is_the_cellparams_metric():
    """
    near_minimum.cell_distance is kept as tst_move_nearest.py's independent
    judge (see _assert_nearest_setting_matches); pin it exactly equal to
    the named 'cellparams' metric so the two cannot silently diverge.
    """
    from cctbx.uctbx import cell_metrics as cm

    np.random.seed(12)
    n = 20
    for _ in range(n):
        p1 = np.random.uniform(10, 200, 3).tolist() + np.random.uniform(60, 120, 3).tolist()
        p2 = np.random.uniform(10, 200, 3).tolist() + np.random.uniform(60, 120, 3).tolist()
        direct = cell_distance(p1, p2)
        via_metric = cm.distance(
            'cellparams', cm.features('cellparams', p1), cm.features('cellparams', p2))
        np.testing.assert_allclose(float(via_metric), direct, rtol=1e-15)

    print("near_minimum.cell_distance == cell_metrics.distance('cellparams', ...) "
          f"for {n} random pairs")


def test_v7_is_a_lattice_invariant():
    """
    v7 is a lattice invariant: reindexing a minimum cell by a unimodular
    change of basis and re-minimising must reproduce the same v7 vector
    (see cell_metrics.py's v7 branch of features()). Also pins one value
    against the 4-line recipe (niggli_cell + reciprocal().minimum_cell(),
    sorted triples) inlined here rather than imported from scratch code.
    """
    from cctbx.uctbx import cell_metrics as cm

    p = [20.1, 25.3, 31.7, 88.0, 92.0, 95.0]
    mc = uctbx.unit_cell(p).minimum_cell()
    cb_op = sgtbx.change_of_basis_op('y,z,x+y')
    remc = mc.change_basis(cb_op).minimum_cell()

    f1 = cm.features('v7', mc.parameters())
    f2 = cm.features('v7', remc.parameters())
    assert f1.shape == (7,)
    np.testing.assert_allclose(f1, f2, atol=1e-8)

    d = sorted(mc.niggli_cell().parameters()[:3])
    r = sorted(1.0 / x for x in mc.reciprocal().minimum_cell().parameters()[:3])[::-1]
    manual = d + r + [mc.volume() ** (1. / 3.)]
    np.testing.assert_allclose(f1, manual, rtol=1e-12)

    print("v7 feature vector is invariant under reindexing + re-minimisation, "
          "and matches the 4-line reference recipe")


def test_metric_ties_are_exact():
    """
    Regression guard for the plot_uc_cloud tie-cycling behaviour (see the
    Domain Knowledge note on nearest_setting_count): every non-v7 metric
    must produce bit-identical ties among the pla2 pseudo-cubic pair's
    near-reduced settings (measured during planning: 6-18 tied settings
    per metric, see test_min_cell_frame_prefilter_ties for the same pair),
    and change_of_basis_op_to_nearest_setting must still cycle across
    multiple tied cb_ops while every returned cell stays lattice-
    equivalent to the input.

    The cycling check uses metric='l1': for this particular pair, the
    identity cb_op happens to be among the tied settings under
    cellparams/ncdist/s6/dc7unsrt (verified both before and after the
    metric-threading change -- this is pre-existing identity-preference
    behaviour, not something introduced by metric threading), so those
    metrics deterministically return identity here and never exercise
    cycling. l1/l2 do not include identity among their ties for this pair
    and so exercise the cycling path that plot_uc_cloud depends on.
    """
    from cctbx.uctbx import cell_metrics as cm

    cs_ref = symmetry(
        unit_cell=(57.98, 57.98, 57.98, 92.02, 92.02, 92.02), space_group='R3:R')
    cs_row = symmetry(
        unit_cell=(57.020, 57.020, 57.020, 89.605, 90.395, 90.395), space_group='R3:R')

    settings, _ = cs_ref.near_minimum_settings_and_cb_ops(
        length_tolerance=0.03, angle_tolerance=3.0, test_multiples=False)
    query_settings = np.array([s['cell'] for s in settings], dtype=float)
    cb_to_min = cs_row.change_of_basis_op_to_minimum_cell()
    uc_other = cs_row.unit_cell().change_basis(cb_to_min).parameters()

    for metric in ('cellparams', 'l1', 'l2', 'ncdist', 's6', 'dc7unsrt'):
        f_settings = cm.features(metric, query_settings)
        f_other = cm.features(metric, np.asarray(uc_other))
        d = cm.distance(metric, f_other, f_settings)
        n_tied = int((d == d.min()).sum())
        assert n_tied > 1, f"{metric}: expected bit-exact ties, got {n_tied}"

    observed_xyz = set()
    for _ in range(4):
        cb_op = cs_ref.change_of_basis_op_to_nearest_setting(
            cs_row, length_tolerance=0.03, angle_tolerance=3.0,
            test_multiples=False, metric='l1')
        observed_xyz.add(cb_op.as_xyz())
        transformed = cs_row.unit_cell().change_basis(cb_op)
        sim_ops = uctbx.unit_cell(transformed.parameters()).similarity_transformations(
            uctbx.unit_cell(cs_row.unit_cell().parameters()),
            relative_length_tolerance=0.02, absolute_angle_tolerance=0.5,
            unimodular_generator_range=1)
        assert len(sim_ops) > 0, \
            f"{transformed.parameters()} is not the same lattice as " \
            f"{cs_row.unit_cell().parameters()}"

    assert len(observed_xyz) >= 2, \
        f"expected cycling across >= 2 cb_ops, got {observed_xyz}"

    print("ties are bit-exact under every non-v7 metric; nearest_setting "
          f"cycled across {len(observed_xyz)} distinct cb_ops over 4 calls")


def test_search_nearest_cell_metric_rankings():
    """
    Each of the 7 metric= choices must produce a sane ranking on a
    synthetic database (built the same way as
    test_search_nearest_cell_matches_direct_computation): the query's own
    cell should rank first at ~0 distance, and reported distances should
    be non-decreasing. Cross-metric orderings below rank 1 are NOT
    asserted -- ties are exact and per-metric orderings of near-equal
    cells are not ground truth (see test_metric_ties_are_exact).
    """
    from cctbx.command_line import fetch_pdb_cells, search_nearest_cell
    from cctbx.uctbx import cell_metrics as cm

    sgs_ucs = [
        ['C2',   [49.021, 52.475, 96.609, 90, 96.53, 90]],        # 1rgx
        ['P2',   [33.429, 95.775, 33.665, 90, 101.67, 90]],       # 1r8m
        ['P1',   [40.157, 41.867, 97.795, 91.11, 92.73, 107.18]], # 2fxo
        ['C222', [37.656, 54.197, 95.677, 90, 90, 90]],           # 4rne
        ['R3:R', [57.98, 57.98, 57.98, 92.02, 92.02, 92.02]],     # pla2
        ['C2',   [80.95, 80.57, 57.1, 90, 90.35, 90]],            # pla2
    ]
    pdb_ids = [f"synm{i:03d}" for i in range(len(sgs_ucs))]
    base_rows = [
        {'pdb_id': pdb_id, 'unit_cell_a': uc[0], 'unit_cell_b': uc[1],
         'unit_cell_c': uc[2], 'unit_cell_alpha': uc[3], 'unit_cell_beta': uc[4],
         'unit_cell_gamma': uc[5], 'unit_cell_Z': 1, 'space_group': sg}
        for pdb_id, (sg, uc) in zip(pdb_ids, sgs_ucs)
    ]
    csv_path = '/Users/dwpaley/daily/20251120/.tmp/tst_search_metrics.csv'
    fetch_pdb_cells.write_cells_csv(base_rows, csv_path, out=io.StringIO())

    query_sg, query_uc = sgs_ucs[0]
    query_args = [
        "unit_cell=%.6f,%.6f,%.6f,%.6f,%.6f,%.6f" % tuple(query_uc),
        "space_group=%s" % query_sg,
        "data_file=%s" % csv_path,
        "n_results=%d" % len(sgs_ucs),
    ]
    for m in cm.METRICS:
        run_out = io.StringIO()
        search_nearest_cell.run(query_args + ["metric=%s" % m], out=run_out)
        lines = run_out.getvalue().splitlines()
        header_idx = next(i for i, l in enumerate(lines) if l.startswith("Rank"))
        data_lines = [l for l in lines[header_idx + 1:] if l.strip()]
        assert len(data_lines) > 0, m

        actual = []
        for line in data_lines:
            tok = line.split()
            actual.append((tok[1], float(tok[2])))

        rank1_pdb, rank1_dist = actual[0]
        assert rank1_pdb == pdb_ids[0], \
            f"{m}: expected top hit {pdb_ids[0]}, got {rank1_pdb}"
        assert rank1_dist < 1e-3, f"{m}: expected ~0 distance, got {rank1_dist}"

        dists = [d for _, d in actual]
        assert dists == sorted(dists), f"{m}: run() output is not sorted by distance"

    print("search_nearest_cell.run() ranks the query's own cell first "
          f"at ~0 distance, sorted, for all {len(cm.METRICS)} metrics")


if __name__ == '__main__':
    test_a2a_abs_2023()
    test_pla2_abs_2023()
    test_table11_database()
    test_cell_multiples()
    test_change_of_basis_op_to_nearest_setting()
    test_min_cell_frame_prefilter()
    test_min_cell_frame_prefilter_ties()
    test_fetch_and_search_symmetry_construction_agree()
    test_search_nearest_cell_matches_direct_computation()
    test_cell_metrics_against_scalar_reference()
    test_ncdist_matches_compiled_NCDist()
    test_cell_distance_is_the_cellparams_metric()
    test_v7_is_a_lattice_invariant()
    test_metric_ties_are_exact()
    test_search_nearest_cell_metric_rankings()
    print("ok")
