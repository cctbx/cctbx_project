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


if __name__ == '__main__':
    test_a2a_abs_2023()
    test_pla2_abs_2023()
    test_table11_database()
    print("ok")
