from __future__ import division, print_function, absolute_import

import pytest
import numpy as np

from scitbx.array_family.flex import vec3_double

def test_vec3_double_as_numpy_array():
    test_data = [
        [0, 1, 30.5],
        [-4.0, 2.0, 300000],
    ]
    vec3 = vec3_double(test_data)
    np_vec3 = vec3.as_numpy_array()
    assert np.all(np.isclose(np_vec3, np.array(test_data)))
