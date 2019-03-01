from __future__ import absolute_import, division, print_function

import pytest


def test_find_python3_violations():
    import dxtbx
    import libtbx.test_utils.python3_regression as py3test

    result = py3test.find_new_python3_incompatible_code(dxtbx)
    if result is None:
        pytest.skip("No python3 interpreter available")
    elif result:
        pytest.fail(result)
