from __future__ import absolute_import, division, print_function

from libtbx.version import get_version

# =============================================================================
def test_version():

  assert get_version() is not None
  assert get_version(filename='random_nonexisting_filename') is not None
  assert get_version(filename='random_nonexisting_filename', fail_with_none=True) is None

# =============================================================================
if __name__ == '__main__':
  test_version()

# =============================================================================
# end
