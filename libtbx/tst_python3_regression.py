from __future__ import absolute_import, division, print_function

class pytest():
  @classmethod
  def skip(cls, arg):
    print("SKIP: %s" % arg)
  @classmethod
  def fail(cls, arg):
    import sys
    sys.exit("FAIL: %s" % arg)

def test_find_python3_violations():
  import libtbx
  import libtbx.test_utils.python3_regression as py3test
  result = py3test.find_new_python3_incompatible_code(libtbx)
  if result is None:
    pytest.skip('No python3 interpreter available')
  elif result:
    pytest.fail(result)

if __name__ == '__main__':
  test_find_python3_violations()
