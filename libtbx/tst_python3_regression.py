from __future__ import absolute_import, division, print_function

def test_find_python3_violations():
  import libtbx
  import libtbx.test_utils.python3_regression as py3test
  result = py3test.find_new_python3_incompatible_code(libtbx)
  if result is None:
    print("SKIP: No python3 interpreter available")
  elif result:
    import sys
    sys.exit("FAIL: %s" % result)

if __name__ == '__main__':
  test_find_python3_violations()
