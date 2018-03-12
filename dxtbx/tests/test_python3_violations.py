def test_find_python3_violations():
  import pytest
  try:
    import dials.test.python3_regression as py3test
  except ImportError:
    pytest.skip('DIALS required for this test')
  import dxtbx
  result = py3test.find_new_python3_incompatible_code(dxtbx)
  if result is None:
    pytest.skip('No python3 interpreter available')
  elif result:
    pytest.fail(result)
