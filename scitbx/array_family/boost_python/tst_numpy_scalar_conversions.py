from __future__ import absolute_import, division, print_function
import sys

def run(args):
  assert len(args) == 0
  try:
    import numpy as np
  except ImportError:
    print("numpy not available, skipping")
    print("OK")
    return
  from scitbx.array_family import flex

  exercise_original_reproducer(np, flex)
  exercise_element_setitem(np, flex)
  exercise_element_iadd(np, flex)
  exercise_array_iadd(np, flex)
  exercise_construction_from_numpy_scalars(np, flex)
  exercise_value_preservation(np, flex)
  print("OK")

def exercise_original_reproducer(np, flex):
  """Reproducer from https://github.com/cctbx/cctbx_project/issues/1084"""
  arr = flex.double(10)
  arr[0] += np.array([1,2,3], dtype=np.float32)[0]
  assert arr[0] == 1.0

def exercise_element_setitem(np, flex):
  """Test arr[i] = numpy_scalar for various type combinations."""
  # float scalars -> flex.double
  a = flex.double(1)
  for val in [np.float32(3.5), np.float64(3.5)]:
    a[0] = val
    assert a[0] == 3.5, (a[0], type(val))
  # integer scalars -> flex.double (implicit widening)
  for val in [np.int32(7), np.int64(7)]:
    a[0] = val
    assert a[0] == 7.0, (a[0], type(val))
  # float scalars -> flex.float
  b = flex.float(1)
  for val in [np.float32(2.5), np.float64(2.5)]:
    b[0] = val
    assert abs(b[0] - 2.5) < 1e-6, (b[0], type(val))
  # integer scalars -> flex.int
  c = flex.int(1)
  for val in [np.int32(42), np.int64(42)]:
    c[0] = val
    assert c[0] == 42, (c[0], type(val))

def exercise_element_iadd(np, flex):
  """Test arr[i] += numpy_scalar for various type combinations."""
  a = flex.double([10.0])
  a[0] += np.float32(2.5)
  assert a[0] == 12.5
  a[0] += np.float64(1.0)
  assert a[0] == 13.5
  a[0] += np.int32(1)
  assert a[0] == 14.5
  a[0] += np.int64(1)
  assert a[0] == 15.5

  b = flex.int([10])
  b[0] += np.int32(5)
  assert b[0] == 15
  b[0] += np.int64(3)
  assert b[0] == 18

def exercise_array_iadd(np, flex):
  """Test arr += numpy_scalar (whole-array operations)."""
  a = flex.double([1.0, 2.0, 3.0])
  a += np.float32(10.0)
  assert list(a) == [11.0, 12.0, 13.0]
  a += np.float64(1.0)
  assert list(a) == [12.0, 13.0, 14.0]
  a += np.int32(1)
  assert list(a) == [13.0, 14.0, 15.0]

  b = flex.int([1, 2, 3])
  b += np.int32(10)
  assert list(b) == [11, 12, 13]
  b += np.int64(1)
  assert list(b) == [12, 13, 14]

def exercise_construction_from_numpy_scalars(np, flex):
  """Test constructing flex arrays from lists containing numpy scalars."""
  a = flex.double([np.float32(1.0), np.float64(2.0), np.float32(3.0)])
  assert list(a) == [1.0, 2.0, 3.0]
  b = flex.int([np.int32(1), np.int32(2), np.int32(3)])
  assert list(b) == [1, 2, 3]

def exercise_value_preservation(np, flex):
  """Test that values are preserved accurately through conversion."""
  # float32 has ~7 decimal digits of precision
  a = flex.double(1)
  a[0] = np.float32(1.23456789)
  # float32 truncates, so check against float32 precision
  assert abs(a[0] - float(np.float32(1.23456789))) < 1e-10

  # float64 should be exact for representable values
  a[0] = np.float64(1.234567890123456)
  assert a[0] == 1.234567890123456

  # Large integers
  b = flex.int(1)
  b[0] = np.int32(2147483647)  # INT32_MAX
  assert b[0] == 2147483647

if __name__ == "__main__":
  run(args=sys.argv[1:])
