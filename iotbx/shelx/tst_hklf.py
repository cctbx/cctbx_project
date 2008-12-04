from libtbx.test_utils import approx_equal, Exception_expected
from iotbx.shelx import hklf
import tempfile

def exercise_fast_hkl_reading():
  filename = tempfile.mktemp()
  f = open(filename, 'w')
  f.write('   1   2  -1   23.34    4.56   1\n'
          '   2  -3   9   12.45    6.12   2\n'
          '99999999999999999.9999999.999999\n'
          '-999-999-999-9999.99-9999.99-999\n'
          '   0   0   0    0.00    0.00   0\n')
  f.close()
  fast = hklf.reader(filename=filename)
  slow = hklf.python_reader(filename=filename)
  for r in (fast, slow):
    assert list(r.indices()) == [ (1, 2, -1), (2, -3, 9), (9999, 9999, 9999), (-999, -999, -999), ]
    assert approx_equal(r.data(), [23.34, 12.45, 99999.99, -9999.99, ])
    assert approx_equal(r.sigmas(), [4.56, 6.12, 99999.99, -9999.99, ])
    assert approx_equal(r.batch_numbers(), [1, 2, 9999, -999, ])
    assert approx_equal(r.alphas(), [1, 2, 9999, -999, ])
    for ma in r.as_miller_arrays():
      assert ma.indices().all_eq(r.indices())
    ma = r.as_miller_arrays()[0]
    assert ma.data().all_approx_equal(r.data())
    assert ma.sigmas().all_approx_equal(r.sigmas())
    ma = r.as_miller_arrays()[1]
    assert ma.data().all_eq(r.alphas())
    assert ma.sigmas() is None

  f = open(filename, 'w')
  f.write('   1   2  -1   23.34    4.56\n'
          '   2  -3   9   12.45    6.12\n'
          '99999999999999999.9999999.99\n'
          '-999-999-999-9999.99-9999.99\n'
          '   0   0   0    0.00    0.00\n')
  f.close()
  r = hklf.reader(filename=filename)
  assert approx_equal(r.sigmas(), [4.56, 6.12, 99999.99, -9999.99, ])

  f = open(filename, 'w')
  f.close()
  try: fast = hklf.reader(filename=filename)
  except RuntimeError: pass
  else: raise Exception_expected

  f = open(filename, 'w')
  f.write('   1   2  -1   23.34    4.56   1     45.36\n'
          '   2  -3   9   12.45    6.12   2     45.36\n'
          '   0   0   0    0.00    0.00   0         0\n')
  f.close()
  try: fast = hklf.reader(filename=filename)
  except RuntimeError: pass
  else: raise Exception_expected

  r = hklf.reader(filename=filename, strict=False)
  assert list(r.indices()) == [ (1, 2, -1), (2, -3, 9), ]
  assert approx_equal(r.data(), [23.34, 12.45, ])
  assert approx_equal(r.sigmas(), [4.56, 6.12, ])
  assert approx_equal(r.batch_numbers(), [1, 2, ])

  f = open(filename, 'w')
  f.write("""King Arthur: [after Arthur's cut off both of the Black Knight's arms]
  Look, you stupid Bastard. You've got no arms left.
Black Knight: Yes I have.
King Arthur: *Look*!
Black Knight: It's just a flesh wound.""")
  f.close()
  try: fast = hklf.reader(filename=filename)
  except RuntimeError: pass
  else: raise Exception_expected

  try: fast = hklf.reader(filename=filename, strict=False)
  except RuntimeError: pass
  else: raise Exception_expected

def run():
  exercise_fast_hkl_reading()
  print "OK"

if __name__ == '__main__':
  run()
