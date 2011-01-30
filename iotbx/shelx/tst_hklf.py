from libtbx.test_utils import approx_equal
from libtbx.test_utils import Exception_expected
from iotbx.shelx import hklf
from cStringIO import StringIO

def exercise_hklf_reader():
  s = ('   1   2  -1  -23.34    4.56   1\n'
       '   2  -3   9   12.45    6.12   2\r\n'
       '99999999999999999.9999999.999999\n'
       '-999-999-999-9999.99-9999.99-999\r\n'
       '   0   0   0    0.00    0.00   0\n')
  r = hklf.reader(file_object=StringIO(s))
  assert list(r.indices()) == [
    (1, 2, -1), (2, -3, 9), (9999, 9999, 9999), (-999, -999, -999), ]
  assert approx_equal(r.data(), [-23.34, 12.45, 99999.99, -9999.99, ])
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

  s = ('   0   2   3 1816.00   20.00\n'
       '   0   2   415508.00  138.00\n'
       '   0   2   5 4776.00   40.00\n')
  r = hklf.reader(file_object=StringIO(s))
  assert list(r.indices()) == [ (0,2,3), (0,2,4), (0,2,5) ]
  assert r.batch_numbers() is None
  assert r.alphas() is None
  assert r.wavelengths() is None

  for end_line in (True, False):
    for something_after_end_line in (True, False):
      s = ('   1   2  -1   23.34    4.56\n'
           '   2  -3   9   12.45    6.12\r\n'
           '99999999999999999.9999999.99\n'
           '-999-999-999-9999.99-9999.99\n')
      if end_line:
        s += '   0   0   0    0.00    0.00\n'
        if something_after_end_line:
          s += '  -5   1   0  123.45   66.12\n'
      s = (s)
      r = hklf.reader(file_object=StringIO(s))
      assert approx_equal(r.sigmas(), [4.56, 6.12, 99999.99, -9999.99, ])
      assert r.batch_numbers() is None
      assert r.alphas() is None
      assert r.wavelengths() is None

  s = ''
  try: r = hklf.reader(file_object=StringIO(s))
  except RuntimeError: pass
  else: raise Exception_expected

  s = ('   1   2  -1   23.34    4.56   1\n'
       '   2  -3  a9   12.45    6.12   2\n'
       '   0   0   0    0.00    0.00   0\n')
  try: r = hklf.reader(file_object=StringIO(s))
  except Exception: pass
  else: raise Exception_expected

  s = ('   1   2  -1   23.34    4.56   1   45.36\n'
       '   2  -3   9  -12.45   -6.12   2   54.63\n'
       '   0   0   0    0.00    0.00   0       0\n')
  r = hklf.reader(file_object=StringIO(s))
  assert list(r.indices()) == [ (1, 2, -1), (2, -3, 9)]
  assert approx_equal(r.data(), [23.34, -12.45])
  assert approx_equal(r.sigmas(), [4.56, -6.12])
  assert approx_equal(r.batch_numbers(), [1, 2])
  assert approx_equal(r.wavelengths(), [45.36, 54.63])

  s = ('   1   2  -1     23.      4.   2\n'
       '  -2   1   3     -1.      3.   1\n'
       '   0   0   0      0.      0.   0\n')
  r = hklf.reader(file_object=StringIO(s))
  assert list(r.indices()) == [(1, 2, -1), (-2, 1, 3)]
  assert approx_equal(r.data(), [23, -1])
  assert approx_equal(r.sigmas(), [4, 3])
  assert approx_equal(r.alphas(), [2, 1])
  assert r.wavelengths() is None

  s = ('   3   2  -1     32.      5.\n'
       '   0   0   0\n')
  r = hklf.reader(file_object=StringIO(s))
  assert list(r.indices()) == [(3, 2, -1)]
  assert approx_equal(r.data(), [32])
  assert approx_equal(r.sigmas(), [5])
  assert r.alphas() is None
  assert r.wavelengths() is None

  s = (
"""King Arthur: [after Arthur's cut off both of the Black Knight's arms]
  Look, you stupid Bastard. You've got no arms left.
Black Knight: Yes I have.
King Arthur: *Look*!
Black Knight: It's just a flesh wound.""")
  try: r = hklf.reader(file_object=StringIO(s))
  except RuntimeError: pass
  else: raise Exception_expected

def exercise_miller_export_as_shelx_hklf():
  s = """\
   1   2  -1   23.34    4.56
   2  -3   9   12.45    6.12
99999999999999999.9999999.99
-999-999-999-9999.99-9999.99
"""
  ma = hklf.reader(file_object=StringIO(s)).as_miller_arrays()[0]
  so = StringIO()
  ma.export_as_shelx_hklf(file_object=so)
  ma2 = hklf.reader(file_object=StringIO(so.getvalue())).as_miller_arrays()[0]
  assert approx_equal(ma.indices(), ma2.indices())
  assert approx_equal(ma.data(), ma2.data())
  assert approx_equal(ma.sigmas(), ma2.sigmas())

def run():
  exercise_hklf_reader()
  exercise_miller_export_as_shelx_hklf()
  print "OK"

if __name__ == '__main__':
  run()
