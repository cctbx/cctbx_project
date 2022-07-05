from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file

def exercise_1(prefix='tst_wavelength_units'):
  easy_run.call("cctbx.wavelength_units &> %s" % prefix)
  assert_lines_in_file('%s' % prefix, "Usage: cctbx.wavelength_units 1A|1keV [...]")

  easy_run.call("cctbx.wavelength_units 2.3a &> %s" % prefix)
  assert_lines_in_file('%s' % prefix, "2.30000 A = 5.39062 keV")

  easy_run.call("cctbx.wavelength_units 9.8kev &> %s" % prefix)
  assert_lines_in_file('%s' % prefix, "1.26514 A = 9.80000 keV")

if (__name__ == "__main__"):
  exercise_1()
  print("OK")
