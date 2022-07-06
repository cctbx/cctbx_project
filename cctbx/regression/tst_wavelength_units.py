from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_text

def exercise_1():
  fb = easy_run.fully_buffered("cctbx.wavelength_units")
  assert_lines_in_text('\n'.join(fb.stderr_lines), "Usage: cctbx.wavelength_units 1A|1keV [...]")

  fb = easy_run.fully_buffered("cctbx.wavelength_units 2.3a")
  assert_lines_in_text('\n'.join(fb.stdout_lines), "2.30000 A = 5.39062 keV")

  fb = easy_run.fully_buffered("cctbx.wavelength_units 9.8kev")
  assert_lines_in_text('\n'.join(fb.stdout_lines), "1.26514 A = 9.80000 keV")

if (__name__ == "__main__"):
  exercise_1()
  print("OK")
