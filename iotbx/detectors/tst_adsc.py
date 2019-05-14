from __future__ import absolute_import, division, print_function
from iotbx.detectors import adsc
from libtbx.test_utils import approx_equal
import os
import bz2

adsc_file = 'adsc.img'
adsc_file_bz2 = 'adsc.img.bz2'
# file originally from
# http://cci.lbl.gov/cctbx_downloads/regression/iotbx/adsc.img

adsc_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), adsc_file)
adsc_file_bz2 = os.path.join(os.path.dirname(os.path.realpath(__file__)), adsc_file_bz2)

def exercise_adscread():
  if (not os.path.isfile(adsc_file)):
    if (not os.path.isfile(adsc_file_bz2)):
      print("Skipping exercise_adscread(): input file not available")
      return
    with open(adsc_file, 'wb') as orig, bz2.BZ2File(adsc_file_bz2, 'rb') as comp:
      for data in iter(lambda : comp.read(100 * 1024), b''):
        orig.write(data)

  a = adsc.ADSCImage(adsc_file)
  a.read()
  assert a.size1 == 2304
  assert a.size2 == 2304
  assert a.npixels == 5308416
  assert approx_equal(a.pixel_size, 0.0816)
  assert a.saturation == 65535
  assert approx_equal(a.osc_start, 0)

def run():
  exercise_adscread()
  print("OK")

if __name__=="__main__":
  run()
