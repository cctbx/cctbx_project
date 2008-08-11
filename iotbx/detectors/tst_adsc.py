from iotbx.detectors import adsc
from libtbx.test_utils import approx_equal
import urllib
import os

def get_test_files():
  try:
    urllib.urlretrieve(
      'http://cci.lbl.gov/cctbx_downloads/regression/iotbx/adsc.img',
      'adsc.img')
  except IOError:
    pass

def exercise_adscread():
  if (not os.path.isfile("adsc.img")):
    print "Skipping exercise_adscread(): input file not available"
    return
  a = adsc.ADSCImage('adsc.img')
  a.read()
  assert a.size1 == 2304
  assert a.size2 == 2304
  assert a.npixels == 5308416
  assert approx_equal(a.pixel_size, 0.0816)
  assert a.saturation == 65535
  assert approx_equal(a.osc_start, 0)

def run():
  get_test_files()
  exercise_adscread()
  print "OK"

if __name__=="__main__":
  run()
