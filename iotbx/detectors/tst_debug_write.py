from __future__ import absolute_import, division, print_function
from six.moves import range
import os,math
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from iotbx.detectors import adsc
from iotbx.detectors.detectorbase import DetectorImageBase

filesize=64
filename='debug_wtest_001.img'

def sinfunc(x):
  return math.sin(math.pi*x/32)

def exercise_debug_write():
  #---------Write out a test file
  D = DetectorImageBase('no_file')
  D.parameters = {'SIZE1':filesize,
                  'SIZE2':filesize,
                  'PIXEL_SIZE':0.1,
                  'DISTANCE':100.0,
                  'TWOTHETA':0.0,
                  'OSC_START':0.0,
                  'OSC_RANGE':1.0,
                  'WAVELENGTH':1.0,
                  'BEAM_CENTER_X':12.5,
                  'BEAM_CENTER_Y':12.5,
  }
  def getEndian():
    return 0
  D.getEndian = getEndian
  sindata = flex.int()
  for x in range(filesize):
    sindata.append(min(255,abs(int(256*sinfunc(x)))))
  moddata = flex.int()
  accum = 0
  for x in range(filesize):
    for y in range(filesize):
      value = sindata[x]*sindata[y]
      moddata.append(value)
      accum+=value
  D.debug_write(filename,mod_data=moddata)

  #---------Read back the test file
  a = adsc.ADSCImage(filename)
  a.read()
  assert a.size1 == filesize
  assert a.size2 == filesize
  assert a.npixels == filesize*filesize
  assert approx_equal(a.pixel_size, 0.1)
  assert approx_equal(a.osc_start, 0)
  checkaccum = 0
  for x in range(a.npixels):
    if moddata[x]!=a.linearintdata[x]:
      print(x,moddata[x],a.linearintdata[x])
    checkaccum+=a.linearintdata[x]
  assert accum == checkaccum
  os.remove(filename)

def run():
  exercise_debug_write()
  print("OK")

if __name__=="__main__":
  run()
