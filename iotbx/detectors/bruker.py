from __future__ import absolute_import, division, print_function
import math
from iotbx.detectors import Bruker_base
from iotbx.detectors.detectorbase import DetectorImageBase

class BrukerImage(DetectorImageBase):
  def __init__(self,filename):
    self.filename=filename
    self.bruker = Bruker_base(filename)
    self.bin=1
    self.bin_safe_set_data(self.bruker.linearintdata())
    self.vendortype = "Bruker Proteus CCD"
    self.readHeader()

  def readHeader(self):
    self.parameters={
      # image saturation not available.  take the maximum pixel instead
      'CCD_IMAGE_SATURATION': self.bruker.ccd_image_saturation,
      'SIZE1':1024,
      'SIZE2':1024,
      'PIXEL_SIZE':self.bruker.pixel_size,
      'OSC_START':self.bruker.osc_start,
      'DISTANCE':10.0*(self.bruker.distance_cm+self.bruker.distance_delta),
      'WAVELENGTH':self.bruker.wavelength,

      #first attempt formula to support two theta offset.
      # But the resulting distl_overlay did not put the beam in exactly
      # the right position.
      'BEAM_CENTER_X':self.bruker.pixel_size*self.bruker.centerx_pix +
                      10.*self.bruker.distance_cm *
                      math.tan (math.pi*self.bruker.twotheta/180.),
      'BEAM_CENTER_Y':self.bruker.pixel_size*self.bruker.centery_pix,
      'OSC_RANGE':self.bruker.osc_range,
      'TWOTHETA':self.bruker.twotheta,
      'DETECTOR_SN':0,
    }

  def read(self):
    return

if __name__=="__main__":
  import sys
  if len(sys.argv)<2:
    pass#file = "/net/cci/dials/from_adder/sauter/rawdata/mckee/bruker/lyziph6p5_01_0001.sfrm"
  else:
    file = sys.argv[1]
  B = BrukerImage(file)
  print("B.vendortype,B.saturation",B.vendortype,B.saturation)
  print("pixel size",B.pixel_size)
  print("osc start",B.osc_start, end=' ')
  print("distance",B.distance)
  print("wavelength",B.wavelength)
  print("beamx",B.beamx)
  print("beamy",B.beamy)
  print("deltaphi",B.deltaphi)
  print("twotheta",B.twotheta)
  #for item in B.linearintdata:
  #  print item
