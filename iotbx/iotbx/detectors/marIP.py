import re,struct
from iotbx.detectors.detectorbase import DetectorImageBase

class NullAdaptor:
  def size1(self): return 0
  def size2(self): return 0
  def overload(self): return 0
  def pixel_size(self): return 0.0
  def osc_start(self): return 0.0
  def osc_range(self): return 0.0
  def distance(self): return 0.0
  def wavelength(self): return 0.0
  def twotheta(self): return 0.0
  def rawdata(self): return [0.0]

class MARIPImage(DetectorImageBase):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    try:
      from cbflib_ext import Mar345Adaptor # optional package
      self.adaptor = Mar345Adaptor(filename)
    except:
      self.adaptor = NullAdaptor()
    self.vendortype = "MARIP"

  def readHeader(self):
    self.parameters = {'SIZE1':self.adaptor.size1(),
                       'SIZE2':self.adaptor.size2(),
                       'CCD_IMAGE_SATURATION':self.adaptor.overload(),
                       'PIXEL_SIZE':self.adaptor.pixel_size(),
                       'OSC_START':self.adaptor.osc_start(),
                       'DISTANCE':self.adaptor.distance(),
                       'WAVELENGTH':self.adaptor.wavelength(),
        'BEAM_CENTER_X':self.adaptor.size1()*self.adaptor.pixel_size()/2.0,
        'BEAM_CENTER_Y':self.adaptor.size2()*self.adaptor.pixel_size()/2.0,
                       'OSC_RANGE':self.adaptor.osc_range(),
                       'TWOTHETA':self.adaptor.twotheta(),
                       'DETECTOR_SN':0
                       }

  def fileLength(self):
    return 0

  def getEndian(self):
    return 0

  def read(self):
    self.linearintdata = self.adaptor.rawdata()
    if self.bin==2:
      from iotbx.detectors import Bin2_by_2
      self.linearintdata = Bin2_by_2(self.linearintdata)

  def dataoffset(self):
    return 0

  def integerdepth(self):
    return 0

if __name__=='__main__':
  pass
