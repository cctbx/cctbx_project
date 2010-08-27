import re,struct
from iotbx.detectors.detectorbase import DetectorImageBase
from iotbx.detectors.marIP import MARIPImage

class CBFImage(MARIPImage):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    try:
      from cbflib_adaptbx import CBFAdaptor # optional package
      self.adaptor = CBFAdaptor(filename)
      # for testing only:
      '''
      self.adaptor.read_header()
      print "     size1:", self.adaptor.size1()
      print "     size2:", self.adaptor.size2()
      print "  overload:", self.adaptor.overload()
      print "wavelength:", self.adaptor.wavelength()
      print "  distance:", self.adaptor.distance()
      print "pixel size:", self.adaptor.pixel_size()
      print "     beamx:", self.adaptor.beam_index_slow*self.adaptor.pixel_size()
      print "     beamy:", self.adaptor.beam_index_fast*self.adaptor.pixel_size()
      print " osc_start:", self.adaptor.osc_start()
      print " osc_range:", self.adaptor.osc_range()
      self.adaptor.read_data()
      '''
    except:
      from iotbx.detectors.marIP import NullAdaptor
      self.adaptor = NullAdaptor()
    self.vendortype = "CBF"
    self.readHeader()

  def beam_center_slow(self):
    return self.adaptor.beam_index_slow*self.adaptor.pixel_size()

  def beam_center_fast(self):
    return self.adaptor.beam_index_fast*self.adaptor.pixel_size()

if __name__=='__main__':
  import sys
  C = CBFImage(sys.argv[1])
