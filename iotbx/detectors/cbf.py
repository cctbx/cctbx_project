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
      print self.adaptor.raster_description()
      flags = self.adaptor.transform_flags()
      print flags.transpose, flags.reverse_slow, flags.reverse_fast
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

  def read(self):
    T = self.adaptor.transform_flags()

    #future: use these tests to trigger in-place column & row swaps
    assert T.reverse_fast == False
    assert T.reverse_slow == False

    from cbflib_adaptbx import cbf_binary_adaptor
    M = cbf_binary_adaptor( self.filename )
    data = M.uncompress_implementation("buffer_based").uncompress_data()
    if T.transpose==True:
      #very inefficient; 0.09 sec for 3K x 3K uncompression
      #              yet 0.54 sec for in-place transpose
      # other options not tried: a) alloc & set new matrix; b) use data as is
      data.matrix_transpose_in_place()

    self.bin_safe_set_data(data)

if __name__=='__main__':
  import sys
  C = CBFImage(sys.argv[1])
