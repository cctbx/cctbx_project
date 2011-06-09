from iotbx.detectors.detectorbase import DetectorImageBase
from iotbx.detectors.raxisbase import Raxis
from iotbx.detectors import ReadRAXIS

class RAXISImage(DetectorImageBase,Raxis):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    Raxis.__init__(self,filename)
    self.vendortype = "RAXIS"

  def readHeader(self,):
    if not self.parameters:
      Raxis.readHeader(self)
      self.generic_param_from_vendor_head()

  def generic_param_from_vendor_head(self):
      self.parameters={}
      self.parameters['SIZE1'] = self.head['nSlow']
      self.parameters['SIZE2'] = self.head['nFast']
      self.parameters['CCD_IMAGE_SATURATION'] = 32767 * self.head['Ratio']
      self.parameters['PIXEL_SIZE'] = self.head['sizeFast']
      self.parameters['OSC_START'] = self.head['phistart']
      self.parameters['OSC_RANGE'] = self.head['phiend']-self.head['phistart']
      self.parameters['DISTANCE'] = self.head['distance']
      self.parameters['WAVELENGTH'] = self.head['wavelength']
      self.parameters['BEAM_CENTER_X'] = self.head['beampixels_x']*self.head[
                                                              'sizeFast']
      self.parameters['BEAM_CENTER_Y'] = self.head['beampixels_y']*self.head[
                                                              'sizeSlow']
      self.parameters['TWOTHETA'] = self.head['twotheta']
      self.parameters['DETECTOR_SN'] = self.head['operatorname']

  def dataoffset(self):
    return self.head['record_length']

  def integerdepth(self):
    return self.head['record_length']//self.head['nFast']

  def getEndian(self):
    return 1

  def read(self):
    self.fileLength()
    self.data()
    self.bin_safe_set_data( ReadRAXIS(self.CharTemp,self.integerdepth(),
         self.size1*self.bin,self.size2*self.bin,self.endian_swap_required())
    )

if __name__=='__main__':
  import sys
  i = sys.argv[1]
  a = RAXISImage(i)
  a.read()
  print a.linearintdata
  print a.linearintdata.size()
  print a.linearintdata.accessor()
  #
  #from labelit.dptbx.graphics_support import GenericImageWorker
  #W = GenericImageWorker(i)
  #W.output(sys.argv[2])
