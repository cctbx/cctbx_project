from iotbx.detectors.raxisbase import Raxis
from iotbx.detectors.raxis import RAXISImage
from iotbx.detectors import ReadRAXIS,MakeSquareRAXIS

#Special class to handle Raxis-II images, containing rectangular pixels
# Array is resized so that the internal representation has square pixels
class NonSquareRAXISImage(RAXISImage):
  def __init__(self,filename):
    RAXISImage.__init__(self,filename)

  def readHeader(self,):
    if not self.parameters:
      Raxis.readHeader(self)

      self.np = int(round(self.head['sizeSlow']/self.head['sizeFast']*self.head['nFast']))

      # assume that the slow dimension has longer pixels, as for the Raxis-II
      assert self.head['sizeSlow'] > self.head['sizeFast']

      # assume that the original pixel array is square
      assert  self.head['nFast'] == self.head['nSlow']

      # assume that the extra pixels can be evenly divided by two
      assert (self.np-self.head['nSlow'])%2==0

      self.extra = (self.np-self.head['nSlow'])/2

      self.generic_param_from_vendor_head()
      self.generic_param_from_adapt_head()

  def generic_param_from_adapt_head(self):

      self.parameters['SIZE1'] = self.np
      self.parameters['SIZE2'] = self.np
      self.parameters['PIXEL_SIZE'] = self.head['sizeFast']
      self.parameters['BEAM_CENTER_X'] = (self.extra+self.head['beampixels_x']
                                         )*self.head['sizeFast']
      self.parameters['BEAM_CENTER_Y'] = self.head['beampixels_y']*self.head[
                                                              'sizeSlow']

  def read(self):
    self.fileLength()
    self.data()
    self.rawlinearintdata = ReadRAXIS(self.CharTemp,self.integerdepth(),
         self.head['nSlow']*self.bin,
         self.head['nFast']*self.bin,
         bool(self.getEndian()))

    self.bin_safe_set_data( MakeSquareRAXIS(self.np,self.extra,
                                         self.head['nSlow'],
                                         self.rawlinearintdata)
    )
    del self.rawlinearintdata
