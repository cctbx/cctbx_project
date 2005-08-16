from iotbx.detectors.raxisbase import Raxis
from iotbx.detectors.raxis import RAXISImage
from iotbx.detectors import ReadRAXIS,MakeSquareRAXIS
from scitbx.array_family import flex

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

      print "readHeader of non-square image"
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
         self.getEndian())

    self.linearintdata = MakeSquareRAXIS(self.np,self.extra,
                                         self.head['nSlow'],
                                         self.rawlinearintdata)
    '''Reference Python implementation for naive array resizing
    self.linearintdata = flex.int(flex.grid(self.np,self.np))

    for slow in xrange(2*self.extra):  #padding at beginning
      for fast in xrange(self.np):
        self.linearintdata[slow*self.np+fast]=fast
    for slow in xrange(self.head['nSlow']):
      for fast in xrange(self.extra):
        self.linearintdata[(2*self.extra+slow)*self.np+fast]=0
      for fast in xrange(self.head['nFast']):
        self.linearintdata[(2*self.extra+slow)*self.np+(self.extra+fast)]=(self.rawlinearintdata[
          slow*self.head['nFast']+fast])
      for fast in xrange(self.extra):
        self.linearintdata[(2*self.extra+slow)*self.np+(self.extra+self.head['nFast']+fast)]=0

    #Now rearrange the data in-place
    for newslow in xrange(self.np):
      print newslow,
      float_begin = (float(newslow)/self.np)*float(self.head['nSlow'])+(2*self.extra)
      float_end = (float(newslow+1)/self.np)*float(self.head['nSlow'])+(2*self.extra)
      signatures = []
      for xf in xrange(int(float_begin), 1+int(float_end)):
        #xf_min
        if float_begin<float(xf): xf_min=float(xf)
        else: xf_min=float_begin

        #xf_max
        if float_end>float(xf+1): xf_max=float(xf+1)
        else: xf_max=float_end

        if xf_max-xf_min>1.e-7: signatures.append((xf,xf_max-xf_min))

      print float_begin,float_end,signatures
      row_temp = flex.int(flex.grid(self.np))
      for item in signatures:
        for fast in xrange(self.np):
          row_temp[fast]+=int(round(self.linearintdata[item[0]*self.np+fast]*item[1]))
      for fast in xrange(self.np):
        self.linearintdata[newslow*self.np+fast]=row_temp[fast]
    '''
    if self.bin==2:
      from iotbx.detectors import Bin2_by_2
      self.linearintdata = Bin2_by_2(self.linearintdata)
    del self.rawlinearintdata
