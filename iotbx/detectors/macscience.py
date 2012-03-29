import struct,os
from iotbx.detectors.detectorbase import DetectorImageBase
from iotbx.detectors import ReadDIP

class DIPImage(DetectorImageBase):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.filename = filename
    self.vendortype = "MacScience"
    self.getImageSize()

  def getImageSize(self):
    self.filesize = os.stat(self.filename)[6]
    # MacScience DIP 2030b
    assert self.filesize == 18001024 #necessary filelength for 3000x3000 image
    # other file sizes can be added here (for compatibility with
    # MacScience DIP 2020; 2000 x 2000)
    self.deduce_size = 3000

  def fileLength(self): return self.filesize

  def readHeader(self,):
    headerstart = self.deduce_size * self.deduce_size * 2 #short integers
    F = open(self.filename,'rb')
    F.seek(headerstart)
    self.rawheader = F.read(1024)
    F.close()
    assert self.rawheader[0:3] # Image file self-identifies as a DIP
    #self.getEndian() not cached anyway

    if not self.parameters:
      self.parameters={'CCD_IMAGE_SATURATION':1048576} #according to Stan Svenson
      self.parameters['SIZE1']=struct.unpack(
             self.intformat(),self.rawheader[28:32])[0]
      self.parameters['SIZE2']=struct.unpack(
             self.intformat(),self.rawheader[32:36])[0]
      assert self.size1==self.deduce_size
      assert self.size2==self.deduce_size
      self.parameters['PIXEL_SIZE'] = 0.1 # mm, just a guess
      self.parameters['OSC_START']=struct.unpack(
             self.floatformat(),self.rawheader[192:196])[0]
      self.parameters['DISTANCE']=struct.unpack(
             self.floatformat(),self.rawheader[124:128])[0]
      self.parameters['WAVELENGTH']=struct.unpack(
             self.floatformat(),self.rawheader[120:124])[0]
      self.parameters['BEAM_CENTER_X']=struct.unpack(
             self.intformat(),self.rawheader[164:168])[0] * self.pixel_size
      self.parameters['BEAM_CENTER_Y']=struct.unpack(
             self.intformat(),self.rawheader[168:172])[0] * self.pixel_size
      self.parameters['OSC_RANGE']=struct.unpack(
             self.floatformat(),self.rawheader[196:200])[0] - self.osc_start
      self.parameters['TWOTHETA']=0.0

  def integerdepth(self):
    return 2

  def floatformat(self):
    if self.getEndian(): return '>f'
    else: return '<f'

  def intformat(self):
    if self.getEndian(): return '>i'
    else: return '<i'

  def getEndian(self):
    return struct.unpack('>i',self.rawheader[4:8])[0]==self.deduce_size
    # True for big_endian; False for little_endian

  def endian_swap_required(self):
    return struct.unpack('i',self.rawheader[4:8])[0]!=self.deduce_size

  def read(self):
    self.bin_safe_set_data( ReadDIP(self.filename,
         self.size1*self.bin,self.size2*self.bin,self.endian_swap_required())
         )
