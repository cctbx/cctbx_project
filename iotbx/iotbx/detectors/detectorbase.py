from iotbx.detectors import ReadADSC

class DetectorImageBase(object):
  def __init__(self,filename):
    self.filename=filename
    self.parameters=None
    self.linearintdata=None
    self.bin=1
    self.vendortype = "baseclass"

  def setBin(self,bin): #software binning.
                        # the only bin values supported are 1 & 2
    if self.bin!=1 or bin!=2: return
    if self.size1%bin!=0: return
    self.parameters['SIZE1']=self.parameters['SIZE1']/bin
    self.parameters['SIZE2']=self.parameters['SIZE2']/bin
    if self.parameters.has_key('CCD_IMAGE_SATURATION'):
      self.parameters['CCD_IMAGE_SATURATION']=self.parameters['CCD_IMAGE_SATURATION']*bin*bin
    self.parameters['PIXEL_SIZE']=self.parameters['PIXEL_SIZE']*bin
    self.bin = bin

  def fileLength(self):
    self.readHeader()
    return self.dataoffset()+self.size1*self.size2*self.integerdepth()
    # dataoffset() and integerdepth() must be defined in derived class
    # pure supposition:
    #  size1 corresponds to number of rows.  Columns are slow.
    #  size2 corresponds to number of columns.  Rows are fast.

  def getEndian(self): pass
    # must be defined in derived class

  def endian_swap_required(self):
    data_is_big_endian = self.getEndian()
    import struct
    platform_is_big_endian = (
      struct.unpack('i',struct.pack('>i',3000))[0] == 3000
    )
    return data_is_big_endian != platform_is_big_endian

  def read(self):
    self.fileLength()
    self.linearintdata = ReadADSC(self.filename,self.dataoffset(),
         self.size1*self.bin,self.size2*self.bin,self.getEndian())
    if self.bin==2:
      from iotbx.detectors import Bin2_by_2
      self.linearintdata = Bin2_by_2(self.linearintdata)

  def debug_write(self,fileout,mod_data=None):
    if not self.parameters.has_key("TWOTHETA"):
      self.parameters["TWOTHETA"]=0.0
    if self.getEndian()==1:
      self.parameters["BYTE_ORDER"]="big_endian"
    else:
      self.parameters["BYTE_ORDER"]="little_endian"
    info = """{
HEADER_BYTES= 1024;
DIM=2;
BYTE_ORDER=%(BYTE_ORDER)s;
TYPE=unsigned_short;
SIZE1=%(SIZE1)4d;
SIZE2=%(SIZE1)4d;
PIXEL_SIZE=%(PIXEL_SIZE)8.6f;
TIME=0.000000;
DISTANCE=%(DISTANCE).2f;
TWOTHETA=%(TWOTHETA).2f;
PHI=%(OSC_START).3f;
OSC_START=%(OSC_START).3f;
OSC_RANGE=%(OSC_RANGE).3f;
WAVELENGTH=%(WAVELENGTH).6f;
BEAM_CENTER_X=%(BEAM_CENTER_X).2f;
BEAM_CENTER_Y=%(BEAM_CENTER_Y).2f;
CCD_IMAGE_SATURATION=65535;
}\f"""%self.parameters
    F = open(fileout,"wb")
    F.write(info)
    len_null=1024-len(info)
    F.write('\0'*len_null)
    F.close()
    from iotbx.detectors import WriteADSC
    if mod_data==None: mod_data=self.linearintdata
    WriteADSC(fileout,mod_data,self.size1,self.size2,self.getEndian())

  def __getattr__(self, attr):
    if   attr=='size1' : return self.parameters['SIZE1']
    elif attr=='size2' : return self.parameters['SIZE2']
    elif attr=='npixels' : return self.parameters['SIZE1'] * self.parameters['SIZE2']
    elif attr=='saturation' : return self.parameters['CCD_IMAGE_SATURATION']
    elif attr=='rawdata' : return self.linearintdata
    elif attr=='pixel_size' : return self.parameters['PIXEL_SIZE']
    elif attr=='osc_start' : return self.parameters['OSC_START']
    elif attr=='distance' : return self.parameters['DISTANCE']
    elif attr=='wavelength' : return self.parameters['WAVELENGTH']
    elif attr=='beamx' : return self.parameters['BEAM_CENTER_X']
    elif attr=='beamy' : return self.parameters['BEAM_CENTER_Y']
    elif attr=='deltaphi' : return self.parameters['OSC_RANGE']
    elif attr=='twotheta' : return self.parameters['TWOTHETA']
    elif attr=='serial_number' : return self.parameters['DETECTOR_SN']
