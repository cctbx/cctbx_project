import copy,re
from iotbx.detectors.detectorbase import DetectorImageBase
from iotbx.detectors import ImageException


class PilatusImage(DetectorImageBase):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "Pilatus"

  mandatory_keys = ['PIXEL_SIZE_UNITS', 'DISTANCE', 'PHI', 'WAVELENGTH', 'SIZE1',
    'SIZE2', 'TWOTHETA', 'DISTANCE_UNITS', 'OSC_RANGE',
    'BEAM_CENTER_X', 'BEAM_CENTER_Y',
    'CCD_IMAGE_SATURATION', 'OSC_START', 'DETECTOR_SN', 'PIXEL_SIZE',
    'AXIS']

  def fileLength(self):
    raise ImageException("file length not computed for miniCBF")

  def getEndian(self):
    raise ImageException("endian-ness not computed for miniCBF")

  def endian_swap_required(self):
    return False

  def read(self,algorithm="cbflib_adaptbx_optimized"):
    self.readHeader()
    if self.linearintdata != None and\
      self.linearintdata.size()==self.size1*self.size2:
      #data has already been read
      return
    if self.bin==2:
      raise ImageException("2-by-2 binning not supported for miniCBF")
    try:
      from cbflib_adaptbx import MiniCBFAdaptor # optional package
      self.adaptor = MiniCBFAdaptor(self.filename)
      assert algorithm in ["cbflib","cbflib_adaptbx_optimized"]
      if algorithm=="cbflib_adaptbx_optimized":
        self.bin_safe_set_data( self.adaptor.optimized_read_data(self.size1,self.size2) )
      else:
        self.bin_safe_set_data( self.adaptor.read_data(self.size1,self.size2) )
    except:
      raise ImageException("unable to read miniCBF data; contact authors")

  def readHeader(self,maxlength=12288): # usually 1024 is OK; require 12288 for ID19
    if not self.parameters:
      rawdata = open(self.filename,"rb").read(maxlength)
      # The tag _array_data.header_convention "SLS_1.0" could be with/without quotes "..."
      pattern = re.compile(r'''_array_data.header_convention[ "]*SLS''')
      match = pattern.findall(rawdata)
      assert len(match)>=1

      # read SLS header
      headeropen = rawdata.index("_array_data.header_contents")
      headerclose= rawdata.index("_array_data.data")
      self.header = rawdata[headeropen+1:headerclose]
      self.headerlines = [x.strip() for x in self.header.split("#")]
      for idx in xrange(len(self.headerlines)):
        for character in '\r\n,();':
          self.headerlines[idx] = self.headerlines[idx].replace(character,'')

      self.parameters={'CCD_IMAGE_SATURATION':65535}
      for tag,search,idx,datatype in [
          ('CCD_IMAGE_SATURATION','Count_cutoff',1,int),
          ('DETECTOR_SN','Detector:',-1,str),
          ('PIXEL_SIZE','Pixel_size',1,float),
          ('PIXEL_SIZE_UNITS','Pixel_size',2,str),
          ('OSC_START','Start_angle',1,float),
          ('DISTANCE','Detector_distance',1,float),
          ('DISTANCE_UNITS','Detector_distance',2,str),
          ('WAVELENGTH',r'Wavelength',1,float),
          ('BEAM_CENTER_X',r'Beam_xy',1,float),
          ('BEAM_CENTER_Y',r'Beam_xy',2,float),
          ('OSC_RANGE','Angle_increment',1,float),
          ('TWOTHETA','Detector_2theta',1,float),
          ('AXIS','Oscillation_axis',1,str),
          ('PHI','Phi',1,float),
          ('OMEGA','OMEGA',1,float),
          ('DATE','DATE',1,str),
          ]:
          for line in self.headerlines:
            if line.find(search)==0:
              if idx==-1:
                tokens=line.split(" ")
                self.parameters[tag] = " ".join(tokens[1:len(tokens)])
                break
              self.parameters[tag] = datatype(line.split(" ")[idx])
              break
      #unit fixes
      self.parameters['DISTANCE']*={
                  'mm':1,'m':1000}[self.parameters['DISTANCE_UNITS']]
      self.parameters['PIXEL_SIZE']*={
                  'mm':1,'m':1000}[self.parameters['PIXEL_SIZE_UNITS']]
      self.parameters['BEAM_CENTER_X']*=self.parameters['PIXEL_SIZE']
      self.parameters['BEAM_CENTER_Y']*=self.parameters['PIXEL_SIZE']
      # x,y beam center swap; do not know why
      swp = copy.copy(self.parameters['BEAM_CENTER_X'])
      self.parameters['BEAM_CENTER_X']=copy.copy(self.parameters['BEAM_CENTER_Y'])
      self.parameters['BEAM_CENTER_Y']=copy.copy(swp)

      # read array size
      headeropen = rawdata.index("_array_data.data")
      headerclose= rawdata.index("X-Binary-Size-Padding")
      self.header = rawdata[headeropen+1:headerclose]
      self.headerlines = [x.strip() for x in self.header.split("\n")]
      for idx in xrange(len(self.headerlines)):
        for character in '\r\n,();':
          self.headerlines[idx] = self.headerlines[idx].replace(character,'')

      for tag,search,idx,datatype in [
          ('SIZE1','X-Binary-Size-Second-Dimension',1,int),
          ('SIZE2','X-Binary-Size-Fastest-Dimension',1,int),
          ]:
          for line in self.headerlines:
            if line.find(search)==0:
              self.parameters[tag] = datatype(line.split(" ")[idx])
              break
      if self.size1==2527 and self.size2==2463:
        self.vendortype="Pilatus-6M"
      elif self.size1==1679 and self.size2==1475:
        self.vendortype="Pilatus-2M"

if __name__=='__main__':
  import sys
  i = sys.argv[1]
  a = PilatusImage(i)
  a.read()
  print a
  print a.parameters
  print a.rawdata, len(a.rawdata), a.size1*a.size2
