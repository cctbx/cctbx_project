from __future__ import absolute_import, division, print_function
import copy
import re
from iotbx.detectors.detectorbase import DetectorImageBase
from iotbx.detectors import ImageException

try:
  import bz2
except ImportError:
  bz2 = None

try:
  import gzip
except ImportError:
  gzip = None

class PilatusImage(DetectorImageBase):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "Pilatus"
    self.vendor_specific_null_value = -1

  mandatory_keys = ['PIXEL_SIZE_UNITS', 'DISTANCE', 'PHI', 'WAVELENGTH', 'SIZE1',
    'SIZE2', 'TWOTHETA', 'DISTANCE_UNITS', 'OSC_RANGE',
    'BEAM_CENTER_X', 'BEAM_CENTER_Y',
    'CCD_IMAGE_SATURATION', 'OSC_START', 'DETECTOR_SN', 'PIXEL_SIZE',
    'AXIS']

  def fileLength(self):
    raise ImageException("file length not computed for miniCBF")

  def getEndian(self):
    raise ImageException("endian-ness not computed for miniCBF")

  @staticmethod
  def is_bz2(filename):
    '''Check if a file pointed at by filename is bzip2 format.'''

    if not filename.endswith('.bz2'):
      return False

    with open(filename, 'rb') as fh:
      return b'BZh' == fh.read(3)

  @staticmethod
  def is_gzip(filename):
    '''Check if a file pointed at by filename is gzip compressed.'''

    if not filename.endswith('.gz'):
      return False

    with open(filename, 'rb') as fh:
      return fh.read(2) == b"\x1f\x8b"

  @staticmethod
  def open_file(filename, mode='rb'):
    '''Open file for reading, decompressing silently if necessary,
       caching transparently if possible.'''

    if PilatusImage.is_bz2(filename):
      if bz2 is None:
        raise RuntimeError('bz2 file provided without bz2 module')
      fh_func = lambda: bz2.BZ2File(filename, mode)

    elif PilatusImage.is_gzip(filename):
      if gzip is None:
        raise RuntimeError('gz file provided without gzip module')
      fh_func = lambda: gzip.GzipFile(filename, mode)

    else:
      fh_func = lambda: open(filename, mode)

    return fh_func()

  def endian_swap_required(self):
    return False

  def read(self,algorithm="buffer_based"):
    self.readHeader()
    if self.linearintdata != None and\
      self.linearintdata.size()==self.size1*self.size2:
      #data has already been read
      return
    if self.bin==2:
      raise ImageException("2-by-2 binning not supported for miniCBF")
    try:
      from cbflib_adaptbx import cbf_binary_adaptor # optional package
      self.adaptor = cbf_binary_adaptor(self.filename)

      # assert algorithm in ["cbflib","cbflib_optimized","buffer_based"]

      data = self.adaptor.uncompress_implementation( algorithm
             ).uncompress_data(self.size1,self.size2)
      self.bin_safe_set_data( data )

    except Exception as e:
      raise ImageException(
          "unable to read miniCBF data; contact authors; error=\"%s\"" % \
          str(e).strip())

  def readHeader(self,maxlength=12288): # usually 1024 is OK; require 12288 for ID19
    if not self.parameters:
      with self.open_file(self.filename,"rb") as fh:
        rawdata = fh.read(maxlength)

      # The tag _array_data.header_convention "SLS_1.0" could be with/without quotes "..."
      # SLS_match = re.findall(b'_array_data.header_convention[ "]*SLS', rawdata)
      # PILATUS_match = re.findall(b'_array_data.header_convention[ "]*PILATUS', rawdata)
      #assert len(SLS_match) + len(PILATUS_match)>=1

      # read SLS header
      headeropen = rawdata.index(b"_array_data.header_contents")
      headerclose = rawdata.index(b"_array_data.data")
      self.header = rawdata[headeropen+1:headerclose].decode("latin-1")
      self.headerlines = [x.strip() for x in self.header.split("#")]
      character_filter = re.compile(r"[\r\n,\(\);]")
      self.headerlines = [character_filter.sub("", x) for x in self.headerlines]

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
      header_lines = []
      found_array_data_data = False
      for record in rawdata.decode("latin-1").splitlines():
        if "_array_data.data" in record:
          found_array_data_data = True
        elif not found_array_data_data:
          continue
        elif len(record.strip()) == 0:
          # http://sourceforge.net/apps/trac/cbflib/wiki/ARRAY_DATA%20Category
          #    In an imgCIF file, the encoded binary data begins after
          #    the empty line terminating the header.
          break
        header_lines.append(record)
      self.header = "\n".join(header_lines)
      self.headerlines = [x.strip() for x in self.header.split("\n")]
      self.headerlines = [character_filter.sub("", x) for x in self.headerlines]

      for tag,search,idx,datatype in [
          ('SIZE1','X-Binary-Size-Second-Dimension',-1,int),
          ('SIZE2','X-Binary-Size-Fastest-Dimension',-1,int),
          ]:
          for line in self.headerlines:
            if line.find(search)==0:
              self.parameters[tag] = datatype(line.split(" ")[idx])
              break
      if self.size1==2527 and self.size2==2463:
        self.vendortype="Pilatus-6M"
      elif self.size1==1679 and self.size2==1475:
        self.vendortype="Pilatus-2M"
      elif self.size1==619 and self.size2==487:
        self.vendortype="Pilatus-300K"


if __name__=='__main__':
  import sys
  i = sys.argv[1]
  a = PilatusImage(i)
  a.read()
  print(a)
  print(a.parameters)
  print(a.rawdata, len(a.rawdata), a.size1*a.size2)
