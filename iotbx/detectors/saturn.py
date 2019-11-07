from __future__ import absolute_import, division, print_function
import re
from iotbx.detectors.adsc import ADSCImage
from iotbx.detectors.detectorbase import DetectorImageBase

INT = (int,)
FLOAT = (float,)
STR = (str,)

class SaturnImage(ADSCImage):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "RigakuSaturn"

  def getTupleofType(self,inputstr,typefunc):
    parsed = inputstr.split(' ')
    return [typefunc(I) for I in parsed if I != '']

  def readHeader(self,maxlength=6144):
    if not self.parameters:
      with open(self.filename,"rb") as fh:
        rawdata = fh.read(maxlength)
      headeropen = rawdata.index(b"{")
      headerclose = rawdata.index(b"}")
      self.header = rawdata[headeropen+1:headerclose-headeropen].decode("latin-1")

      self.parameters={}
      for tag,search,datatype in [
          ('CCD_IMAGE_SATURATION','SATURATED_VALUE',float),
          ('HEADER_BYTES','HEADER_BYTES',int),
          ('BitmapSize','BitmapSize',int),
          ('SIZE1','SIZE1',int),
          ('SIZE2','SIZE2',int),
          ('CCD_DETECTOR_DESCRIPTION','CCD_DETECTOR_DESCRIPTION',str),
          ('CCD_DETECTOR_DIMENSIONS','CCD_DETECTOR_DIMENSIONS',INT),
          ('CCD_DETECTOR_SIZE','CCD_DETECTOR_SIZE',FLOAT),
          ('CCD_GONIO_DESCRIPTION','CCD_GONIO_DESCRIPTION',str),
          ('CCD_GONIO_NAMES','CCD_GONIO_NAMES',STR),
          ('CCD_GONIO_NUM_VALUES','CCD_GONIO_NUM_VALUES',int),
          ('CCD_GONIO_UNITS','CCD_GONIO_UNITS',str),
          ('CCD_GONIO_VALUES','CCD_GONIO_VALUES',FLOAT),
          ('CCD_GONIO_VALUES_MAX','CCD_GONIO_VALUES_MAX',FLOAT),
          ('CCD_GONIO_VALUES_MIN','CCD_GONIO_VALUES_MIN',FLOAT),
#          ('PIXEL_SIZE','PIXEL_SIZE',float),
#          ('OSC_START','OSC_START',float),
#          ('DISTANCE','DISTANCE',float),
          ('WAVELENGTH','SCAN_WAVELENGTH',float),
          ('CCD_SPATIAL_BEAM_POSITION','CCD_SPATIAL_BEAM_POSITION',FLOAT),
#          ('BEAM_CENTER_X',r'\nBEAM_CENTER_X',float),
#          ('BEAM_CENTER_Y',r'\nBEAM_CENTER_Y',float),
#          ('OSC_RANGE','OSC_RANGE',float),
#          ('TWOTHETA','TWOTHETA',float),
          ('BYTE_ORDER','BYTE_ORDER',str),
          ('AXIS','ROTATION_AXIS_NAME',str),
#          ('PHI','PHI',float),
#          ('OMEGA','OMEGA',float),
          ('DATE','DTREK_DATE_TIME',str),
          ('ROTATION',r'\nROTATION',FLOAT),
          ]:
          matches = re.findall(search+'='+r'(.*);', self.header)
          if len(matches)>0:
            if type(datatype) == type((0,1)):
              self.parameters[tag] = self.getTupleofType(
                matches[-1],datatype[0])
            else:
              self.parameters[tag] = datatype(matches[-1])
      assert self.parameters['CCD_DETECTOR_DESCRIPTION'].find('Saturn')>=0
      assert self.parameters['CCD_DETECTOR_DIMENSIONS'][0]==self.size1
      self.parameters['PIXEL_SIZE'] = self.parameters['CCD_DETECTOR_SIZE'
        ][0] / self.size1
      # rounding to hundreth of degree since encoder reports six (too many) decimal places
      self.parameters['OSC_START'] = round(self.parameters['ROTATION'][0],2)
      assert self.parameters['CCD_GONIO_NAMES'][5]=='Distance'
      self.parameters['DISTANCE'] = self.parameters['CCD_GONIO_VALUES'][5]
      self.parameters['BEAM_CENTER_X'] = self.parameters[
        'CCD_SPATIAL_BEAM_POSITION'][0] * self.pixel_size
      self.parameters['BEAM_CENTER_Y'] = self.parameters[
        'CCD_SPATIAL_BEAM_POSITION'][1] * self.pixel_size
      self.parameters['OSC_RANGE'] = round(self.parameters[
        'ROTATION'][1] - self.parameters['ROTATION'][0],2)
      assert self.parameters['CCD_GONIO_NAMES'][1]=='2Theta'
      self.parameters['TWOTHETA'] = self.parameters['CCD_GONIO_VALUES'][2]

  def read(self):
    from iotbx.detectors import ReadRAXIS
    with open(self.filename,'rb') as fh:
      fh.seek(self.dataoffset())
      chardata = fh.read(self.size1 * self.size2 * self.integerdepth() )
    self.bin_safe_set_data( ReadRAXIS(chardata,self.dataoffset(),
         self.size1*self.bin,self.size2*self.bin,
         self.endian_swap_required())
    )

if __name__=='__main__':
  import sys
  i = sys.argv[1]
  a = SaturnImage(i)
  a.readHeader()
  a.read()
