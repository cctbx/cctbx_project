from __future__ import absolute_import, division, print_function
import re
from iotbx.detectors.adsc import ADSCImage
from iotbx.detectors.detectorbase import DetectorImageBase

INT = (int,)
FLOAT = (float,)
STR = (str,)

class NoirImage(ADSCImage):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "RAXIS"

  def getTupleofType(self,inputstr,typefunc):
    parsed = inputstr.split(' ')
    return [typefunc(I) for I in parsed if I != '']

  def readHeader(self,maxlength=6144):
    if not self.parameters:
      with open(self.filename,"rb") as fh:
        rawdata = fh.read(maxlength)
      headeropen = rawdata.index(b"{")
      headerclose= rawdata.index(b"}")
      self.header = rawdata[headeropen+1:headerclose-headeropen].decode("latin-1")

      self.parameters={}
      for tag,search,datatype in [
          ('SATURATED_VALUE','SATURATED_VALUE',float),
          ('HEADER_BYTES','HEADER_BYTES',int),
          #('BitmapSize','BitmapSize',int),
          ('SIZE1','SIZE1',int),
          ('SIZE2','SIZE2',int),
          ('NOIR1_DETECTOR_DESCRIPTION','NOIR1_DETECTOR_DESCRIPTION',str),
          ('NOIR1_DETECTOR_DIMENSIONS','NOIR1_DETECTOR_DIMENSIONS',INT),
          ('NOIR1_DETECTOR_SIZE','NOIR1_DETECTOR_SIZE',FLOAT),
          ('NOIR1_GONIO_DESCRIPTION','NOIR1_GONIO_DESCRIPTION',str),
          ('NOIR1_GONIO_NAMES','NOIR1_GONIO_NAMES',STR),
          ('NOIR1_GONIO_NUM_VALUES','NOIR1_GONIO_NUM_VALUES',int),
          ('NOIR1_GONIO_UNITS','NOIR1_GONIO_UNITS',str),
          ('NOIR1_GONIO_VALUES','NOIR1_GONIO_VALUES',FLOAT),
          #('NOIR1_GONIO_VALUES_MAX','NOIR1_GONIO_VALUES_MAX',FLOAT),
          #('NOIR1_GONIO_VALUES_MIN','NOIR1_GONIO_VALUES_MIN',FLOAT),
#          ('PIXEL_SIZE','PIXEL_SIZE',float),
#          ('OSC_START','OSC_START',float),
          ('DISTANCE','NOIR1_ACTUAL_DISTANCE',float),
          ('WAVELENGTH','NOIR1_ACTUAL_ENERGY',float),
          ('NOIR1_SPATIAL_BEAM_POSITION','NOIR1_SPATIAL_BEAM_POSITION',FLOAT),
#          ('BEAM_CENTER_X',r'\nBEAM_CENTER_X',float),
#          ('BEAM_CENTER_Y',r'\nBEAM_CENTER_Y',float),
#          ('OSC_RANGE','OSC_RANGE',float),
          ('TWOTHETA','NOIR1_ACTUAL_THETA',float),
          ('BYTE_ORDER','BYTE_ORDER',str),
          ('AXIS','ROTATION_AXIS_NAME',str),
#          ('PHI','PHI',float),
#          ('OMEGA','OMEGA',float),
          #('DATE','DTREK_DATE_TIME',str),
          ('ROTATION',r'\nROTATION',FLOAT),
          ]:
          matches = re.findall(search+'='+r'(.*);', self.header)
          if len(matches)>0:
            if type(datatype) == type((0,1)):
              self.parameters[tag] = self.getTupleofType(
                matches[-1],datatype[0])
            else:
              self.parameters[tag] = datatype(matches[-1])
      assert self.parameters['NOIR1_DETECTOR_DESCRIPTION'].find('NOIR')>=0
      assert self.parameters['NOIR1_DETECTOR_DIMENSIONS'][0]==self.size1
      self.parameters['PIXEL_SIZE'] = self.parameters['NOIR1_DETECTOR_SIZE'][0] / self.size1
      # rounding to hundreth of degree since encoder reports six (too many) decimal places
      self.parameters['OSC_START'] = round(self.parameters['ROTATION'][0],2)
      #assert self.parameters['NOIR1_GONIO_NAMES'][5]=='Distance'
      #self.parameters['DISTANCE'] = self.parameters['NOIR1_GONIO_VALUES'][5]
      self.parameters['BEAM_CENTER_X'] = self.parameters[
        'NOIR1_SPATIAL_BEAM_POSITION'][0] * self.pixel_size
      self.parameters['BEAM_CENTER_Y'] = self.parameters[
        'NOIR1_SPATIAL_BEAM_POSITION'][1] * self.pixel_size
      self.parameters['OSC_RANGE'] = round(self.parameters[
        'ROTATION'][1] - self.parameters['ROTATION'][0],2)
      #assert self.parameters['NOIR1_GONIO_NAMES'][1]=='2Theta'
      #self.parameters['TWOTHETA'] = self.parameters['NOIR1_GONIO_VALUES'][2]
      hc = 12398
      self.parameters['WAVELENGTH'] = hc / self.parameters['WAVELENGTH']

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
  a = NoirImage(i)
  a.readHeader()
  a.read()
