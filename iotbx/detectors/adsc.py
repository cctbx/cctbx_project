import re
from iotbx.detectors.detectorbase import DetectorImageBase

class ADSCImage(DetectorImageBase):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "ADSC"

  def readHeader(self,maxlength=12288, external_keys=None): # usually 1024 is OK; require 12288 for ID19
    if not self.parameters:
      rawdata = open(self.filename,"rb").read(maxlength)
      headeropen = rawdata.index("{")
      headerclose= rawdata.index("}")
      self.header = rawdata[headeropen+1:headerclose-headeropen]

      self.parameters={'CCD_IMAGE_SATURATION':65535}

      library = [
          ('HEADER_BYTES','HEADER_BYTES',int),
          ('SIZE1','SIZE1',int),
          ('SIZE2','SIZE2',int),
          ('CCD_IMAGE_SATURATION','CCD_IMAGE_SATURATION',int),
          ('DETECTOR_SN','DETECTOR_SN',int),
          ('PIXEL_SIZE','PIXEL_SIZE',float),
          ('OSC_START','OSC_START',float),
          ('DISTANCE','DISTANCE',float),
          ('WAVELENGTH',r'\nWAVELENGTH',float),
          ('BEAM_CENTER_X',r'\nBEAM_CENTER_X',float),
          ('BEAM_CENTER_Y',r'\nBEAM_CENTER_Y',float),
          ('OSC_RANGE','OSC_RANGE',float),
          ('TWOTHETA','TWOTHETA',float),
          ('BYTE_ORDER','BYTE_ORDER',str),
          ('AXIS','AXIS',str),
          ('PHI','PHI',float),
          ('OMEGA','OMEGA',float),
          ('DATE','DATE',str),
          ]
      if external_keys is not None:
          assert len(external_keys[0]) == 3
          library = library + external_keys

      for tag,search,datatype in library:
          pattern = re.compile(search+'='+r'(.*);')
          matches = pattern.findall(self.header)
          if len(matches)>0:
            if matches[-1] not in [None,"None","unknown"]:
              self.parameters[tag] = datatype(matches[-1])
      if not self.parameters.has_key("TWOTHETA"):
        self.parameters["TWOTHETA"]=0.0


  def dataoffset(self):
    return self.parameters['HEADER_BYTES']

  def integerdepth(self):
    return 2

  #ADSC Quantum 210, ALS beamline 5.0.2; SUN: unsigned short big endian
  #ADSC Quantum 4R, ALS beamline 5.0.3; WINDOWS: unsigned short little endian
  def getEndian(self):
    if self.parameters['BYTE_ORDER'].lower().find('big')>=0:
      return 1 #big_endian
    else:
      return 0 #little_endian

if __name__=='__main__':
  import sys
  i = sys.argv[1]
  a = ADSCImage(i)
  a.read()
  print a.linearintdata
  print a.linearintdata.size()
  print a.linearintdata.accessor().grid()
  from labelit.detectors.jpeg import JPEGImage
  j = JPEGImage(a)
  j.calcimage()
  j.write(sys.argv[2])
