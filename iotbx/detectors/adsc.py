from __future__ import absolute_import, division, print_function
import re
from iotbx.detectors.detectorbase import DetectorImageBase

class ADSCImage(DetectorImageBase):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "ADSC"
    self.open_file = open # default: open files with built-in, unless otherwise instructed by dxtbx format

  def readHeader(self,maxlength=12288, external_keys=None): # usually 1024 is OK; require 12288 for ID19
    if not self.parameters:
      MAGIC_NUMBER = b'{\nHEADER_BYTES='
      stream = self.open_file(self.filename, 'rb')

      # Check the magic number and get the size of the header before
      # the start of the image data
      # (http://strucbio.biologie.uni-konstanz.de/ccp4wiki/index.php/SMV_file_format).
      # This obsoletes the maxlength parameter.
      assert stream.read(len(MAGIC_NUMBER)) == MAGIC_NUMBER
      while True:
        c = stream.read(1)
        if not c.isspace():
          break
      header_bytes = c
      while True:
        c = stream.read(1)
        if not c.isdigit():
          break
        header_bytes += c

      # Reread the entire header.  Only storing the remaining parts of
      # the header in self.header may break the API.
      stream.seek(0)
      self.header = stream.read(int(header_bytes))
      stream.close()

      self.parameters={'CCD_IMAGE_SATURATION':65535}

      library = [
          ('HEADER_BYTES','HEADER_BYTES',int),
          ('SIZE1','SIZE1',int),
          ('SIZE2','SIZE2',int),
          ('CCD_IMAGE_SATURATION','CCD_IMAGE_SATURATION',int),
          ('DETECTOR_SN','DETECTOR_SN',int),
          ('PIXEL_SIZE','PIXEL_SIZE',float),
          ('OSC_START','OSC_START',float),
          ('DISTANCE',r'\nDISTANCE',float),
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
          matches = pattern.findall(self.header.decode("ascii"))
          if len(matches)>0:
            if matches[-1] not in [None,"None","unknown"]:
              self.parameters[tag] = datatype(matches[-1])
      if "TWOTHETA" not in self.parameters:
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
  print(a.linearintdata)
  print(a.linearintdata.size())
  print(a.linearintdata.accessor().grid())
  from labelit.detectors.jpeg import JPEGImage
  j = JPEGImage(a)
  j.calcimage()
  j.write(sys.argv[2])
