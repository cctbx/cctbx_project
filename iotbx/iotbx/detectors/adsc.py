import re
from iotbx.detectors import ReadADSC

class ADSCImage:
  def __init__(self,filename):
    self.filename=filename
    self.parameters=None
    self.linearintdata=None

  def readHeader(self,maxlength=1024):
    if not self.parameters:
      rawdata = open(self.filename,"rb").read(maxlength)
      headeropen = rawdata.index("{")
      headerclose= rawdata.index("}")
      header = rawdata[headeropen+1:headerclose-headeropen]

      parameters={}
      for item in ['HEADER_BYTES','SIZE1','SIZE2','CCD_IMAGE_SATURATION']:
        pattern = re.compile(item+'='+r'(.*);')
        matches = pattern.findall(header)
        parameters[item] = int(matches[-1])
      for item in ['PIXEL_SIZE','OSC_START',]:
        pattern = re.compile(item+'='+r'(.*);')
        matches = pattern.findall(header)
        parameters[item] = float(matches[-1])
      self.parameters=parameters

  def fileLength(self):
    self.readHeader()
    self.ptr = self.parameters['HEADER_BYTES']
    self.size1 = self.parameters['SIZE1']
    self.size2 = self.parameters['SIZE2']
    self.file_length = self.ptr+2*self.size1*self.size2
    return self.file_length

  def read(self):
    self.fileLength()
    self.linearintdata = ReadADSC(self.filename,self.ptr,self.size1,self.size2)

  def __getattr__(self, attr):
    if   attr=='size1' : return self.parameters['SIZE1']
    elif attr=='size2' : return self.parameters['SIZE2']
    elif attr=='npixels' : return self.parameters['SIZE1'] * self.parameters['SIZE2']
    elif attr=='saturation' : return self.parameters['CCD_IMAGE_SATURATION']
    elif attr=='rawdata' : return self.linearintdata
    elif attr=='pixel_size' : return self.parameters['PIXEL_SIZE']
    elif attr=='osc_start' : return self.parameters['OSC_START']

if __name__=='__main__':
  i = "./procrun0000035903/run35903_1_001.img"
  i = "/net/boa/scratch1/sauter/lyso1128_4_001.img"
  i = "/net/boa/scratch1/sauter/19-july-02/ProjectNorth/DD5257/52813A12a_1_001.img"
  i = "./procrun0000035905/run35905_1_001.img"
  a = ADSCImage(i)
  a.read()
  print a.linearintdata
  print a.linearintdata.size()
  print a.linearintdata.accessor().grid()
  from iotbx.detectors.jpeg import JPEGImage
  j = JPEGImage(a)
  j.calcimage()
  j.write("/net/cci/sauter/public_html/004.jpg")
