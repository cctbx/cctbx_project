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
      for item in ['HEADER_BYTES','SIZE1','SIZE2','CCD_IMAGE_SATURATION',
                   'DETECTOR_SN']:
        pattern = re.compile(item+'='+r'(.*);')
        matches = pattern.findall(header)
        parameters[item] = int(matches[-1])
      for item in ['PIXEL_SIZE','OSC_START','DISTANCE','WAVELENGTH',
                   'BEAM_CENTER_X','BEAM_CENTER_Y','OSC_RANGE',
                   'TWOTHETA']:
        pattern = re.compile(item+'='+r'(.*);')
        matches = pattern.findall(header)
        parameters[item] = float(matches[-1])
      for item in ['BYTE_ORDER']:
        pattern = re.compile(item+'='+r'(.*);')
        matches = pattern.findall(header)
        parameters[item] = matches[-1]
      self.parameters=parameters

  def fileLength(self):
    self.readHeader()
    self.ptr = self.parameters['HEADER_BYTES']
    self.size1 = self.parameters['SIZE1']
    self.size2 = self.parameters['SIZE2']
    self.file_length = self.ptr+2*self.size1*self.size2
    return self.file_length
    # pure supposition:
    #  size1 corresponds to number of rows.  Columns are slow. 
    #  size2 corresponds to number of columns.  Rows are fast.

  def read(self):
    self.fileLength()
    #ADSC Quantum 210, ALS beamline 5.0.2; SUN: unsigned short little endian 
    #ADSC Quantum 4R, ALS beamline 5.0.3; WINDOWS: unsigned short big endian
    if self.parameters['BYTE_ORDER'].lower().find('big')>=0:
      self.linearintdata = ReadADSC(self.filename,self.ptr,
                                    self.size1,self.size2,1) #big_endian
    else:
      self.linearintdata = ReadADSC(self.filename,self.ptr,
                                    self.size1,self.size2,0) #little_endian
      

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
