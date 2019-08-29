from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from iotbx.detectors.detectorbase import DetectorImageBase

class GenericDetector(DetectorImageBase):
  """Generic image data, for example from electron crystallography. User
     should supply their own data here, and use the
     phenix.example_viewer command in rstbx/command_line/example_viewer.py
  """
  def __init__(self,filename):
    self.filename = filename
    self.size2 = 200
    self.size1 = 250
    self.pixel_size = 0.1
    #Fictitious vendor type to work around the beam cener convention code
    self.vendortype = "npy_raw"
    self.beamx = 10.0
    self.beamy = 12.5

  def readHeader(self):
    self.distance = 100.
    self.twotheta = 0.0
    self.wavelength = 1.0
    self.saturation = 255
    return

  def show_header(self):
    return "Generic detector with nothing in it"

  def read(self):
    # it is intended that the filename should be used to read in the raw
    #  data; but in this example just use random numbers:
    rawdata = 256*flex.random_double(self.size1*self.size2)
    rawdata.reshape(flex.grid(self.size1,self.size2))
    # this could equally well have been a imported numpy array
    #  import numpy
    #  rawdata2 = 256*numpy.random.rand(self.size1,self.size2)
    #  rawdata = flex.double(rawdata2) #conversion of numpy to cctbx-flex type

    self.linearintdata = rawdata.iround()
