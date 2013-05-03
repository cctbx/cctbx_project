from __future__ import division
from iotbx.detectors.detectorbase import DetectorImageBase

class EIGERImage(DetectorImageBase):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "EIGER"
    self.supports_multiple_images = True
    self.img_number = 0 # 0-indexed to clients, internally 1-indexed

  def readHeader(self,maxlength=12288): # XXX change maxlength!!!
    if not self.parameters:
      import h5py

      masterFilename = self.filename
      self.hdf5File = h5py.File(masterFilename, 'r')
      params = self.hdf5File['entry']['instrument']['detector']
      dparams = self.hdf5File['entry']['instrument']['detector']['detectorSpecific']
      mparams = self.hdf5File['entry']['instrument']['monochromator']
      sparams = self.hdf5File['entry']['sample']

      self.LUT = {} #look up table
      entry = self.hdf5File['entry']

      for datalink in list(entry):
          if not(datalink[0:4] == 'data'):
              continue

          ### open the link ###
          try:
              data = entry[datalink]
          except KeyError, exception: ### cannot open link, probably file does not exist
              continue


          ### read the image_nr_low and image_nr_high attributes ###
          image_nr_low  = data.attrs['image_nr_low']
          image_nr_high = data.attrs['image_nr_high']

          for imgNr in range(image_nr_low, image_nr_high+1):
              self.LUT[imgNr] = (datalink, imgNr-image_nr_low)


      self.parameters={'CCD_IMAGE_SATURATION':65535}

      # XXX paramters listed as ? still need to be matched to metadata in the HDF5 file
      library = [
          ('HEADER_BYTES','HEADER_BYTES',int),                   # ?
          ('SIZE1','y_pixels_in_detector',int),
          ('SIZE2','x_pixels_in_detector',int),
          ('CCD_IMAGE_SATURATION','CCD_IMAGE_SATURATION',int),   # ?
          ('DETECTOR_SN','detector_number',str),
          ('PIXEL_SIZE','x_pixel_size',float),                   # should account for x or y size
          ('OSC_START','OSC_START',float),                       # ?
          ('DISTANCE','detector_distance',float),
          ('WAVELENGTH','wavelength',float),
          ('BEAM_CENTER_X','beam_center_x',float),
          ('BEAM_CENTER_Y','beam_center_y',float),
          ('OSC_RANGE','rotation_angle_step',float),             # note, we return only the first of a table of values here, one for each image
          ('TWOTHETA','TWOTHETA',float),                         # N/A, uses it's own geometery
          ('BYTE_ORDER','BYTE_ORDER',str),                       # ?
          ('AXIS','AXIS',str),                                   # ?
          ('PHI','PHI',float),                                   # ?
          ('OMEGA','OMEGA',float),                               # ?
          ('DATE','data_collection_date',str),
          ]

      for tag,search,datatype in library:
        try:
            self.parameters[tag] = datatype(params[search][0])
        except KeyError:
          try:
            self.parameters[tag] = datatype(dparams[search][0])
          except KeyError:
            try:
              self.parameters[tag] = datatype(mparams[search][0])
            except KeyError:
              try:
                self.parameters[tag] = datatype(sparams[search][0])
              except KeyError: pass

      #these parameters have to be set here due to the call to image_coords_as_detector_cords below.  they are normally set later.
      if self.parameters.has_key("SIZE2"):
        self.image_size_fast = self.size2
      if self.parameters.has_key("SIZE1"):
        self.image_size_slow = self.size1
      if self.parameters.has_key("PIXEL_SIZE"):
        self.pixel_resolution = self.pixel_size

      if not self.parameters.has_key("TWOTHETA"):
        self.parameters["TWOTHETA"]=0.0
      if not self.parameters.has_key("OSC_START"):
        self.parameters["OSC_START"]=0.0
      if self.parameters.has_key("BEAM_CENTER_X") and self.parameters.has_key("BEAM_CENTER_Y"):
        self.parameters["BEAM_CENTER_X"], self.parameters["BEAM_CENTER_Y"] = self.image_coords_as_detector_coords(
          self.parameters["BEAM_CENTER_X"], self.parameters["BEAM_CENTER_Y"])

  def read(self):
      datalink = ''
      try:
          (datalink,imageNrOffset) = self.LUT[self.img_number+1]
      except KeyError, e:
          raise ImageReadException('imgNr out of range')
      data = self.hdf5File['entry'][datalink]
      # use slicing access to get images with image number self.img_number
      image = data[imageNrOffset, : , : ] ## z / y / x

      if len(image) <=0 or len(image[0]) <=0:
        raise ImageReadException('an image dimension is zero')

      from scitbx.array_family import flex
      import numpy
      tmp = image.reshape(len(image), len(image[0]))
      tmp = image.astype(numpy.int32)
      tmp = flex.int(tmp)
      self.bin_safe_set_data(tmp)

  def image_count(self):
    return len(self.LUT)

  def integerdepth(self):
    return 2

  def dataoffset(self):
    return 0

  def get_data_link(self,index=None):
    if index == None:
      return self.LUT[self.img_number+1][0]
    else:
      return self.LUT[index+1][0]

if __name__=="__main__":
  import sys
  E = EIGERImage(sys.argv[1])
  E.readHeader()
  E.show_header()
