import re,struct
from iotbx.detectors.detectorbase import DetectorImageBase

class MARImage(DetectorImageBase):
  def __init__(self,filename):
    DetectorImageBase.__init__(self,filename)
    self.vendortype = "MARCCD"

    byte_order = str(open(self.filename,"rb").read(2))
    if byte_order == 'II':
      self.endian = 0
    else:
      self.endian = 1

    assert not self.isCompressed()

  def isCompressed(self):
    if self.getEndian(): format = '>'
    else: format = '<'

    f = open(self.filename,"rb")
    try:
      f.seek(4)
      rawdata = f.read(4)
      ifd = struct.unpack(format+'i',rawdata)[0]

      # if ifd is 0 then there are no more Image File Directories
      while ifd:
        f.seek(ifd)

        rawdata = f.read(2)
        # get the number of directory entries in the IFD
        numentries = struct.unpack(format+'h',rawdata)[0]

        # search for compression tag
        for x in range(numentries):
          f.seek(ifd+x*12+2)
          rawdata = f.read(2)
          tag = struct.unpack(format+'h',rawdata)[0]
          if tag == 259: # value of compression tag
            f.seek(ifd+x*12+2+8) # seek to value offset
            # value is left justified in 4 byte field
            # read two bytes so can unpack as short
            rawdata = f.read(2)
            value = struct.unpack(format+'h',rawdata)[0]
            #print value
            if value == 1: # no compression
              return 0
            else:
              return 1

        f.seek(ifd+numentries*12+2)
        rawdata = f.read(4)
        ifd = struct.unpack(format+'i',rawdata)[0]

    finally:
      f.close()

    # control should never reach this point
    assert 1==0

  # returns 0 for little endian 'II'
  # returns 1 for big endian 'MM'
  def getEndian(self):
    return self.endian

  def readHeader(self,offset=1024):
    if not self.parameters:
      if self.getEndian(): format = '>'
      else: format = '<'

      f = open(self.filename,"rb")
      try:
        f.seek(2464) # seek to file_comments
        file_comments = f.read(512) # read file_comments

        parameters={}
        for item in [('Detector Serial Number','DETECTOR_SN')]: #expected integers
          pattern = re.compile(item[0]+' = '+r'(.*)')
          matches = pattern.findall(file_comments)
          if len(matches) > 0:
            parameters[item[1]] = int(matches[-1])
          else:
            parameters[item[1]] = 0

        f.seek(offset+28)
        rawdata = f.read(8)
        header_byte_order,data_byte_order = struct.unpack(format+'ii',rawdata)

        f.seek(offset+80)
        rawdata = f.read(8)
        parameters['SIZE1'],parameters['SIZE2'] = struct.unpack(format+'ii',rawdata)
        assert parameters['SIZE1'] == parameters['SIZE2']
        f.seek(offset+88)
        rawdata = f.read(4)
        self.depth = struct.unpack(format+'i',rawdata)[0]

        f.seek(offset+104)
        rawdata = f.read(4)
        parameters['CCD_IMAGE_SATURATION'] = struct.unpack(format+'i',rawdata)[0]

        f.seek(offset+116)
        rawdata = f.read(12)
        origin,orientation,view_direction = struct.unpack(format+'iii',rawdata)

        f.seek(offset+696)
        rawdata = f.read(4)
        start_xtal_to_detector = struct.unpack(format+'i',rawdata)[0]/1000.
        f.seek(offset+728)
        rawdata = f.read(4)
        end_xtal_to_detector = struct.unpack(format+'i',rawdata)[0]/1000.
        #assert start_xtal_to_detector == end_xtal_to_detector
        #that assertion would've been nice but ESRF BM14 frames fail; instead:
        assert start_xtal_to_detector>0.
        parameters['DISTANCE'] = start_xtal_to_detector

        f.seek(offset+772)
        rawdata = f.read(8)
        pixelsize_x,pixelsize_y = struct.unpack(format+'ii',rawdata)
        assert pixelsize_x == pixelsize_y
        parameters['PIXEL_SIZE'] = pixelsize_x*1.0e-6 # convert from nano to milli


        f.seek(offset+644)
        rawdata = f.read(8)
        beam_center_x,beam_center_y = struct.unpack(format+'ii',rawdata)
        parameters['BEAM_CENTER_X'] = beam_center_x/1000.*parameters['PIXEL_SIZE']
        parameters['BEAM_CENTER_Y'] = beam_center_y/1000.*parameters['PIXEL_SIZE']

        # ----- phi analysis
        f.seek(offset+684)
        rawdata = f.read(4)
        parameters['OSC_START'] = struct.unpack(format+'i',rawdata)[0]/1000.

        f.seek(offset+716)
        rawdata = f.read(4)
        end_phi = struct.unpack(format+'i',rawdata)[0]/1000.

        #parameters['OSC_RANGE'] = end_phi - parameters['OSC_START']
        #would have thought this would work; but turns out unreliable because
        # software doesn't always fill in the end_phi

        # ----- rotation analysis
        f.seek(offset+736)
        rawdata = f.read(4)
        rotation_range = struct.unpack(format+'i',rawdata)[0]/1000.
        parameters['OSC_RANGE'] = rotation_range

        f.seek(offset+732)
        rawdata = f.read(4)
        rotation_axis = struct.unpack(format+'i',rawdata)[0]
        #assert rotation_axis == 4 # if it isn't phi; go back and recode to cover all cases

        # ----- omega analysis
        f.seek(offset+672)
        rawdata = f.read(4)
        parameters['OMEGA_START'] = struct.unpack(format+'i',rawdata)[0]/1000.

        f.seek(offset+704)
        rawdata = f.read(4)
        parameters['OMEGA_END'] = struct.unpack(format+'i',rawdata)[0]/1000.

        if rotation_axis == 4: # rotation axis is phi
          pass
        elif rotation_axis == 1: # rotation about omega
          parameters['OSC_START'] = parameters['OMEGA_START']

        f.seek(offset+668)
        rawdata = f.read(4)
        start_twotheta = struct.unpack(format+'i',rawdata)[0]/1000.
        f.seek(offset+700)
        rawdata = f.read(4)
        end_twotheta = struct.unpack(format+'i',rawdata)[0]/1000.
        assert start_twotheta == end_twotheta
        parameters['TWOTHETA'] = start_twotheta

        f.seek(offset+908)
        rawdata = f.read(4)
        parameters['WAVELENGTH'] = struct.unpack(format+'i',rawdata)[0]*1.0e-5 # convert from femto to angstrom

      finally:
        f.close()

      self.parameters=parameters

  def dataoffset(self):
    return 4096

  def integerdepth(self):
    return self.depth

if __name__=='__main__':
  i = "/net/racer/scratch1/ttleese/lyso201.0002"
  #i = "/net/racer/scratch1/ttleese/oxford.tif"
  m = MARImage(i)
  print m.isCompressed()
  #m.read()
  #print 'endian:',m.getEndian()
  #print 'serial number:',m.serial_number
  #print 'size 1:',m.size1
  #print 'size 2:',m.size2
  #print 'npixels:',m.npixels
  #print 'saturation:',m.saturation
  #print 'beamx:',m.beamx
  #print 'beamy:',m.beamy
  #print 'pixel size:',m.pixel_size
  #print 'osc start:',m.osc_start
  #print 'delta phi:',m.deltaphi
  #print 'two theta:',m.twotheta
  #print 'wav:',m.wavelength
  #print 'distance:',m.distance
  #print 'file length:',m.fileLength()
