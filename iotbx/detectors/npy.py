# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8; -*-

#import re
#from iotbx.detectors.adsc         import ADSCImage
from iotbx.detectors.detectorbase import DetectorImageBase
from scitbx.array_family          import flex
import cPickle as pickle

#INT   = (int,)
#FLOAT = (float,)
#STR   = (str,)

class NpyImage(DetectorImageBase):
  def __init__(self, filename):
    DetectorImageBase.__init__(self, filename)
    self.vendortype = "npy_raw"

#  def getTupleofType(self,inputstr,typefunc):
#    parsed = inputstr.split(' ')
#    while '' in parsed:
#      parsed.remove('')
#    return [typefunc(I) for I in parsed]

  def readHeader(self):
    import numpy

    stream      = open(self.filename, "rb")
    cspad_data  = pickle.load(stream)
    stream.close()

    # XXX assert that cspad_data['image'].ndim is 2?

    # From Philipp et al. (2007): pixel size 110 um by 110 um, 14-bit
    # counters.  XXX The beamEnrg thing is still horribly wrong!
    self.parameters                         = {}
    self.parameters['SIZE1']                = cspad_data['image'].shape[0] # XXX order?
    self.parameters['SIZE2']                = cspad_data['image'].shape[1] # XXX order?
    self.parameters['PIXEL_SIZE']           = 110e-3 # XXX fiction
    self.parameters['CCD_IMAGE_SATURATION'] = 2**14 - 1
    self.parameters['SATURATED_VALUE']      = 2**14 - 1
    self.parameters['DISTANCE']             = 93   # XXX fiction
    self.parameters['TWOTHETA']             = 0    # XXX fiction
    self.parameters['OSC_START']            = 0    # XXX fiction
    self.parameters['OSC_RANGE']            = 0    # XXX fiction
    self.parameters['BEAM_CENTER_X']        = 0.5 * self.parameters['SIZE1'] * self.parameters['PIXEL_SIZE']  # XXX order?
    self.parameters['BEAM_CENTER_Y']        = 0.5 * self.parameters['SIZE2'] * self.parameters['PIXEL_SIZE']  # XXX order?
    self.parameters['WAVELENGTH']           = 12398.0 / cspad_data['beamEnrg'] # XXX correct?

    SI = cspad_data['image'].astype(numpy.int32)
    print SI.dtype

#    x        = 626
#    y        = 458
#    x_size   = 185
#    y_size   = 392
#    x_off    = 2
#    y_off    = 3
#
#    section1                           = SI[y:(y + y_size), x:(x + x_size)].cop#y()
#    SI[y:(y + y_size), x:(x + x_size)] = numpy.ones((y_size, x_size), dtype="in#t32")
#    SI[(y + y_off):(y + y_off + y_size),
#       (x + x_off):(x + x_off + x_size)] = section1

    SI = flex.int(SI)

    self.bin_safe_set_data(SI)

#    if not self.parameters:
#      rawdata = open(self.filename,"rb").read(maxlength)
#      headeropen = rawdata.index("{")
#      headerclose= rawdata.index("}")
#      self.header = rawdata[headeropen+1:headerclose-headeropen]
#
#      self.parameters={}
#      for tag,search,datatype in [
#          ('CCD_IMAGE_SATURATION','SATURATED_VALUE',float),
#          ('HEADER_BYTES','HEADER_BYTES',int),
#          ('BitmapSize','BitmapSize',int),
#          ('SIZE1','SIZE1',int),
#          ('SIZE2','SIZE2',int),
#          ('CCD_DETECTOR_DESCRIPTION','CCD_DETECTOR_DESCRIPTION',str),
#          ('CCD_DETECTOR_DIMENSIONS','CCD_DETECTOR_DIMENSIONS',INT),
#          ('CCD_DETECTOR_SIZE','CCD_DETECTOR_SIZE',FLOAT),
#          ('CCD_GONIO_DESCRIPTION','CCD_GONIO_DESCRIPTION',str),
#          ('CCD_GONIO_NAMES','CCD_GONIO_NAMES',STR),
#          ('CCD_GONIO_NUM_VALUES','CCD_GONIO_NUM_VALUES',int),
#          ('CCD_GONIO_UNITS','CCD_GONIO_UNITS',str),
#          ('CCD_GONIO_VALUES','CCD_GONIO_VALUES',FLOAT),
#          ('CCD_GONIO_VALUES_MAX','CCD_GONIO_VALUES_MAX',FLOAT),
#          ('CCD_GONIO_VALUES_MIN','CCD_GONIO_VALUES_MIN',FLOAT),
##          ('PIXEL_SIZE','PIXEL_SIZE',float),
##          ('OSC_START','OSC_START',float),
##          ('DISTANCE','DISTANCE',float),
#          ('WAVELENGTH','SCAN_WAVELENGTH',float),
#          ('CCD_SPATIAL_BEAM_POSITION','CCD_SPATIAL_BEAM_POSITION',FLOAT),
##          ('BEAM_CENTER_X',r'\nBEAM_CENTER_X',float),
##          ('BEAM_CENTER_Y',r'\nBEAM_CENTER_Y',float),
##          ('OSC_RANGE','OSC_RANGE',float),
##          ('TWOTHETA','TWOTHETA',float),
#          ('BYTE_ORDER','BYTE_ORDER',str),
#          ('AXIS','ROTATION_AXIS_NAME',str),
##          ('PHI','PHI',float),
##          ('OMEGA','OMEGA',float),
#          ('DATE','DTREK_DATE_TIME',str),
#          ('ROTATION',r'\nROTATION',FLOAT),
#          ]:
#          pattern = re.compile(search+'='+r'(.*);')
#          matches = pattern.findall(self.header)
#          if len(matches)>0:
#            if type(datatype) == type((0,1)):
#              self.parameters[tag] = self.getTupleofType(
#                matches[-1],datatype[0])
#            else:
#              self.parameters[tag] = datatype(matches[-1])
#      assert self.parameters['CCD_DETECTOR_DESCRIPTION'].find('Saturn')>=0
#      assert self.parameters['CCD_DETECTOR_DIMENSIONS'][0]==self.size1
#      self.parameters['PIXEL_SIZE'] = self.parameters['CCD_DETECTOR_SIZE'
#        ][0] / self.size1
#      # rounding to hundreth of degree since encoder reports six (too many) decimal places
#      self.parameters['OSC_START'] = round(self.parameters['ROTATION'][0],2)
#      assert self.parameters['CCD_GONIO_NAMES'][5]=='Distance'
#      self.parameters['DISTANCE'] = self.parameters['CCD_GONIO_VALUES'][5]
#      self.parameters['BEAM_CENTER_X'] = self.parameters[
#        'CCD_SPATIAL_BEAM_POSITION'][0] * self.pixel_size
#      self.parameters['BEAM_CENTER_Y'] = self.parameters[
#        'CCD_SPATIAL_BEAM_POSITION'][1] * self.pixel_size
#      self.parameters['OSC_RANGE'] = round(self.parameters[
#        'ROTATION'][1] - self.parameters['ROTATION'][0],2)
#      assert self.parameters['CCD_GONIO_NAMES'][1]=='2Theta'
#      self.parameters['TWOTHETA'] = self.parameters['CCD_GONIO_VALUES'][2]

  # This is nop, because all the data has been read by readHeader().
  # The header information and the data are all contained in the same
  # pickled object.
  def read(self):
    pass


  def translate_tiles(self, phil):

    assert 2 * len(phil.distl.tile_translations) == len(phil.distl.detector_tiling)

    SI_old = self.__getattr__('rawdata') # XXX Why are these called SI?
    SI_new = flex.int(flex.grid(SI_old.focus()))

    # XXX The // operator for integer division was introduced when?
    for i in xrange(len(phil.distl.tile_translations) // 2):
      shift_slow = phil.distl.tile_translations[2 * i + 0]
      shift_fast = phil.distl.tile_translations[2 * i + 1]

      ur_slow = phil.distl.detector_tiling[4 * i + 0]
      ur_fast = phil.distl.detector_tiling[4 * i + 1]
      ll_slow = phil.distl.detector_tiling[4 * i + 2]
      ll_fast = phil.distl.detector_tiling[4 * i + 3]

      #print "Shifting tile at (%d, %d) by (%d, %d)" % (ur_slow, ur_fast, shift_slow, shift_fast)

      SI_new.matrix_paste_block_in_place(
        block = SI_old.matrix_copy_block(
          i_row=ur_slow,i_column=ur_fast,
          n_rows=ll_slow-ur_slow, n_columns=ll_fast-ur_fast),
        i_row = ur_slow + shift_slow,
        i_column = ur_fast + shift_fast
      )

    self.bin_safe_set_data(SI_new)

#if __name__=='__main__':
#  import sys
#  i = sys.argv[1]
#  a = SaturnImage(i)
#  a.readHeader()
#  a.read()
