# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8; -*-
#
# $Id$

from iotbx.detectors.detectorbase import DetectorImageBase
from scitbx.array_family          import flex
import cPickle as pickle

#INT   = (int,)
#FLOAT = (float,)
#STR   = (str,)

class NpyImage(DetectorImageBase):
  def __init__(self, filename, source_data = None):
    DetectorImageBase.__init__(self, filename)
    self.vendortype = "npy_raw"
    self.source_data = source_data

#  def getTupleofType(self,inputstr,typefunc):
#    parsed = inputstr.split(' ')
#    while '' in parsed:
#      parsed.remove('')
#    return [typefunc(I) for I in parsed]

  def readHeader(self, version_control = "CXI 3.2"):
    import numpy

    if self.source_data == None:
      stream      = open(self.filename, "rb")
      cspad_data  = pickle.load(stream)
      stream.close()
    else:
      cspad_data  = self.source_data

    # XXX assert that cspad_data['image'].ndim is 2?

    self.parameters                         = {}

    if version_control == "CXI 3.1":
      self.parameters['SIZE1']                = cspad_data['image'].shape[0] # XXX order?
      self.parameters['SIZE2']                = cspad_data['image'].shape[1] # XXX order?
      self.parameters['PIXEL_SIZE']           = 110e-3 # XXX fiction
      self.parameters['BEAM_CENTER_X']        = 0.5 * self.parameters['SIZE1'] * self.parameters['PIXEL_SIZE']  # XXX order?
      self.parameters['BEAM_CENTER_Y']        = 0.5 * self.parameters['SIZE2'] * self.parameters['PIXEL_SIZE']  # XXX order?
      self.parameters['CCD_IMAGE_SATURATION'] = 2**14 - 1
      self.parameters['DISTANCE']             = 93   # XXX fiction
      self.parameters['OSC_START']            = 0    # XXX fiction
      self.parameters['OSC_RANGE']            = 0    # XXX fiction
      self.parameters['SATURATED_VALUE']      = 2**14 - 1
      self.parameters['TWOTHETA']             = 0    # XXX fiction
      # From Margaritondo & Rebernik Ribic (2011): the dimensionless
      # relativistic gamma-factor is derived from beam energy in MeV and
      # the electron rest mass, K is a dimensionless "undulator
      # parameter", and L is the macroscopic undulator period in
      # Aangstroem (XXX).  See also
      # http://ast.coe.berkeley.edu/srms/2007/Lec10.pdf.  XXX This
      # should really move into the pyana code, since the parameters are
      # SLAC-specific.
      gamma                         = cspad_data['beamEnrg'] / 0.510998910
      K                             = 3.5
      L                             = 3.0e8
      self.parameters['WAVELENGTH'] = L / (2 * gamma**2) * (1 + K**2 / 2)
      SI = cspad_data['image'].astype(numpy.int32)
      SI = flex.int(SI)
      self.bin_safe_set_data(SI)
    elif version_control == "CXI 3.2":
      self.parameters['BEAM_CENTER_X']        = cspad_data['BEAM_CENTER_X']
      self.parameters['BEAM_CENTER_Y']        = cspad_data['BEAM_CENTER_Y']
      self.parameters['CCD_IMAGE_SATURATION'] = cspad_data['CCD_IMAGE_SATURATION']
      self.parameters['DISTANCE']             = cspad_data['DISTANCE']
      self.parameters['OSC_RANGE']            = 0 # XXX fiction
      self.parameters['OSC_START']            = 0 # XXX fiction
      self.parameters['PIXEL_SIZE']           = cspad_data['PIXEL_SIZE']
      self.parameters['SATURATED_VALUE']      = cspad_data['SATURATED_VALUE']
      self.parameters['SIZE1']                = cspad_data['SIZE1']
      self.parameters['SIZE2']                = cspad_data['SIZE2']
      self.parameters['TWOTHETA']             = 0 # XXX fiction
      self.parameters['WAVELENGTH']           = cspad_data['WAVELENGTH']
      self.bin_safe_set_data(cspad_data['DATA'])

  # This is nop, because all the data has been read by readHeader().
  # The header information and the data are all contained in the same
  # pickled object.
  def read(self):
    pass


  def translate_tiles(self, phil):
    if phil.distl.tile_translations==None: return
    assert 2 * len(phil.distl.tile_translations) == len(phil.distl.detector_tiling)

    SI_old = self.__getattr__('rawdata') # XXX Why are these called SI?
    SI_new = flex.int(flex.grid(SI_old.focus()))

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
