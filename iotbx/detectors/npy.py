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

  def readHeader(self, horizons_phil):
    import numpy

    version_control = horizons_phil.distl.detector_format_version

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
    elif version_control in ["CXI 3.2","CXI 4.1"]:
      self.parameters['ACTIVE_AREAS']         = cspad_data.get('ACTIVE_AREAS', None)
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

      if (self.parameters['ACTIVE_AREAS'] != None):
        horizons_phil.distl.detector_tiling = self.parameters['ACTIVE_AREAS']

    if version_control == "CXI 4.1":
      if horizons_phil.distl.tile_translations==None and \
         horizons_phil.distl.detector_tiling is not None:
          horizons_phil.distl.tile_translations = [0]*(len(horizons_phil.distl.detector_tiling)/2)


  # This is nop, because all the data has been read by readHeader().
  # The header information and the data are all contained in the same
  # pickled object.
  def read(self):
    pass

  def translate_tiles(self, phil):
    if phil.distl.detector_tiling==None: return
    if phil.distl.tile_translations==None: return

    if len(phil.distl.detector_tiling) <= 16:
      # assume this is the 2x2 CS Pad for spectroscopy; do not use tile translations
      if phil.distl.detector_format_version in ["CXI 4.1"]:
        # For the Run 4 CXI detector, the first sensor is inactive and pegged high(16K).
        # For calculating display contrast it is better to eliminate the sensor.
        if self.size1 == 370: #there are two sensors; we should eliminate the first
          self.parameters['SIZE1'] = 185
          self.linearintdata = self.linearintdata[int(len(self.linearintdata)/2):]
          self.linearintdata.reshape(flex.grid(self.size1,self.size2))
        print "CXI 2x2 size",self.size1,self.size2, self.linearintdata.focus()
      return

    assert 2 * len(phil.distl.tile_translations) == len(phil.distl.detector_tiling)

    shifted_int_data_old = self.__getattr__('rawdata')
    shifted_int_data_new = flex.int(flex.grid(shifted_int_data_old.focus()))

    manager = self.get_tile_manager(phil)

    for i,shift in enumerate(manager.effective_translations()):
      shift_slow = shift[0]
      shift_fast = shift[1]

      ur_slow = phil.distl.detector_tiling[4 * i + 0]
      ur_fast = phil.distl.detector_tiling[4 * i + 1]
      ll_slow = phil.distl.detector_tiling[4 * i + 2]
      ll_fast = phil.distl.detector_tiling[4 * i + 3]

      #print "Shifting tile at (%d, %d) by (%d, %d)" % (ur_slow, ur_fast, shift_slow, shift_fast)

      shifted_int_data_new.matrix_paste_block_in_place(
        block = shifted_int_data_old.matrix_copy_block(
          i_row=ur_slow,i_column=ur_fast,
          n_rows=ll_slow-ur_slow, n_columns=ll_fast-ur_fast),
        i_row = ur_slow + shift_slow,
        i_column = ur_fast + shift_fast
      )

    self.bin_safe_set_data(shifted_int_data_new)

  def get_tile_manager(self, phil):
    return tile_manager(phil,beam=(int(self.beamx/self.pixel_size),
                                   int(self.beamy/self.pixel_size)))

class tile_manager:
  def __init__(self,working_params,beam=None):
    self.working_params = working_params
    self.beam = beam # direct beam position supplied as slow,fast pixels

  def effective_translations(self):

    # if there are quadrant translations, do some extra work to apply them
    if self.working_params.distl.quad_translations != None:
      from scitbx.matrix import col
      beam = col(self.beam)
      for itile in xrange(len(self.working_params.distl.detector_tiling) // 4):
        tile_center = (
          col(self.working_params.distl.detector_tiling[4*itile:4*itile+2]) +
          col(self.working_params.distl.detector_tiling[4*itile+2:4*itile+4]))/2
        delta = tile_center-beam
        iquad = [(True,True),(True,False),(False,True),(False,False)
                ].index((delta[0]<0, delta[1]<0)) # UL,UR,LL,LR
        yield (self.working_params.distl.tile_translations[2 * itile + 0] +
               self.working_params.distl.quad_translations[2 * iquad + 0],
               self.working_params.distl.tile_translations[2 * itile + 1] +
               self.working_params.distl.quad_translations[2 * iquad + 1])
      return

    for i in xrange(len(self.working_params.distl.tile_translations) // 2):
       yield (self.working_params.distl.tile_translations[2 * i + 0],
              self.working_params.distl.tile_translations[2 * i + 1])

  def effective_tiling_as_flex_int(self,reapply_peripheral_margin=False,**kwargs):
    import copy
    IT = flex.int(copy.copy(self.working_params.distl.detector_tiling))

    assert len(IT)%4==0 # only meaningful for groups of 4
    for itl in xrange(0,len(IT),4): # validate upper-left/ lower-right ordering
      assert IT[itl] < IT[itl+2]; assert IT[itl+1] < IT[itl+3]

    if self.working_params.distl.tile_translations!=None and \
      2*len(self.working_params.distl.tile_translations) == len(IT):
      #assume that the tile translations have already been applied at the time
      #the file is read; now they need to be applied to the spotfinder tile boundaries

      #check if beam position has been supplied
      self.beam = kwargs.get("beam",self.beam)

      for i,shift in enumerate(self.effective_translations()):
        shift_slow = shift[0]
        shift_fast = shift[1]
        IT[4 * i + 0] += shift_slow
        IT[4 * i + 1] += shift_fast
        IT[4 * i + 2] += shift_slow
        IT[4 * i + 3] += shift_fast

    if reapply_peripheral_margin:
      try:    peripheral_margin = self.working_params.distl.peripheral_margin
      except Exception: peripheral_margin = 0
      for i in xrange(len(self.working_params.distl.detector_tiling) // 4):
          IT[4 * i + 0] += peripheral_margin
          IT[4 * i + 1] += peripheral_margin
          IT[4 * i + 2] -= peripheral_margin
          IT[4 * i + 3] -= peripheral_margin

    return IT

#if __name__=='__main__':
#  import sys
#  i = sys.argv[1]
#  a = SaturnImage(i)
#  a.readHeader()
#  a.read()
