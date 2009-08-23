"""
The position of the direct beam center is assumed to be stated in mm throughout.
The beam center is defined as the position of the beam when the twotheta setting is 0.0 degrees.

We consider two reference frames here:
  1) Instrument reference frame: Coordinate system tied to the face of the detector instrument.
    It is assumed that a) the origin of this coordinate system is located at one of the four corners of
                          the image block.
                       b) the origin is in the center of the corner pixel.
  2) Imageblock reference frame: Coordinate system tied to the slow (x) and fast (y) directions in the
     pixel data.

There are eight possible relationships between the instrument and imageblock frames.  The instrument origin
  can be in any of the four corners of the image; and for each of these, the x and y axes can be either the
  fast or slow dimensions of the image.

It is assumed that the image header gives the beam position in the instrument reference frame.
Diffraction is obviously observed and recorded in the image reference frame.
This module provides a function, convert_beam_instrument_to_imageblock() to give the corresponding
  beam position in the imageblock reference frame.  The Image Object attribute "beam_center_reference_frame"
  defines in which reference frame the beam is expressed.
"""

class beam_center_convention_definitions:
  def __init__(self,beam_center_convention):
    #axis reverse flag:
    #     False:  slow position determines BEAM_CENTER_X. fast position determines BEAM_CENTER_Y
    #     True:   fast position determines BEAM_CENTER_X. slow position determines BEAM_CENTER_Y
    self.ar_flag = bool(beam_center_convention & 1)

    #y-direction reverse flag:
    #     False:  BEAM_CENTER_Y increases in same direction as data pixel index
    #     True:   BEAM_CENTER_Y increases in opposite direction as data pixel index
    self.yr_flag = bool(beam_center_convention & 2)

    #x-direction reverse flag:
    #     False:  BEAM_CENTER_X increases in same direction as data pixel index
    #     True:   BEAM_CENTER_X increases in opposite direction as data pixel index
    self.xr_flag = bool(beam_center_convention & 4)

  def explicit_formulae(self):
    if beam_center_convention==0:  beam_center_in_pixels = slow,fast
    if beam_center_convention==1:  beam_center_in_pixels = fast,slow
    if beam_center_convention==2:  beam_center_in_pixels = slow,fastwidth - fast
    if beam_center_convention==3:  beam_center_in_pixels = fast,slowwidth - slow
    if beam_center_convention==4:  beam_center_in_pixels = slowwidth - slow,fast
    if beam_center_convention==5:  beam_center_in_pixels = fastwidth - fast,slow
    if beam_center_convention==6:  beam_center_in_pixels = slowwidth - slow,fastwidth - fast
    if beam_center_convention==7:  beam_center_in_pixels = fastwidth - fast,slowwidth - slow

class instrument_to_imageblock_relation:
  #while intended for conversion from instrument frame to imageblock frame, formulae will probably work
  # in both directions.
  def __init__(self,imageobject):
    self.p = imageobject.pixel_size
    self.input_beam_mm = (imageobject.beamx,imageobject.beamy)
    width_in_pixels = (imageobject.size1,imageobject.size2)
    self.width_in_mm = (width_in_pixels[0]*self.p,width_in_pixels[1]*self.p)

  def select(self,beam_center_convention):
    C = beam_center_convention_definitions(beam_center_convention)
    output_beam_mm = []
    for outidx in [0,1]:
      srcidx = int((not bool(outidx)) != (not C.ar_flag)) # is xor(outidx, ar_flag)
      direction_reverse = [C.xr_flag,C.yr_flag][outidx]
      imageblock_mm = self.input_beam_mm[srcidx]
      if direction_reverse:
        imageblock_mm = self.width_in_mm[srcidx] - imageblock_mm
      output_beam_mm.append(imageblock_mm)
    return tuple(output_beam_mm)

def convert_beam_instrument_to_imageblock(imageobject,beam_center_convention,force=False):
  if not force and imageobject.beam_center_reference_frame != "instrument": return
  converter = instrument_to_imageblock_relation(imageobject)
  imageobject.parameters['BEAM_CENTER_X'],\
  imageobject.parameters['BEAM_CENTER_Y']= converter.select(beam_center_convention)
  imageobject.beam_center_reference_frame = "imageblock"
  imageobject.beam_center_convention = beam_center_convention

def convert_beam_instrument_to_module(input_parameters,image_divider,moduleindex,beam_center_convention):
    C = beam_center_convention_definitions(beam_center_convention)

    input_beam_mm = (input_parameters['BEAM_CENTER_X'],input_parameters['BEAM_CENTER_Y'])

    width_in_pixels = (input_parameters['SIZE1'],input_parameters['SIZE2'])

    output_beam_mm = []
    for outidx in [0,1]:
      srcidx = int((not bool(outidx)) != (not C.ar_flag)) # is xor(outidx, ar_flag)
      direction_reverse = [C.xr_flag,C.yr_flag][outidx]

      #convert instrument beam_mm back to instrument pixels
      instrument_pixels = input_beam_mm[srcidx]/input_parameters['PIXEL_SIZE']
      if direction_reverse:
        instrument_pixels = width_in_pixels[srcidx] - instrument_pixels

      #translate to module_pixels
      module_tile = [image_divider.tile_slow_interval(moduleindex),
                     image_divider.tile_fast_interval(moduleindex)][srcidx]
      module_pixels = instrument_pixels - module_tile.first

      #convert module pixels to output_beam_mm
      if direction_reverse:
        module_pixels = module_tile.size() - module_pixels
      beam_mm = module_pixels * input_parameters['PIXEL_SIZE']
      output_beam_mm.append(beam_mm)
    return tuple(output_beam_mm)

if __name__=="__main__":
  from libtbx import adopt_init_args
  from libtbx.test_utils import show_diff
  import StringIO,sys
  class test_tile:
    def __init__(self,first,last):
      adopt_init_args(self, locals())
    def size(self): return self.last-self.first+1
  class test_divider:
    def tile_slow_interval(self,idx):
      if idx in [0,1,2]:  return test_tile(4,2043)
      if idx in [3,4,5]:  return test_tile(2052,4091)
      if idx in [6,7,8]:  return test_tile(4100,6139)
    def tile_fast_interval(self,idx):
      if idx in [0,3,6]:  return test_tile(4,2043)
      if idx in [1,4,7]:  return test_tile(2052,4091)
      if idx in [2,5,8]:  return test_tile(4100,6139)
  input_parameters = {'BEAM_CENTER_X':154.9,'BEAM_CENTER_Y':148.7,
                      'SIZE1':6144,'SIZE2':6144,'PIXEL_SIZE':0.051294}
  ID = test_divider()
  S = StringIO.StringIO()
  for convention in xrange(8):
    for moduleidx in xrange(9):
      B = convert_beam_instrument_to_module(input_parameters,ID,moduleidx,convention)
      print >>S,"(%.1f,%.1f)"%(B[0],B[1]),
    print >>S
  assert not show_diff(S.getvalue(),
"""(154.7,148.5) (154.7,43.4) (154.7,-61.6) (49.6,148.5) (49.6,43.4) (49.6,-61.6) (-55.4,148.5) (-55.4,43.4) (-55.4,-61.6)
(148.5,154.7) (43.4,154.7) (-61.6,154.7) (148.5,49.6) (43.4,49.6) (-61.6,49.6) (148.5,-55.4) (43.4,-55.4) (-61.6,-55.4)
(154.7,-61.6) (154.7,43.4) (154.7,148.5) (49.6,-61.6) (49.6,43.4) (49.6,148.5) (-55.4,-61.6) (-55.4,43.4) (-55.4,148.5)
(148.5,-55.4) (43.4,-55.4) (-61.6,-55.4) (148.5,49.6) (43.4,49.6) (-61.6,49.6) (148.5,154.7) (43.4,154.7) (-61.6,154.7)
(-55.4,148.5) (-55.4,43.4) (-55.4,-61.6) (49.6,148.5) (49.6,43.4) (49.6,-61.6) (154.7,148.5) (154.7,43.4) (154.7,-61.6)
(-61.6,154.7) (43.4,154.7) (148.5,154.7) (-61.6,49.6) (43.4,49.6) (148.5,49.6) (-61.6,-55.4) (43.4,-55.4) (148.5,-55.4)
(-55.4,-61.6) (-55.4,43.4) (-55.4,148.5) (49.6,-61.6) (49.6,43.4) (49.6,148.5) (154.7,-61.6) (154.7,43.4) (154.7,148.5)
(-61.6,-55.4) (43.4,-55.4) (148.5,-55.4) (-61.6,49.6) (43.4,49.6) (148.5,49.6) (-61.6,154.7) (43.4,154.7) (148.5,154.7)
""")
  print "OK"
