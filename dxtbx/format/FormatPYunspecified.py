from __future__ import division

from dxtbx.format.FormatPY import FormatPY

class FormatPYunspecified(FormatPY):

  @staticmethod
  def understand(image_file):
    """Seems like the static method wastes a lot of effort here; it's not possible to
    just read the first few bytes; instead understand() reads the entire first data
    item in the file; an entire binary image.  This data is then read again in the
    _start() method and again in the detectorbase constructor."""
    try:
      stream = FormatPYunspecified.open_file(image_file, 'rb')
      import cPickle as pickle
      data = pickle.load(stream)
    except IOError,e:
      return False

    wanted_header_items = ['SIZE1','SIZE2','TIMESTAMP']

    for header_item in wanted_header_items:
      if not header_item in data:
        return False

    return True

  def __init__(self, image_file):
    '''Initialise the image structure from the given file.'''

    assert(self.understand(image_file))

    FormatPY.__init__(self, image_file)

  def _start(self):
    import os
    if isinstance(self._image_file, basestring) and os.path.isfile(self._image_file):
      stream = FormatPYunspecified.open_file(self._image_file, 'rb')

      import cPickle as pickle
      data = pickle.load(stream)
    else:
      data = self._image_file

    if not "DETECTOR_ADDRESS" in data:
      # legacy format; try to guess the address
      self.LCLS_detector_address = 'CxiDs1-0|Cspad-0'
      if "DISTANCE" in data and data["DISTANCE"] > 1000:
      # downstream CS-PAD detector station of CXI instrument
        self.LCLS_detector_address = 'CxiDsd-0|Cspad-0'
    else:
      self.LCLS_detector_address = data["DETECTOR_ADDRESS"]
    from xfel.detector_formats import reverse_timestamp
    self._timesec = reverse_timestamp( data["TIMESTAMP"] )[0]
    from xfel.detector_formats import detector_format_version as detector_format_function
    version_lookup = detector_format_function(self.LCLS_detector_address,self._timesec)
    self.start_helper(version_token="distl.detector_format_version=%s"%version_lookup)

  def start_helper(self, version_token):

    from spotfinder.applications.xfel import cxi_phil
    from iotbx.detectors.npy import NpyImage
    import os,copy

    is_file = isinstance(self._image_file, basestring) and os.path.isfile(self._image_file)

    if is_file:
      file_name = self._image_file
    else:
      file_name = "inmem"

    args = [file_name,
            version_token,
            "viewer.powder_arcs.show=False",
            "viewer.powder_arcs.code=3n9c",
           ]

    params = cxi_phil.cxi_versioned_extract(args)
    horizons_phil = params.persist.commands

    if is_file:
      I = NpyImage(file_name)
    else:
      print "This is not a file; assume the data are in the defined dictionary format"
      I = NpyImage(file_name, source_data=self._image_file)
    I.readHeader(horizons_phil)
    I.translate_tiles(horizons_phil)
    # necessary to keep the phil parameters for subsequent calls to get_tile_manager()
    I.horizons_phil_cache = copy.deepcopy(horizons_phil)
    self.detectorbase = I

  def _goniometer(self):

    return self._goniometer_factory.single_axis()

  def _detector(self):
    '''Return a model for a simple detector'''
    trusted_range = (self.detectorbase.parameters.get('MIN_TRUSTED_VALUE',0),
      self.detectorbase.saturation)

    return self._detector_factory.simple(
        sensor = 'PAD',
        distance = self.detectorbase.distance,
        beam_centre = (self.detectorbase.beamx,
                       self.detectorbase.beamy),
        fast_direction = '+x',
        slow_direction = '-y',
        pixel_size = (self.detectorbase.pixel_size,
                      self.detectorbase.pixel_size),
        image_size = (self.detectorbase.size2,
                      self.detectorbase.size1),
        trusted_range = trusted_range,
        mask = [])  # a list of dead rectangles

  def _beam(self):
    '''Return a simple model for the beam.'''

    return self._beam_factory.simple(self.detectorbase.wavelength)

  def _scan(self):
    '''Return the scan information for this image.'''

    if self.detectorbase.osc_start is not None and self.detectorbase.deltaphi is not None \
        and self.detectorbase.deltaphi > 0:
      format = ''

      exposure_time = self.detectorbase.parameters.get('TIME', 1)
      osc_start = self.detectorbase.osc_start
      osc_range = self.detectorbase.deltaphi
      timestamp = self._timesec

      return self._scan_factory.single(
          self._image_file, format, exposure_time,
          osc_start, osc_range, timestamp)

    else:
      return self._scan_factory.make_scan(
        image_range = (1,1),
        # femtosecond X-ray pulse, set this to a dummy argument
        exposure_times = [1.],
        oscillation = (0.0,0.0),
        epochs = {1:self._timesec}
        )

  def get_mask(self):
    '''Creates a mask merging untrusted pixels with active areas.'''
    from scitbx.array_family import flex
    detector_base = self.detectorbase
    # get effective active area coordinates
    tile_manager = detector_base.get_tile_manager(detector_base.horizons_phil_cache)
    tiling = tile_manager.effective_tiling_as_flex_int(reapply_peripheral_margin = True)
    # get the raw data to get the size of the mask
    data = self.get_raw_data()
    if tiling is None or len(tiling) == 0:
      return None

    # set the mask to the same dimensions as the data
    mask = flex.bool(flex.grid(data.focus()))

    # set active areas to True so they are not masked
    for i in xrange(len(tiling)//4):
      x1,y1,x2,y2=tiling[4*i:(4*i)+4]
      sub_array = flex.bool(flex.grid(x2-x1,y2-y1),True)
      mask.matrix_paste_block_in_place(sub_array,x1,y1)

    # create untrusted pixel mask
    detector = self.get_detector()
    assert len(detector) == 1
    trusted_mask = detector[0].get_trusted_range_mask(data)

    # returns merged untrusted pixels and active areas using bitwise AND (pixels are accepted
    # if they are inside of the active areas AND inside of the trusted range)
    return (mask & trusted_mask,)

class FormatPYunspecifiedInMemory(FormatPYunspecified):
  """ Overrides the Format object's init method to accept an image dictionary
      instead of a file name. Used with XFELs when it is desirable to never write
      a file to disk, but to process it only in memory.
  """

  @staticmethod
  def understand(image_file):
    """ If it's an image dictionary, we understand this """
    try:
      wanted_header_items = ['SIZE1','SIZE2','TIMESTAMP']

      for header_item in wanted_header_items:
        if not header_item in image_file.keys():
          return False

      return True

    except AttributeError:
      return False

  def __init__(self, data):
    """ @param data In memory image dictionary, alredy initialized """
    FormatPYunspecified.__init__(self, data)

    import copy
    self._image_file = copy.deepcopy(data)


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatPYunspecified.understand(arg)
