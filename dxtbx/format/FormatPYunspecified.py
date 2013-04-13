from __future__ import division

from dxtbx.format.FormatPY import FormatPY

class FormatPYunspecified(FormatPY):

    @staticmethod
    def understand(image_file):
        try:
          stream = FormatPYunspecified.open_file(image_file, 'rb')
          import cPickle as pickle
          data = pickle.load(stream)
        except IOError,e:
          return False

        wanted_header_items = ['SIZE1','SIZE2']

        for header_item in wanted_header_items:
            if not header_item in data:
                return False

        unwanted_header_items = ['DETECTOR_FORMAT_VERSION']

        for header_item in unwanted_header_items:
            if header_item in data:
                return False

        return True

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        assert(self.understand(image_file))

        FormatPY.__init__(self, image_file)

    def start_helper(self, version_token):

        from spotfinder.applications.xfel import cxi_phil
        from iotbx.detectors.npy import NpyImage
        import os
        args = [self._image_file,
                version_token,
                "viewer.powder_arcs.show=False",
                "viewer.powder_arcs.code=3n9c",
               ]

        params = cxi_phil.cxi_versioned_extract(args)
        horizons_phil = params.persist.commands
        if isinstance(self._image_file, basestring) and os.path.isfile(self._image_file):
          I = NpyImage(self._image_file)
        else:
          print "This is not a file; assume the data are in the defined dictionary format"
          I = NpyImage(self._image_file, source_data=params.indexing.data)
        I.readHeader(horizons_phil)
        I.translate_tiles(horizons_phil)
        self.detectorbase = I

    def _goniometer(self):

        return self._goniometer_factory.single_axis()

    def _detector(self):
        '''Return a model for a simple detector'''

        return self._detector_factory.simple(
            sensor = 'PAD',
            distance = self.detectorbase.distance,
            beam_centre = (self.detectorbase.beamx,
                           self.detectorbase.beamy),
            fast_direction = '+x',
            slow_direction = '-y',
            pixel_size = (self.detectorbase.pixel_size,
                          self.detectorbase.pixel_size),
            image_size = (self.detectorbase.size1,
                          self.detectorbase.size2),
            trusted_range = (0, self.detectorbase.saturation),
            mask = [])  # a list of dead rectangles

    def _beam(self):
        '''Return a simple model for the beam.'''

        return self._beam_factory.simple(self.detectorbase.wavelength)

    def _scan(self):
        '''Return the scan information for this image.'''

        return self._scan_factory.make_scan(
          image_range = (1,1),
          # femtosecond X-ray pulse, set this to a dummy argument
          exposure_time = 1.,
          oscillation = (0.0,0.0),
          epochs = {1:0.}  # Later see if we can actually get the millisecond time stamp in here
          )


if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatPYunspecified.understand(arg)
