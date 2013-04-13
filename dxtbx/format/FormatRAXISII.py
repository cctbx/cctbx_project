from __future__ import division

from dxtbx.format.Format import Format

class FormatRAXISII(Format):

    @staticmethod
    def understand(image_file):
        try:
          tag = FormatRAXISII.open_file(image_file, 'rb').read(7)
        except IOError,e:
          return False

        return tag == "R-AXIS2"

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        assert(self.understand(image_file))

        Format.__init__(self, image_file)

    def _start(self):
        from iotbx.detectors.raxis_nonsquare import NonSquareRAXISImage
        self.detectorbase = NonSquareRAXISImage(self._image_file)
        self.detectorbase.readHeader()

    def _goniometer(self):

        return self._goniometer_factory.single_axis()

    def _detector(self):
        '''Return a model for a simple detector'''

        return self._detector_factory.simple(
            sensor = 'IMAGE_PLATE',
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

        return self._scan_factory.single(
          filename = self._image_file,
          format = "Raxis2 image plate",
          exposure_time = 1,
          osc_start = self.detectorbase.parameters["OSC_START"],
          osc_width = self.detectorbase.parameters["OSC_RANGE"],
          epoch = None)

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatRAXISII.understand(arg)
