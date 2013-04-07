# Implementation of an ImageFormat class to read MarIP-format image

from __future__ import division

from dxtbx.format.Format import Format

class FormatMarIP(Format):
    '''An image reading class for MarIP-format images
    Positive identification:  first 140 bytes contain the string "mar research"
    '''

    @staticmethod
    def understand(image_file):
        try:
          tag = FormatMarIP.open_file(image_file, 'rb').read(140)
        except IOError:
          return False
        return tag.find("mar research")>0

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        assert(self.understand(image_file))

        Format.__init__(self, image_file)

    def _start(self):
        '''Open the image file, read the image header, copy the key / value
        pairs into an internal dictionary self._header_dictionary along with
        the length of the header in bytes self._header_size.'''
        from iotbx.detectors.marIP import MARIPImage
        self.maripimage = MARIPImage(self._image_file)
        self.maripimage.readHeader()

    def _goniometer(self):

        return self._goniometer_factory.single_axis()

    def _detector(self):
        '''Return a model for a simple detector, which at the moment insists
        that the offsets and rotations are all 0.0.'''

        assert self.maripimage.parameters["TWOTHETA"] == 0.0

        return self._detector_factory.simple(
            sensor = 'IMAGE_PLATE',
            distance = self.maripimage.parameters["DISTANCE"],
            beam_centre = (self.maripimage.parameters["BEAM_CENTER_X"],
                           self.maripimage.parameters["BEAM_CENTER_Y"]),
            fast_direction = '+x',
            slow_direction = '-y',
            pixel_size = (self.maripimage.parameters["PIXEL_SIZE"],
                          self.maripimage.parameters["PIXEL_SIZE"]),
            image_size = (self.maripimage.parameters["SIZE1"],
                          self.maripimage.parameters["SIZE2"]),
            trusted_range = (0, self.maripimage.parameters["CCD_IMAGE_SATURATION"]),
            mask = [])

    def _beam(self):
        '''Return a simple model for the beam.'''

        return self._beam_factory.simple(self.maripimage.parameters["WAVELENGTH"])

    def _scan(self):
        '''Return the scan information for this image.'''

        return self._scan_factory.single(
          filename = self._image_file,
          format = "MARIP",
          exposure_time = self.maripimage.adaptor.exposure_time(),
          osc_start = self.maripimage.parameters["OSC_START"],
          osc_width = self.maripimage.parameters["OSC_RANGE"],
          epoch = None)

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatMarIP.understand(arg)
