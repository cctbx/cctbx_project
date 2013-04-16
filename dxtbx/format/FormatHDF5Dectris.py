from __future__ import division

from dxtbx.format.Format import Format

class FormatHDF5Dectris(Format):
    """ This is a placeholder implementation only.  Open example dataset
        provided by Dectris Ltd, Jan 2013.  Read the first image only.
        Later replace this with a sweep-oriented implementation that
        reads the whole dataset."""
    @staticmethod
    def understand(image_file):
        try:
          tag = FormatHDF5Dectris.open_file(image_file, 'rb').read(8)
        except IOError,e:
          return False

        return tag == "\211HDF\r\n\032\n"

    def __init__(self, image_file):

        assert(self.understand(image_file))

        Format.__init__(self, image_file)

    def _start(self):
        from iotbx.detectors.eiger import EIGERImage
        self.detectorbase = EIGERImage(self._image_file)
        self.detectorbase.readHeader()

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
        '''Replace this with true values later
           when HDF5 sweep support is implemented.'''

        return self._scan_factory.make_scan(
          image_range = (1,1),
          #dummy value--return to this later please
          exposure_time = 1,
          oscillation = (self.detectorbase.osc_start,
            self.detectorbase.osc_start + self.detectorbase.deltaphi),
          epochs = {1:0.} # Later get the true time values from HDF5 file
          )

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatHDF5Dectris.understand(arg)
