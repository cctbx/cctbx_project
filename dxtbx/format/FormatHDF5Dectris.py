from __future__ import absolute_import, division, print_function

from dxtbx.format.Format import Format
from dxtbx.format.FormatHDF5 import FormatHDF5


class FormatHDF5Dectris(FormatHDF5):
    """ This is a placeholder implementation only.  Open example dataset
      provided by Dectris Ltd, Jan 2013.  Read the first image only.
      Later replace this with a sweep-oriented implementation that
      reads the whole dataset."""

    @staticmethod
    def understand(image_file):
        return False
        try:
            tag = FormatHDF5Dectris.open_file(image_file, "rb").read(8)
        except IOError:
            return False

        return tag == "\211HDF\r\n\032\n"

    def __init__(self, image_file, **kwargs):

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatHDF5.__init__(self, image_file, **kwargs)

    def _start(self):
        from iotbx.detectors.eiger import EIGERImage

        self.detectorbase = EIGERImage(self._image_file)
        self.detectorbase.readHeader()

    def _goniometer(self):

        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector"""

        return self._detector_factory.simple(
            sensor="PAD",
            distance=self.detectorbase.distance,
            beam_centre=(self.detectorbase.beamx, self.detectorbase.beamy),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(self.detectorbase.pixel_size, self.detectorbase.pixel_size),
            image_size=(self.detectorbase.size1, self.detectorbase.size2),
            trusted_range=(-1, self.detectorbase.saturation),
            mask=[],
        )  # a list of dead rectangles

    def _beam(self):
        """Return a simple model for the beam."""

        return self._beam_factory.simple(self.detectorbase.wavelength)

    def _scan(self):
        """Replace this with true values later
       when HDF5 sweep support is implemented."""

        return self._scan_factory.make_scan(
            image_range=(1, 1),
            # dummy value--return to this later please
            exposure_time=1,
            oscillation=(
                self.detectorbase.osc_start,
                self.detectorbase.osc_start + self.detectorbase.deltaphi,
            ),
            epochs={1: 0.0},  # Later get the true time values from HDF5 file
        )

    def get_num_images(self):
        return self.detectorbase.image_count()

    def get_goniometer(self, index=None):
        return Format.get_goniometer(self)

    def get_detector(self, index=None):
        return Format.get_detector(self)

    def get_beam(self, index=None):
        return Format.get_beam(self)

    def get_scan(self, index=None):
        return Format.get_scan(self)

    def get_raw_data(self, index=None):
        return Format.get_raw_data(self)

    def get_detectorbase(self, index=None):
        self.detectorbase.img_number = index
        return self.detectorbase

    def get_image_file(self, index=None):
        return self.detectorbase.get_data_link(index)


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatHDF5Dectris.understand(arg))
