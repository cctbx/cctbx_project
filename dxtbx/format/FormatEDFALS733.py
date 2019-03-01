# Implementation of an ImageFormat class to read MarIP-format image

from __future__ import absolute_import, division, print_function

from dxtbx.format.Format import Format


class FormatEDFALS733(Format):
    """EDF is the ESRF data format.
  In present context, ALS 7.3.3 uses this format for its single-tile Pilatus
  """

    @staticmethod
    def understand(image_file):
        try:
            tag = FormatEDFALS733.open_file(image_file, "rb").read(10)
        except IOError:
            return False

        return tag == "{\nHeaderID"

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        Format.__init__(self, image_file, **kwargs)

    def detectorbase_start(self):
        pass

    def _start(self):
        """Open the image file, read the image header, copy the key / value
    pairs into an internal dictionary self._header_dictionary along with
    the length of the header in bytes self._header_size."""
        from iotbx.detectors import EDFWrapper

        self.detectorbase = EDFWrapper(self._image_file)
        self.detectorbase.readHeader()

    def _goniometer(self):

        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector"""

        twotheta = 0.0

        return self._detector_factory.simple(
            sensor="PAD",
            distance=self.detectorbase.parameters["DISTANCE"],
            beam_centre=(
                self.detectorbase.parameters["BEAM_CENTER_X"],
                self.detectorbase.parameters["BEAM_CENTER_Y"],
            ),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(
                self.detectorbase.parameters["PIXEL_SIZE"],
                self.detectorbase.parameters["PIXEL_SIZE"],
            ),
            image_size=(
                self.detectorbase.parameters["SIZE1"],
                self.detectorbase.parameters["SIZE2"],
            ),
            trusted_range=(0, self.detectorbase.saturation),
            mask=[],
        )  # a list of dead rectangles

    def _beam(self):
        """Return a simple model for the beam."""

        return self._beam_factory.simple(wavelength=1.0)  # dummy argument

    def _scan(self):
        """Return the scan information for this image."""

        return self._scan_factory.single(
            filename=self._image_file,
            format="EDF",
            exposure_times=self.detectorbase.parameters["count_time"],
            osc_start=0.0,
            osc_width=0.0,
            epoch=None,
        )


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatEDFALS733.understand(arg))
