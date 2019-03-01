#!/usr/bin/env python
# FormatSMVADSCSNID19.py
#   Copyright (C) 2015 Diamond Light Source, Richard Gildea
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for ADSC images. Inherits from
# FormatSMVADSC, customised for example on APS ID19 SN 458 and 914
# which have reversed phi.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN


class FormatSMVADSCSNAPSID19(FormatSMVADSCSN):
    """A class for reading SMV format ADSC images, and correctly constructing
  a model for the experiment from this, for instrument numbers 458 and 914
  from the APS ID19 beamline."""

    @staticmethod
    def understand(image_file):
        """Check to see if this is ADSC SN 458 or 914."""

        # check this is detector serial number 458 or 914

        size, header = FormatSMVADSCSN.get_smv_header(image_file)

        if int(header["DETECTOR_SN"]) not in (458, 914):
            return False

        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
    proper model of the experiment."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatSMVADSCSN.__init__(self, image_file, **kwargs)

    def _detector(self):
        """Return a model for a simple detector, presuming no one has
    one of these on a two-theta stage. Assert that the beam centre is
    provided in the Mosflm coordinate frame."""

        distance = float(self._header_dictionary["DISTANCE"])
        pixel_size = float(self._header_dictionary["PIXEL_SIZE"])
        image_size = (
            float(self._header_dictionary["SIZE1"]),
            float(self._header_dictionary["SIZE2"]),
        )
        overload = 65535
        underload = 1

        key = [
            s for s in self._header_dictionary if s.endswith("_SPATIAL_BEAM_POSITION")
        ][0]

        beam_x, beam_y = [
            float(f) * pixel_size for f in self._header_dictionary[key].split()
        ]

        return self._detector_factory.simple(
            "CCD",
            distance,
            (beam_y, (image_size[1] * pixel_size) - beam_x),
            "+x",
            "-y",
            (pixel_size, pixel_size),
            image_size,
            (underload, overload),
            [],
            gain=self._adsc_module_gain(),
        )

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header."""

        return self._goniometer_factory.single_axis_reverse()

    def _scan(self):
        """Return the scan information for this image. There may be
    no timestamps in there..."""

        format = self._scan_factory.format("SMV")
        exposure_time = float(self._header_dictionary["TIME"])
        epoch = 0
        osc_start = float(self._header_dictionary["OSC_START"])
        osc_range = float(self._header_dictionary["OSC_RANGE"])

        return self._scan_factory.single(
            self._image_file, format, exposure_time, osc_start, osc_range, epoch
        )


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatSMVADSCSN.understand(arg))
