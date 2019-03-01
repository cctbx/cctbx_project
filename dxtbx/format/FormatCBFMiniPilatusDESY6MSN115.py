#!/usr/bin/env python
# FormatCBFMiniPilatusDESY6MSN115.py
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the FormatCBFMiniPilatus image reader for the P6M
# detector at PETRA III beamline P14, which has a vertical goniometer.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus


class FormatCBFMiniPilatusDESY6MSN115(FormatCBFMiniPilatus):
    """A class for reading mini CBF format Pilatus images for 6M SN 115 @ DESY,
    which has a vertical goniometer axis."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an Pilatus mini CBF format image,
        i.e. we can make sense of it."""

        header = FormatCBFMiniPilatus.get_cbf_header(image_file)

        for record in header.split("\n"):
            if "Detector: PILATUS 6M-F, S/N 60-0115-F" in record:
                return True

        return False

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header, though for miniCBF
        there are limited options for this."""

        return self._goniometer_factory.known_axis((0, 1, 0))


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatCBFMiniPilatusDESY6MSN115.understand(arg))
