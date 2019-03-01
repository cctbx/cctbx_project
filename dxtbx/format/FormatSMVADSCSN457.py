#!/usr/bin/env python
# FormatSMVADSCSN457.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for ADSC images. Inherits from
# FormatSMVADSC, customised for example on Australian Synchrotron SN 457
# which has reversed phi.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN


class FormatSMVADSCSNSN457(FormatSMVADSCSN):
    """A class for reading SMV format ADSC images, and correctly constructing
  a model for the experiment from this, for instrument number 457."""

    @staticmethod
    def understand(image_file):
        """Check to see if this is ADSC SN 457."""

        # check this is detector serial number 457

        size, header = FormatSMVADSCSN.get_smv_header(image_file)

        if int(header["DETECTOR_SN"]) != 457:
            return False

        return True

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
    proper model of the experiment."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatSMVADSCSN.__init__(self, image_file, **kwargs)

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header."""

        return self._goniometer_factory.single_axis_reverse()


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatSMVADSC.understand(arg))
