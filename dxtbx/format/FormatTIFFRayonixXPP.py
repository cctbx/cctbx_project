#!/usr/bin/env python
# FormatTIFFRayonixXPP.py
# Sub class of FormatTIFFRayonix specialized for the XPP Rayonix dectector at LCLS
#
# Images from the XPP Rayonix detector have several unitialized values, such as
# distance, wavelength, etc.  Set these values to zero so the images can be at
# least viewed.
#

from __future__ import absolute_import, division, print_function

import re
import struct

from dxtbx.format.FormatTIFFRayonix import FormatTIFFRayonix


def check(l):
    """ Sets l or values in l that are less than zero to zero """
    if not isinstance(l, list) and not isinstance(l, tuple):
        if l < 0:
            return 0
        else:
            return l
    ret = []
    for val in l:
        if isinstance(val, list):
            ret.append(check(val))
        else:
            if val < 0:
                val = 0
            ret.append(val)

    if isinstance(l, tuple):
        return tuple(ret)
    return ret


class FormatTIFFRayonixXPP(FormatTIFFRayonix):
    """A class for reading TIFF format Rayonix images, and correctly
  constructing a model for the experiment from this."""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an XPP Rayonix TIFF """

        def get(ftype, index, typelen, offset=1024):
            f = open(image_file, "rb")
            si_format = "<"
            f.seek(offset + index)
            rawdata = f.read(typelen)
            return struct.unpack(si_format + ftype, rawdata)[0]

        data = [get("c", 1152 + i, 1) for i in xrange(128)]
        filepath = "".join(
            [c for c in data if c != " " and ord(c) < 127 and ord(c) > 32]
        )

        pattern = re.compile(".*xpp[a-zA-Z][0-9][0-9][0-9][0-9].*")

        return pattern.match(filepath)

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
    proper model of the experiment."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        from iotbx.detectors.mar import MARImage

        MARImage._read_header_asserts = lambda self: None
        FormatTIFFRayonix.__init__(self, image_file, **kwargs)

        return

    ####################################################################
    #                                                                  #
    # Helper methods to get all of the values out of the TIFF header   #
    # - separated out to assist with code clarity                      #
    #                                                                  #
    ####################################################################

    def _get_rayonix_beam_xy(self):
        """Get the beam x, y positions which are defined in the standard
    to be in pixels. X and Y are not defined by the documentation, so
    taking as a given that these are horizontal and vertical. N.B.
    the documentation states that the horizontal direction is fast."""

        beam_x, beam_y = struct.unpack(self._ii, self._tiff_header_bytes[1668:1676])[:2]
        pixel_x, pixel_y = struct.unpack(self._ii, self._tiff_header_bytes[1796:1804])[
            :2
        ]

        return beam_x * 1000 / pixel_x, beam_y * 1000 / pixel_y

    def _get_rayonix_detector_rotations(self):
        return check(FormatTIFFRayonix._get_rayonix_detector_rotations(self))

    def _get_rayonix_distance(self):
        return check(FormatTIFFRayonix._get_rayonix_distance(self))

    def _get_rayonix_pixel_size(self):
        return check(FormatTIFFRayonix._get_rayonix_pixel_size(self))

    def _get_rayonix_times(self):
        return check(FormatTIFFRayonix._get_rayonix_times(self))

    def _get_rayonix_timestamp(self):
        return check(FormatTIFFRayonix._get_rayonix_timestamp(self))

    def _get_rayonix_scan_angles(self):
        return check(FormatTIFFRayonix._get_rayonix_scan_angles(self))

    def _get_rayonix_detector_rotations(self):
        return check(FormatTIFFRayonix._get_rayonix_detector_rotations(self))


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatTIFFRayonixXPP.understand(arg))
