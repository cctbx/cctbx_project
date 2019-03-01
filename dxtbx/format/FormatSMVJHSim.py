#!/usr/bin/env python
# FormatSMVJHSim.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for JHSim images. Inherits from
# FormatSMV.

from __future__ import division, print_function

import time

from dxtbx.format.FormatSMV import FormatSMV


class FormatSMVJHSim(FormatSMV):
    """A class for reading SMV format JHSim images, and correctly constructing
  a model for the experiment from this."""

    # all ADSC detectors generate images with an ADC offset of 40
    # for Mar/Rayonix it is 10
    # Rigaku SMV uses 20, and 5 for image plate formats
    # for one particular simulation, I used 1
    ADC_OFFSET = 1
    image_pedestal = 1

    @staticmethod
    def understand(image_file):
        '''Check to see if this looks like an JHSim SMV format image, i.e. we can
    make sense of it. From JH: "The best way to identify images from any of my
    simulators is to look for BEAMLINE=fake in the header."'''

        size, header = FormatSMV.get_smv_header(image_file)

        if header.get("BEAMLINE") == "fake":
            return True
        else:
            return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file, including a
    proper model of the experiment."""

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatSMV.__init__(self, image_file, **kwargs)

        return

    def _start(self):

        FormatSMV._start(self)

    def detectorbase_start(self):
        if not hasattr(self, "detectorbase") or self.detectorbase is None:
            from iotbx.detectors import SMVImage

            self.detectorbase = SMVImage(self._image_file)
            self.detectorbase.open_file = self.open_file
            self.detectorbase.readHeader()

    def _goniometer(self):
        """Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header."""

        return self._goniometer_factory.single_axis()

    def _detector(self):
        """Return a model for a simple detector, presuming no one has
    one of these on a two-theta stage. Assert that the beam centre is
    provided in the Mosflm coordinate frame."""

        distance = float(self._header_dictionary["DISTANCE"])
        beam_x = float(self._header_dictionary["BEAM_CENTER_X"])
        beam_y = float(self._header_dictionary["BEAM_CENTER_Y"])
        pixel_size = float(self._header_dictionary["PIXEL_SIZE"])
        image_size = (
            float(self._header_dictionary["SIZE1"]),
            float(self._header_dictionary["SIZE2"]),
        )
        image_pedestal = 1
        try:
            image_pedestal = float(self._header_dictionary["ADC_OFFSET"])
        except (KeyError):
            pass
        overload = 65535 - image_pedestal
        underload = 1 - image_pedestal

        # interpret beam center conventions
        image_width_mm = pixel_size * image_size[0]
        image_height_mm = pixel_size * image_size[1]
        adxv_beam_center = (beam_x, beam_y)
        mosflm_beam_center = (
            image_height_mm - adxv_beam_center[1],
            adxv_beam_center[0],
        )
        xds_beam_center = (
            adxv_beam_center[0] / pixel_size,
            (image_height_mm - adxv_beam_center[1]) / pixel_size,
        )
        cctbx_beam_center = (
            adxv_beam_center[0] + pixel_size,
            image_height_mm - adxv_beam_center[1] + pixel_size,
        )

        # Guess whether this is mimicking a Pilatus, if so set detector type so
        # that spot-finding parameters are appropriate
        if pixel_size == 0.172:
            stype = "SENSOR_PAD"
        else:
            stype = "CCD"

        return self._detector_factory.simple(
            stype,
            distance,
            cctbx_beam_center,
            "+x",
            "-y",
            (pixel_size, pixel_size),
            image_size,
            (underload, overload),
            [],
        )

    def _beam(self):
        """Return a simple model for the beam."""

        wavelength = float(self._header_dictionary["WAVELENGTH"])

        return self._beam_factory.simple(wavelength)

    def _scan(self):
        """Return the scan information for this image."""
        import calendar

        format = self._scan_factory.format("SMV")
        exposure_time = 1
        epoch = None

        # PST, PDT timezones not recognised by default...

        epoch = 0
        try:
            date_str = self._header_dictionary["DATE"]
            date_str = date_str.replace("PST", "").replace("PDT", "")
        except KeyError:
            date_str = ""
        for format_string in ["%a %b %d %H:%M:%S %Y", "%a %b %d %H:%M:%S %Z %Y"]:
            try:
                epoch = calendar.timegm(time.strptime(date_str, format_string))
                break
            except ValueError:
                pass

        # assert(epoch)
        osc_start = float(self._header_dictionary["OSC_START"])
        osc_range = float(self._header_dictionary["OSC_RANGE"])

        return self._scan_factory.single(
            self._image_file, format, exposure_time, osc_start, osc_range, epoch
        )

    def get_raw_data(self):
        """Get the pixel intensities (i.e. read the image and return as a
    flex array of integers.)"""

        from boost.python import streambuf
        from dxtbx import read_uint16, read_uint16_bs, is_big_endian
        from scitbx.array_family import flex

        assert len(self.get_detector()) == 1
        image_pedestal = 1
        try:
            image_pedestal = float(self._header_dictionary["ADC_OFFSET"])
        except KeyError:
            pass
        panel = self.get_detector()[0]
        size = panel.get_image_size()
        f = FormatSMVJHSim.open_file(self._image_file, "rb")
        f.read(self._header_size)

        if self._header_dictionary["BYTE_ORDER"] == "big_endian":
            big_endian = True
        else:
            big_endian = False

        if big_endian == is_big_endian():
            raw_data = read_uint16(streambuf(f), int(size[0] * size[1]))
        else:
            raw_data = read_uint16_bs(streambuf(f), int(size[0] * size[1]))

        # apply image pedestal, will result in *negative pixel values*
        # this is a horrible idea since raw_data is unsigned
        # see all instances of image_pedestal in dxtbx

        raw_data -= int(image_pedestal)

        image_size = panel.get_image_size()
        raw_data.reshape(flex.grid(image_size[1], image_size[0]))

        return raw_data


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatSMVJHSim.understand(arg))
