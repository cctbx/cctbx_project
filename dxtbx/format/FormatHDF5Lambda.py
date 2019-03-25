#!/usr/bin/env python
# FormatHDF5Lambda.py
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Experimental format for the X-Spectrum LAMBDA detector
# http://www.x-spectrum.de/

from __future__ import absolute_import, division, print_function

from dxtbx.format.Format import Format
from dxtbx.format.FormatHDF5 import FormatHDF5


class FormatHDF5Lambda(FormatHDF5):
    @staticmethod
    def understand(image_file):
        try:
            tag = FormatHDF5.open_file(image_file, "rb").read(8)
        except IOError:
            return False

        # check that this is a HDF5 file (should not have got here if not
        # anyway...)

        if tag != "\211HDF\r\n\032\n":
            return False

        import h5py

        h5_handle = h5py.File(image_file, "r")

        try:
            desc = h5_handle["entry/instrument/detector/description"]
        except KeyError:
            h5_handle.close()
            return False

        if "Lambda" in desc[()][0]:
            h5_handle.close()
            return True

        return False

    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatHDF5.__init__(self, image_file, **kwargs)

    def _start(self):
        import h5py

        self._h5_handle = h5py.File(self.get_image_file(), "r")

    def _goniometer(self):
        """Dummy goniometer - assume vertical (EMBl-HH P14)"""

        return self._goniometer_factory.known_axis((0, 1, 0))

    def _detector(self):
        """Partly dummy detector"""
        from scitbx import matrix

        # Get the detector geometry
        entry = self._h5_handle["entry"]
        instrument = entry["instrument"]
        detector = instrument["detector"]

        # Initialise detector frame - origin at 0,0,0
        fast = matrix.col((1.0, 0.0, 0.0))
        slow = matrix.col((0.0, 1.0, 0.0))
        orig = matrix.col((0.0, 0.0, 0.0))

        # Get the pixel and image size
        pixel_size = (
            1.0e-3 * detector["x_pixel_size"].value,
            1.0e-3 * detector["y_pixel_size"].value,
        )
        layout = detector["layout"].value[0].split("X")
        image_size = int(layout[0]), int(layout[1])
        trusted_range = (-1, detector["saturation_value"][0])
        thickness = float(detector["sensor_thickness"].value) / 1000.0
        material = str(detector["sensor_material"].value[0])

        # Make the detector
        detector = self._detector_factory.make_detector(
            "PAD",
            fast,
            slow,
            orig,
            pixel_size,
            image_size,
            trusted_range,
            name="Panel",
            thickness=thickness,
            material=material,
        )

        # At the moment, beam is a dummy object because wavelength is not set in
        # the header. Therefore, the px<-->mm strategy will generally be
        # incorrect. Set it anyway, to override later.
        beam = self._beam()
        wavelength = beam.get_wavelength()

        from cctbx.eltbx import attenuation_coefficient
        from dxtbx.model import ParallaxCorrectedPxMmStrategy

        # this will fail for undefined composite materials
        table = attenuation_coefficient.get_table(material)

        # mu_at_angstrom returns cm^-1, but need mu in mm^-1
        mu = table.mu_at_angstrom(wavelength) / 10.0

        for panel in detector:
            panel.set_mu(mu)
            panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, thickness))

        return detector

    def _beam(self):
        """Dummy beam"""

        wavelength = 1.0
        return self._beam_factory.simple(wavelength)

    def _scan(self):
        """Dummy scan"""

        entry = self._h5_handle["entry"]
        nframes = int(entry["instrument/detector/collection/frame_numbers"].value)
        image_range = (1, nframes)
        exposure_times = 0.0
        oscillation = (0, 1)
        epochs = [0] * nframes

        return self._scan_factory.make_scan(
            image_range, exposure_times, oscillation, epochs, deg=True
        )

    def get_num_images(self):
        detector = self._h5_handle["entry/instrument/detector"]
        data = detector["data"]
        return data.shape[0]

    def get_goniometer(self, index=None):
        return Format.get_goniometer(self)

    def get_detector(self, index=None):
        return Format.get_detector(self)

    def get_beam(self, index=None):
        return Format.get_beam(self)

    def get_scan(self, index=None):
        if index == None:
            return Format.get_scan(self)
        else:
            scan = Format.get_scan(self)
            return scan[index]

    def get_raw_data(self, index):
        from scitbx.array_family import flex

        detector = self._h5_handle["entry/instrument/detector"]
        data = detector["data"]
        im = data[index, :, :].astype("int32")  # convert from int16
        return flex.int(im)

    def get_image_file(self, index=None):
        return Format.get_image_file(self)


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        print(FormatHDF5Lambda.understand(arg))
