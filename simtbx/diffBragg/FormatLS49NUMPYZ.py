from __future__ import absolute_import, division

import numpy as np

from dxtbx.format.Format import Format
from dxtbx.format.FormatStill import FormatStill
from dials.array_family import flex
from dxtbx.model.detector import DetectorFactory
from dxtbx.model.beam import BeamFactory


class FormatLS49NUMPYZ(FormatStill):
    """
    Class for reading D9114 simulated monolithic cspad data
    """
    @staticmethod
    def understand(image_file):
        return image_file.endswith(".npz")

    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError
        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatStill.__init__(self, image_file, **kwargs)

        self._handle = np.load(image_file)
        self._geometry_define()

    def _geometry_define(self):
        self._cctbx_detector = self._detector_factory.from_dict(self._handle["det"][()])
        self._cctbx_beam = self._beam_factory.from_dict(self._handle["beam"][()])

    def get_num_images(self):
        return 1

    def load_panel_img(self):
        self.panel_img = self._handle["img"]
        if not self.panel_img.dtype == np.float64:
            self.panel_img = self.panel_img.astype(np.float64)

    def get_raw_data(self, index=0):
        self.load_panel_img()
        return flex.double(self.panel_img)

    def get_detectorbase(self, index=None):
        raise NotImplementedError

    def get_image_file(self, index=None):
        return Format.get_image_file(self)

    def get_detector(self, index=None):
        return self._cctbx_detector

    def get_beam(self, index=None):
        return self._cctbx_beam

    @staticmethod
    def get_instrument_name(handle):
        if "short_name" in handle["/entry/instrument"].attrs:
            name = handle["/entry/instrument"].attrs["short_name"]
        elif "/entry/instrument/name" in handle:
            if "short_name" in handle["/entry/instrument/name"].attrs:
                name = handle["/entry/instrument/name"].attrs["short_name"]
            else:
                name = handle["/entry/instrument/name"][()]
        else:
            name = None
        return name


if __name__ == '__main__':
    import sys
    for arg in sys.argv[1:]:
        print(FormatLS49NUMPYZ.understand(arg))
