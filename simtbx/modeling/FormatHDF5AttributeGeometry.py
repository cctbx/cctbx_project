from __future__ import absolute_import, division

import numpy as np
import h5py
from copy import deepcopy
import ast

from dxtbx.format.FormatHDF5 import FormatHDF5
from dials.array_family import flex
from dxtbx.format.FormatStill import FormatStill
from dxtbx.model import Beam


class FormatHDF5AttributeGeometry(FormatHDF5, FormatStill):
    """
    Class for reading HDF5 files for arbitrary geometries
    focused on performance
    """
    @staticmethod
    def understand(image_file):
        try:
            img_handle = h5py.File(image_file, "r")
            keys = img_handle.keys()
        except (IOError, AttributeError) as err:
            return False
        if "images" not in keys:
            return False
        images = img_handle["images"]
        if "dxtbx_detector_string" not in images.attrs:
            return False
        if "dxtbx_beam_string" not in images.attrs:
            return False
        return True

    def _start(self):
        self._handle = h5py.File(self._image_file, "r")
        self.HAS_SPECTRUM_BEAM = False
        if hasattr(Beam, "set_spectrum"):
            self.HAS_SPECTRUM_BEAM = True
        self._image_dset = self._handle["images"]
        self._geometry_define()
        self._has_spectra = False
        self._has_central_wavelengths = False
        self._energies = None
        self._weights = None
        self._central_wavelengths = None
        self._check_per_shot_spectra()
        self._ENERGY_CONV = 12398.419739640716

    def _geometry_define(self):
        det_str = self._image_dset.attrs["dxtbx_detector_string"]
        beam_str = self._image_dset.attrs["dxtbx_beam_string"]
        try:
            det_str = det_str.decode()
            beam_str = beam_str.decode()
        except AttributeError:
            pass
        det_dict = ast.literal_eval(det_str)
        beam_dict = ast.literal_eval(beam_str)
        self._cctbx_detector = self._detector_factory.from_dict(det_dict)
        self._cctbx_beam = self._beam_factory.from_dict(beam_dict)

    def _check_per_shot_spectra(self):
        keys = list(self._handle.keys())
        has_energies = "spectrum_energies" in keys
        has_weights = "spectrum_weights" in keys
        has_central_wavelengths = "central_wavelengths" in keys
        if has_energies and has_weights:
            self._has_spectra = True
            self._energies = self._handle["spectrum_energies"]
            self._weights = self._handle["spectrum_weights"]
        elif has_central_wavelengths:
            self._has_central_wavelengths = True
            self._central_wavelengths = self._handle["central_wavelengths"]

    def get_num_images(self):
        return self._image_dset.shape[0]

    def get_raw_data(self, index=0):
        self.panels = self._image_dset[index]
        if self.panels.dtype == np.float64:
            flex_data = [flex.double(p) for p in self._image_dset[index]]
        else:
            flex_data = [flex.double(p.astype(np.float64)) for p in self._image_dset[index]]
        return tuple(flex_data)

    def get_detectorbase(self, index=None):
        raise NotImplementedError

    def get_detector(self, index=None):
        return self._cctbx_detector

    def _get_wavelength(self, index):
        if self._has_spectra:
            w = self._weights[index]
            E = self._energies[index]
            ave_E = (w*E).sum() / (w.sum())
            wavelength = self._ENERGY_CONV / ave_E
            if self.HAS_SPECTRUM_BEAM:
                self._w = w
                self._E = E
        elif self._has_central_wavelengths:
            wavelength = self._central_wavelengths[index]
        else:
            wavelength = None
        return wavelength

    def get_beam(self, index=0):
        beam = self._cctbx_beam
        wavelength = self._get_wavelength(index)
        if wavelength is not None:
            beam = deepcopy(self._cctbx_beam)
            beam.set_wavelength(wavelength)
        return beam


if __name__ == '__main__':
    import sys
    for arg in sys.argv[1:]:
        print(FormatHDF5AttributeGeometry.understand(arg))
