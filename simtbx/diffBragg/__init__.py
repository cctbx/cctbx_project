from __future__ import absolute_import, division, print_function
from boost_adaptbx import boost
import boost_adaptbx.boost.python as bp
import cctbx.uctbx # possibly implicit
from simtbx import nanoBragg # implicit import
ext = boost.python.import_ext("simtbx_diffBragg_ext")
import numpy as np
from simtbx_diffBragg_ext import *

@bp.inject_into(ext.diffBragg)
class _():
    def get_derivative_pixels(self, refine_id):
        return self.__get_derivative_pixels(refine_id)[:self.__number_of_pixels_modeled_using_diffBragg]

    def update_Fhkl_scale_factors(self, scales):
        assert isinstance(scales, np.ndarray)
        scales = self._check_contig(scales)
        if scales.dtype != np.float64:
            print("Warning, converting scale factors to double!")
            scales = scales.astype(np.float64)
        # TODO add num_asu property, then do the assertion
        #assert len(scales) == self.num_asu
        self.__update_Fhkl_scale_factors(scales)

    def _check_contig(self, arr):
        if not arr.flags.c_contiguous:
            print("WARNINGS:Making contiguous")
            arr = np.ascontiguousarray(arr)
        return arr

    def add_Fhkl_gradients(self, psf, residuals, variance, trusted):

        assert isinstance(residuals, np.ndarray)
        assert isinstance(variance, np.ndarray)
        assert isinstance(trusted, np.ndarray)
        residuals = self._check_contig(residuals)
        variance = self._check_contig(variance)
        trusted = self._check_contig(trusted)
        Npix = len(psf)/3
        assert Npix == len(residuals)
        assert Npix == len(variance)
        assert Npix == len(trusted)
        # TODO check or contiguous arrays..

        assert trusted.dtype==bool
        assert residuals.dtype==np.float64
        assert variance.dtype==np.float64
        return self.__add_Fhkl_gradients(psf, residuals, variance, trusted)
