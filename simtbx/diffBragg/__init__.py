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

    def update_Fhkl_channels(self, channel_ids):
        """channel_ids is a  numpy int array the same length as the number of energy sources,
        that specifies the mapping of structure factor to energy channel,
        allowing one to refine multiple energy-dependent structure factors
        for example in a two-color experiment
        """
        assert isinstance(channel_ids, np.ndarray)
        channel_ids = self._check_contig(channel_ids)
        if channel_ids.dtype != int:
            channel_ids = channel_ids.astype(int)
        # assert len(channel_ids)== number_of_sources
        self.__update_Fhkl_channels(channel_ids)

    def update_Fhkl_scale_factors(self, scales, num_Fhkl_channels):
        assert isinstance(scales, np.ndarray)
        scales = self._check_contig(scales)
        if scales.dtype != np.float64:
            print("Warning, converting scale factors to double!")
            scales = scales.astype(np.float64)
        # TODO add num_asu property, then do the assertion
        #assert len(scales) == self.num_asu*num_Fhkl_channels
        self.__update_Fhkl_scale_factors(scales, num_Fhkl_channels)

    def _check_contig(self, arr):
        if not arr.flags.c_contiguous:
            print("WARNINGS:Making contiguous")
            arr = np.ascontiguousarray(arr)
        return arr

    def add_Fhkl_gradients(self, psf, residuals, variance, trusted, freq, num_Fhkl_channels, spot_scale, track=False):
        """

        :param psf:
        :param residuals:
        :param variance:
        :param trusted:
        :param freq:
        :param num_Fhkl_channels:
        :param spot_scale:
        :return:
        """

        assert isinstance(residuals, np.ndarray)
        assert isinstance(variance, np.ndarray)
        assert isinstance(trusted, np.ndarray)
        assert isinstance(freq, np.ndarray)
        residuals = self._check_contig(residuals)
        variance = self._check_contig(variance)
        trusted = self._check_contig(trusted)
        Npix = len(psf)/3
        assert Npix == len(residuals)
        assert Npix == len(variance)
        assert Npix == len(trusted)
        assert Npix == len(freq)
        # TODO check or contiguous arrays..

        assert trusted.dtype==bool
        assert freq.dtype==np.uintc
        assert residuals.dtype==np.float64
        assert variance.dtype==np.float64
        return self.__add_Fhkl_gradients(psf, residuals, variance, trusted, freq, num_Fhkl_channels, spot_scale, track)
