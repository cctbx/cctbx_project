from __future__ import absolute_import, division, print_function
from boost_adaptbx import boost
import boost_adaptbx.boost.python as bp
import cctbx.uctbx # possibly implicit
from simtbx import nanoBragg # implicit import
ext = boost.python.import_ext("simtbx_diffBragg_ext")
from simtbx_diffBragg_ext import *

@bp.inject_into(ext.diffBragg)
class _():
    def get_derivative_pixels(self, refine_id):
        return self.__get_derivative_pixels(refine_id)[:self.__number_of_pixels_modeled_using_diffBragg]
