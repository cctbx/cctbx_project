from __future__ import absolute_import, division, print_function
from six.moves import range
import spotfinder.array_family.flex # implicit import

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("spotfinder_dxtbx_ext")
from spotfinder_dxtbx_ext import *

from libtbx import adopt_init_args

class Distl:

  def __init__(self, params, detector, beam, data):
    adopt_init_args(self, locals())

    npanels = len(detector)
    self.finderlist = []
    for n in range(npanels):
      SF = w_Distl(optionstring="",report_overloads=False)

      panel = detector[n]

      saturation = panel.get_trusted_range()[1]

      SF.setspotimg(panel = panel, beam = beam,
                      rawdata = data,
                      peripheral_margin = params.distl.peripheral_margin,
                      saturation = saturation)

      special_vendortype = {
          (2463,2527):"Pilatus-6M",
          (1475,1679):"Pilatus-2M",
          (487,195):"Pilatus-300K",
        }.get(panel.get_image_size(),"")
      print(special_vendortype)
      SF.set_tiling(special_vendortype)

      if params.distl.minimum_spot_area != None:
        SF.set_minimum_spot_area(params.distl.minimum_spot_area)
      if params.distl.minimum_signal_height != None:
        SF.set_minimum_signal_height(params.distl.minimum_signal_height)
      if params.distl.minimum_spot_height != None:
        SF.set_minimum_spot_height(params.distl.minimum_spot_height)
      if params.distl.spot_area_maximum_factor != None:
        SF.set_spot_area_maximum_factor(params.distl.spot_area_maximum_factor)
      if params.distl.peak_intensity_maximum_factor != None:
        SF.set_peak_intensity_maximum_factor(params.distl.peak_intensity_maximum_factor)

      SF.set_scanbox_windows(params.distl.scanbox_windows)

      #SF.parameter_guarantees() # XXX return to this.  As implemented the guarantees ruin the results
      SF.get_underload()

      SF.pxlclassify()
# Functions having resolution and thus need overrides:
# r2_to_resol and resol_to_r2 (search_icerings)
# xy2resol (pixelclassify_scanbox, search_spots, pixelisonice) DONE
      #SF.search_icerings()
      SF.search_maximas()
      SF.search_spots()
      SF.search_overloadpatches()
      SF.finish_analysis()
      from matplotlib import pyplot as plt
      plt.plot( [f.ctr_mass_x() for f in SF.spots],
                [f.ctr_mass_y() for f in SF.spots],"b.")
      plt.show()

      self.finderlist.append(SF)
