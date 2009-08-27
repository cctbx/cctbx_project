import boost.python
boost.python.import_ext("spotfinder_distltbx_ext")
from spotfinder_distltbx_ext import *
import spotfinder_distltbx_ext as ext
boost.python.import_ext("spotfinder_hough_ext")
from spotfinder_hough_ext import *

class Distl(w_Distl):

  def __init__(self,options,image,pd,report_overloads=False):
    w_Distl.__init__(self,options,report_overloads)
    try:    saturation = image.saturation
    except: saturation = 65535
    self.setspotimg(image.pixel_size, image.distance, image.wavelength,
                    float(pd['xbeam']),float(pd['ybeam']),image.rawdata,
                    saturation)
    #Fixes a longstanding gremlin:  my corrected xy must be propagated
    # to zepu's code; or else ice rings are treated incorrectly.

    #Setup tiling, if any.
    self.set_tiling(image.vendortype)

    self.get_underload()
    self.pxlclassify()
    self.search_icerings()
    self.search_maximas()
    self.search_spots()
    self.search_overloadpatches()
    self.finish_analysis()

class _SpotFilterAgent(boost.python.injector, SpotFilterAgent):
  def __getinitargs__(self):
    return (self.pixel_size, self.xbeam, self.ybeam, self.distance,
            self.wavelength, self.icerings)
