import spotfinder.array_family.flex # implicit import

import boost.python
ext = boost.python.import_ext("spotfinder_distltbx_ext")
from spotfinder_distltbx_ext import *
boost.python.import_ext("spotfinder_hough_ext")
from spotfinder_hough_ext import *
from libtbx import adopt_init_args
from libtbx.utils import Sorry

class Distl(w_Distl):

  def __init__(self,options,image,pd,report_overloads=False,params=None):
    w_Distl.__init__(self,options,report_overloads)
    adopt_init_args(self, locals())

    try:    saturation = image.saturation
    except Exception: saturation = 65535
    try:    peripheral_margin = params.distl.peripheral_margin
    except Exception: peripheral_margin = 20
    self.setspotimg(pixel_size = image.pixel_size, distance = image.distance,
                    wavelength = image.wavelength, xbeam = float(pd['xbeam']),
                    ybeam = float(pd['ybeam']), rawdata = image.rawdata,
                    peripheral_margin = peripheral_margin,
                    saturation = saturation)
    #Fixes a longstanding gremlin:  my corrected xy must be propagated
    # to zepu's code; or else ice rings are treated incorrectly.

    #Setup tiling, if any.
    self.set_tiling(image.vendortype)

    self.deprecation_warnings()
    if params!=None:
        if params.distl.minimum_spot_area != None:
          self.set_minimum_spot_area(params.distl.minimum_spot_area)
        if params.distl.minimum_signal_height != None:
          self.set_minimum_signal_height(params.distl.minimum_signal_height)
        if params.distl.minimum_spot_height != None:
          self.set_minimum_spot_height(params.distl.minimum_spot_height)
        if params.distl.spot_area_maximum_factor != None:
          self.set_spot_area_maximum_factor(params.distl.spot_area_maximum_factor)
        if params.distl.peak_intensity_maximum_factor != None:
          self.set_peak_intensity_maximum_factor(params.distl.peak_intensity_maximum_factor)
        self.set_scanbox_windows(params.distl.scanbox_windows)
        if params.distl.detector_tiling != None:
          IT = image.get_tile_manager(params
               ).effective_tiling_as_flex_int()

          self.set_tiling(detector_tiling = IT,
                          peripheral_margin = peripheral_margin)

    self.parameter_guarantees()

    self.get_underload()
    try:
      self.pxlclassify()
    except Exception, e:
      if str(e).find("cannot distinguish signal")>0: print e
      else: raise e
    self.search_icerings()
    self.search_maximas()
    self.search_spots()
    self.search_overloadpatches()
    self.finish_analysis()

  def deprecation_warnings(self):
    """Eventually migrate away from dataset_preferences.py mechanism, toward
    100% use of phil for specifying parameters.  For now, simply guard against
    specifying a given parameter by both mechanisms."""

    template = "%s on the command line and %s parameter (%s) of dataset_preferences.py file specify the same thing."

    if self.params==None: return

    # spotarealowcut <==> -s2 <==> minimum_spot_area
    if self.params.distl.minimum_spot_area != None:
      if self.options.find("-s2") >= 0:
        raise Sorry( (template%("minimum_spot_area (%d)","-s2",self.options)%(
        self.params.distl.minimum_spot_area,)) )

    if self.params.distl.minimum_signal_height != None:
      template1 = template%("minimum_signal_height (%.2f)","-bg%1d",self.options)

      # bgupperint <==> -bg0 <==> minimum_signal_height
      if self.options.find("-bg0") >= 0:
        raise Sorry( (template1%(self.params.distl.minimum_signal_height,0)) )

      # bgupperint <==> -bg1 <==> minimum_signal_height
      if self.options.find("-bg1") >= 0:
        raise Sorry( (template1%(self.params.distl.minimum_signal_height,1)) )

      # bgupperint <==> -bg2 <==> minimum_signal_height
      if self.options.find("-bg2") >= 0:
        raise Sorry( (template1%(self.params.distl.minimum_signal_height,2)) )

    # difflowerint <==> -d1 <==> minimum_spot_height
    if self.params.distl.minimum_spot_height != None:
      if self.options.find("-d1") >= 0:
        raise Sorry( (template%("minimum_spot_height (%.2f)","-d1",self.options)%(
        self.params.distl.minimum_spot_height,)) )

    # spotareamaxfactor <==> -s7 <==> spot_area_maximum_factor
    if self.params.distl.spot_area_maximum_factor != None:
      if self.options.find("-s7") >= 0:
        raise Sorry( (template%("spot_area_maximum_factor (%.2f)","-s7",self.options)%(
        self.params.distl.spot_area_maximum_factor,)) )

    # spotpeakintmaxfactor <==> -s8 <==> peak_intensity_maximum_factor
    if self.params.distl.peak_intensity_maximum_factor != None:
      if self.options.find("-s8") >= 0:
        raise Sorry( (template%("peak_intensity_maximum_factor (%.2f)","-s8",self.options)%(
        self.params.distl.peak_intensity_maximum_factor,)) )

    # scanbox_window <==> -bx0,1,2 <==> spot_area_maximum_factor
    if self.options.find("-bx0") >= 0 or self.options.find("-bx1") >= 0 or self.options.find("-bx2") >= 0:
        raise Sorry( """dataset_preferences.py parameters -bx0 -bx1 -bx2 are no longer allowed; found %s.
    Use command line option distl.scanbox_windows=bx0,bx1,bx2 instead with three integer arguments."""%self.options)

class _(boost.python.injector, SpotFilterAgent):
  def __getinitargs__(self):
    return (self.pixel_size, self.xbeam, self.ybeam, self.distance,
            self.wavelength, self.icerings)
