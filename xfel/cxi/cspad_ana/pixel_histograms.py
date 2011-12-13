import cspad_tbx
import common_mode
import math
import os

from scitbx.array_family import flex
import scitbx.math
from libtbx import easy_pickle


class pixel_histograms(common_mode.common_mode_correction):

  def __init__(self,
               address,
               pickle_dirname=".",
               pickle_basename="hist",
               roi=None,
               **kwds):
    """

    @param address         Address string XXX Que?!
    @param pickle_dirname     Directory portion of output pickle file
                           XXX mean, mu?
    @param pickle_basename    Filename prefix of output pickle file
                           image XXX mean, mu?
    @param calib_dir       Directory with calibration information
    @param dark_path       Path to input dark image
    @param dark_stddev     Path to input dark standard deviation
    """
    super(pixel_histograms, self).__init__(
      address=address,
      **kwds
    )
    self.pickle_dirname = cspad_tbx.getOptString(pickle_dirname)
    self.pickle_basename = cspad_tbx.getOptString(pickle_basename)
    self.histograms = {}
    self.dimensions = None
    self.roi = cspad_tbx.getOptROI(roi)
    self.values = flex.long()

    self.sigma_scaling = True
    if self.sigma_scaling:
      self.hist_min = -4
      self.hist_max = 8
      self.n_slots = 80
    else:
      self.hist_min = -40
      self.hist_max = 80
      self.n_slots = 80


  def beginjob(self, evt, env):
    print "beginjob!"
    self.hist_mgr = env.hmgr()
    super(pixel_histograms, self).beginjob(evt, env)


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.
    Once self.nshots shots are accumulated, this function turns into
    a nop.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(pixel_histograms, self).event(evt, env)

    if (evt.get("skip_event")):
      return

    if self.sigma_scaling:
      flex_cspad_img = self.cspad_img.as_double()
      flex_cspad_img_sel = flex_cspad_img.as_1d().select(self.dark_mask.as_1d())
      flex_dark_stddev = self.dark_stddev.select(self.dark_mask.as_1d()).as_double()
      assert flex_dark_stddev.count(0) == 0
      flex_cspad_img_sel /= flex_dark_stddev
      flex_cspad_img.as_1d().set_selected(self.dark_mask.as_1d().iselection(), flex_cspad_img_sel)
      self.cspad_img = flex_cspad_img.iround()

    pixels = self.cspad_img.deep_copy()
    dimensions = pixels.all()
    if self.roi is None:
      self.roi = (0, dimensions[1], 0, dimensions[0])

    for i in range(self.roi[2], self.roi[3]):
      for j in range(self.roi[0], self.roi[1]):
        if (i,j) not in self.histograms:
          self.histograms[(i,j)] = flex.histogram(flex.double(), self.hist_min, self.hist_max, self.n_slots)
        self.histograms[(i,j)].update(pixels[i,j])

    self.nmemb += 1
    if math.log(self.nmemb, 2) % 1 == 0:
      self.endjob(env)
    print self.nmemb


  def endjob(self, env):
    """The endjob() function finalises the mean and standard deviation
    images and writes them to disk.

    @param env Environment object
    """

    super(pixel_histograms, self).endjob(env)

    d = {
      "nmemb": self.nmemb,
      "histogram": self.histograms,
    }

    pickle_path = os.path.join(self.pickle_dirname,
                               self.pickle_basename+str(env.subprocess())+".pickle")
    easy_pickle.dump(pickle_path, d)
    self.logger.info(
      "Pickle written to %s" % self.pickle_dirname)

    if (self.nfail == 0):
      self.logger.info(
        "%d images processed" % self.nmemb)
    else:
      self.logger.warn(
        "%d images processed, %d failed" % (self.nmemb, self.nfail))
