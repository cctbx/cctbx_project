from __future__ import division
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import common_mode
import math
from scitbx.array_family import flex
from xfel import compute_projection

class mod_spectra(common_mode.common_mode_correction):

  def __init__(self,
               address,
               angle=0,
               common_mode_correction = "none",
               **kwds):

    """

    @param address         Address string XXX Que?!
    @param dirname         Directory portion of output pickle file
    @param basename        Filename prefix of output pickle file
    @param dark_path       Path to input dark image
    @param dark_stddev     Path to input dark standard deviation
    """
    super(mod_spectra, self).__init__(
      address=address,
      common_mode_correction=common_mode_correction,
      **kwds
    )

    self.angle=cspad_tbx.getOptFloat(angle)
    self.nv = 0

  def beginjob(self, evt, env):
    super(mod_spectra, self).beginjob(evt, env)
    if self.angle != 0:
       self.cos0=math.cos(math.radians(self.angle))   # calculate the sin and cos in the begin job to save computer time
       self.sin0=math.sin(math.radians(self.angle))
       self.limit=int(self.cos0*1024+self.sin0*1024) # get the dimesntion of the array (because the signal is rotated you need larger array)

  def project(self, pixels):
    sum=compute_projection(pixels, self.limit, self.sin0, self.cos0)
    return sum

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.
    Once self.nshots shots are accumulated, this function turns into
    a nop.
    @param evt Event data object, a configure object
    @param env Environment object
    """
    super(mod_spectra, self).event(evt, env)

    if (evt.get("skip_event")):
      return
    self.nv=self.nv+1
    self.logger.info("processing event Number: " + str(self.nv))
    pixels = self.cspad_img.deep_copy()
    if self.angle==0:  #if the signal is not rotated sum each raw and add it to oneD array
      oneD=[]
      for r in flex.rows(pixels):
        oneD.append(flex.sum(r))
    else:              #if the signal is rotated do projection and sum the signal based on the projection line
      oneD=self.project(pixels)
      count=0
    evt.put(oneD, "cctbx_spectra")


  def endjob(self, env):
    """The endjob() function terminates the viewer process by sending
    it a @c None object, and waiting for it to finish.

    @param env Environment object
    """

    super(mod_spectra, self).endjob(env)
    self.logger.info("endjob(): end of stream")

