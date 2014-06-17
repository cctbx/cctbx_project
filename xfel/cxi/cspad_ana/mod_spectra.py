from __future__ import division
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import common_mode
import math
from scitbx.array_family import flex
from xfel import compute_projection
from pyana.event import Event

class mod_spectra(common_mode.common_mode_correction):

  def __init__(self,
               address,
               angle=0,
               n_collate = None,
               n_update    = 1,
               store = "spectrum",
               threshold= 0.4,
               mode="E1",
               filter="False",
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
    self.threshold=cspad_tbx.getOptFloat(threshold)
    self.angle=cspad_tbx.getOptFloat(angle)
    self.nv = 0
    self.collate=None
    self.data=None
    self.store=store
    self.ncollate = cspad_tbx.getOptInteger(n_collate)
    self.nvalid   = 0
    self.n_shots  = 0
    self.nupdate  = cspad_tbx.getOptInteger(n_update)
    self.mode=mode
    self.filter=filter
    if (self.ncollate is None):
      self.ncollate = self.nupdate
    if (self.ncollate > self.nupdate):
      self.ncollate = self.nupdate
      self.logger.warning("n_collate capped to %d" % self.nupdate)

  def beginjob(self, evt, env):
    super(mod_spectra, self).beginjob(evt, env)
    if self.angle != 0:
      self.cos0=math.cos(math.radians(self.angle))   # calculate the sin and cos in the begin job to save computer time
      self.sin0=math.sin(math.radians(self.angle))
      self.limit=int(self.cos0*1024+self.sin0*1024) # get the dimesntion of the array (because the signal is rotated you need larger array)

  def project(self, pixels):
    sum=compute_projection(pixels, self.limit, self.sin0, self.cos0)
    return sum
  def Gauss(self, oneD):
        from scipy.signal import gaussian
        from scipy.ndimage import filters
        b = gaussian(10, 5, sym=False)
        ga = filters.convolve1d(oneD, b/b.sum())
        return ga

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
    self.n_shots += 1
    self.logger.info("processing event Number: " + str(self.nv))
    next_update = (self.nupdate - 1) - (self.nshots - 1) % self.nupdate
    if (self.ncollate > 0 and next_update >= self.ncollate):
      return
    count=0
    file=open("spectra.txt","w+")
    if (self.nvalid == 0 or self.ncollate > 0 and self.nvalid >= self.ncollate):
      self.img_sum = self.cspad_img
      self.nvalid  = 1
    else:
      self.img_sum += self.cspad_img
      self.nvalid  += 1
    if next_update==0:
      pixels=self.img_sum/self.nvalid
      if self.angle==0:  #if the signal is not rotated sum each raw and add it to oneD array
        oneD=[]
        pixels=pixels.matrix_transpose()
        for r in flex.rows(pixels):
          oneD.append(flex.sum(r))
      else:              #if the signal is rotated do projection and sum the signal based on the projection line
       oneD=self.project(pixels)
       oneD=self.Gauss(oneD)
      if self.filter=="False":
#         import numpy
#         print numpy.std(oneD[350:650])
         evt.put(oneD, "cctbx_spectra")
      if self.filter=="True":
        filter=0
        location=0
        track=0  # track if the event should be skipped
        flag=0   # flag to tell the viewer if the shot is E1,E2 or 2 Color
        for i in range(50,250):
          if (oneD[i]>self.threshold):
            filter=filter+1
            location=i
            break
        for i in range(675,875):
          if (oneD[i]>self.threshold):
            filter=filter+1
            location=i
            break
        if filter>1 and (self.mode=="2C" or self.mode=="all"):
          evt.put(oneD, "cctbx_spectra")
          track=1
          flag=1
          self.logger.warning("event(): Two color shot")
        if filter==1:
          if location<500 and (self.mode=="E1" or self.mode=="all"):
            evt.put(oneD, "cctbx_spectra")
            track=1
            flag=2
            self.logger.warning("event(): E1 shot")
          if location>500 and (self.mode=="E2" or self.mode=="all"):
            evt.put(oneD, "cctbx_spectra")
            track=1
            flag=3
            self.logger.warning("event(): E2 shot")
        if track!=0: evt.put(flag,"flag")
        if track==0:
            self.logger.warning("event(): Event Skipped")
            evt.put(True, "skip_event")


  def endjob(self, env):
    """The endjob() function terminates the viewer process by sending
    it a @c None object, and waiting for it to finish.

    @param env Environment object
    """

    super(mod_spectra, self).endjob(env)
    self.logger.info("endjob(): end of stream")
