from __future__ import division
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import common_mode
import math
from scitbx.array_family import flex
from xfel import compute_projection
from pyana.event import Event
import numpy as np

class mod_spectra(common_mode.common_mode_correction):

  def __init__(self,
               address,
               target_energy,
               angle=0,
               n_collate = None,
               n_update    = 1,
               threshold= 0.4,
               filter="False",
               clean="False",
               mode="E1",
               peak_ratio=0.17,
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
    self.peak_ratio=cspad_tbx.getOptFloat(peak_ratio)
    self.angle=cspad_tbx.getOptFloat(angle)
    self.nv = 0
    self.collate=None
    self.data=None
    self.ncollate = cspad_tbx.getOptInteger(n_collate)
    self.nvalid   = 0
    self.n_shots  = 0
    self.nupdate  = cspad_tbx.getOptInteger(n_update)
    self.filter=filter
    self.clean=clean
    self.target=map(float,target_energy.split(","))
    self.mode=mode
    self.all=[]
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

  def Filterstdev(self, oneD, target, t):
    import numpy
    stddev=numpy.std(oneD)
    filter=0
    location=0
    for i in target:
      for j in range(len(oneD)):
        if j==i:
          print oneD[j]/stddev
          if (oneD[j]/stddev)>t:
            filter=filter+1
            location=i
    if filter==2: return 1
    if filter==1:
       if location==target[0]: return 2
       if location==target[1]: return 3
    return 0

  def Filter_prev(self, oneD, target, t):
    filter=0
    var=15
    signal=[]
    count=0
    max=0
    location=0
    if len(target)>0:
      for i in target:
        sumS=0
        sumN=0
        for j in range(len(oneD)):
          if j<i-var or j>i+var:
            sumN=sumN+math.fabs(oneD[j]) # sum the noise only
          else:                          # sum signal only
            sumS=sumS+oneD[j]
        if max<sumS:
          max=sumS
          location=i
        signal.append((location,sumS))
    if (max/sumN)>t: filter=filter+1  # check if the maximum peak above the threshold
    if filter==1:              # if the maximum peak above the threshold check if the second peak 25% of the first peak
      for i in signal:
        if i[1]!=max and (i[1]/max)>self.peak_ratio: filter=filter+1
    return filter, location


  def Filter(self, oneD, target, t):
    filter=0
    var=20
    signal=[]
    count=0
    Max=0
    location=0
    if len(target)>0:
        S1=[]
        S2=[]
        N=[]
        for j in range(len(oneD)):
          if j<target[0]-var or (j>target[0]+var and j<target[1]-var) or j>target[1]+var:
            N.append(oneD[j]) # Noise only
          elif j>target[0]-var and j<target[0]+var:
            S1.append(oneD[j]) # Signal first peak only
            if Max<oneD[j]:
              Max=oneD[j]
              location=target[0]
          elif j>target[1]-var and j<target[1]+var:
            S2.append(oneD[j]) # Signal second peak only
            if Max<oneD[j]:
              Max=oneD[j]
              location=target[1]
    import numpy
    stddev=numpy.std(N)
    if (Max/stddev)>t:
       filter=filter+1  # check if the maximum peak above the threshold
#       print Max/stddev
    if filter==1:              # if the maximum peak above the threshold check if the second peak 25% of the first peak
      print Max/max(S1), Max/max(S2)
      if location==target[1]:
        if (max(S1)/Max)>self.peak_ratio: filter=filter+1;print Max/max(S1)
      if location==target[0]:
        if (max(S2)/Max)>self.peak_ratio: filter=filter+1;print Max/max(S2)
#    print filter, location
    return filter, location


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
    self.nv=self.nv+1
    super(mod_spectra, self).event(evt, env)
    if (evt.get("skip_event")):
      return
#    self.nv=self.nv+1
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
        if self.clean=="True":
          oneD=self.Gauss(oneD)
      else:              #if the signal is rotated do projection and sum the signal based on the projection line
       oneD=self.project(pixels)
      if self.clean=="True":
         oneD=self.Gauss(oneD)
      (filter,location)=self.Filter(oneD, self.target, self.threshold)
      if filter==2:
         flag='2C'
      elif filter==1:
        if location==self.target[0]:
           flag='E1'
        if location==self.target[1]:
           flag='E2'
      else:
        flag='bad'
      if self.filter=='False':
          evt.put(flag, "flag")
          evt.put(oneD, "cctbx_spectra")
      if self.filter=="True":
        if self.mode==flag:
          evt.put(flag, "flag")
          evt.put(oneD, "cctbx_spectra")
        else:
          self.logger.warning("event(): Event Skipped")
          evt.put(True, "skip_event")

#         gas=evt.getFeeGasDet()
#         if gas is not None:
#           self.all.append((np.average(gas),oneD))

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """The endjob() function terminates the viewer process by sending
    it a @c None object, and waiting for it to finish.

    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    import pickle
    output = open('dark_fee.pkl', 'wb')
    pickle.dump(self.dark_img, output)
    super(mod_spectra, self).endjob(env)
    self.logger.info("endjob(): end of stream")
