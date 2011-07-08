import math
from spotfinder.array_family import flex
from spotfinder.applications.heuristic_tbx.method2_resolution\
          import ResolutionShells
from libtbx.development.timers import Timer

# for identifying rings, the ratio of ring bin population to
#  adjacent non-ring bin population
ring_threshold = 2.0

class RingFinder:
  def __init__(self,spot_resolutions,firstBinCount,fractionCalculator):
    self.ref_spots = spot_resolutions
    self.firstBinCount = firstBinCount
    self.targetResolution = self.ref_spots[self.firstBinCount-1]
    self.fractionCalculator = fractionCalculator
    self.intervals = []
    self.force_finer_bins = False
    self.Shell = ResolutionShells(self.ref_spots,self.targetResolution,
                             self.fractionCalculator)
    while self.finer_bins_go():
      self.targetResolution *= math.pow(0.5,-1./3.) # double the reciprocal volume
      #T = Timer("instance of Resolution shells")
      #  this is an important performance bottleneck for virus work
      #  still needs to be optimized, probably with C++ code

      self.Shell = ResolutionShells(self.ref_spots,self.targetResolution,
                             self.fractionCalculator)
      #self.Shell.show()
      #del T
  def finer_bins_go(self):
    self.find_target_intervals()

    #signal must be high enough to detect
    if max(self.Shell.adjustPop) < ring_threshold:
      return False

    #the first bin should never have less than 1 spot
    elif self.ref_spots[0] <= self.targetResolution: return False

    #return True now if all spots are clustered at low-res end
    elif self.force_finer_bins:
      self.force_finer_bins = False
      return True

    #want at least an average of 4 spots per bin
    elif self.Shell.rows() >= 0.25 * len(self.ref_spots): return False

    else: return True

  def find_target_intervals(self):
    self.intervals = []
    if self.Shell.rows()<10: return # not enough bins
    self.Shell.addColumn('interval','---')
    self.Shell.interval.format = '%20s'

    # Get the signal from the central three bins:
    # ...___***___...
    # where .=ignored, _=background, *=signal, with x centered on middle *
    for x in xrange(2,self.Shell.rows()-1):
      first_bk = x - 4
      last_bk = x + 4
      if first_bk < 0:
        adjust = -first_bk
        first_bk += adjust
        last_bk += adjust
      if last_bk >= self.Shell.rows():
        adjust = last_bk - (self.Shell.rows() - 1)
        first_bk -= adjust
        last_bk -= adjust
      signal = self.Shell.adjustPop[x] + self.Shell.adjustPop[x-1] +\
               self.Shell.adjustPop[x+1]
      nsignal = 3.0

      background = 0.0
      for y in xrange(first_bk,last_bk+1):
        background += self.Shell.adjustPop[y]
      nbackground = 6
      background-=signal

      if (signal/nsignal) > ring_threshold * (background/nbackground):
        # bin population seems to be higher than surroundings, but
        # make sure this is significant based on counting statistics
        # use toy formula background + bkgd_error < signal - signal_error
        # signal minimum of 6 is still a little arbitrary
        if ( signal > 6.0 and
             ((signal-math.sqrt(signal))/nsignal) > ring_threshold *
             ((background + math.sqrt(background))/nbackground) ):
          if x < 4: self.force_finer_bins = True
          self.Shell.interval[x]='***'
          self.add_interval((self.Shell.Limit[x-2],self.Shell.Limit[x+1]))
    #self.Shell.show(['Limit','Population','Fract','interval'])

  def add_interval(self,interval):
    assert interval[0]>interval[1]
    eps = 0.0#0.05 # angstroms
    #consolidate so that end product is a set of disjoint resolution rings
    #self.intervals is a list of rings of form (low res end,high res end)
    for x in xrange(len(self.intervals)):
      existing = self.intervals[x]

      # check for disjointness
      if interval[1] > existing[0]+eps or interval[0] < existing[1]-eps:
        continue

      # check for containment of new interval
      if interval[0] < existing[0]+eps and interval[1] > existing[1]-eps:
        return

      # check for containment of old interval
      if interval[0] > existing[0]+eps and interval[1] < existing[1]-eps:
        self.intervals[x] = interval
        return

      # check for new interval overlap at low resolution end
      if (interval[1] <= existing[0]+eps) and (interval[0] >existing[0]+eps):
        #print existing,interval,"==>",(interval[0],existing[1])
        self.intervals.pop(x)
        self.add_interval( (interval[0],existing[1]) )
        return
      # check for new interval overlap at high resolution end
      if interval[0] >= existing[1]-eps:
        #print existing,interval,"-->",(existing[0],interval[1])
        self.intervals.pop(x)
        self.add_interval( (existing[0],interval[1]) )
        return

    self.intervals.append(interval)

  def filtered(self,sorted_order):
    # recoded in C++
    from spotfinder.core_toolbox import bin_exclusion_and_impact
    assert len(sorted_order)==len(self.ref_spots)

    self.impact = flex.int([0]*len(self.intervals))

    limits_low_hi_res = flex.double();
    for x in xrange(len(self.intervals)):
      limits_low_hi_res.append(self.intervals[x][0]);
      limits_low_hi_res.append(self.intervals[x][1]);

    N = bin_exclusion_and_impact(self.ref_spots,
                                    sorted_order,
                                    limits_low_hi_res,
                                    self.impact)
    return N

  def ice_ring_impact(self):
    #don't want to bother the user with information if the number
    # of filtered spots in the ring is less than 10.
    #
    return len( [x for x in self.impact if x>10] )

  def ice_ring_bounds(self):
    bounds = []
    for i,interval in enumerate(self.intervals):
      if self.impact[i] > 10:
        bounds.append( self.intervals[i] )
    return bounds
