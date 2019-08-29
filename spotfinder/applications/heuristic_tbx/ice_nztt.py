from __future__ import absolute_import, division, print_function
from six.moves import range
from spotfinder.array_family import flex

from spotfinder.applications.heuristic_tbx.ice2 import RingFinder
from spotfinder.core_toolbox import hough
from spotfinder.diffraction.geometry import radius_to_resol
#not sure why he uses radius_to_resol if the equi-resolution curves are ellipses

class RingFinder_nztt(RingFinder):
  def __init__(self,spot_resolutions,firstBinCount,fractionCalculator,
               image):
    RingFinder.__init__(self,spot_resolutions,firstBinCount,fractionCalculator)
    # use Hough transform to locate ellipses
    frame = hough()
    frame.importData(image.linearintdata,image.pixel_size)
    frame.setGeometry( float(fractionCalculator.pd['xbeam']),
                       float(fractionCalculator.pd['ybeam']),
                       float(fractionCalculator.pd['distance']),
                       float(fractionCalculator.pd['twotheta']) )
    frame.cannyEdge(2,0.95,0.97)
    frame.findEllipse(2,750.0,3.5)
    self.rings = frame.getRings()
    self.col = 7
    self.nRings = len(self.rings)//self.col

    # add ellipses to existing list of rings
    for i in range(self.nRings):
      res = radius_to_resol(self.rings[i*self.col+6]*image.pixel_size,
                            fractionCalculator.pd)
      RingFinder.add_interval(self,(res+0.015,res-0.015))

    # ensure that rings from hough are definitely ice rings
    self.impact = flex.int([0]*len(self.intervals))
    for i in range(len(self.intervals)):
      for j in range(self.nRings):
        res = radius_to_resol(self.rings[j*self.col+6]*image.pixel_size,
                              fractionCalculator.pd)
        if (res <= self.intervals[i][0] and
            res >= self.intervals[i][1]):
          self.impact[i] = 15

  def filtered(self,sorted_order):
    # recoded in C++
    from spotfinder.core_toolbox import bin_exclusion_and_impact
    assert len(sorted_order)==len(self.ref_spots)

    limits_low_hi_res = flex.double();
    for x in range(len(self.intervals)):
      limits_low_hi_res.append(self.intervals[x][0]);
      limits_low_hi_res.append(self.intervals[x][1]);

    N = bin_exclusion_and_impact(self.ref_spots,
                                    sorted_order,
                                    limits_low_hi_res,
                                    self.impact)
    return N
