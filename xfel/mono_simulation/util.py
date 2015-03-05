from __future__ import division
from scitbx.array_family import flex

def green_curve_area(twotheta, deltaphi):
  # the area under the green curve, in AD14 Figure 7.
  # some caveats:
  #   changing the cell parameters and wavelength will change two theta slightly
  #   changing the spot selection will change the area, since the area
  #   is sampled by the spotfinder spots that are actually used for the fit.
  #   So this is approximate only!

  order = flex.sort_permutation(twotheta)

  ordered_two_theta = twotheta.select(order)
  ordered_deltaphi = deltaphi.select(order)

  area = 0.
  for idx in xrange(len(ordered_deltaphi)-1):
    deltaX = ordered_two_theta[idx+1] - ordered_two_theta[idx]

    averageY = ordered_deltaphi[idx+1] + ordered_deltaphi[idx]

    area_of_trapezoid = deltaX * averageY
    area += area_of_trapezoid
  return area
