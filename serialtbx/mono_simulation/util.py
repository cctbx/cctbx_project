from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex
import math

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
  for idx in range(len(ordered_deltaphi)-1):
    deltaX = ordered_two_theta[idx+1] - ordered_two_theta[idx]

    averageY = ordered_deltaphi[idx+1] + ordered_deltaphi[idx]

    area_of_trapezoid = deltaX * averageY
    area += area_of_trapezoid
  return area

def ewald_proximal_volume(wavelength_ang,resolution_cutoff_ang,domain_size_ang,full_mosaicity_rad):
    """computes the volume of reciprocal space (actually, half the volume, in this implementation) in which
    reciprocal lattice centroids will fall under the green curve.  In other words, this is proportional to the
    number of predicted reflections."""

    R_L = 1./wavelength_ang # radius of Ewald sphere

    # TT is the outermost two-theta angle to perform the volume integration (hi-resolution cutoff)
    TT = 2. * math.asin( wavelength_ang /
                         (2. * resolution_cutoff_ang) )

    part_vol = math.pi * (2./3.) * (1. - math.cos(TT))
    Ewald_sphere_volume = part_vol * math.pow(R_L, 3.) # base volume of Ewald sphere segment
    R_prime = R_L + 1./domain_size_ang
    domain_size_volume = part_vol * math.pow(R_prime, 3.) # expanded volume accomodating spot size

    # compicated integral for mosaic spread volume, must be calculated numerically
    summation = 0.
    N_terms = 100
    for x in range(N_terms):
      phi = (x/N_terms) * TT
      # inner integral over radius r
      integral = math.pow( R_prime + (full_mosaicity_rad * R_L * math.sin(phi)/2.), 3.) - \
                 math.pow( R_prime, 3. )
      summation += (integral * math.sin(phi)) * (TT/N_terms)
    mosaicity_volume = (2./3.) * math.pi * summation

    return (domain_size_volume - Ewald_sphere_volume) + mosaicity_volume
