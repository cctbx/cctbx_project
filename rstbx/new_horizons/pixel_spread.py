import math
from scitbx.array_family import flex
from scitbx.matrix import col

switch = 0.12366
upper_cutoff = 0.5 - switch
xupper_cutoff = 0.5 + switch
lower_cutoff = -0.5 + switch # units of pixels
xlower_cutoff = -0.5 - switch # units of pixels
"""Above formulae reflect the line-spread results given in Koerner et al (2009)
Journal of Instrumentation 4: P03001.  Therefore this code is valid for the
CS-PAD, only."""

def line_spread_response_1d(x):
    if x > xupper_cutoff:
      return 0.
    elif x > upper_cutoff:
      return 1.+ math.cos(((x - upper_cutoff)/switch)*math.pi/2.)
    elif x < xlower_cutoff:
      return 0.
    elif x < lower_cutoff:
      return 1.+ math.cos(((x - lower_cutoff)/switch)*math.pi/2.)
    else:
      return 2.
def plot_line_spread_response_1d():
  xs = flex.double()
  ys = flex.double()
  for a in xrange(-15,16,1):
    x = 0.1 * a
    y = line_spread_response_1d(x)
    xs.append(x)
    ys.append(y)
  ys /= flex.sum(ys)
  from matplotlib import pyplot as plt
  plt.plot(xs,ys,"r.")
  plt.show()
def project_2d_response_onto_line(vector):
  xs = flex.double()
  ys = flex.double()
  axis2 = vector.rotate_2d(angle=90.,deg=True)
  #print vector.elems,axis2.elems
  for a in xrange(-15,16,1):
    ax = 0.1 * a
    xs.append(ax)
    sumx = 0.0
    for b in xrange(-15,16,1):
      bx = 0.1 * b
      #convert these rotating-frame coords to stationary frame
      stat = ax*vector + bx*axis2
      #print stat.elems
      sumx += line_spread_response_1d(stat[0]) * line_spread_response_1d(stat[1])
    ys.append(sumx)
  ys /= flex.sum(ys)
  return xs,ys

"""basic idea:  have a collection of bodypixels and intensities. Project each bodypixel
onto the projection direction of interest, either radial or azimuthal.  Convolute
each pixel value intensity with the projection of the line-spread-response function,
to get an observed spot projection.  Take the full-width-half-max of this,
then deconvolute with the point spread, essentially subtract one,
thus giving the FWHM of the diffracted rays, in units of pixels."""

class fwhm_2d_response:

  def __init__(self,rawdata,projection_vector,spotfinder_spot,verbose=False):
      # projection vector is either the radial or azimuthal unit vector
      #   at a specific Bragg spot position
      model_center = col((spotfinder_spot.ctr_mass_x(),spotfinder_spot.ctr_mass_y()))

      px_x,px_y = project_2d_response_onto_line(projection_vector)

      point_projections = flex.double()
      pixel_values = flex.double()
      for point in spotfinder_spot.bodypixels:
        point_projection = (col((point.x,point.y)) - model_center).dot( projection_vector )
        point_projections.append(point_projection)
        pxval = rawdata[(point.x,point.y)]
        if verbose:
          print "point_projection",point_projection,
          print "signal",pxval
        pixel_values.append(  pxval  )
      Lmin = flex.min(point_projections)
      Lmax = flex.max(point_projections)
      #print "Range %6.2f"%(Lmax-Lmin)
      Rmin = round(Lmin-2.0,1)
      Rmax = round(Lmax+2.0,1)
      #print "Range %6.2f"%(Rmax-Rmin)
      def histogram_bin (j) : return int(10.*(j-Rmin)) # bin units of 1/10 pixel

      histo_x = flex.double((int(10*(Rmax-Rmin))))
      histo_y = flex.double(len(histo_x))
      for ihis in xrange(len(histo_x)): histo_x[ihis] = Rmin + 0.1*ihis
      for ipp, point_projection in enumerate(point_projections):
        value = pixel_values[ipp]
        for isample in xrange(len(px_x)):
          histo_y[int(10*(point_projection + px_x[isample] - Rmin))] += value * px_y[isample]
      self.histo_x = histo_x
      self.histo_y = histo_y

  def show_plot(self):
      from matplotlib import pyplot as plt
      plt.plot(self.histo_x,self.histo_y,"r.")
      plt.show()
      del plt

  def fwhm_pix(self):
      half_max = flex.max(self.histo_y) / 2.

      min_idx = 0
      for x in xrange(len(self.histo_y)):
        if self.histo_y[x] > half_max: min_idx = x; break
      max_idx = len(self.histo_y)
      for x in xrange(len(self.histo_y)-1,0,-1):
        if self.histo_y[x] > half_max: max_idx = x; break
      #min_idx and max_idx give FWHM to 0.1 pixel as already constructed
      #but use simple linear interpolation to get FWHM to better than 0.1 pixel.
      lower_bound = min_idx - (self.histo_y[min_idx]-half_max)/(
                          self.histo_y[min_idx] - self.histo_y[min_idx-1])
      upper_bound = max_idx + (self.histo_y[max_idx]-half_max)/(
                          self.histo_y[max_idx] - self.histo_y[max_idx+1])
      potential_fwhm = 0.1 * (upper_bound - lower_bound) - 1.0
      return max(0.,potential_fwhm)

if __name__=="__main__":
  from scitbx.matrix import col
  starting_vec = col((1.0,0.,))
  for angle in xrange(0,360,5):
    f_angle = float(angle)
    projection_vec = starting_vec.rotate_2d(angle=f_angle,deg=True)
    print f_angle,projection_vec.elems
    xs, ys = project_2d_response_onto_line(projection_vec)
    from matplotlib import pyplot as plt
    plt.plot(xs,ys,"r.")
    plt.show()
