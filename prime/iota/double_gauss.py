from __future__ import division
from libtbx import easy_pickle as ep
import numpy as np
from scipy.misc import toimage
import scipy.ndimage as ndi
import matplotlib.pyplot as plt

def get_img_data(img_path):
  d = ep.load(img_path)
  return  d['DATA'].as_numpy_array()

def to_im(im):
    return toimage(im,high=20763)

def get_spots(blured_img, background, thresh):
    return np.where(blured_img > thresh * background)

def double_gauss_spotfinder(img_data, g1=2, g2=10, thresh=1.5):
    """ spotfinder prototype using gaussian local background estimation, and gaussian smoothing.
    :param img_data: scipy 2d array of pixel intensities
    :param g1: sigma of backgorund for smoothing
    :param g2: sigma of gaussian for background estimation.
    :param thresh: pixel in g2-convoluted image must be at least thresh times pixel in g1-convoluted image
    :return: splot coordinates, the covoluted smoothed image, and the convoluted background image.
    """
    blured_img = ndi.filters.gaussian_filter(img_data, g1)
    background = ndi.filters.gaussian_filter(img_data, g2)
    spots = get_spots(blured_img, background, thresh)
    print "{} spots".format(len(spots[0]))
    return spots, blured_img, background


class SpotManager:

  def __init__(self, img_data, init_thresh=1.5, init_bg=10, init_smooth=1.5,
               thresh_range=(1, 4), bg_range=(0, 100), smooth_range=(0, 15)):

    self.img_data = img_data

    # Thresholds
    self.spot_thresh = init_thresh
    self.bg_thresh = init_bg
    self.smooth_thresh = init_smooth

    # Slider ranges
    self.thresh_range = thresh_range
    self.bg_range = bg_range
    self.smooth_range = smooth_range

    # Convoluted images
    self.bg_img = None  # Create
    self.update_bg(self.bg_thresh)  # Update
    self.smooth_img = None  # Create
    self.update_smoothing(self.smooth_thresh)  # Update


  def update_bg(self, bg_thresh):
    """ update the state of `self.bg_thresh` and `self.bg_img` """
    self.bg_thresh = bg_thresh
    self.bg_img = ndi.filters.gaussian_filter(self.img_data, self.bg_thresh)

  def update_smoothing(self, smooth_thresh):
    """ update the state of `self.smooth_thresh` and `self.smooth_img` """
    self.smooth_thresh = smooth_thresh
    self.smooth_img = ndi.filters.gaussian_filter(self.img_data,
                                                  self.smooth_thresh)

  def calc_spots(self):
    """ Return spots as a list of coordinates using current state """
    return get_spots(self.smooth_img, self.bg_img, self.spot_thresh)

  def plot_spots(self):
      """ Create an interactive plot using the data produced by `double_gauss_spotfinder`"""

      from matplotlib.widgets import Slider, Button

      def update_ax_prop(ax):
          """ Convenience function to make axes full and square"""
          ax.set_xlim(0, len(self.img_data[1]))
          ax.set_ylim(0, len(self.img_data[0]))
          ax.set_aspect('equal')


      # 1. Plot image to 1st set of axes
      fig = plt.figure(figsize=(8, 8))
      ax = fig.add_subplot(1,1,1)
      update_ax_prop(ax)
      plt.subplots_adjust(bottom=0.25, left=0.2)
      ax.imshow(to_im(self.img_data), origin='lower')

      # 2. Create new axes on top of old ones for plotting to
      ax2 = fig.add_axes(ax.get_position(), frameon=False)

      # 3. Create axes for the slidersi and button:
      axthresh = plt.axes([0.25, 0.05, 0.6, 0.03])
      axsmooth = plt.axes([0.25, 0.10, 0.6, 0.03])
      axbg     = plt.axes([0.25, 0.15, 0.6, 0.03])
      axbut    = plt.axes([0.05, 0.45, 0.1, 0.1])

      # Create sliders and button
      sthresh = Slider(axthresh, 'Threshold',
                       self.thresh_range[0], self.thresh_range[1],
                       valinit=self.spot_thresh, dragging=True)
      ssmooth = Slider(axsmooth, 'Peak smoothing',
                       self.smooth_range[0], self.smooth_range[1],
                       valinit=self.smooth_thresh, dragging=False)
      sbg = Slider(axbg, 'Background smoothing',
                   self.bg_range[0], self.bg_range[1],
                   valinit=self.bg_thresh, dragging=False)
      button = Button(axbut, 'Write\nparams', color='0.95', hovercolor='0.4')


      # Define behaviour when they are changed
      def update_thresh(val):
          """ Update on move of slider """
          ax2.clear()
          update_ax_prop(ax2)
          self.spot_thresh = sthresh.val
          spots = self.calc_spots()
          ax2.scatter(spots[1], spots[0], s=10, c='r', marker='o', linewidths=0, alpha=1)
          plt.draw()
      sthresh.on_changed(update_thresh)

      def update_smooth(val):
        self.update_smoothing(val)
        update_thresh(self.spot_thresh)
      ssmooth.on_changed(update_smooth)

      def update_bg(val):
        self.update_bg(val)
        update_thresh(self.spot_thresh)
      sbg.on_changed(update_bg)

      def write_params(event):
        print "Background smoothing: {:.2f}\nSpot smoothing: {:.2f}\nSpot threshold: {:.2f}" \
              .format(self.bg_thresh, self.smooth_thresh, self.spot_thresh)
      button.on_clicked(write_params)


      # Manually plot starting values
      update_thresh(self.spot_thresh)
      plt.show()

if __name__ == '__main__':
  import sys
  img_data = get_img_data(sys.argv[1])
  spot_man  = SpotManager(img_data)
  spot_man.plot_spots()

################################################################################
################################### Deprecated #################################
################################################################################


def plot_spots_simple(img_data, blured_img, background, init_thresh=1.5,
    thresh_range=(1, 4)):
    """ Create an interactive plot using the data produced by `double_gauss_spotfinder`"""

    print "warning, deprecated!!"

    from matplotlib.widgets import Slider

    def update_ax_prop(ax):
        """ Convenience function to make axes full and square"""
        ax.set_xlim(0, len(img_data[1]))
        ax.set_ylim(0, len(img_data[0]))
        ax.set_aspect('equal')

    def update(val):
        """ Update on move of slider """
        ax2.clear()
        update_ax_prop(ax2)
        thresh = sthresh.val
        spots = get_spots(blured_img, background, thresh)
        ax2.scatter(spots[1], spots[0], s=10, c='r', marker='o', linewidths=0, alpha=1)
        plt.draw()

    # 1. Plot image to 1st set of axes
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    update_ax_prop(ax)
    plt.subplots_adjust(bottom=0.15)
    ax.imshow(to_im(img_data), origin='lower')

    # 2. Create new axes on top of old ones for plotting to
    ax2 = fig.add_axes(ax.get_position(), frameon=False)

    # 3. Create axes for the sliders:
    axthresh = plt.axes([0.25, 0.05, 0.65, 0.03])
    sthresh = Slider(axthresh, 'Threshold', thresh_range[0], thresh_range[1], valinit=init_thresh)
    sthresh.on_changed(update)

    # Manually plot starting values
    update(init_thresh)
    plt.show()

def analyze_threshold(img, bg):
  thresholds = np.arange(3, 0.5, -0.05)
  n_spots = [len(get_spots(img, bg, thresh)[0]) / (len(img) ** 2) for thresh in thresholds]
  plt.plot(thresholds, n_spots)
  plt.semilogy()
  plt.grid()
  plt.ylabel("fraction of pixels counted as spots")
  plt.xlabel("Spot threshold")
  plt.show()

