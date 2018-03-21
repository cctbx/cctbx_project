
from __future__ import absolute_import, division, print_function

if __name__ == '__main__':
  import sys
  from dxtbx.datablock import DataBlockFactory
  from dxtbx.model import ParallaxCorrectedPxMmStrategy
  datablocks = DataBlockFactory.from_args(sys.argv[1:])
  assert(len(datablocks) == 1)
  detectors = datablocks[0].unique_detectors()
  assert(len(detectors) == 1)
  detector = detectors[0]
  assert(len(detector) == 1)
  px_mm = detector[0].get_px_mm_strategy()
  assert(isinstance(px_mm, ParallaxCorrectedPxMmStrategy))
  print("Mu: %f mm^-1 " % px_mm.mu())
  print("t0: %f mm" % px_mm.t0())
  from matplotlib import pylab
  from scitbx.array_family import flex
  image_size = detector[0].get_image_size()[::-1]
  xcorr = flex.double(flex.grid(image_size))
  ycorr = flex.double(flex.grid(image_size))
  pixel_size = detector[0].get_pixel_size()
  for j in range(xcorr.all()[0]):
    for i in range(xcorr.all()[1]):
      x1, y1 = detector[0].pixel_to_millimeter((i,j))
      x0, y0 = i * pixel_size[0], j * pixel_size[1]
      xcorr[j,i] = x1 - x0
      ycorr[j,i] = y1 - y0
  vmin = min([flex.min(xcorr), flex.min(ycorr)])
  vmax = max([flex.max(xcorr), flex.max(ycorr)])
  fig, ax = pylab.subplots()
  pylab.subplot(121)
  pylab.imshow(xcorr.as_numpy_array(), interpolation='none',
               vmin=vmin, vmax=vmax)
  pylab.subplot(122)
  im = pylab.imshow(ycorr.as_numpy_array(), interpolation='none',
               vmin=vmin, vmax=vmax)
  fig.subplots_adjust(right=0.8)
  cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
  cbar = fig.colorbar(im, cax=cax)
  pylab.show()
