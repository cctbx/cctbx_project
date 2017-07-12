from __future__ import absolute_import, division
# LIBTBX_SET_DISPATCHER_NAME dxtbx.radial_average
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import libtbx.phil
from libtbx.utils import Usage, Sorry
from libtbx import easy_pickle
import sys
import os
import math

master_phil = libtbx.phil.parse("""
  file_path = None
    .type = str
    .multiple = True
  beam_x = None
    .type = float
  beam_y = None
    .type = float
  n_bins = 0
    .type = int
  verbose = True
    .type = bool
  output_bins = True
    .type = bool
  output_file = None
    .type = str
  plot_x_max = None
    .type = int
  plot_y_max = None
    .type = int
  low_max_two_theta_limit = None
    .type = float
  normalize = False
    .type = bool
  show_plots = True
    .type = bool
  mask = None
    .type = str
""")

def distance (a,b): return math.sqrt((math.pow(b[0]-a[0],2)+math.pow(b[1]-a[1],2)))

def run (args, image = None):
  from xfel import radial_average
  from scitbx.array_family import flex
  import os, sys
  import dxtbx

  # Parse input
  try:
    n = len(args)
  except Exception:
    params = args
  else:
    user_phil = []
    for arg in args :
      if (not "=" in arg) :
        try :
          user_phil.append(libtbx.phil.parse("""file_path=%s""" % arg))
        except ValueError, e :
          raise Sorry("Unrecognized argument '%s'" % arg)
      else :
        try :
          user_phil.append(libtbx.phil.parse(arg))
        except RuntimeError, e :
          raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
    params = master_phil.fetch(sources=user_phil).extract()
  if image is None:
    if params.file_path is None or len(params.file_path) == 0 or not all([os.path.isfile(f) for f in params.file_path]):
      master_phil.show()
      raise Usage("file_path must be defined (either file_path=XXX, or the path alone).")
  assert params.n_bins is not None
  assert params.verbose is not None
  assert params.output_bins is not None

  # Allow writing to a file instead of stdout
  if params.output_file is None:
    logger = sys.stdout
  else:
    logger = open(params.output_file, 'w')
    logger.write("%s "%params.output_file)

  if params.show_plots:
    from matplotlib import pyplot as plt
    import numpy as np
    colormap = plt.cm.gist_ncar
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(params.file_path))])

  if params.mask is not None:
    params.mask = easy_pickle.load(params.mask)

  if image is None:
    iterable = params.file_path
    load_func = lambda x: dxtbx.load(x)
  else:
    iterable = [image]
    load_func = lambda x: x

  # Iterate over each file provided
  for item in iterable:
    img = load_func(item)
    beam = img.get_beam()
    detector = img.get_detector()

    # Search the detector for the panel farthest from the beam. The number of bins in the radial average will be
    # equal to the farthest point from the beam on the detector, in pixels, unless overridden at the command line
    extent = 0
    extent_two_theta = 0
    for panel in detector:
      size2, size1 = panel.get_image_size()
      bc = panel.get_beam_centre_px(beam.get_s0())
      p_extent = int(math.ceil(max(distance((0,0),bc),
                                   distance((size1,0),bc),
                                   distance((0,size2),bc),
                                   distance((size1,size2),bc))))
      p_extent_two_theta = 2*math.asin(beam.get_wavelength()/(2*panel.get_max_resolution_at_corners(beam.get_s0())))*180/math.pi
      if p_extent > extent:
        extent = p_extent
        assert p_extent_two_theta >= extent_two_theta
        extent_two_theta = p_extent_two_theta

    if params.n_bins < extent:
      params.n_bins = extent

    # These arrays will store the radial average info
    sums    = flex.double(params.n_bins) * 0
    sums_sq = flex.double(params.n_bins) * 0
    counts  = flex.int(params.n_bins) * 0

    all_data    = img.get_raw_data()
    if not isinstance(all_data, tuple):
      all_data = (all_data,)

    for tile, (panel, data) in enumerate(zip(detector, all_data)):
      if params.mask is None:
        mask = flex.bool(flex.grid(data.focus()), True)
      else:
        mask = params.mask[tile]

      if hasattr(data,"as_double"):
        data = data.as_double()

      logger.flush()
      if params.verbose:
        logger.write("Average intensity tile %d: %9.3f\n"%(tile, flex.mean(data)))
        logger.write("N bins: %d\n"%params.n_bins)
        logger.flush()

      x1,y1,x2,y2 = 0,0,panel.get_image_size()[1],panel.get_image_size()[0]
      bc = panel.get_beam_centre_px(beam.get_s0())
      bc = int(round(bc[1])), int(round(bc[0]))

      # compute the average
      radial_average(data,mask,bc,sums,sums_sq,counts,panel.get_pixel_size()[0],panel.get_distance(),
                     (x1,y1),(x2,y2))

    # average the results, avoiding division by zero
    results = sums.set_selected(counts <= 0, 0)
    results /= counts.set_selected(counts <= 0, 1).as_double()

    # calculate standard devations
    std_devs = [math.sqrt((sums_sq[i]-sums[i]*results[i])/counts[i])
                if counts[i] > 0 else 0 for i in xrange(len(sums))]

    xvals = flex.double(len(results))
    max_twotheta = float('-inf')
    max_result   = float('-inf')

    for i in xrange(len(results)):
      twotheta = i * extent_two_theta/params.n_bins
      xvals[i] = twotheta

      if params.output_bins and "%.3f"%results[i] != "nan":
       #logger.write("%9.3f %9.3f\n"%     (twotheta,results[i]))        #.xy  format for Rex.cell.
        logger.write("%9.3f %9.3f %9.3f\n"%(twotheta,results[i],std_devs[i])) #.xye format for GSASII
       #logger.write("%.3f %.3f %.3f\n"%(twotheta,results[i],ds[i]))  # include calculated d spacings
      if params.low_max_two_theta_limit is None or twotheta >= params.low_max_two_theta_limit:
        if results[i] > max_result:
          max_twotheta = twotheta
          max_result = results[i]

    #logger.write("Maximum 2theta for %s: %f, value: %f\n"%(file_path, max_twotheta, max_result))

    if params.show_plots:
      if params.plot_x_max is not None:
        results = results.select(xvals <= params.plot_x_max)
        xvals = xvals.select(xvals <= params.plot_x_max)
      if params.normalize:
        plt.plot(xvals.as_numpy_array(),(results/flex.max(results)).as_numpy_array(),'-')
      else:
        plt.plot(xvals.as_numpy_array(),results.as_numpy_array(),'-')
      plt.xlabel("2 theta")
      plt.ylabel("Avg ADUs")
      if params.plot_y_max is not None:
        plt.ylim(0, params.plot_y_max)

  if params.show_plots:
    #plt.legend([os.path.basename(os.path.splitext(f)[0]) for f in params.file_path], ncol=2)
    plt.show()

  return xvals, results

if (__name__ == "__main__") :
  run(sys.argv[1:])
