from __future__ import division
from __future__ import print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME xpp.beamcenter
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import sys,os
from scitbx.array_family import flex
from scitbx.matrix import sqr, col
from math import sin, cos, pi
import dxtbx
from xfel.metrology.legacy_scale import quadrant_self_correlation
import libtbx.phil

master_phil = libtbx.phil.parse("""
  beam_center_fast = None
    .type = int
    .help = Initial estimate of beam center (fast coordinate)
  beam_center_slow = None
    .type = int
    .help = Initial estimate of beam center (slow coordinate)
  px_max = None
    .type = int
    .help = Only test pixels this distance or less from the initial beam center estimate
  px_min = None
    .type = int
    .help = Only test pixels this distance or more from the initial beam center estimate
  show_plots = False
    .type = bool
    .help = If True, show the pixels that will be tested and the rotational self-correlations from the grid search
""")

if (__name__ == "__main__"):
  files = [arg for arg in sys.argv[1:] if os.path.isfile(arg)]
  arguments = [libtbx.phil.parse(arg) for arg in sys.argv[1:] if not os.path.isfile(arg)]
  params = master_phil.fetch(sources=arguments).extract()

  rot45 = sqr((sin(pi/4.),-cos(pi/4.),cos(pi/4.),sin(pi/4.)))

  for file in files:
    message="""Based on the file %s, this program will compute incremental translations to circularize
    powder rings.  The algorithm  scores based on self-correlation upon 45-degree rotation.
    Increments are determined ON TOP OF beam center value in image header.  Output is given in the form of
    delta to that value."""%file
    print(message)

    img = dxtbx.load(file)
    beam = img.get_beam()
    s0 = beam.get_s0()

    raw_data = img.get_raw_data()
    if not isinstance(raw_data, tuple):
      raw_data = (raw_data,)

    for panel_id, panel in enumerate(img.get_detector()):
      beam_center = col(panel.get_beam_centre_px(s0))
      data = raw_data[panel_id]

      print("Assembling mask...", end=' '); sys.stdout.flush()
      mask = panel.get_trusted_range_mask(data)
      trusted_min = panel.get_trusted_range()[0]

      mask_center = col((params.beam_center_slow,params.beam_center_fast))
      px_max = params.px_max
      px_min = params.px_min
      data = data[mask_center[0] - px_max:mask_center[0] + px_max, mask_center[1] - px_max:mask_center[1] + px_max]
      mask = mask[mask_center[0] - px_max:mask_center[0] + px_max, mask_center[1] - px_max:mask_center[1] + px_max]
      panel_origin = col((mask_center[0] - px_max,mask_center[1] - px_max))

      for y in range(mask.focus()[1]):
        for x in range(mask.focus()[0]):
          l = (col((x-px_max+mask_center[0],y-px_max+mask_center[1])) - mask_center).length()
          if l < px_min or l > px_max:
            mask[x,y] = False

      data.set_selected(~mask, trusted_min-1)
      print("done")
      if params.show_plots:
        from matplotlib import pyplot as plt
        plt.imshow(data.as_numpy_array())
        plt.show()

      grid_radius = 20
      mapp = flex.double(flex.grid(2*grid_radius+1, 2*grid_radius+1))
      print(mapp.focus())

      gmax = 0.0
      coordmax = (0,0)
      for xi in range(-grid_radius, grid_radius+1):
        for yi in range(-grid_radius, grid_radius+1):
          test_bc = beam_center + col((xi,yi))
          print("Testing beam center", test_bc.elems, end=' ')
          REF,ROT = quadrant_self_correlation(data,panel_origin,test_bc,rot45,trusted_min)
          CCRR = flex.linear_correlation(REF,ROT)
          VV = CCRR.coefficient()
          if VV>gmax:
            gmax = VV
            coordmax = col((xi,yi))
          mapp[(xi+grid_radius,yi+grid_radius)]=VV
          print(VV)

      print("max cc %7.4F is at "%gmax, (beam_center + coordmax).elems, "(slow, fast). Delta:", coordmax.elems)
      if params.show_plots:
        npy = mapp.as_numpy_array()
        from matplotlib import pyplot as plt
        plt.imshow(npy, cmap="hot")
        plt.plot([coordmax[1]+grid_radius],[coordmax[0]+grid_radius],"k.")
        plt.show()
