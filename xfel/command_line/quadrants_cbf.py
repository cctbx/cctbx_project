from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cspad.quadrants_cbf
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import sys,os
from scitbx.matrix import col
from xfel.cftbx.detector.cspad_cbf_tbx import center
import dxtbx
import libtbx.load_env
from libtbx.utils import Usage


if (__name__ == "__main__"):
  if len(sys.argv) == 1 or '-h' in sys.argv or '--help' in sys.argv or '-c' in sys.argv:
    raise Usage("%s files"%libtbx.env.dispatcher_name)

  files = [arg for arg in sys.argv[1:] if os.path.isfile(arg)]
  arguments = [arg for arg in sys.argv[1:] if not os.path.isfile(arg)]

  for file in files:
    message="""Based on the file %s,
    this program will compute incremental quadrant translations to circularize
    powder rings on the inner four sensors.  The algorithm treats each quadrant
    independently and scores based on self-correlation upon 45-degree rotation.
    """%file
    print message

    image = dxtbx.load(file)
    detector = image.get_detector()
    beam = image.get_beam()

    from xfel.metrology.quadrant import one_panel
    for i_quad, quad in enumerate(detector.hierarchy()):
      # find panel closest to the beam center
      panels = []
      def recursive_get_panels(group):
        if hasattr(group, "children"):
          for child in group:
            recursive_get_panels(child)
        else:
          panels.append(group)
      recursive_get_panels(quad)

      smallest_dist = float("inf")

      for panel in panels:
        p_w, p_h = panel.get_image_size()
        c = center([col(panel.get_pixel_lab_coord((0    ,0  ))),
                    col(panel.get_pixel_lab_coord((p_w-1,0  ))),
                    col(panel.get_pixel_lab_coord((p_w-1,p_h-1))),
                    col(panel.get_pixel_lab_coord((0    ,p_h-1)))])
        beam_center = col(panel.get_beam_centre_lab(beam.get_s0()))

        dist = (c-beam_center).length()
        if dist < smallest_dist:
          smallest_dist = dist
          key_panel = panel

      print "Doing cross-correlation on panel", key_panel.get_name()
      Q = one_panel(image,key_panel,i_quad,quad)
      print Q.coordmax[0], Q.coordmax[1]
      delta = panel.pixel_to_millimeter((Q.coordmax[0], Q.coordmax[1]))

      quad.set_frame(quad.get_fast_axis(),
                     quad.get_slow_axis(),
                     col(quad.get_origin())-col((delta[0],delta[1],0)))

    import pycbf
    image.sync_detector_to_cbf()
    dest_path = os.path.splitext(file)[0]+"_cc.cbf"
    print "Saving result to", dest_path
    image._cbf_handle.write_widefile(dest_path,pycbf.CBF,\
        pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K,0)
