from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME cspad.quadrants_cbf
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import sys,os
from scitbx.matrix import col
from xfel.cftbx.detector.cspad_cbf_tbx import center
import dxtbx
import libtbx.load_env
from libtbx.utils import Usage
from scitbx.array_family import flex
from libtbx.phil import parse
from libtbx.utils import Sorry

phil_scope = parse("""
  show_plots = False
    .type = bool
    .help = Show CC gridmap plots
  multi_angle = True
    .type = bool
    .help = If true, compute CC over many angles at each gridpoint (20-70 degrees \
            in increments of 2.5 degrees). Otherwise, just compute CC at 45 \
            degrees at each grid point (much faster but not as robust)
  plot_range = None
    .type = floats(size=2)
    .help = Min and max CC values for gridmap plots
  pdf_file = None
    .type = path
    .help = If not None and show_plots is True, then save CC plots as multi- \
            page cbf file
  save_cbf = True
    .type = bool
    .help = If True, write cbf with best corrections applied
""")

def run(args):
  if len(args) == 0 or '-h' in args or '--help' in args or '-c' in args:
    print("Usage: %s [-p] files"%libtbx.env.dispatcher_name)
    phil_scope.show(attributes_level = 2)
    return

  files = [arg for arg in args if os.path.isfile(arg)]
  arguments = [arg for arg in args if not os.path.isfile(arg)]
  user_phil = []
  for arg in arguments:
    if arg == '-p':
      user_phil.append(parse("show_plots=True"))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)
  params = phil_scope.fetch(sources = user_phil).extract()

  for file in files:
    message="""Based on the file %s,
    this program will compute incremental quadrant translations to circularize
    powder rings on the inner four sensors.  The algorithm treats each quadrant
    independently and scores based on self-correlation upon 45-degree rotation.
    """%file
    print(message)

    image = dxtbx.load(file)
    detector = image.get_detector()
    beam = image.get_beam()

    from xfel.metrology.quadrant import one_panel
    ccs = flex.double()
    for i_quad, quad in enumerate(detector.hierarchy()):
      # find panel closest to the beam center
      panels = []
      def recursive_get_panels(group):
        if group.is_group():
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

      print("Doing cross-correlation on panel", key_panel.get_name())
      Q = one_panel(image,key_panel,i_quad,quad,params.show_plots,params.multi_angle,params.plot_range,params.pdf_file is None)
      delta = panel.pixel_to_millimeter((Q.coordmax[0], Q.coordmax[1]))

      quad.set_frame(quad.get_fast_axis(),
                     quad.get_slow_axis(),
                     col(quad.get_origin())-col((delta[0],delta[1],0)))
      ccs.append(Q.ccmax)
    print("Average CC: %7.4f"%flex.mean(ccs))
    if params.pdf_file is not None:
      print("Saving plots to", params.pdf_file)
      from matplotlib import pyplot as plt
      from matplotlib.backends.backend_pdf import PdfPages
      pp = PdfPages(params.pdf_file)
      for i in plt.get_fignums():
        pp.savefig(plt.figure(i), dpi=300)
      pp.close()

    if params.save_cbf:
      import pycbf
      image.sync_detector_to_cbf()
      dest_path = os.path.splitext(file)[0]+"_cc.cbf"
      print("Saving result to", dest_path)
      image._cbf_handle.write_widefile(dest_path,pycbf.CBF,\
          pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K,0)

if (__name__ == "__main__"):
  run(sys.argv[1:])
