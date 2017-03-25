from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.calibdir2cbfheader
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage, Sorry
from xfel.cftbx.detector.metrology2phil import metrology2phil, sections2phil
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf
from scitbx.matrix import col

master_phil = libtbx.phil.parse("""
metrology_dir = None
  .type = str
  .help = Directory with optical metrology information posistioning quadrants and sensors, corrected for
  .help = rectangularity
  .optional = False
corrections_phil = None
  .type = str
  .help = Phil file with quad/unit pixel translations and optionally subpixel metrology
  .optional = True
out = None
  .type = str
  .help = Output file name
  .optional = False
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isdir(arg)) :
      user_phil.append(libtbx.phil.parse("""metrology_dir=\"%s\" """ % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  if (params.metrology_dir is None) :
    master_phil.show()
    raise Usage("metrology_dir must be defined (either metrology_dir=XXX, or the directory path alone).")
  assert params.out is not None

  if params.corrections_phil is None:
    metro_phil = metrology2phil(params.metrology_dir, verbose=True)
  else:
    from xfel.cxi.cspad_ana.parse_calib import calib2sections
    sections = calib2sections(params.metrology_dir)

    from spotfinder.applications.xfel import cxi_phil
    horizons_phil = cxi_phil.cxi_versioned_extract()
    horizons_phil = horizons_phil.persist.phil_scope.fetch(source=libtbx.phil.parse(file_name=params.corrections_phil))
    corrections_params = horizons_phil.extract()

    bc = [0,0]
    for q in xrange(len(sections)):
      corner = sections[q][1].corners(True)[0]
      bc     = [bc[0] + corner[1] / len(sections), bc[1] + corner[0] / len(sections)]
    beam = col(bc)

    for itile in xrange(len(corrections_params.distl.tile_translations) // 4): # 128 tile_translations/4 = 32 sensors
      # order of quads in sections (UL,UR,LR,LL) is not the same as in quad_translations (UL,UR,LL,LR)
      sections_quad = itile//8
      phil_quad = [0,1,3,2].index(sections_quad)

      s = sections[sections_quad][itile%8] # not the same as iquad!
      c = list(s.center)
      c[0] += corrections_params.distl.quad_translations[2 * phil_quad + 0] + \
              corrections_params.distl.tile_translations[4 * itile + 0]
      c[1] += corrections_params.distl.quad_translations[2 * phil_quad + 1] + \
              corrections_params.distl.tile_translations[4 * itile + 1]
      s.center = tuple(c)
    metro_phil = sections2phil(sections, verbose=True)

  write_cspad_cbf(None, metro_phil, 'calibdir', None, params.out, None, 0, header_only=True)
