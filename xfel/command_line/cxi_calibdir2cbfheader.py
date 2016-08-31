from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.calibdir2cbfheader
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage, Sorry
from xfel.cftbx.detector.metrology2phil import metrology2phil
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf

master_phil = libtbx.phil.parse("""
metrology_dir = None
  .type = str
  .help = Directory with optical metrology information posistioning quadrants and sensors, corrected for
  .help = rectangularity
  .optional = False
out = None
  .type = str
  .help = Output file name
  .optional = False
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isdir(arg)) :
      user_phil.append(libtbx.phil.parse("""metrology_dir=\"%s\"""" % arg))
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

  metro_phil = metrology2phil(params.metrology_dir, verbose=True)

  write_cspad_cbf(None, metro_phil, 'calibdir', None, params.out, None, 0, header_only=True)

