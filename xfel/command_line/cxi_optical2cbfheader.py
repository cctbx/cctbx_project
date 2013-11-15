from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.optical2cbfheader
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage
from xfel.cftbx.detector.cspad_cbf_tbx import read_optical_metrology_from_flat_file

master_phil = libtbx.phil.parse("""
metrology_file = None
  .type = str
  .help = File with optical metrology information posistioning quadrants and sensors
detector = *CxiDs1 XppDs1
  .type = choice
  .optional = False
  .help = Specifiy CxiDs1 for the CXI Ds1 or DsD detectors which have relative coordinates for each quadrant,
  .help = or XppDs1 for XPP Ds1 detector which specifies absolute positions for each quadrant
plot = False
  .type = bool
  .help = show plots during processing
old_style_diff_path = None
  .type = str
  .help = path to old style metrology directory to compare with the given metrology file
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isfile(arg)) :
      user_phil.append(libtbx.phil.parse("""metrology_file=\"%s\"""" % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  if (params.metrology_file is None) :
    master_phil.show()
    raise Usage("metrology_file must be defined (either metrology_file=XXX, or the file path alone).")
  assert params.detector is not None
  assert params.plot is not None

  print params.metrology_file, params.detector

  try:
    from xfel.cxi.cspad_ana.cspad_tbx import pixel_size as ps
    from xfel.cftbx.detector.cspad_cbf_tbx import sensor_dimension as dim
    pixel_size = ps
    sensor_dimension = dim
  except ImportError:
    pixel_size = 110e-3
    sensor_dimension = ((194*2)+3,185)

  quadrants = read_optical_metrology_from_flat_file(params.metrology_file, params.detector,
                                                    pixel_size,sensor_dimension,
                                                    plot=params.plot,old_style_diff_path=params.old_style_diff_path)
