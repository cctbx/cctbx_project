from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.optical2cbfheader
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage, Sorry
from xfel.cftbx.detector.cspad_cbf_tbx import read_optical_metrology_from_flat_file
from xfel.cftbx.detector.cspad_cbf_tbx import asic_dimension, asic_gap, write_cspad_cbf

master_phil = libtbx.phil.parse("""
metrology_file = None
  .type = str
  .help = File with optical metrology information posistioning quadrants and sensors
  .optional = False
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
out = None
  .type = str
  .help = Output file name
  .optional = False
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isfile(arg)) :
      user_phil.append(libtbx.phil.parse("""metrology_file=\"%s\"""" % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError as e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  if (params.metrology_file is None) :
    master_phil.show()
    raise Usage("metrology_file must be defined (either metrology_file=XXX, or the file path alone).")
  assert params.detector is not None
  assert params.plot is not None
  assert params.out is not None

  print(params.metrology_file, params.detector)

  from xfel.cftbx.detector.cspad_cbf_tbx import pixel_size

  metro = read_optical_metrology_from_flat_file(params.metrology_file, params.detector,
                                                pixel_size,asic_dimension,asic_gap,
                                                plot=params.plot,old_style_diff_path=params.old_style_diff_path)

  write_cspad_cbf(None, metro, 'cbf', None, params.out, None, 0, header_only=True)

