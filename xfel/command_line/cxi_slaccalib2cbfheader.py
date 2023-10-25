from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.slaccalib2cbfheader
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage, Sorry
from serialtbx.detector.cspad import read_slac_metrology
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf

master_phil = libtbx.phil.parse("""
metrology_file = None
  .type = str
  .help = File with optical metrology information posistioning quadrants and sensors.
  .help = Usually in the calib/geometry folder of the experiment, in the form of N-end.data
  .optional = False
plot = False
  .type = bool
  .help = show plots during processing
out = None
  .type = str
  .help = Output file name
  .optional = False
""")

def run(args):
  if ("--help" in args or "-h" in args) :
    print("Write a CBF header from a SLAC metrology file. Parameters:")
    master_phil.show(attributes_level=2)
    return

  user_phil = []
  for arg in args:
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
  assert params.plot is not None
  assert params.out is not None

  print(params.metrology_file)

  metro = read_slac_metrology(params.metrology_file, plot=params.plot)

  write_cspad_cbf(None, metro, 'cbf', None, params.out, None, 0, header_only=True)

if (__name__ == "__main__") :
  args = sys.argv[1:]
  run(args)
