# LIBTBX_SET_DISPATCHER_NAME cxi.finalise
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import sys

from libtbx.option_parser import option_parser
from xfel.cxi.cspad_ana.xes_finalise import xes_finalise

if (__name__ == "__main__"):
  args = sys.argv[1:]
  assert len(args) > 0
  command_line = (option_parser()
                  .option("--roi",
                          type="string",
                          help="Region of interest for summing up spectrum.")
                  .option("--output_dirname", "-o",
                          type="string",
                          default=".",
                          help="Directory for output files.")
                  ).process(args=args)
  roi = command_line.options.roi
  output_dirname = command_line.options.output_dirname
  runs = command_line.args
  if not len(runs) > 0:
    print "Usage: cxi.finalise [-o result directory] [data directories as r0xxx/nnn ...]"
  xes_finalise(runs, output_dirname=output_dirname, roi=roi)
