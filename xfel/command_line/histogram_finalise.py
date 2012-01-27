# LIBTBX_SET_DISPATCHER_NAME cxi.hist_finalise
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import os
import sys

from libtbx.option_parser import option_parser
from xfel.cxi.cspad_ana.histogram_finalise import histogram_finalise

if __name__ == '__main__':
  import sys
  args = sys.argv[1:]
  assert len(args) > 0
  command_line = (option_parser()
                  .option("--output_dirname", "-o",
                          type="string",
                          help="Directory for output files.")
                  .option("--pickle_pattern",
                          type="string",
                          help="regex for matching pickle files.")
                  ).process(args=args)
  output_dirname = command_line.options.output_dirname
  pickle_pattern = command_line.options.pickle_pattern
  runs = command_line.args
  if output_dirname is None:
    output_dirname = runs[0]
  print "Output directory: %s" %output_dirname
  histogram_finalise(output_dirname, runs, pickle_pattern=pickle_pattern)
