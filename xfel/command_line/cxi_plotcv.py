#
# LIBTBX_SET_DISPATCHER_NAME cxi.plotcv
#
"""Main idea:  go through the integration log files and grep out the difference
   vectors between observed and predicted spot positions.  Plot these "correction
   vectors", but do so separately for each ASIC tile, giving an independent
   check on whether a tile is positioned properly by the tile_translation
   parameters."""

import iotbx.phil
import sys

master_phil="""
run_numbers = None
  .type = ints
  .help = List of run numbers to be aggregated together to make the plots.
outdir_template = None
  .type = str
  .help = Full path directory containing the stdout logs, with %%04d tag for run number
bravais_setting_id = None
  .type = int
  .help = ID number for the Bravais setting of interest (Labelit format).  eg, 1=triclinic, 12=hexagonal
show_plots = False
  .type = bool
  .help = Show graphical plots using matplotlib
"""

#-----------------------------------------------------------------------
def run(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  if ("--help" in args) :
    libtbx.phil.parse(master_phil.show())
    return

  if ((work_params.run_numbers is None) or
      (work_params.outdir_template is None) or
      (work_params.bravais_setting_id is None)) :
    raise Usage("cxi.plotcv "
                "run_numbers=16,17,18,19,20,21,22,23,24,25,26,27,71,72,73 "
                "outdir_template=/reg/data/ana11/cxi/cxi49812/scratch/april_2012/r%%04d/042/stdout "
                "bravais_setting_id=5")
  if work_params.show_plots is True:
    from matplotlib import pyplot as plt # special import

  from xfel.cxi.correction_vector_plot import run_correction_vector_plot
  run_correction_vector_plot(work_params)

  return None

if (__name__ == "__main__"):

  result = run(args=sys.argv[1:])
