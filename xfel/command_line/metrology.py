#
# LIBTBX_SET_DISPATCHER_NAME cspad.metrology
#

from __future__ import division
from scipy.optimize import leastsq # special import
from xfel.metrology.mark10 import fit_translation4
from libtbx.utils import Sorry
import os.path

master_phil="""
data = None
  .type = path
  .multiple = True
  .help = Directory containing integrated data in pickle format.
  .help = Repeat to specify additional directories, as many as desired.
  .help = Single entry can be a valid unix glob like /pathto/run000[3-8]/integration
bravais_setting_id = None
  .type = int
  .help = ID number for the Bravais setting of interest (Labelit format).
  .help = triclinic=1 monoclinic=2 orthorhombic=5 rhombohedral=5
  .help = tetragonal=9 hexagonal=12 cubic=22
max_frames = None
  .type = int(value_min=2)
  .help = From all the input images, only use the first max_frames for metrology analysis.
show_plots = False
  .type = bool
  .help = Show graphical plots using matplotlib (expert only)
show_consistency = False
  .type = bool
  .help = Run the consistency controls (expert only)
effective_tile_boundaries = None
  .type = ints
  .help = effective integer tile boundaries applied to convert xtc stream to pickled image files. Must be 64 * 4 integers.
"""

def get_phil(args):
  print args

  import iotbx.phil
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()

  if work_params.show_plots is True:
    from matplotlib import pyplot as plt # special import
  return work_params

def validate_phil(wp):
  if wp.data == []:
    raise Sorry("need one or more data directories")

  #handle globs at this level:
  import glob
  provisional = []
  for glob_item in wp.data:
    provisional = provisional + glob.glob(glob_item)
  wp.data = provisional

  for item in wp.data:
    if not os.path.isdir(item):
      raise Sorry("can't stat directory",item)

  if wp.bravais_setting_id is None or \
     wp.bravais_setting_id not in [1,2,5,9,12,22]:
    raise Sorry("""must set the bravais_setting_id
  triclinic=1 monoclinic=2 orthorhombic=5 rhombohedral=5
  tetragonal=9 hexagonal=12 cubic=22""")


def run(args):

  work_params = get_phil(args)
  try:
    validate_phil(work_params)
  except Sorry, s:
    raise s
  except Exception, s:
    raise s

  C = fit_translation4(work_params)
  print "done constructor"

  print
  C.run_cycle_a()
  C.run_cycle_b(0)
  print
  C.run_cycle_a()
  C.run_cycle_b(1)
  print
  C.run_cycle_a()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

#cspad.metrology data=/reg/d/psdm/CXI/cxid9114/ftc/sauter/results/r010[2-9]/018/integration max_frames=300 bravais_setting_id=9
#[nksauter@psexport02][~/phenix/src/cctbx_project/xfel/metrology]$ diff ../../../cctbx_project/xfel/mono_simulation/parameters.h ../../../scale/parameters.h
#12c12
#< #include <xfel/mono_simulation/vector_collection.h>
#---
#> #include <scale/vector_collection.h>
#suggests that the code I need may already be in mono_simulation
