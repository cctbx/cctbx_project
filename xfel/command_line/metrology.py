#
# LIBTBX_SET_DISPATCHER_NAME cspad.metrology
#

from __future__ import absolute_import, division, print_function
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
min_count = 25
  .type = int(value_min=3)
  .help = On each sensor, require this minimum number of spots before proposing unit pixel translations.
  .help = On the CSPAD, 1 sensor is a pair of 2 ASICs.
detector_format_version = None
  .type = str
  .help = CSPAD format version, as defined in spotfinder/applications/xfel/cxi_phil.py
  .help = Case I: If given, program will print out old and new unit pixel translations, to be pasted back into code.
  .help = The new replaces the old.
  .help = Case II: If not given, program will print increments to the whatever translations were used for integration,
  .help = also to be pasted back into the code.  The new increments the old, and this must be explicitly coded.
show_plots = False
  .type = bool
  .help = Show graphical plots using matplotlib (expert only)
show_consistency = False
  .type = bool
  .help = Run the consistency controls (expert only)
effective_tile_boundaries = None
  .type = ints
  .help = Effective integer tile boundaries applied to convert xtc stream to pickled image files. Must be 64 * 4 integers.
  .help = Boundaries should not be under user control, they are provided by the input data (pickle files).
diff_cutoff = 5
  .type = float
  .help = throw out image if any correction vector exceeds this length in pixels
"""

def get_phil(args):
  print(args)

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
  except Sorry as s:
    raise s
  except Exception as s:
    raise s

  C = fit_translation4(work_params)
  print("done constructor")

  #parse a phil file from the args
  for item in args:
    if os.path.isfile(item):
      from xfel.phil_preferences import load_cxi_phil
      C.optional_params = load_cxi_phil(item,[])
      break

  print()
  C.run_cycle_a()
  C.run_cycle_b(0)
  print()
  C.run_cycle_a()
  C.run_cycle_b(1)
  print()
  C.run_cycle_a()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

#cspad.metrology data=/reg/d/psdm/CXI/cxid9114/ftc/sauter/results/r010[2-9]/018/integration max_frames=300 bravais_setting_id=9
#cspad.metrology data=/reg/d/psdm/CXI/cxid9114/ftc/brewster/results/r00[3-4]*/003/integration max_frames=300 bravais_setting_id=9
#cspad.metrology data=/reg/d/psdm/CXI/cxid9114/ftc/brewster/results/r00[3-4]*/003/integration max_frames=300 bravais_setting_id=9 min_count=20
