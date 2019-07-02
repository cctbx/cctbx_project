# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME xpp.isoform
#

from __future__ import absolute_import, division, print_function
import iotbx.phil
import sys

from xfel.command_line.experiment_manager import master_phil as ex_master_phil
master_phil=ex_master_phil + """

isoform
  .help=Constrain the unit cell to specific values during refinement
  .help=As presently implemented, applies only to dials_refinement_preceding_integration
  .help=and applies only to the higher-symmetry integration trial, not the initial triclinic
  .multiple=False
  {
    name=None
      .type=str
    cell=None
      .type=unit_cell
    lookup_symbol=None
      .type=str
      .help=The sgtbx lookup symbol of the reflections pointgroup
    resolution_limit=10.
      .type = float
      .help = an outer resolution limit considerably beyond anything possible in this experiment
  }
"""

#-----------------------------------------------------------------------
def run(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.xpp.isoform import phil_validation
  phil_validation(work_params)
  if ("--help" in args) :
    libtbx.phil.parse(master_phil.show())
    return

  from xfel.xpp.isoform import application
  application(work_params)

if (__name__ == "__main__"):
  result = run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)

# typical usage for experiment LI61:
# xpp.isoform experiment=xppi6115 experiment_tag=xppi6115 db.user=xxxx isoform.name=A isoform.cell=1,2,3,90,90,90 isoform.lookup_symbol="P m m m"
