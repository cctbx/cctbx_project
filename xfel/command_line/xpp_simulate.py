# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME xpp.simulate
#
# $Id: cxi_merge.py 22906 2015-09-15 22:32:08Z phyy-nx $

from __future__ import absolute_import, division, print_function
import iotbx.phil
import sys

master_phil="""
data = None
  .type = str
  .multiple = True
  .help = LCLS experiment, like xppi6113
web {
  user = None
    .type = str
  password = None
    .type = str
}
speedup {
  factor = 100
    .type = float
}
runlimits = None
  .type = ints(2)
enforce80 = False
  .type = bool
  .help = report only on stream 80, FEE spectrometer
enforce81 = False
  .type = bool
  .help = report only on stream 81, FEE spectrometer
"""

#-----------------------------------------------------------------------
def run(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.xpp.simulate import phil_validation
  phil_validation(work_params)
  if ("--help" in args) :
    libtbx.phil.parse(master_phil.show())
    return

  from xfel.xpp.simulate import application
  application(work_params)

if (__name__ == "__main__"):
  result = run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)

# typical usage for experiment LI61:
# xpp.simulate web.user=<LCLSuser> web.password=XXXXXXXX data=xppi6115 runlimits=86,144
