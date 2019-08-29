# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME xpp.progress_detail
#
# $Id: cxi_merge.py 22906 2015-09-15 22:32:08Z phyy-nx $

from __future__ import absolute_import, division, print_function
import iotbx.phil
import sys

master_phil="""
experiment = None
  .type = str
experiment_tag = None
  .type = str
db {
  host = None
    .type=str
  name = None
    .type = str
  user = None
    .type=str
  password = None
    .type = str
}
trial = None
  .type = int
resolution = 2.5
  .type = float
n_bins = 15
  .type = int
run_tags = None
  .type = str
include_negatives = False
  .type = bool
"""
#xpp.progress_detail db.host= db.user= data=xppi6115 trial=12
#-----------------------------------------------------------------------
def run(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.xpp.progress_utils import phil_validation
  phil_validation(work_params)
  if ("--help" in args) :
    libtbx.phil.parse(master_phil.show())
    return

  from xfel.xpp.progress_utils import application
  application(work_params)

if (__name__ == "__main__"):
  result = run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)

# typical usage for experiment LI61:
# xpp.progress_detail password=XXXXXXXX db.host=psdb db.user=xppi6115 data=xppi6115 trial=12
