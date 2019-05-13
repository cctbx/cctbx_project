#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# detector_congruence.py
#
#  Copyright (C) 2016 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_SET_DISPATCHER_NAME cspad.detector_statistics
#
from __future__ import division
from __future__ import print_function
from six.moves import range
from libtbx.phil import parse
import libtbx.load_env
from libtbx.utils import Usage
from libtbx import easy_run

help_message = '''

This helper program is for looking at all hierarchy levels of the cspad detector after joint hierarchical
refinement using cspad.cbf_metrology

Example:

  %s tag=v1metrology
''' % libtbx.env.dispatcher_name

# Create the phil parameters
phil_scope = parse('''
tag = None
  .type = str
  .help = Used in the plot titles
''')

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser

    # Create the option parser
    usage = "usage: %s tag=tagname" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params

    if params.tag is None:
      raise Usage(self.parser.usage)

    level_json = "%s_%d_refined_experiments_level%d.json"
    level_pickle = "%s_%d_refined_reflections_level%d.pickle"

    command = "cspad.detector_congruence %s %s %s %s hierarchy_level=%d show_plots=False"

    help_strs = [
      "detector as a whole block",
      "quadrants",
      "sensors, I.E. 2x1s",
      "ASICs, I.E. individual tiles"]

    for i in range(3):
      c = command%(level_json%(params.tag, 1, i),
                   level_pickle%(params.tag, 1, i),
                   level_json%(params.tag, 2, i),
                   level_pickle%(params.tag, 2, i),
                   i)

      print("*"*80)
      print("Showing statistics for detector at level %d (%s)"%(i, help_strs[i]))
      print("*"*80)
      print(c)
      result = easy_run.fully_buffered(c).raise_if_errors()
      result.show_stdout()


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
