#!/usr/bin/env python
#
# cspad_detector_shifts.py
#
#  Copyright (C) 2016 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_SET_DISPATCHER_NAME dxtbx.detector_superpose
#
from __future__ import division
from scitbx.array_family import flex
from scitbx.matrix import col
from libtbx.phil import parse
from libtbx.utils import Sorry
import libtbx.load_env
import math
from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx.array_family import flex
from scitbx.math.superpose import least_squares_fit
from xfel.command_line.cspad_detector_congruence import iterate_detector_at_level

help_message = '''
This program is used to superpose a moving detector onto a reference detector
  %s reference.json moving.json
''' % libtbx.env.dispatcher_name

# Create the phil parameters
phil_scope = parse('''
reference_experiments = None
  .type = str
  .help = Experiments list with reference detector
moving_experiments = None
  .type = str
  .help = Experiments list with moving detector
panel_list = None
  .type = ints
  .help = List of panels to use as reference. Use all if set to None.
apply_at_hierarchy_level = None
  .type = int
  .help = If None, apply shift at panel level.  If not None, apply shift \
          at specified hierarchy level.
''', process_includes=True)

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s reference.json moving.json " % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      check_format=False,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)

    reference_experiments = ExperimentListFactory.from_json_file(params.reference_experiments)
    if len(reference_experiments.detectors()) != 1:
      raise Sorry("Please ensure reference has only 1 detector model")
    reference = reference_experiments.detectors()[0]

    moving_experiments = ExperimentListFactory.from_json_file(params.moving_experiments)
    if len(moving_experiments.detectors()) != 1:
      raise Sorry("Please ensure moving has only 1 detector model")
    moving = moving_experiments.detectors()[0]

    reference_sites = flex.vec3_double()
    moving_sites = flex.vec3_double()

    # Get list of panels to compare
    if params.panel_list is None or len(params.panel_list) == 0:
      assert len(reference) == len(moving), "Detectors not same length"
      panel_ids = range(len(reference))
    else:
      max_p_id = max(params.panel_list)
      assert max_p_id < len(reference), "Reference detector must be at least %d panels long given the panel list"%(max_p_id+1)
      assert max_p_id < len(moving), "Moving detector must be at least %d panels long given the panel list"%(max_p_id+1)
      panel_ids = params.panel_list

    # Treat panels a a list of 4 sites for use with lsq superpose
    for panel_id in panel_ids:
      for detector, sites in zip([reference, moving], [reference_sites, moving_sites]):
        panel = detector[panel_id]
        size = panel.get_image_size()
        for point in [(0,0),(0,size[1]),(size[0],0),(size[0],size[1])]:
          sites.append(panel.get_pixel_lab_coord(point))

    # Compute super position
    lsq = least_squares_fit(reference_sites, moving_sites)
    rmsd = 1000*math.sqrt((reference_sites-lsq.other_sites_best_fit()).sum_sq()/len(reference_sites))
    print "RMSD of fit: %.1f microns"%rmsd

    # Apply the shifts
    if params.apply_at_hierarchy_level == None:
      iterable = moving
    else:
      iterable = iterate_detector_at_level(moving.hierarchy(), level = params.apply_at_hierarchy_level)

    for group in iterable:
      fast = col(group.get_fast_axis())
      slow = col(group.get_slow_axis())
      ori = col(group.get_origin())

      group.set_frame(lsq.r * fast, lsq.r * slow, ori + lsq.t)

    from dxtbx.serialize import dump
    dump.experiment_list(moving_experiments, "superposed_experiments.json")

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
