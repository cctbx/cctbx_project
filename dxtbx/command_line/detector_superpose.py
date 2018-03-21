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
from __future__ import print_function
from __future__ import absolute_import
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
from libtbx.test_utils import approx_equal

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
output_experiments = superposed_experiments.json
  .type = str
  .help = Output file name with superposed detector
panel_list = None
  .type = ints
  .help = List of panels to use as reference. Use all if set to None.
apply_at_hierarchy_level = None
  .type = int
  .help = If None, apply shift at panel level.  If not None, apply shift \
          at specified hierarchy level.
fit_target = *corners centers
  .type = choice
  .help = Corners: perform superpose using corners of panels. Centers:\
          perform superpose using centers of panels (requires at least 3\
          panels
repeat_until_converged = True
  .type = bool
  .help = If True, do rounds of fitting until the angle of rotation and \
          magnitude of translation stop changing.
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

    reference_experiments = ExperimentListFactory.from_json_file(params.reference_experiments, check_format=False)
    if len(reference_experiments.detectors()) != 1:
      raise Sorry("Please ensure reference has only 1 detector model")
    reference = reference_experiments.detectors()[0]

    moving_experiments = ExperimentListFactory.from_json_file(params.moving_experiments, check_format=False)
    if len(moving_experiments.detectors()) != 1:
      raise Sorry("Please ensure moving has only 1 detector model")
    moving = moving_experiments.detectors()[0]

    # Get list of panels to compare
    if params.panel_list is None or len(params.panel_list) == 0:
      assert len(reference) == len(moving), "Detectors not same length"
      panel_ids = range(len(reference))
    else:
      max_p_id = max(params.panel_list)
      assert max_p_id < len(reference), "Reference detector must be at least %d panels long given the panel list"%(max_p_id+1)
      assert max_p_id < len(moving), "Moving detector must be at least %d panels long given the panel list"%(max_p_id+1)
      panel_ids = params.panel_list

    if params.fit_target == "centers":
      assert len(panel_ids) >= 3, "When using centers as target for superpose, detector needs at least 3 panels"

    def rmsd_from_centers(a, b):
      assert len(a) == len(b)
      assert len(a)%4 == len(b)%4 == 0
      ca = flex.vec3_double()
      cb = flex.vec3_double()
      for i in xrange(len(a)//4):
        ca.append(a[i:i+4].mean())
        cb.append(b[i:i+4].mean())
      return 1000*math.sqrt((ca-cb).sum_sq()/len(ca))

    cycles = 0
    while True:
      cycles += 1

      # Treat panels as a list of 4 sites (corners) or 1 site (centers) for use with lsq superpose
      reference_sites = flex.vec3_double()
      moving_sites = flex.vec3_double()
      for panel_id in panel_ids:
        for detector, sites in zip([reference, moving], [reference_sites, moving_sites]):
          panel = detector[panel_id]
          size = panel.get_image_size()
          corners = flex.vec3_double([panel.get_pixel_lab_coord(point) for point in [(0,0),(0,size[1]-1),(size[0]-1,size[1]-1),(size[0]-1,0)]])
          if params.fit_target == "corners":
            sites.extend(corners)
          elif params.fit_target == "centers":
            sites.append(corners.mean())

      # Compute super position
      rmsd = 1000*math.sqrt((reference_sites-moving_sites).sum_sq()/len(reference_sites))
      print("RMSD before fit: %.1f microns"%rmsd)
      if params.fit_target =="corners":
        rmsd = rmsd_from_centers(reference_sites, moving_sites)
        print("RMSD of centers before fit: %.1f microns"%rmsd)
      lsq = least_squares_fit(reference_sites, moving_sites)
      rmsd = 1000*math.sqrt((reference_sites-lsq.other_sites_best_fit()).sum_sq()/len(reference_sites))
      print("RMSD of fit: %.1f microns"%rmsd)
      if params.fit_target =="corners":
        rmsd = rmsd_from_centers(reference_sites, lsq.other_sites_best_fit())
        print("RMSD of fit of centers: %.1f microns"%rmsd)
      angle, axis = lsq.r.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)
      print("Axis and angle of rotation: (%.3f, %.3f, %.3f), %.2f degrees"%(axis[0], axis[1], axis[2], angle))
      print("Translation (x, y, z, in microns): (%.3f, %.3f, %.3f)"% (1000 * lsq.t).elems)

      # Apply the shifts
      if params.apply_at_hierarchy_level == None:
        iterable = moving
      else:
        iterable = iterate_detector_at_level(moving.hierarchy(), level = params.apply_at_hierarchy_level)

      for group in iterable:
        fast = col(group.get_fast_axis())
        slow = col(group.get_slow_axis())
        ori = col(group.get_origin())

        group.set_frame(lsq.r * fast, lsq.r * slow, (lsq.r*ori) + lsq.t)

        fast = col(group.get_fast_axis())
        slow = col(group.get_slow_axis())
        ori = col(group.get_origin())

      if not params.repeat_until_converged:
        break

      if approx_equal(angle, 0.0, out=None) and approx_equal((1000*lsq.t).length(), 0.0, out=None):
        print("Converged after", cycles, "cycles")
        break
      else:
        print("Movement not close to zero, repeating fit")
        print()

    from dxtbx.serialize import dump
    dump.experiment_list(moving_experiments, params.output_experiments)

    moved_sites = flex.vec3_double()
    for panel_id in panel_ids:
      panel = moving[panel_id]
      size = panel.get_image_size()
      corners = flex.vec3_double([panel.get_pixel_lab_coord(point) for point in [(0,0),(0,size[1]-1),(size[0]-1,size[1]-1),(size[0]-1,0)]])
      if params.fit_target == "corners":
        moved_sites.extend(corners)
      elif params.fit_target == "centers":
        moved_sites.append(corners.mean())


    # Re-compute RMSD after moving detector components
    rmsd = 1000*math.sqrt((reference_sites-moved_sites).sum_sq()/len(reference_sites))
    print("RMSD of fit after movement: %.1f microns"%rmsd)
    if params.fit_target =="corners":
      rmsd = rmsd_from_centers(reference_sites, moved_sites)
      print("RMSD of fit of centers after movement: %.1f microns"%rmsd)

    if params.panel_list is not None:
      reference_sites = flex.vec3_double()
      moved_sites = flex.vec3_double()
      for panel_id in xrange(len(reference)):
        for detector, sites in zip([reference, moving], [reference_sites, moved_sites]):
          panel = detector[panel_id]
          size = panel.get_image_size()
          corners = flex.vec3_double([panel.get_pixel_lab_coord(point) for point in [(0,0),(0,size[1]-1),(size[0]-1,size[1]-1),(size[0]-1,0)]])
          if params.fit_target == "corners":
            sites.extend(corners)
          elif params.fit_target == "centers":
            sites.append(corners.mean())
      # Re-compute RMSD for full detector after moving detector components
      rmsd = 1000*math.sqrt((reference_sites-moved_sites).sum_sq()/len(reference_sites))
      print("RMSD of whole detector fit after movement: %.1f microns"%rmsd)
      if params.fit_target =="corners":
        rmsd = rmsd_from_centers(reference_sites, moved_sites)
        print("RMSD of whole detector fit of centers after movement: %.1f microns"%rmsd)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
