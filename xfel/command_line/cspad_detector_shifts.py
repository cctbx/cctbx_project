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
# LIBTBX_SET_DISPATCHER_NAME cspad.detector_shifts
#
from __future__ import division
from scitbx.array_family import flex
from scitbx.matrix import col
from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.command_line.cspad_detector_congruence import get_center
import libtbx.load_env

help_message = '''

This program is used to show differences between a reference and a moving set of detectors
Example:

  %s experiment1.json experiment2.json reflections1.pickle reflections2.pickle
''' % libtbx.env.dispatcher_name

# Create the phil parameters
phil_scope = parse('''
max_hierarchy_level=Auto
  .type = int
  .help = Maximum hierarchy level to compute shifts to
''', process_includes=True)

from dials_scratch.asb.detector_congruence import iterate_detector_at_level, iterate_panels, id_from_name
from dials_scratch.asb.detector_congruence import Script as ParentScript

class Script(ParentScript):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s experiment1.json experiment2.json reflections1.pickle reflections2.pickle" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_datablocks=True,
      read_reflections=True,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_datablocks, flatten_reflections
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params
    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)
    reflections = flatten_reflections(params.input.reflections)

    # Find all detector objects
    detectors = []
    detectors.extend(experiments.detectors())
    dbs = []
    for datablock in datablocks:
      dbs.extend(datablock.unique_detectors())
    detectors.extend(dbs)

    # Verify inputs
    if len(detectors) != 2:
      raise Sorry("Please provide a reference and a moving set of experiments and or datablocks")

    reflections = reflections[1]
    detector = detectors[1]

    if not hasattr(detector, 'hierarchy'):
      raise Sorry("Script intended for hierarchical detectors")

    if params.max_hierarchy_level is None or str(params.max_hierarchy_level).lower() == 'auto':
      params.max_hierarchy_level = 0
      root = detector.hierarchy()
      while hasattr(root, 'children'):
        root = root[0]
        params.max_hierarchy_level += 1
      print "Found", params.max_hierarchy_level+1, "hierarchy levels"

    reference_root = detectors[0].hierarchy()
    moving_root = detector.hierarchy()
    rori = get_center(reference_root)
    rf = col(reference_root.get_fast_axis())
    rs = col(reference_root.get_slow_axis())
    r_norm = col(reference_root.get_normal())

    table_header = ["Hierarchy","Delta XY","Delta XY","R Offsets","R Offsets","T Offsets","T Offsets","Z Offsets","Z Offsets","dR Norm","dR Norm","dT Norm","dT Norm","Local dNorm", "Local dNorm", "Rot Z","Rot Z"]
    table_header2 = ["Level","","Sigma","","Sigma","","Sigma","","Sigma","","Sigma","","Sigma","","Sigma","","Sigma"]
    table_header3 = ["","(microns)","(microns)","(microns)","(microns)","(microns)","(microns)","(microns)","(microns)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)"]
    table_data = []
    table_data.append(table_header)
    table_data.append(table_header2)
    table_data.append(table_header3)

    def get_parent(pg):
      if hasattr(pg, 'children'):
        return pg.parent()
      else:
        return pg.parent

    # Iterate through the hierarchy levels
    for level in xrange(params.max_hierarchy_level+1):
      delta_xy = flex.double()
      r_offsets = flex.double()
      t_offsets = flex.double()
      z_offsets = flex.double()
      rot_z = flex.double()
      delta_r_norm = flex.double()
      delta_t_norm = flex.double()
      local_dnorm = flex.double()
      weights = flex.double()

      for pg_id, (pg1, pg2) in enumerate(zip(iterate_detector_at_level(reference_root, 0, level),
                                             iterate_detector_at_level(moving_root, 0, level))):
        weight = 0
        for panel_id, p in enumerate(iterate_panels(pg2)):
          weight += len(reflections.select(reflections['panel'] == id_from_name(detector, p.get_name())))
        weights.append(weight)

        z_dists = []
        ori_xy = []
        for pg in [pg1,pg2]:
          ori = pg.get_local_origin()
          ori_xy.append(col((ori[0], ori[1])))
          z_dists.append(ori[2]*1000)
        delta_xy.append((ori_xy[1]-ori_xy[0]).length()*1000)
        z_offsets.append(z_dists[1]-z_dists[0])

        pgo1 = get_center(pg1)
        pgo2 = get_center(pg2)
        ro_pgo = pgo2 - rori # vector from the detector origin to the panel group origin
        if ro_pgo.length() == 0:
          radial = col((0,0,0))
          transverse = col((0,0,0))
        else:
          radial = ((rf.dot(ro_pgo) * rf) + (rs.dot(ro_pgo) * rs)).normalize() # component of ro_pgo in rf rs plane
          transverse = r_norm.cross(radial).normalize()
        # now radial and transverse are vectors othogonal to each other and the detector normal, such that
        # radial points at the panel group origin

        # After the radial and transverse vectors are determined, adjust the origins to be relative to their parent objects' origins.
        parent1 = get_parent(pg1)
        if parent1 is not None:
          parent2 = get_parent(pg2)
          pgo1 = pgo1 - get_center(parent1)
          pgo2 = pgo2 - get_center(parent2)

        # v is the component of delta_pgo along the radial vector
        delta_pgo = pgo2-pgo1
        v = (radial.dot(delta_pgo) * radial)
        r = v.length() * 1000
        angle = r_norm.angle(v, deg=True)
        if r_norm.cross(v).dot(transverse) < 0:
          r = -r
        r_offsets.append(r)
        # v is the component of delta_pgo along the transverse vector
        v = (transverse.dot(delta_pgo) * transverse)
        t = v.length() * 1000
        angle = r_norm.angle(v, deg=True)
        if r_norm.cross(v).dot(radial) < 0:
          t = -t
        t_offsets.append(t)

        pgn1 = col(pg1.get_normal())
        pgf1 = col(pg1.get_fast_axis())
        pgs1 = col(pg1.get_slow_axis())
        pgn2 = col(pg2.get_normal())
        pgf2 = col(pg2.get_fast_axis())

        # v1 and v2 are the component of pgf1 and pgf2 in the rf rs plane
        v1 = (rf.dot(pgf1) * rf) + (rs.dot(pgf1) * rs)
        v2 = (rf.dot(pgf2) * rf) + (rs.dot(pgf2) * rs)
        angle = v1.angle(v2, deg=True)
        rot_z.append(angle)

        # v1 and v2 are the components of pgn1 and pgn2 in the r_norm radial plane
        v1 = (r_norm.dot(pgn1) * r_norm) + (radial.dot(pgn1) * radial)
        v2 = (r_norm.dot(pgn2) * r_norm) + (radial.dot(pgn2) * radial)
        angle = v1.angle(v2, deg=True)
        if v2.cross(v1).dot(transverse) < 0:
          angle = -angle
        delta_r_norm.append(angle)

        # v1 and v2 are the components of pgn1 and pgn2 in the r_norm transverse plane
        v1 = (r_norm.dot(pgn1) * r_norm) + (transverse.dot(pgn1) * transverse)
        v2 = (r_norm.dot(pgn2) * r_norm) + (transverse.dot(pgn2) * transverse)
        angle = v1.angle(v2, deg=True)
        if v2.cross(v1).dot(radial) < 0:
          angle = -angle
        delta_t_norm.append(angle)

        # Determine angle between normals in local space
        lpgf1 = col(pg1.get_local_fast_axis())
        lpgs1 = col(pg1.get_local_slow_axis())
        lpgn1 = lpgf1.cross(lpgs1)
        lpgf2 = col(pg2.get_local_fast_axis())
        lpgs2 = col(pg2.get_local_slow_axis())
        lpgn2 = lpgf2.cross(lpgs2)
        local_dnorm.append(lpgn1.angle(lpgn2, deg=True))

      row = ["%d"%level]
      iterable = zip([delta_xy, r_offsets, t_offsets, z_offsets, delta_r_norm, delta_t_norm, local_dnorm, rot_z],
                     ["%6.1f","%6.1f","%6.1f","%6.1f","%.4f","%.4f","%.4f","%.4f"])
      if len(z_offsets) == 0:
        row.extend(["%6.1f"%0]*6)
      elif len(z_offsets) == 1:
        for data, fmt in iterable:
          row.append(fmt%data[0])
          row.append(fmt%0)
      else:
        for data, fmt in iterable:
          stats = flex.mean_and_variance(data, weights)
          row.append(fmt%stats.mean())
          row.append(fmt%stats.gsl_stats_wsd())
      table_data.append(row)

    from libtbx import table_utils
    print "Detector shifts"
    print table_utils.format(table_data,has_header=3,justify='center',delim=" ")

    print
    print """
For each hierarchy level, the average shifts in are computed among objects at that level, weighted by the number of reflections recorded on each object. For example, for a four quadrant detector, the average Z shift will be the average of the four quadrant Z values, each weighted by the number of reflections on that quadrant.

Individual columns:
Delta XY: magnitude of the shift in the local XY frame.
R, T offsets: shifts relative to the parent object's location in the radial and transverse directions (relative to the detector center).
Z offsets: relative shifts in the local frame in the local Z direction.
R, T Norm: angle between normal vectors in lab space, projected onto the radial or transverse plane.
Local dNorm: local relative angle between normal vectors.
Rot Z: rotation around detector normal in lab space
"""

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
