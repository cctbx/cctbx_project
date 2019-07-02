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
from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex
from scitbx.matrix import col
from libtbx.phil import parse
from libtbx.utils import Sorry
from xfel.command_line.cspad_detector_congruence import get_center
import libtbx.load_env
from six.moves import zip

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

from xfel.command_line.cspad_detector_congruence import iterate_detector_at_level, iterate_panels, id_from_name
from xfel.command_line.cspad_detector_congruence import Script as ParentScript

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
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_reflections
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    # Find all detector objects
    detectors = []
    detectors.extend(experiments.detectors())

    # Verify inputs
    if len(detectors) != 2:
      raise Sorry("Please provide a reference and a moving set of experiments")

    reflections = reflections[1]
    detector = detectors[1]

    if not hasattr(detector, 'hierarchy'):
      raise Sorry("Script intended for hierarchical detectors")

    if params.max_hierarchy_level is None or str(params.max_hierarchy_level).lower() == 'auto':
      params.max_hierarchy_level = 0
      root = detector.hierarchy()
      while root.is_group():
        root = root[0]
        params.max_hierarchy_level += 1
      print("Found", params.max_hierarchy_level+1, "hierarchy levels")

    reference_root = detectors[0].hierarchy()
    moving_root = detector.hierarchy()
    rori = get_center(reference_root)
    rf = col(reference_root.get_fast_axis())
    rs = col(reference_root.get_slow_axis())
    r_norm = col(reference_root.get_normal())
    s0 = col(flex.vec3_double([col(b.get_s0()) for b in experiments.beams()]).mean())

    summary_table_header = ["Hierarchy","Delta XY","Delta XY","R Offsets","R Offsets","T Offsets","T Offsets","Z Offsets","Z Offsets","dR Norm","dR Norm","dT Norm","dT Norm","Local dNorm", "Local dNorm", "Rot Z","Rot Z"]
    summary_table_header2 = ["Level","","Sigma","","Sigma","","Sigma","","Sigma","","Sigma","","Sigma","","Sigma","","Sigma"]
    summary_table_header3 = ["","(microns)","(microns)","(microns)","(microns)","(microns)","(microns)","(microns)","(microns)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)"]
    summary_table_data = []
    summary_table_data.append(summary_table_header)
    summary_table_data.append(summary_table_header2)
    summary_table_data.append(summary_table_header3)

    table_header = ["PanelG","BC dist","Delta XY","R Offsets","T Offsets","Z Offsets","dR Norm","dT Norm","Local dNorm","Rot Z","N Refls"]
    table_header2 = ["ID","(mm)","(microns)","(microns)","(microns)","(microns)","(deg)","(deg)","(deg)","(deg)",""]

    from xfel.cftbx.detector.cspad_cbf_tbx import basis
    def get_full_basis_shift(pg):
      """Compute basis shift from pg to lab space"""
      shift = basis(panelgroup=pg)
      while True:
        parent = pg.parent()
        if parent is None:
          break
        shift = basis(panelgroup=parent) * shift
        pg = parent
      return shift

    # Iterate through the hierarchy levels
    for level in range(params.max_hierarchy_level+1):
      delta_xy = flex.double()
      r_offsets = flex.double()
      t_offsets = flex.double()
      z_offsets = flex.double()
      rot_z = flex.double()
      delta_r_norm = flex.double()
      delta_t_norm = flex.double()
      local_dnorm = flex.double()
      bc_dists = flex.double()
      weights = flex.double()

      rows = []

      for pg_id, (pg1, pg2) in enumerate(zip(iterate_detector_at_level(reference_root, 0, level),
                                             iterate_detector_at_level(moving_root, 0, level))):
        weight = 0
        for panel_id, p in enumerate(iterate_panels(pg2)):
          weight += len(reflections.select(reflections['panel'] == id_from_name(detector, p.get_name())))
        weights.append(weight)

        bc = col(pg1.get_beam_centre_lab(s0))
        ori = get_center(pg1)
        bc_dist = (ori-bc).length()
        bc_dists.append(bc_dist)

        z_dists = []
        ori_xy = []
        for pg in [pg1,pg2]:
          ori = pg.get_local_origin()
          ori_xy.append(col((ori[0], ori[1])))
          z_dists.append(ori[2]*1000)
        dxy = (ori_xy[1]-ori_xy[0]).length()*1000
        delta_xy.append(dxy)

        z_off = z_dists[1]-z_dists[0]
        z_offsets.append(z_off)

        pgo1 = col(pg1.get_origin())
        ro_pgo = pgo1 - rori # vector from the detector origin to the panel group origin
        if ro_pgo.length() == 0:
          radial = col((0,0,0))
          transverse = col((0,0,0))
        else:
          radial = ((rf.dot(ro_pgo) * rf) + (rs.dot(ro_pgo) * rs)).normalize() # component of ro_pgo in rf rs plane
          transverse = r_norm.cross(radial).normalize()
        # now radial and transverse are vectors othogonal to each other and the detector normal, such that
        # radial points at the panel group origin

        # compute shift in local frame, then convert that shift to lab space, then make it relative to the reference's origin, in lab space
        lpgo1 = col(pg1.get_local_origin())
        lpgo2 = col(pg2.get_local_origin())
        delta_pgo = (get_full_basis_shift(pg1) * (lpgo2-lpgo1)) - pgo1

        # v is the component of delta_pgo along the radial vector
        v = (radial.dot(delta_pgo) * radial)
        r_offset = v.length() * 1000
        angle = r_norm.angle(v, deg=True)
        if r_norm.cross(v).dot(transverse) < 0:
          r_offset = -r_offset
        r_offsets.append(r_offset)
        # v is the component of delta_pgo along the transverse vector
        v = (transverse.dot(delta_pgo) * transverse)
        t_offset = v.length() * 1000
        angle = r_norm.angle(v, deg=True)
        if r_norm.cross(v).dot(radial) < 0:
          t_offset = -t_offset
        t_offsets.append(t_offset)

        pgn1 = col(pg1.get_normal())
        pgf1 = col(pg1.get_fast_axis())
        pgs1 = col(pg1.get_slow_axis())
        pgn2 = col(pg2.get_normal())
        pgf2 = col(pg2.get_fast_axis())

        # v1 and v2 are the component of pgf1 and pgf2 in the rf rs plane
        v1 = (rf.dot(pgf1) * rf) + (rs.dot(pgf1) * rs)
        v2 = (rf.dot(pgf2) * rf) + (rs.dot(pgf2) * rs)
        rz = v1.angle(v2, deg=True)
        rot_z.append(rz)

        # v1 and v2 are the components of pgn1 and pgn2 in the r_norm radial plane
        v1 = (r_norm.dot(pgn1) * r_norm) + (radial.dot(pgn1) * radial)
        v2 = (r_norm.dot(pgn2) * r_norm) + (radial.dot(pgn2) * radial)
        drn = v1.angle(v2, deg=True)
        if v2.cross(v1).dot(transverse) < 0:
          drn = -drn
        delta_r_norm.append(drn)

        # v1 and v2 are the components of pgn1 and pgn2 in the r_norm transverse plane
        v1 = (r_norm.dot(pgn1) * r_norm) + (transverse.dot(pgn1) * transverse)
        v2 = (r_norm.dot(pgn2) * r_norm) + (transverse.dot(pgn2) * transverse)
        dtn = v1.angle(v2, deg=True)
        if v2.cross(v1).dot(radial) < 0:
          dtn = -dtn
        delta_t_norm.append(dtn)

        # Determine angle between normals in local space
        lpgf1 = col(pg1.get_local_fast_axis())
        lpgs1 = col(pg1.get_local_slow_axis())
        lpgn1 = lpgf1.cross(lpgs1)
        lpgf2 = col(pg2.get_local_fast_axis())
        lpgs2 = col(pg2.get_local_slow_axis())
        lpgn2 = lpgf2.cross(lpgs2)
        ldn = lpgn1.angle(lpgn2, deg=True)
        local_dnorm.append(ldn)

        row = ["%3d"%pg_id, "%6.1f"%bc_dist, "%6.1f"%dxy,
               "%6.1f"%r_offset, "%6.1f"%t_offset, "%6.1f"%z_off,
               "%.4f"%drn, "%.4f"%dtn, "%.4f"%ldn, "%.4f"%rz, "%8d"%weight]
        rows.append(row)

      wm_row = ["Weighted mean", ""]
      ws_row = ["Weighted stddev", ""]
      s_row = ["%d"%level]
      iterable = zip([delta_xy, r_offsets, t_offsets, z_offsets, delta_r_norm, delta_t_norm, local_dnorm, rot_z],
                     ["%6.1f","%6.1f","%6.1f","%6.1f","%.4f","%.4f","%.4f","%.4f"])
      if len(z_offsets) == 0:
        wm_row.extend(["%6.1f"%0]*8)
        ws_row.extend(["%6.1f"%0]*8)
        s_row.extend(["%6.1f"%0]*8)
      elif len(z_offsets) == 1:
        for data, fmt in iterable:
          wm_row.append(fmt%data[0])
          ws_row.append(fmt%0)
          s_row.append(fmt%data[0])
          s_row.append(fmt%0)
      else:
        for data, fmt in iterable:
          stats = flex.mean_and_variance(data, weights)
          wm_row.append(fmt%stats.mean())
          ws_row.append(fmt%stats.gsl_stats_wsd())
          s_row.append(fmt%stats.mean())
          s_row.append(fmt%stats.gsl_stats_wsd())
      wm_row.append("")
      ws_row.append("")
      summary_table_data.append(s_row)

      table_data = [table_header, table_header2]
      table_d = {d:row for d, row in zip(bc_dists, rows)}
      table_data.extend([table_d[key] for key in sorted(table_d)])
      table_data.append(wm_row)
      table_data.append(ws_row)

      from libtbx import table_utils
      print("Hierarchy level %d Detector shifts"%level)
      print(table_utils.format(table_data,has_header=2,justify='center',delim=" "))

    print("Detector shifts summary")
    print(table_utils.format(summary_table_data,has_header=3,justify='center',delim=" "))

    print()
    print("""
For each hierarchy level, the average shifts in are computed among objects at that level, weighted by the number of reflections recorded on each object. For example, for a four quadrant detector, the average Z shift will be the average of the four quadrant Z values, each weighted by the number of reflections on that quadrant.

-------------------
Column descriptions
-------------------

Individual hierarchy level tables only:
PanelG id: ID of the panel group.
BC dist: distance of the panel group from the beam center.
N Refls: number of reflections on this panel group

All tables:
Delta XY: magnitude of the shift in the local XY frame.
R, T offsets: shifts relative to the parent object's location in the radial and transverse directions (relative to the detector center).
Z offsets: relative shifts in the local frame in the local Z direction.
R, T Norm: angle between normal vectors in lab space, projected onto the radial or transverse plane.
Local dNorm: local relative angle between normal vectors.
Rot Z: rotation around detector normal in lab space
""")

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
