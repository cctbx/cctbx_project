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
# LIBTBX_SET_DISPATCHER_NAME cspad.detector_congruence
#
from __future__ import absolute_import, division, print_function
from six.moves import range
from dials.array_family import flex
from dials.util import show_mail_on_error
from scitbx.matrix import col
import numpy as np
from libtbx.phil import parse
import libtbx.load_env
import math
from six.moves import zip
from serialtbx.detector import iterate_detector_at_level, iterate_panels, id_from_name, get_center

try:
  from matplotlib import pyplot as plt
  from matplotlib.patches import Polygon
  from matplotlib.colors import Normalize
  from matplotlib import cm
except ImportError:
  pass # Can happen on non-GUI nodes

help_message = '''

This program is used to calculate statisical measurements of consistency
between two detectors.

Example:

  %s experiment1.expt experiment2.expt reflections1.refl reflections2.refl
''' % libtbx.env.dispatcher_name

# Create the phil parameters
phil_scope = parse('''
tag = None
  .type = str
  .help = Used in the plot titles
hierarchy_level=0
  .type=int
  .help=Provide congruence statistics for detector modules at the given hierarchy level.
colormap=RdYlGn_r
  .type=str
  .help=matplotlib color map. See e.g.: \
        http://matplotlib.org/examples/color/colormaps_reference.html
show_plots=True
  .type=bool
  .help=Whether to show congruence plots
draw_normal_arrows=False
  .type=bool
  .help=Whether to draw the XY components of each panel group's normal vector. Useful \
        for visualizing tilt.
''')

def get_bounds(root, pg):
  """ Find the max extent of the panel group pg, projected onto the fast/slow plane of root """
  def panel_bounds(root, panel):
    size = panel.get_image_size()
    p0 = col(panel.get_pixel_lab_coord((0,0)))
    p1 = col(panel.get_pixel_lab_coord((size[0]-1,0)))
    p2 = col(panel.get_pixel_lab_coord((size[0]-1,size[1]-1)))
    p3 = col(panel.get_pixel_lab_coord((0,size[1]-1)))

    rn = col(root.get_normal())
    rf = col(root.get_fast_axis())
    rs = col(root.get_slow_axis())

    return [col((p.dot(rf), + p.dot(rs),0)) for p in [p0, p1, p2, p3]]

  if pg.is_group():
    minx = miny = float('inf')
    maxx = maxy = float('-inf')
    for panel in iterate_panels(pg):
      bounds = panel_bounds(root, panel)
      for v in bounds:
        if v[0] < minx:
          minx = v[0]
        if v[0] > maxx:
          maxx = v[0]
        if v[1] < miny:
          miny = v[1]
        if v[1] > maxy:
          maxy = v[1]
    return [col((minx, miny, 0)),
            col((maxx, miny, 0)),
            col((maxx, maxy, 0)),
            col((minx, maxy, 0))]

  else:
    return panel_bounds(root, pg)

def detector_plot_dict(params, detector, data, title, units_str, show=True, reverse_colormap=False):
  """
  Use matplotlib to plot a detector, color coding panels according to data
  @param detector detector reference detector object
  @param data python dictionary of panel names as keys and numbers as values
  @param title title string for plot
  @param units_str string with a formatting statment for units on each panel
  """
  # initialize the color map
  values = flex.double(list(data.values()))
  norm = Normalize(vmin=flex.min(values), vmax=flex.max(values))
  cmap = plt.get_cmap(params.colormap + ("_r" if reverse_colormap else ''))
  sm = cm.ScalarMappable(norm=norm, cmap=cmap)
  if len(values) == 0:
    print("no values")
    return
  elif len(values) == 1:
    sm.set_array(np.arange(values[0], values[0], 1)) # needed for colorbar
  else:
    sm.set_array(np.arange(flex.min(values), flex.max(values), (flex.max(values)-flex.min(values))/20)) # needed for colorbar

  fig = plt.figure()
  ax = fig.add_subplot(111, aspect='equal')
  max_dim = 0
  root = detector.hierarchy()
  rf = col(root.get_fast_axis())
  rs = col(root.get_slow_axis())
  for pg_id, pg in enumerate(iterate_detector_at_level(root, 0, params.hierarchy_level)):
    if pg.get_name() not in data:
      continue
    # get panel coordinates
    p0, p1, p2, p3 = get_bounds(root, pg)

    v1 = p1-p0
    v2 = p3-p0
    vcen = ((v2/2) + (v1/2)) + p0

    # add the panel to the plot
    ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color=sm.to_rgba(data[pg.get_name()]), fill=True))
    ax.annotate("%d %s"%(pg_id, units_str%data[pg.get_name()]), vcen[0:2], ha='center')

    if params.draw_normal_arrows:
      pgn = col(pg.get_normal())
      v = col((rf.dot(pgn), rs.dot(pgn), 0))
      v *= 10000
      ax.arrow(vcen[0], vcen[1], v[0], v[1], head_width=5.0, head_length=10.0, fc='k', ec='k')

    # find the plot maximum dimensions
    for p in [p0, p1, p2, p3]:
      for c in p[0:2]:
        if abs(c) > max_dim:
          max_dim = abs(c)

  # plot the results
  ax.set_xlim((-max_dim,max_dim))
  ax.set_ylim((-max_dim,max_dim))
  ax.set_xlabel("mm")
  ax.set_ylabel("mm")
  fig.colorbar(sm, ax=ax)
  plt.title(title)
  if show:
    plt.show()

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import ArgumentParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s experiment1.expt experiment2.expt reflections1.refl reflections2.refl" % libtbx.env.dispatcher_name
    self.parser = ArgumentParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
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
      print("Please provide two experiments for comparison")
      return

    # These lines exercise the iterate_detector_at_level and iterate_panels functions
    # for a detector with 4 hierarchy levels
    """
    print "Testing iterate_detector_at_level"
    for level in range(4):
      print "iterating at level", level
      for panelg in iterate_detector_at_level(detectors[0].hierarchy(), 0, level):
        print panelg.get_name()

    print "Testing iterate_panels"
    for level in range(4):
      print "iterating at level", level
      for panelg in iterate_detector_at_level(detectors[0].hierarchy(), 0, level):
        for panel in iterate_panels(panelg):
          print panel.get_name()
    """
    tmp = []
    for refls in reflections:
      print("N reflections total:", len(refls))
      refls = refls.select(refls.get_flags(refls.flags.used_in_refinement))
      print("N reflections used in refinement", len(refls))
      print("Reporting only on those reflections used in refinement")

      refls['difference_vector_norms'] = (refls['xyzcal.mm']-refls['xyzobs.mm.value']).norms()
      tmp.append(refls)
    reflections = tmp

    # Iterate through the detectors, computing the congruence statistics
    delta_normals = {}
    z_angles = {}
    f_deltas = {}
    s_deltas = {}
    z_deltas = {}
    o_deltas = {} # overall
    z_offsets_d = {}
    refl_counts = {}
    all_delta_normals = flex.double()
    all_rdelta_normals = flex.double()
    all_tdelta_normals = flex.double()
    all_z_angles = flex.double()
    all_f_deltas = flex.double()
    all_s_deltas = flex.double()
    all_z_deltas = flex.double()
    all_deltas = flex.double()
    all_refls_count = flex.int()

    all_normal_angles = flex.double()
    all_rnormal_angles = flex.double()
    all_tnormal_angles = flex.double()
    pg_normal_angle_sigmas = flex.double()
    pg_rnormal_angle_sigmas = flex.double()
    pg_tnormal_angle_sigmas = flex.double()
    all_rot_z = flex.double()
    pg_rot_z_sigmas = flex.double()
    pg_bc_dists = flex.double()
    all_bc_dist = flex.double()
    all_f_offsets = flex.double()
    all_s_offsets = flex.double()
    all_z_offsets = flex.double()
    pg_f_offset_sigmas = flex.double()
    pg_s_offset_sigmas = flex.double()
    pg_z_offset_sigmas = flex.double()
    pg_offset_sigmas = flex.double()
    all_weights = flex.double()

    congruence_table_data = []
    detector_table_data = []
    rmsds_table_data = []
    root1 = detectors[0].hierarchy()
    root2 = detectors[1].hierarchy()

    s0 = col(flex.vec3_double([col(b.get_s0()) for b in experiments.beams()]).mean())

    # Compute a set of radial and transverse displacements for each reflection
    print("Setting up stats...")
    tmp_refls = []
    for refls, expts in zip(reflections, [wrapper.data for wrapper in params.input.experiments]):
      tmp = flex.reflection_table()
      assert len(expts.detectors()) == 1
      dect = expts.detectors()[0]
      # Need to construct a variety of vectors
      for panel_id, panel in enumerate(dect):
        panel_refls = refls.select(refls['panel'] == panel_id)
        bcl = flex.vec3_double()
        # Compute the beam center in lab space (a vector pointing from the origin to where the beam would intersect
        # the panel, if it did intersect the panel)
        for expt_id in set(panel_refls['id']):
          beam = expts[expt_id].beam
          s0 = beam.get_s0()
          expt_refls = panel_refls.select(panel_refls['id'] == expt_id)
          beam_centre = panel.get_beam_centre_lab(s0)
          bcl.extend(flex.vec3_double(len(expt_refls), beam_centre))
        panel_refls['beam_centre_lab'] = bcl

        # Compute obs in lab space
        x, y, _ = panel_refls['xyzobs.mm.value'].parts()
        c = flex.vec2_double(x, y)
        panel_refls['obs_lab_coords'] = panel.get_lab_coord(c)
        # Compute deltaXY in panel space. This vector is relative to the panel origin
        x, y, _ = (panel_refls['xyzcal.mm'] - panel_refls['xyzobs.mm.value']).parts()
        # Convert deltaXY to lab space, subtracting off of the panel origin
        panel_refls['delta_lab_coords'] = panel.get_lab_coord(flex.vec2_double(x,y)) - panel.get_origin()
        tmp.extend(panel_refls)
      refls = tmp
      # The radial vector points from the center of the reflection to the beam center
      radial_vectors = (refls['obs_lab_coords'] - refls['beam_centre_lab']).each_normalize()
      # The transverse vector is orthogonal to the radial vector and the beam vector
      transverse_vectors = radial_vectors.cross(refls['beam_centre_lab']).each_normalize()
      # Compute the raidal and transverse components of each deltaXY
      refls['radial_displacements']     = refls['delta_lab_coords'].dot(radial_vectors)
      refls['transverse_displacements'] = refls['delta_lab_coords'].dot(transverse_vectors)

      tmp_refls.append(refls)
    reflections = tmp_refls

    for pg_id, (pg1, pg2) in enumerate(zip(iterate_detector_at_level(root1, 0, params.hierarchy_level),
                                           iterate_detector_at_level(root2, 0, params.hierarchy_level))):
      """ First compute statistics for detector congruence """
      # Count up the number of reflections in this panel group pair for use as a weighting scheme
      total_refls = 0
      pg1_refls = 0
      pg2_refls = 0
      for p1, p2 in zip(iterate_panels(pg1), iterate_panels(pg2)):
        r1 = len(reflections[0].select(reflections[0]['panel'] == id_from_name(detectors[0], p1.get_name())))
        r2 = len(reflections[1].select(reflections[1]['panel'] == id_from_name(detectors[1], p2.get_name())))
        total_refls += r1 + r2
        pg1_refls += r1
        pg2_refls += r2
      if pg1_refls == 0 and pg2_refls == 0:
        print("No reflections on panel group", pg_id)
        continue

      assert pg1.get_name() == pg2.get_name()
      refl_counts[pg1.get_name()] = total_refls

      row = ["%d"%pg_id]
      for pg, refls, det in zip([pg1, pg2], reflections, detectors):
        pg_refls = flex.reflection_table()
        for p in iterate_panels(pg):
          pg_refls.extend(refls.select(refls['panel'] == id_from_name(det, p.get_name())))
        if len(pg_refls) == 0:
          rmsd = r_rmsd = t_rmsd = 0
        else:
          rmsd = math.sqrt(flex.sum_sq(pg_refls['difference_vector_norms'])/len(pg_refls))*1000
          r_rmsd = math.sqrt(flex.sum_sq(pg_refls['radial_displacements'])/len(pg_refls))*1000
          t_rmsd = math.sqrt(flex.sum_sq(pg_refls['transverse_displacements'])/len(pg_refls))*1000

        row.extend(["%6.1f"%rmsd, "%6.1f"%r_rmsd, "%6.1f"%t_rmsd, "%8d"%len(pg_refls)])
      rmsds_table_data.append(row)

      # Angle between normals of pg1 and pg2
      delta_norm_angle = col(pg1.get_normal()).angle(col(pg2.get_normal()), deg=True)
      all_delta_normals.append(delta_norm_angle)

      # compute radial and transverse components of the delta between normal angles
      pgo = (get_center(pg1)+get_center(pg2))/2
      ro = (get_center(root1)+get_center(root2))/2
      rn = (col(root1.get_normal())+col(root2.get_normal()))/2
      rf = (col(root1.get_fast_axis())+col(root2.get_fast_axis()))/2
      rs = (col(root1.get_slow_axis())+col(root2.get_slow_axis()))/2

      ro_pgo = pgo - ro # vector from the detector origin to the average panel group origin
      if ro_pgo.length() == 0:
        radial = col((0,0,0))
        transverse = col((0,0,0))
      else:
        radial = ((rf.dot(ro_pgo) * rf) + (rs.dot(ro_pgo) * rs)).normalize() # component of ro_pgo in rf rs plane
        transverse = rn.cross(radial).normalize()
      # now radial and transverse are vectors othogonal to each other and the detector normal, such that
      # radial points at the panel group origin
      # v1 and v2 are the components of pg 1 and 2 normals in the rn radial plane
      v1 = (radial.dot(col(pg1.get_normal())) * radial) + (rn.dot(col(pg1.get_normal())) * rn)
      v2 = (radial.dot(col(pg2.get_normal())) * radial) + (rn.dot(col(pg2.get_normal())) * rn)
      rdelta_norm_angle = v1.angle(v2, deg=True)
      if v1.cross(v2).dot(transverse) < 0:
        rdelta_norm_angle = -rdelta_norm_angle
      all_rdelta_normals.append(rdelta_norm_angle)
      # v1 and v2 are the components of pg 1 and 2 normals in the rn transverse plane
      v1 = (transverse.dot(col(pg1.get_normal())) * transverse) + (rn.dot(col(pg1.get_normal())) * rn)
      v2 = (transverse.dot(col(pg2.get_normal())) * transverse) + (rn.dot(col(pg2.get_normal())) * rn)
      tdelta_norm_angle = v1.angle(v2, deg=True)
      if v1.cross(v2).dot(radial) < 0:
        tdelta_norm_angle = -tdelta_norm_angle
      all_tdelta_normals.append(tdelta_norm_angle)

      # compute the angle between fast axes of these panel groups
      z_angle = col(pg1.get_fast_axis()[0:2]).angle(col(pg2.get_fast_axis()[0:2]), deg=True)
      all_z_angles.append(z_angle)
      z_angles[pg1.get_name()] = z_angle

      all_refls_count.append(total_refls)
      all_weights.append(pg1_refls)
      all_weights.append(pg2_refls)


      """ Now compute statistics measuring the reality of the detector. For example, instead of the distance between two things,
      we are concerned with the location of those things relative to laboratory space """
      # Compute distances between panel groups and beam center
      # Also compute offset along Z axis
      dists = flex.double()
      f_offsets = flex.double()
      s_offsets = flex.double()
      z_offsets = flex.double()
      for pg, r in zip([pg1, pg2], [root1, root2]):
        bc = col(pg.get_beam_centre_lab(s0))
        ori = get_center(pg)

        dists.append((ori-bc).length())

        rori = col(r.get_origin())
        delta_ori = ori-rori
        r_norm = col(r.get_normal())
        r_fast = col(r.get_fast_axis())
        r_slow = col(r.get_slow_axis())
        f_offsets.append(r_fast.dot(delta_ori)*1000)
        s_offsets.append(r_slow.dot(delta_ori)*1000)
        z_offsets.append(r_norm.dot(delta_ori)*1000)

      fd = abs(f_offsets[0]-f_offsets[1])
      sd = abs(s_offsets[0]-s_offsets[1])
      zd = abs(z_offsets[0]-z_offsets[1])
      od = math.sqrt(fd**2+sd**2+zd**2)
      f_deltas[pg1.get_name()] = fd
      s_deltas[pg1.get_name()] = sd
      z_deltas[pg1.get_name()] = zd
      o_deltas[pg1.get_name()] = od
      all_f_deltas.append(fd)
      all_s_deltas.append(sd)
      all_z_deltas.append(zd)
      all_deltas.append(od)

      all_f_offsets.extend(f_offsets)
      all_s_offsets.extend(s_offsets)
      all_z_offsets.extend(z_offsets)

      # Compute angle between detector normal and panel group normal
      # Compute rotation of panel group around detector normal
      pg_rotz = flex.double()
      norm_angles = flex.double()
      rnorm_angles = flex.double()
      tnorm_angles = flex.double()
      for pg, r in zip([pg1, pg2], [root1, root2]):

        pgo = get_center(pg)
        pgn = col(pg.get_normal())
        pgf = col(pg.get_fast_axis())

        ro = get_center(r)
        rn = col(r.get_normal())
        rf = col(r.get_fast_axis())
        rs = col(r.get_slow_axis())

        norm_angle = rn.angle(pgn, deg=True)
        norm_angles.append(norm_angle)
        all_normal_angles.append(norm_angle)

        ro_pgo = pgo - ro # vector from the detector origin to the panel group origin
        if ro_pgo.length() == 0:
          radial = col((0,0,0))
          transverse = col((0,0,0))
        else:
          radial = ((rf.dot(ro_pgo) * rf) + (rs.dot(ro_pgo) * rs)).normalize() # component of ro_pgo in rf rs plane
          transverse = rn.cross(radial).normalize()
        # now radial and transverse are vectors othogonal to each other and the detector normal, such that
        # radial points at the panel group origin
        # v is the component of pgn in the rn radial plane
        v = (radial.dot(pgn) * radial) + (rn.dot(pgn) * rn)
        angle = rn.angle(v, deg=True)
        if rn.cross(v).dot(transverse) < 0:
          angle = -angle
        rnorm_angles.append(angle)
        all_rnormal_angles.append(angle)
        # v is the component of pgn in the rn transverse plane
        v = (transverse.dot(pgn) * transverse) + (rn.dot(pgn) * rn)
        angle = rn.angle(v, deg=True)
        if rn.cross(v).dot(radial) < 0:
          angle = -angle
        tnorm_angles.append(angle)
        all_tnormal_angles.append(angle)

        # v is the component of pgf in the rf rs plane
        v = (rf.dot(pgf) * rf) + (rs.dot(pgf) * rs)
        angle = rf.angle(v, deg=True)
        angle = angle-(round(angle/90)*90) # deviation from 90 degrees
        pg_rotz.append(angle)
        all_rot_z.append(angle)

      # Set up table rows using stats aggregated from above
      pg_weights = flex.double([pg1_refls, pg2_refls])
      if 0 in pg_weights:
        dist_m = dist_s = norm_angle_m = norm_angle_s = rnorm_angle_m = rnorm_angle_s = 0
        tnorm_angle_m = tnorm_angle_s = rotz_m = rotz_s = 0
        fo_m = fo_s = so_m = so_s = zo_m = zo_s = o_s = 0

      else:
        stats = flex.mean_and_variance(dists, pg_weights)
        dist_m = stats.mean()
        dist_s = stats.gsl_stats_wsd()

        stats = flex.mean_and_variance(norm_angles, pg_weights)
        norm_angle_m = stats.mean()
        norm_angle_s = stats.gsl_stats_wsd()

        stats = flex.mean_and_variance(rnorm_angles, pg_weights)
        rnorm_angle_m = stats.mean()
        rnorm_angle_s = stats.gsl_stats_wsd()

        stats = flex.mean_and_variance(tnorm_angles, pg_weights)
        tnorm_angle_m = stats.mean()
        tnorm_angle_s = stats.gsl_stats_wsd()

        stats = flex.mean_and_variance(pg_rotz, pg_weights)
        rotz_m = stats.mean()
        rotz_s = stats.gsl_stats_wsd()

        stats = flex.mean_and_variance(f_offsets, pg_weights)
        fo_m = stats.mean()
        fo_s = stats.gsl_stats_wsd()
        stats = flex.mean_and_variance(s_offsets, pg_weights)
        so_m = stats.mean()
        so_s = stats.gsl_stats_wsd()
        stats = flex.mean_and_variance(z_offsets, pg_weights)
        zo_m = stats.mean()
        zo_s = stats.gsl_stats_wsd()

        o_s = math.sqrt(fo_s**2+so_s**2+zo_s**2)

      pg_bc_dists.append(dist_m)
      all_bc_dist.extend(dists)
      pg_normal_angle_sigmas.append(norm_angle_s)
      pg_rnormal_angle_sigmas.append(rnorm_angle_s)
      pg_tnormal_angle_sigmas.append(tnorm_angle_s)
      pg_rot_z_sigmas.append(rotz_s)
      pg_f_offset_sigmas.append(fo_s)
      pg_s_offset_sigmas.append(so_s)
      pg_z_offset_sigmas.append(zo_s)
      pg_offset_sigmas.append(o_s)
      z_offsets_d[pg1.get_name()] = zo_m

      congruence_table_data.append(["%d"%pg_id, "%5.1f"%dist_m, #"%.4f"%dist_s,
                                    "%.4f"%delta_norm_angle, "%.4f"%rdelta_norm_angle,
                                    "%.4f"%tdelta_norm_angle, "%.4f"%z_angle,
                                    "%4.1f"%fd, "%4.1f"%sd, "%4.1f"%zd, "%4.1f"%od, "%6d"%total_refls])
      detector_table_data.append(["%d"%pg_id, "%5.1f"%dist_m, #"%.4f"%dist_s,
                                  "%.4f"%norm_angle_m, "%.4f"%norm_angle_s,
                                  "%.4f"%rnorm_angle_m, "%.4f"%rnorm_angle_s,
                                  "%.4f"%tnorm_angle_m, "%.4f"%tnorm_angle_s,
                                  "%10.6f"%rotz_m, "%.6f"%rotz_s,
                                  #"%9.1f"%fo_m, "%5.3f"%fo_s,
                                  #"%9.1f"%so_m, "%5.3f"%so_s,
                                  "%9.3f"%fo_s,
                                  "%9.3f"%so_s,
                                  "%9.1f"%zo_m, "%9.1f"%zo_s, "%9.3f"%o_s, "%6d"%total_refls])

    # Set up table output
    table_d = {d:row for d, row in zip(pg_bc_dists, congruence_table_data)}
    table_header = ["PanelG","Dist","Normal","RNormal","TNormal","Z rot","Delta","Delta","Delta","Delta","N"]
    table_header2 = ["Id","","Angle","Angle","Angle","Angle","F","S","Z","O","Refls"]
    table_header3 = ["", "(mm)","(mm)","(deg)","(deg)","(microns)","(microns)","(microns)","(microns)","(microns)",""]
    congruence_table_data = [table_header, table_header2, table_header3]
    congruence_table_data.extend([table_d[key] for key in sorted(table_d)])

    table_d = {d:row for d, row in zip(pg_bc_dists, detector_table_data)}
    table_header = ["PanelG","Dist","Normal","Normal","RNormal","RNormal","TNormal","TNormal","RotZ", "RotZ","F Offset","S Offset","Z Offset","Z Offset","Offset","N"]
    table_header2 = ["Id","","","Sigma","","Sigma","","Sigma","","Sigma","Sigma","Sigma","","Sigma","Sigma","Refls"]
    table_header3 = ["", "(mm)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)","(deg)","(microns)","(microns)","(microns)","(microns)","(microns)",""]
    detector_table_data = [table_header, table_header2, table_header3]
    detector_table_data.extend([table_d[key] for key in sorted(table_d)])

    table_d = {d:row for d, row in zip(pg_bc_dists, rmsds_table_data)}
    table_header = ["PanelG"]
    table_header2 = ["Id"]
    table_header3 = [""]
    for i in range(len(detectors)):
      table_header.extend(["D%d"%i]*4)
      table_header2.extend(["RMSD", "rRMSD", "tRMSD", "N refls"])
      table_header3.extend(["(microns)"]*3)
      table_header3.append("")
    rmsds_table_data = [table_header, table_header2, table_header3]
    rmsds_table_data.extend([table_d[key] for key in sorted(table_d)])

    if len(all_refls_count) > 1:
      r1 = ["Weighted mean"]
      r2 = ["Weighted stddev"]
      r1.append("")
      r2.append("")
      #r1.append("")
      #r2.append("")
      stats = flex.mean_and_variance(all_delta_normals, all_refls_count.as_double())
      r1.append("%.4f"%stats.mean())
      r2.append("%.4f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(all_rdelta_normals, all_refls_count.as_double())
      r1.append("%.4f"%stats.mean())
      r2.append("%.4f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(all_tdelta_normals, all_refls_count.as_double())
      r1.append("%.4f"%stats.mean())
      r2.append("%.4f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(all_z_angles, all_refls_count.as_double())
      r1.append("%.4f"%stats.mean())
      r2.append("%.4f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(all_f_deltas, all_refls_count.as_double())
      r1.append("%4.1f"%stats.mean())
      r2.append("%4.1f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(all_s_deltas, all_refls_count.as_double())
      r1.append("%4.1f"%stats.mean())
      r2.append("%4.1f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(all_z_deltas, all_refls_count.as_double())
      r1.append("%4.1f"%stats.mean())
      r2.append("%4.1f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(all_deltas, all_refls_count.as_double())
      r1.append("%4.1f"%stats.mean())
      r2.append("%4.1f"%stats.gsl_stats_wsd())
      r1.append("")
      r2.append("")
      congruence_table_data.append(r1)
      congruence_table_data.append(r2)
      congruence_table_data.append(["Mean", "", "", "","","", "", "", "", "", "", "%6.1f"%flex.mean(all_refls_count.as_double())])

    from libtbx import table_utils
    print("Congruence statistics, I.E. the differences between the input detectors:")
    print(table_utils.format(congruence_table_data,has_header=3,justify='center',delim=" "))

    print("PanelG Id: panel group id or panel id, depending on hierarchy_level. For each panel group, statistics are computed between the matching panel groups between the two input experiments.")
    print("Dist: distance from center of panel group to the beam center")
    print("Dist Sigma: weighted standard deviation of the measurements used to compute Dist")
    print("Normal angle: angle between the normal vectors of matching panel groups.")
    print("RNormal angle: radial component of the angle between the normal vectors of matching panel groups")
    print("TNormal angle: transverse component of the angle between the normal vectors of matching panel groups")
    print("Z rot: angle between the XY components of the fast axes of the panel groups.")
    print("Delta F: shift between matching panel groups along the detector fast axis.")
    print("Delta S: shift between matching panel groups along the detector slow axis.")
    print("Delta Z: Z shift between matching panel groups along the detector normal.")
    print("Delta O: Overall shift between matching panel groups along the detector normal.")
    print("N refls: number of reflections summed between both matching panel groups. This number is used as a weight when computing means and standard deviations.")
    print()
    print()


    if len(all_weights) > 1:
      r1 = ["All"]
      r2 = ["Mean"]
      for data, weights, fmt in [[None,None,None],
                                 #[None,None,None],
                                 [all_normal_angles,       all_weights.as_double(),     "%.4f"],
                                 [pg_normal_angle_sigmas,  all_refls_count.as_double(), "%.4f"],
                                 [all_rnormal_angles,      all_weights.as_double(),     "%.4f"],
                                 [pg_rnormal_angle_sigmas, all_refls_count.as_double(), "%.4f"],
                                 [all_tnormal_angles,      all_weights.as_double(),     "%.4f"],
                                 [pg_tnormal_angle_sigmas, all_refls_count.as_double(), "%.4f"],
                                 [all_rot_z,               all_weights.as_double(),     "%10.6f"],
                                 [pg_rot_z_sigmas,         all_refls_count.as_double(), "%.6f"],
                                 #[all_f_offsets,           all_weights.as_double(),     "%9.1f"],
                                 [pg_f_offset_sigmas,      all_refls_count.as_double(), "%9.3f"],
                                 #[all_s_offsets,           all_weights.as_double(),     "%9.1f"],
                                 [pg_s_offset_sigmas,      all_refls_count.as_double(), "%9.3f"],
                                 [all_z_offsets,           all_weights.as_double(),     "%9.1f"],
                                 [pg_z_offset_sigmas,      all_refls_count.as_double(), "%9.1f"],
                                 [pg_offset_sigmas,        all_refls_count.as_double(), "%9.1f"]]:

        r2.append("")
        if data is None and weights is None:
          r1.append("")
          continue
        stats = flex.mean_and_variance(data, weights)
        r1.append(fmt%stats.mean())

      r1.append("")
      r2.append("%6.1f"%flex.mean(all_refls_count.as_double()))
      detector_table_data.append(r1)
      detector_table_data.append(r2)

    print("Detector statistics, I.E. measurements of parameters relative to the detector plane:")
    print(table_utils.format(detector_table_data,has_header=3,justify='center',delim=" "))

    print("PanelG Id: panel group id or panel id, depending on hierarchy_level. For each panel group, weighted means and weighted standard deviations (Sigmas) for the properties listed below are computed using the matching panel groups between the input experiments.")
    print("Dist: distance from center of panel group to the beam center")
    print("Dist Sigma: weighted standard deviation of the measurements used to compute Dist")
    print("Normal Angle: angle between the normal vector of the detector at its root hierarchy level and the normal of the panel group")
    print("RNormal Angle: radial component of Normal Angle")
    print("TNormal Angle: transverse component of Normal Angle")
    print("RotZ: deviation from 90 degrees of the rotation of each panel group around the detector normal")
    print("F Offset: offset of panel group along the detector's fast axis")
    print("S Offset: offset of panel group along the detector's slow axis")
    print("Z Offset: offset of panel group along the detector normal")
    print("Offset: offset of panel group in F,S,Z space. Sigma is F, S, Z offset sigmas summed in quadrature.")
    print("N refls: number of reflections summed between both matching panel groups. This number is used as a weight when computing means and standard deviations.")
    print("All: weighted mean of the values shown")
    print()
    print("Sigmas in this table are computed using the standard deviation of 2 measurements (I.E. a panel's Z Offset is measured twice, once in each input dataset). This is related by a factor of sqrt(2)/2 to the mean of the Delta Z parameter in the congruence statistics table above, which is the difference between Z parameters.")
    print()

    row = ["Overall"]
    for refls in reflections:
      row.append("%6.1f"%(math.sqrt(flex.sum_sq(refls['difference_vector_norms'])/len(refls))*1000))
      row.append("%6.1f"%(math.sqrt(flex.sum_sq(refls['radial_displacements'])/len(refls))*1000))
      row.append("%6.1f"%(math.sqrt(flex.sum_sq(refls['transverse_displacements'])/len(refls))*1000))
      row.append("%8d"%len(refls))
    rmsds_table_data.append(row)

    print("RMSDs by detector number")
    print(table_utils.format(rmsds_table_data,has_header=3,justify='center',delim=" "))
    print("PanelG Id: panel group id or panel id, depending on hierarchy_level")
    print("RMSD: root mean squared deviation between observed and predicted spot locations")
    print("rRMSD: RMSD of radial components of the observed-predicted vectors")
    print("tRMSD: RMSD of transverse components of the observed-predicted vectors")
    print("N refls: number of reflections")

    # Show stats for detector hierarchy root
    def _print_vector(v):
      for i in v:
        print("%10.5f"%i, end=' ')
      print()
    for d_id, d in enumerate(detectors):
      ori = d.hierarchy().get_origin()
      norm = d.hierarchy().get_normal()
      fast = d.hierarchy().get_fast_axis()
      slow = d.hierarchy().get_slow_axis()
      print("Detector", d_id, "origin:   ", end=' '); _print_vector(ori)
      print("Detector", d_id, "normal:   ", end=' '); _print_vector(norm)
      print("Detector", d_id, "fast axis:", end=' '); _print_vector(fast)
      print("Detector", d_id, "slow axis:", end=' '); _print_vector(slow)

    # Unit cell statstics
    lengths = flex.vec3_double()
    angles = flex.vec3_double()
    weights = flex.double()
    for refls, expts in zip(reflections, [d.data for d in params.input.experiments]):
      for crystal_id, crystal in enumerate(expts.crystals()):
        lengths.append(crystal.get_unit_cell().parameters()[0:3])
        angles.append(crystal.get_unit_cell().parameters()[3:6])
        weights.append(len(refls.select(refls['id'] == crystal_id)))

    print("Unit cell stats (angstroms and degrees), weighted means and standard deviations")
    for subset, tags in zip([lengths, angles], [["Cell a", "Cell b", "Cell c"],["Cell alpha", "Cell beta", "Cell gamma"]]):
      for data, tag in zip(subset.parts(), tags):
        stats = flex.mean_and_variance(data, weights)
        print("%s %5.1f +/- %6.3f"%(tag, stats.mean(), stats.gsl_stats_wsd()))

    if params.tag is None:
      tag = ""
    else:
      tag = "%s "%params.tag

    if params.show_plots:
      # Plot the results
      detector_plot_dict(self.params, detectors[0], refl_counts, u"%sN reflections"%tag, u"%6d", show=False)
      #detector_plot_dict(self.params, detectors[0], delta_normals, u"%sAngle between normal vectors (\N{DEGREE SIGN})"%tag, u"%.2f\N{DEGREE SIGN}", show=False)
      detector_plot_dict(self.params, detectors[0], z_angles, u"%sZ rotation angle between panels (\N{DEGREE SIGN})"%tag, u"%.2f\N{DEGREE SIGN}", show=False)
      detector_plot_dict(self.params, detectors[0], f_deltas, u"%sFast displacements between panels (microns)"%tag, u"%4.1f", show=False)
      detector_plot_dict(self.params, detectors[0], s_deltas, u"%sSlow displacements between panels (microns)"%tag, u"%4.1f", show=False)
      detector_plot_dict(self.params, detectors[0], z_offsets_d, u"%sZ offsets along detector normal (microns)"%tag, u"%4.1f", show=False)
      detector_plot_dict(self.params, detectors[0], z_deltas, u"%sZ displacements between panels (microns)"%tag, u"%4.1f", show=False)
      detector_plot_dict(self.params, detectors[0], o_deltas, u"%sOverall displacements between panels (microns)"%tag, u"%4.1f", show=False)
      plt.show()

if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
