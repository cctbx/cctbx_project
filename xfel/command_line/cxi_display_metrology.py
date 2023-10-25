from __future__ import absolute_import, division, print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cxi.display_metrology
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage, Sorry
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle,Polygon
from xfel.cxi.cspad_ana.cspad_tbx import xpp_active_areas
from scitbx.matrix import col

master_phil = libtbx.phil.parse("""
metrology = None
  .type = str
  .help = File with metrology information or XPP active area name
  .optional = False
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isfile(arg)) or arg in xpp_active_areas or os.path.isdir(arg):
      user_phil.append(libtbx.phil.parse("""metrology=\"%s\" """ % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError as e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  if (params.metrology is None) :
    master_phil.show()
    raise Usage("metrology must be defined (either metrology=XXX, or the name alone).")

  fig = plt.figure()
  ax = fig.add_subplot(111, aspect='equal')

  if os.path.isfile(params.metrology):
    # Try dxtbx first to see if this is a regular diffraction image
    import dxtbx.format.Registry
    try:
      # Read the detector object using dxtbx
      reader = dxtbx.format.Registry.get_format_class_for_file(params.metrology)
      detector = reader(params.metrology).get_detector()
    except (IOError, TypeError):
      # See if it's a json file
      from dxtbx.model.experiment_list import ExperimentListFactory
      try:
        experiments = ExperimentListFactory.from_json_file(params.metrology)
        assert len(experiments) == 1
        detector = experiments[0].detector
      except Exception as e:
        detector = None

    if detector is None:
      # see if it's a SLAC geometry file
      from scitbx import matrix
      try:
        from PSCalib.GeometryAccess import GeometryAccess
        geometry = GeometryAccess(params.metrology)
      except Exception as e:
        geometry = None

      if geometry is None:
        # see if this is a Ginn metrology file (see Helen Ginn et. al. (2015), Acta D Cryst)
        panels = []
        for line in open(params.metrology).readlines():
          if len(line) > 0 and line[0] != '#' and len(line.strip().split()) == 10:
            panels.append(line.strip())

        if len(panels) != 64:
          raise Sorry("Can't parse this metrology file")

        print("Ginn metrology file")
        for panel_id, panel in enumerate(panels):
          PANEL, x1, y1, x2, y2, shiftX, shiftY, tiltX, tiltY, _ = panel.split()
          x1 = float(x1) + float(shiftX)
          x2 = float(x2) + float(shiftX)
          y1 = float(y1) + float(shiftY)
          y2 = float(y2) + float(shiftY)

          p0 = col((x1,y1,0))
          p1 = col((x2,y1,0))
          p2 = col((x2,y2,0))
          p3 = col((x1,y2,0))

          v1 = p1-p0
          v2 = p3-p0
          vcen = ((v2/2) + (v1/2)) + p0

          ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color='green', fill=False, hatch='/'))
          ax.annotate("%d"%(int(PANEL)), vcen[0:2])
        ax.set_xlim((0, 2000))
        ax.set_ylim((0, 2000))
        ax.set_ylim(ax.get_ylim()[::-1])
      else:
        from serialtbx.detector.xtc import basis_from_geo

        root = geometry.get_top_geo()
        root_basis = basis_from_geo(root)

        # get pixel mappings to real space.
        x, y, z = geometry.get_pixel_coords()
        assert x.shape == y.shape == z.shape
        if len(x.shape) == 6:
          x = x.reshape(x.shape[2:6])
          y = y.reshape(y.shape[2:6])
          z = z.reshape(z.shape[2:6])
        elif len(x.shape) == 5:
          x = x.reshape(x.shape[1:5])
          y = y.reshape(y.shape[1:5])
          z = z.reshape(z.shape[1:5])
        else:
          assert len(x.shape) == 4
        sensor_slow = x.shape[2]
        sensor_fast = x.shape[3]
        ori = matrix.col((0,0,0))
        while True:
          if len(root.get_list_of_children()) == 4 or len(root.get_list_of_children()) == 32:
            break
          assert len(root.get_list_of_children()) == 1, len(root.get_list_of_children())
          child = root.get_list_of_children()[0]
          child_basis = root_basis * basis_from_geo(child)

          arrow_end = child_basis * matrix.col((0,0,0))
          dx, dy, _ = arrow_end - ori
          if dx != 0 or dy != 0:
            ax.arrow(ori[0], ori[1], dx, dy, head_width=0.5, head_length=1.0, fc='k', ec='k', length_includes_head=True)
          ori = arrow_end
          root = child
          root_basis = child_basis

        if len(root.get_list_of_children()) == 4:
          for quad_id, quad in enumerate(root.get_list_of_children()):
            quad_basis = root_basis * basis_from_geo(quad)
            arrow_end = quad_basis * matrix.col((0,0,0))
            dx, dy, _ = arrow_end - ori
            if dx != 0 or dy != 0:
              ax.arrow(ori[0], ori[1], dx, dy, head_width=0.5, head_length=1.0, fc='k', ec='k', length_includes_head=True)
            arrow_start = arrow_end
            for sensor_id, sensor in enumerate(quad.get_list_of_children()):
              sensor_x, sensor_y, sensor_z = sensor.get_pixel_coords()
              transformed_x, transformed_y, transformed_z = quad.transform_geo_coord_arrays(sensor_x, sensor_y, sensor_z)

              sensor_basis = quad_basis * basis_from_geo(sensor)
              arrow_end = sensor_basis * matrix.col((0,0,0))
              dx, dy, _ = arrow_end - arrow_start
              ax.arrow(arrow_start[0], arrow_start[1], dx, dy, head_width=0.5, head_length=1.0, fc='k', ec='k', length_includes_head=True)

              p0 = col((x[quad_id,sensor_id,0,0]/1000,
                        y[quad_id,sensor_id,0,0]/1000))
              p1 = col((x[quad_id,sensor_id,sensor_slow-1,0]/1000,
                        y[quad_id,sensor_id,sensor_slow-1,0]/1000))
              p2 = col((x[quad_id,sensor_id,sensor_slow-1,sensor_fast-1]/1000,
                        y[quad_id,sensor_id,sensor_slow-1,sensor_fast-1]/1000))
              p3 = col((x[quad_id,sensor_id,0,sensor_fast-1]/1000,
                        y[quad_id,sensor_id,0,sensor_fast-1]/1000))

              v1 = p1-p0
              v2 = p3-p0
              vcen = ((v2/2) + (v1/2)) + p0

              ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color='green', fill=False, hatch='/'))
              ax.annotate("%d,%d"%(quad_id,sensor_id), vcen[0:2])
          ax.set_xlim((-100, 100))
          ax.set_ylim((-100, 100))
        else:
          for quad_id in range(4):
            for sensor_id, sensor in enumerate(root.get_list_of_children()[quad_id*8:(quad_id+1)*8]):
              sensor_x, sensor_y, sensor_z = sensor.get_pixel_coords()
              #transformed_x, transformed_y, transformed_z = quad.transform_geo_coord_arrays(sensor_x, sensor_y, sensor_z)

              #arrow_start = matrix.col((quad.x0/1000, quad.y0/1000))
              #arrow_end = matrix.col((transformed_x[0,0]/1000, transformed_y[0,0]/1000))
              #dx, dy = arrow_end - arrow_start
              #ax.arrow(quad.x0/1000, quad.y0/1000, dx, dy, head_width=0.05, head_length=0.1, fc='k', ec='k')

              p0 = col((x[0,quad_id*8+sensor_id,0,0]/1000,
                        y[0,quad_id*8+sensor_id,0,0]/1000))
              p1 = col((x[0,quad_id*8+sensor_id,sensor_slow-1,0]/1000,
                        y[0,quad_id*8+sensor_id,sensor_slow-1,0]/1000))
              p2 = col((x[0,quad_id*8+sensor_id,sensor_slow-1,sensor_fast-1]/1000,
                        y[0,quad_id*8+sensor_id,sensor_slow-1,sensor_fast-1]/1000))
              p3 = col((x[0,quad_id*8+sensor_id,0,sensor_fast-1]/1000,
                        y[0,quad_id*8+sensor_id,0,sensor_fast-1]/1000))

              v1 = p1-p0
              v2 = p3-p0
              vcen = ((v2/2) + (v1/2)) + p0

              ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color='green', fill=False, hatch='/'))
              ax.annotate("%d,%d"%(quad_id,sensor_id), vcen[0:2])

            ax.set_xlim((0, 200))
            ax.set_ylim((0, 200))

    else:
      for i, panel in enumerate(detector):
        size = panel.get_image_size()
        p0 = col(panel.get_pixel_lab_coord((0,0)))
        p1 = col(panel.get_pixel_lab_coord((size[0]-1,0)))
        p2 = col(panel.get_pixel_lab_coord((size[0]-1,size[1]-1)))
        p3 = col(panel.get_pixel_lab_coord((0,size[1]-1)))

        v1 = p1-p0
        v2 = p3-p0
        vcen = ((v2/2) + (v1/2)) + p0

        ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color='green', fill=False, hatch='/'))
        ax.annotate(i, vcen[0:2])

        ax.set_xlim((-100, 100))
        ax.set_ylim((-100, 100))

      if hasattr(detector, "hierarchy"):
        def draw_arrow(pg, start):
          o = pg.get_origin()
          if o[0:2] != start[0:2]:
            delta = col(o) - col(start)
            ax.arrow(start[0], start[1], delta[0], delta[1], head_width=0.05, head_length=0.1, fc='k', ec='k', length_includes_head=True)
          if hasattr(pg, 'children'):
            for child in pg:
              draw_arrow(child, o)
        draw_arrow(detector.hierarchy(), (0.0,0.0,0.0))

  elif params.metrology in xpp_active_areas:
    # Read the metrology from the XPP dictionary
    active_areas = xpp_active_areas[params.metrology]['active_areas']

    for asic_number, (y1, x1, y2, x2) in enumerate([(active_areas[(i*4)+0]+1,
                                                     active_areas[(i*4)+1]+1,
                                                     active_areas[(i*4)+2]-1,
                                                     active_areas[(i*4)+3]-1) for i in range(len(active_areas)//4)]):
      ax.add_patch(Rectangle((x1,y1), x2-x1, y2-y1, color="grey"))
      ax.annotate(asic_number, (x1+(x2-x1)/2,y1+(y2-y1)/2))

    ax.set_xlim((0, 2000))
    ax.set_ylim((0, 2000))
    ax.set_ylim(ax.get_ylim()[::-1])
  else:
    # Read the metrology from an LCLS calibration directory
    from xfel.cxi.cspad_ana.parse_calib import calib2sections
    sections = calib2sections(params.metrology)
    for q_id, quad in enumerate(sections):
      for s_id, sensor in enumerate(quad):
        for a_id, asic in enumerate(sensor.corners_asic()):
          y1, x1, y2, x2 = asic
          p0 = col((x1,y1))
          p1 = col((x2,y1))
          p2 = col((x2,y2))
          p3 = col((x1,y2))
          v1 = p1-p0
          v2 = p3-p0
          vcen = ((v2/2) + (v1/2)) + p0

          ax.add_patch(Polygon((p0[0:2],p1[0:2],p2[0:2],p3[0:2]), closed=True, color='green', fill=False, hatch='/'))
          ax.annotate(q_id*16+s_id*2+a_id, vcen[0:2])

    ax.set_xlim((0, 2000))
    ax.set_ylim((0, 2000))
    ax.set_ylim(ax.get_ylim()[::-1])

  plt.show()
