from __future__ import division
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
import sys


master_phil = libtbx.phil.parse("""
metrology = None
  .type = str
  .help = File with metrology information or XPP active area name
  .optional = False
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isfile(arg)) or arg in xpp_active_areas:
      user_phil.append(libtbx.phil.parse("""metrology=\"%s\"""" % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  if (params.metrology is None) :
    master_phil.show()
    raise Usage("metrology must be defined (either metrology=XXX, or the name alone).")

  fig = plt.figure()
  ax = fig.add_subplot(111, aspect='equal')

  if os.path.isfile(params.metrology):
    from dxtbx.format.Registry import Registry
    from scitbx.matrix import col
    try:
      reader = Registry.find(params.metrology)
    except IOError:
      reader = None

    if reader is None:
      # see if it's a SLAC geometry file
      from PSCalib.GeometryAccess import GeometryAccess
      from scitbx import matrix
      try:
        geometry = GeometryAccess(params.metrology)
      except Exception, e:
        raise Sorry("Can't parse this metrology file")

      root = geometry.get_top_geo()

      # get pixel mappings to real space.
      x, y, z = geometry.get_pixel_coords()
      assert x.shape == y.shape == z.shape
      if len(x.shape) == 6:
        x = x.reshape(x.shape[2:6])
        y = y.reshape(y.shape[2:6])
        z = z.reshape(z.shape[2:6])
      else:
        assert len(x.shape) == 4
      sensor_slow = x.shape[2]
      sensor_fast = x.shape[3]
      while True:
        if len(root.get_list_of_children()) == 4:
          break
        assert len(root.get_list_of_children()) == 1
        root = root.get_list_of_children()[0]
      for quad_id, quad in enumerate(root.get_list_of_children()):
        ax.arrow(0, 0, quad.x0/1000, quad.y0/1000, head_width=0.05, head_length=0.1, fc='k', ec='k')
        for sensor_id, sensor in enumerate(quad.get_list_of_children()):
          sensor_x, sensor_y, sensor_z = sensor.get_pixel_coords()
          transformed_x, transformed_y, transformed_z = quad.transform_geo_coord_arrays(sensor_x, sensor_y, sensor_z)

          arrow_start = matrix.col((quad.x0/1000, quad.y0/1000))
          arrow_end = matrix.col((transformed_x[0,0]/1000, transformed_y[0,0]/1000))
          dx, dy = arrow_end - arrow_start
          ax.arrow(quad.x0/1000, quad.y0/1000, dx, dy, head_width=0.05, head_length=0.1, fc='k', ec='k')

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
      img = reader(params.metrology)

      detector = img.get_detector()

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
  else:
    active_areas = xpp_active_areas[params.metrology]['active_areas']

    for asic_number, (y1, x1, y2, x2) in enumerate([(active_areas[(i*4)+0]+1,
                                                     active_areas[(i*4)+1]+1,
                                                     active_areas[(i*4)+2]-1,
                                                     active_areas[(i*4)+3]-1) for i in xrange(len(active_areas)//4)]):
      ax.add_patch(Rectangle((x1,y1), x2-x1, y2-y1, color="grey"))
      ax.annotate(asic_number, (x1+(x2-x1)/2,y1+(y2-y1)/2))

    ax.set_xlim((0, 2000))
    ax.set_ylim((0, 2000))
    ax.set_ylim(ax.get_ylim()[::-1])

  plt.show()
