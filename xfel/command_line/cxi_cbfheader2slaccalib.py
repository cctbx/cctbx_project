from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.cbfheader2slaccalib
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

import sys, os, dxtbx, math
import libtbx.phil
from libtbx.utils import Usage, Sorry
from PSCalib.GeometryAccess import GeometryAccess
from scitbx.matrix import sqr, col

master_phil = libtbx.phil.parse("""
reference_metrology_file = None
  .type = str
  .help = File with optical metrology information posistioning quadrants and sensors.
  .help = Usually in the calib/geometry folder of the experiment, in the form of N-end.data
  .help = Used as a starting point
  .optional = False
cbf_header = None
  .type = str
  .help = Cbf file with metrology to be converted
  .optional = False
out_metrology_file = 0-end.data
  .type = str
  .help = File with optical metrology information posistioning quadrants and sensors.
  .help = Should be placed in the calib/geometry folder of the experiment, in the form of N-end.data
  .optional = False
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    try :
      user_phil.append(libtbx.phil.parse(arg))
    except RuntimeError, e :
      raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()

  if params.reference_metrology_file is None or \
     params.cbf_header is None or \
     params.out_metrology_file is None:
    raise Usage("%s reference_metrology_file=<path> cbf_header=<path> out_metrology_file=<path>"%libtbx.env.dispatcher_name)

  if not os.path.exists(params.reference_metrology_file):
    raise Sorry("File not found: %s"%params.reference_metrology_file)
  if not os.path.exists(params.cbf_header):
    raise Sorry("File not found: %s"%params.cbf_header)

  geometry = GeometryAccess(params.reference_metrology_file)
  top_geo = geometry.get_top_geo()

  img = dxtbx.load(params.cbf_header)
  detector = img.get_detector()

  assert len(top_geo.get_list_of_children()) == len(detector.hierarchy())

  for sq, dq in zip(top_geo.get_list_of_children(), detector.hierarchy()):
    x = col(dq.get_local_fast_axis())
    y = col(dq.get_local_slow_axis())
    z = x.cross(y)
    rot = sqr((x[0],y[0],z[0],
               x[1],y[1],z[1],
               x[2],y[2],z[2]))
    q0, q1, q2, q3 = rot.r3_rotation_matrix_as_unit_quaternion()

    # http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Conversion
    rotx = 180*math.atan2(2*(q0*q1+q2*q3),1-2*(q1**2+q2**2))/math.pi
    roty = 180*math.asin(2*(q0*q2-q3*q1))/math.pi
    rotz = 180*math.atan2(2*(q0*q3+q1*q2),1-2*(q2**2+q3**2))/math.pi
    sq.rot_x = round(rotx/90)*90
    sq.rot_y = round(roty/90)*90
    sq.rot_z = round(rotz/90)*90
    sq.tilt_x = rotx - sq.rot_x
    sq.tilt_y = roty - sq.rot_y
    sq.tilt_z = rotz - sq.rot_z

    ox, oy, oz = dq.get_local_origin()
    dq.x0 = ox * 1000
    dq.y0 = oy * 1000
    dq.z0 = oz * 1000

  geometry.save_pars_in_file(params.out_metrology_file)
