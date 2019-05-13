from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.cbfheader2slaccalib
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

import sys, os, dxtbx, getpass
import libtbx.phil
from libtbx.utils import Usage, Sorry
from PSCalib.GeometryAccess import GeometryAccess, GeometryObject
from scitbx.matrix import sqr, col
from xfel.cxi.cspad_ana.cspad_tbx import evt_timestamp

master_phil = libtbx.phil.parse("""
cbf_header = None
  .type = str
  .help = Cbf file with metrology to be converted
  .optional = False
out_metrology_file = 0-end.data
  .type = str
  .help = File with optical metrology information posistioning quadrants and sensors.
  .help = Should be placed in the calib/<Csapd version>/<detector address/geometry folder
  .help = of the experiment, in the form of N-end.data
  .optional = False
""")

class GeometryAccessFromCspadCBF(GeometryAccess):
  """ Override of GeometryAccess to read the metrology from a CSPAD CBF instead
      of a SLAC metrology file """

  def load_pars_from_file(self, path=None) :
    """ Use the dxtbx object model to build a GeometryObject hierarchy
        @param path Path to the CSPAD CBF file
    """
    if path is not None : self.path = path

    if self.pbits & 32 : print('Load file: %s' % self.path)

    img = dxtbx.load(self.path)
    cbf = img._cbf_handle
    cbf.find_category("diffrn_source")
    cbf.find_column("type")

    self.dict_of_comments = {
      "TITLE"     : "Geometry parameters of CSPAD",
      "DATE_TIME" : evt_timestamp(),
      "AUTHOR"    : getpass.getuser(),
      "EXPERIMENT": cbf.get_value(),
      "DETECTOR"  : "CSPAD",
      "CALIB_TYPE": "geometry",
      "COMMENT:01": "Table contains the list of geometry parameters for alignment of 2x1 sensors, quads, CSPAD, etc",
      "COMMENT:02": " translation and rotation pars of the object are defined w.r.t. parent object Cartesian frame",
      "PARAM:01"  : "PARENT     - name and version of the parent object",
      "PARAM:02"  : "PARENT_IND - index of the parent object",
      "PARAM:03"  : "OBJECT     - name and version of the object",
      "PARAM:04"  : "OBJECT_IND - index of the new object",
      "PARAM:05"  : "X0         - x-coordinate [um] of the object origin in the parent frame",
      "PARAM:06"  : "Y0         - y-coordinate [um] of the object origin in the parent frame",
      "PARAM:07"  : "Z0         - z-coordinate [um] of the object origin in the parent frame",
      "PARAM:08"  : "ROT_Z      - object design rotation angle [deg] around Z axis of the parent frame",
      "PARAM:09"  : "ROT_Y      - object design rotation angle [deg] around Y axis of the parent frame",
      "PARAM:10"  : "ROT_X      - object design rotation angle [deg] around X axis of the parent frame",
      "PARAM:11"  : "TILT_Z     - object tilt angle [deg] around Z axis of the parent frame",
      "PARAM:12"  : "TILT_Y     - object tilt angle [deg] around Y axis of the parent frame",
      "PARAM:13"  : "TILT_X     - object tilt angle [deg] around X axis of the parent frame"
    }

    self.list_of_geos = []

    detector = img.get_detector()
    hierarchy = detector.hierarchy()

    for q, quad in enumerate(hierarchy):
      for s, sensor in enumerate(quad):
        self.list_of_geos.append(self._load_geo(q,"QUAD:V1",s,"SENS2X1:V1",sensor))

    for q, quad in enumerate(hierarchy):
      self.list_of_geos.append(self._load_geo(0,"CSPAD:V1",q,"QUAD:V1",quad))

    # Add placeholder RAIL and IP vectors, including the XY component of the hierarchy's d0 vector
    go = self._load_geo(0,'RAIL',0,'CSPAD:V1',hierarchy)
    go.move_geo(0,0,-go.z0+1000000) # Placeholder
    self.list_of_geos.append(go)
    self.list_of_geos.append(self._null_geo(0,"IP",0,"RAIL"))

    self._set_relations()
    self.valid = True

  def _null_geo(self, pindex, pname, oindex, oname):
    """ Get a GeometryObject whose frameshift is zero
        @param pindex Index of the parent object
        @param pname  Name of the parent object
        @param oindex Index of the current object
        @param oname  Name of the current object

        @return an assembled GeometryObject
    """

    d = {
      'pname' :pname,
      'pindex':pindex,
      'oname' :oname,
      'oindex':oindex,
      'x0':    0,
      'y0':    0,
      'z0':    0,
      'rot_z': 0,
      'rot_y': 0,
      'rot_x': 0,
      'tilt_z':0,
      'tilt_y':0,
      'tilt_x':0
    }

    return GeometryObject(**d)

  def _load_geo(self, pindex, pname, oindex, oname, group):
    """ Given a dxtbx panel group, assemble the appropiate GeometryObject
        @param pindex Index of the parent object
        @param pname  Name of the parent object
        @param oindex Index of the current object
        @param oname  Name of the current object
        @param group  dxtbx panel group

        @return an assembled GeometryObject
    """
    x = col(group.get_local_fast_axis())
    y = col(group.get_local_slow_axis())
    z = x.cross(y)
    rot = sqr((x[0],y[0],z[0],
               x[1],y[1],z[1],
               x[2],y[2],z[2]))
    rotx, roty, rotz = rot.r3_rotation_matrix_as_x_y_z_angles(deg=True)

    # round to the nearest 90 degrees
    rrotx = round(rotx/90)*90
    rroty = round(roty/90)*90
    rrotz = round(rotz/90)*90

    ox, oy, oz = group.get_local_origin()

    d = {
      'pname': pname,
      'pindex':pindex,
      'oname': oname,
      'oindex':oindex,
      'x0':    ox * 1000,
      'y0':    oy * 1000,
      'z0':    oz * 1000,
      'rot_z': rrotz%360, # design
      'rot_y': rroty%360, # design
      'rot_x': rrotx%360, # design
      'tilt_z':rotz - rrotz,
      'tilt_y':roty - rroty,
      'tilt_x':rotx - rrotx
    }

    return GeometryObject(**d)

def run(args):
  if ("--help" in args or "-h" in args) :
    print("Write a SLAC metrology file from a CSPAD CBF. Parameters:")
    master_phil.show(attributes_level=2)
    return

  user_phil = []
  for arg in args:
    try :
      user_phil.append(libtbx.phil.parse(arg))
    except RuntimeError as e :
      raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()

  if params.cbf_header is None or \
     params.out_metrology_file is None:
    raise Usage("%s cbf_header=<path> out_metrology_file=<path>"%libtbx.env.dispatcher_name)

  if not os.path.exists(params.cbf_header):
    raise Sorry("File not found: %s"%params.cbf_header)

  print("Converting", params.cbf_header, "to", params.out_metrology_file)

  geometry = GeometryAccessFromCspadCBF(params.cbf_header)

  geometry.save_pars_in_file(params.out_metrology_file)

if (__name__ == "__main__") :
  run(sys.argv[1:])
