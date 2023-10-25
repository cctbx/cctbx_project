from __future__ import absolute_import, division, print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cxi.calibdir2cbfheader
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage, Sorry
from xfel.cftbx.detector.metrology2phil import metrology2phil, sections2phil
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf
from scitbx.matrix import col

master_phil = libtbx.phil.parse("""
metrology_dir = None
  .type = str
  .help = Directory with optical metrology information posistioning quadrants and sensors, corrected for
  .help = rectangularity
  .optional = False
corrections_phil = None
  .type = str
  .help = Phil file with quad/unit pixel translations and optionally subpixel metrology
  .optional = True
round_and_orthogonalize = True
  .type = str
  .help = Use if supplying a corrections phil. Quad and unit pixel translations assume an orthogonalized \
          and rounded detector, meaning the angles are multiples of 90 and the pixel values are ints. \
          If False, the optical metrology tilts and subpixel measurements are preserved when adding quad \
          and unit pixel corrections
out = None
  .type = str
  .help = Output file name
  .optional = False
""")

def get_asics_center(asics):
  center = col((0,0))
  for asic in asics:
    x1, y1, x2, y2 = asic
    center += col(((x2-x1)/2,(y2-y1)/2)) + col((x1,y1))
  return center / len(asics)

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isdir(arg)) :
      user_phil.append(libtbx.phil.parse("""metrology_dir=\"%s\" """ % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError as e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  if (params.metrology_dir is None) :
    master_phil.show()
    raise Usage("metrology_dir must be defined (either metrology_dir=XXX, or the directory path alone).")
  assert params.out is not None

  if params.corrections_phil is None:
    metro_phil = metrology2phil(params.metrology_dir, verbose=True)
    write_cspad_cbf(None, metro_phil, 'calibdir', None, params.out, None, 0, header_only=True)
  else:
    from xfel.cxi.cspad_ana.parse_calib import calib2sections
    sections = calib2sections(params.metrology_dir)

    from spotfinder.applications.xfel import cxi_phil
    horizons_phil = cxi_phil.cxi_versioned_extract()
    horizons_phil = horizons_phil.persist.phil_scope.fetch(source=libtbx.phil.parse(file_name=params.corrections_phil))
    corrections_params = horizons_phil.extract()

    bc = [0,0]
    for q in range(len(sections)):
      corner = sections[q][1].corners(True)[0]
      bc     = [bc[0] + corner[1] / len(sections), bc[1] + corner[0] / len(sections)]
    beam = col(bc)

    for itile in range(len(corrections_params.distl.tile_translations) // 4): # 128 tile_translations/4 = 32 sensors
      # order of quads in sections (UL,UR,LR,LL) is not the same as in quad_translations (UL,UR,LL,LR)
      sections_quad = itile//8
      phil_quad = [0,1,3,2].index(sections_quad)

      s = sections[sections_quad][itile%8] # not the same as iquad!
      c = list(s.center)
      c[0] += corrections_params.distl.quad_translations[2 * phil_quad + 0] + \
              corrections_params.distl.tile_translations[4 * itile + 0]
      c[1] += corrections_params.distl.quad_translations[2 * phil_quad + 1] + \
              corrections_params.distl.tile_translations[4 * itile + 1]
      s.center = tuple(c)

    if params.round_and_orthogonalize:
      # In order to match the image pickle metrology, we have to discard the tilts and offets in the
      # section objects, and instead use the active areas returned by corners_asic as the tile
      # locations.
      from serialtbx.detector import basis
      from xfel.cftbx.detector.cspad_cbf_tbx import pixel_size, asic_dimension, asic_gap
      null_ori = col((0,0,1)).axis_and_angle_as_unit_quaternion(0, deg=True)

      # the detector is rotated 90 degrees from how corners_asic reports things
      rot_ori = col((0,0,1)).axis_and_angle_as_unit_quaternion(-90, deg=True)
      metro = { (0,): basis(rot_ori, col((0,0,0))) } # basis dictionary used to build cbf header
      for quad_id, quad in enumerate(sections):
        quad_asics = []
        for sensor in quad:
          quad_asics.extend(sensor.corners_asic())
        quad_center = col(get_asics_center(quad_asics))
        v = (quad_center-beam)*pixel_size # vector from beam center to quad center
        metro[(0,quad_id)] = basis(null_ori, col((v[0],v[1],0)))
        for sensor_id, sensor in enumerate(quad):
          sensor_center = get_asics_center(sensor.corners_asic())
          # include sensor rotation, rounded to 90 degrees
          ori = col((0,0,1)).axis_and_angle_as_unit_quaternion(90.0 * round(sensor.angle / 90.0), deg=True)
          v = (sensor_center - quad_center)*pixel_size # vector from quad center to 2x1 center
          metro[(0,quad_id,sensor_id)] = basis(ori, col((v[0],v[1],0)))
          # add the two asics
          w = pixel_size * (asic_dimension[0]/2 + asic_gap/2)
          metro[(0,quad_id,sensor_id,0)] = basis(null_ori,col((-w,0,0)))
          metro[(0,quad_id,sensor_id,1)] = basis(null_ori,col((+w,0,0)))
      write_cspad_cbf(None, metro, 'cbf', None, params.out, None, 0, header_only=True)
    else:
      metro_phil = sections2phil(sections, verbose=True)
      write_cspad_cbf(None, metro_phil, 'calibdir', None, params.out, None, 0, header_only=True)
