from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.pickle2cbf
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf
from libtbx import easy_pickle
from xfel.cxi.cspad_ana.parse_calib import calib2sections

master_phil = libtbx.phil.parse("""
pickle_file = None
  .type = str
  .help = Path to file(s) to convert
  .multiple = True
old_metrology = None
  .type = str
  .help = Path to original metrology this file was created with.  If none, use run 4.
new_metrology = None
  .type = str
  .help = File with optical metrology information posistioning quadrants and sensors or directory with calibration information
  .help = If none, use run 4.
detector = *CxiDs1 XppDs1
  .type = choice
  .optional = False
  .help = Specifiy CxiDs1 for the CXI Ds1 or DsD detectors which have relative coordinates for each quadrant,
  .help = or XppDs1 for XPP Ds1 detector which specifies absolute positions for each quadrant
plot = False
  .type = bool
  .help = show plots during processing
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isfile(arg)) :
      user_phil.append(libtbx.phil.parse("""pickle_file=\"%s\"""" % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  assert params.pickle_file is not None
  if len(params.pickle_file) == 0 :
    master_phil.show()
    raise Usage("pickle_file must be defined (either pickle_file=XXX, or the file path(s) alone).")
  assert params.detector is not None
  assert params.plot is not None

  if params.old_metrology is None:
    params.old_metrology=libtbx.env.find_in_repositories(
      "xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0")
  if os.path.isdir(params.old_metrology):
    sections = calib2sections(params.old_metrology)
  else:
    from xfel.cxi.cspad_ana.cspad_tbx import xpp_active_areas
    assert params.old_metrology in xpp_active_areas

  if params.new_metrology is None:
    params.new_metrology=libtbx.env.find_in_repositories(
      "xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0")
  assert os.path.isdir(params.new_metrology) or os.path.isfile(params.new_metrology)

  # Read the new metrology, either from a calibration directory or an optical metrology
  # flat file
  if os.path.isdir(params.new_metrology):
    metro_style = "calibdir"

    from xfel.cftbx.detector.metrology2phil import metrology2phil
    metro = metrology2phil(params.new_metrology,False)

    args = [
      "beam_center=(%f,%f)"%(img['BEAM_CENTER_X'],img['BEAM_CENTER_Y']),
      "timestamp=%s"%        img['TIMESTAMP'],
      ]

    import iotbx.phil
    for arg in args:
      metro = metro.fetch(sources=[iotbx.phil.parse(arg)])
  else:
    assert os.path.isfile(params.new_metrology)
    ext = os.path.splitext(params.new_metrology)[1].lower()
    if ext in ['.def','.cbf']:
      metro_style = "cbf"

      from xfel.cftbx.detector.cspad_cbf_tbx import cbf_file_to_basis_dict
      metro = cbf_file_to_basis_dict(params.new_metrology)

    else:
      metro_style = "flatfile"

      from xfel.cftbx.detector.cspad_cbf_tbx import read_optical_metrology_from_flat_file, asic_dimension, asic_gap

      metro = read_optical_metrology_from_flat_file(params.new_metrology, params.detector, img['PIXEL_SIZE'],
                                                    asic_dimension, asic_gap, plot=params.plot)

  for filename in params.pickle_file:
    # Read the pickle file and pull the tiles out of it
    img = easy_pickle.load(filename)

    tiles = {}
    data = img['DATA']

    if os.path.isdir(params.old_metrology):

      for p in xrange(len(sections)):
        for s in xrange(len(sections[p])):

          # Pull the sensor block from the image, and rotate it back to
          # the "lying down" convention.
          c = sections[p][s].corners_asic()
          k = (int(round(-sections[p][s].angle / 90.0)) + 1) % 4
          for a in xrange(2):
            asic = data.matrix_copy_block(
              i_row=c[a][0],
              i_column=c[a][1],
              n_rows=c[a][2] - c[a][0],
              n_columns=c[a][3] - c[a][1])
            tiles[(0, p, s, a)] = asic.matrix_rot90(k)

    else:
      active_areas = xpp_active_areas[params.old_metrology]['active_areas']
      rotations    = xpp_active_areas[params.old_metrology]['rotations']
      assert len(active_areas) // 4 == len(rotations) == 64

      active_areas = [(active_areas[(i*4)+0],
                       active_areas[(i*4)+1],
                       active_areas[(i*4)+2],
                       active_areas[(i*4)+3]) for i in xrange(64)]

      tile_id = 0

      for p in xrange(4):
        for s in xrange(8):
          for a in xrange(2):
            x1,y1,x2,y2 = active_areas[tile_id]
            block = data[x1:x2,y1:y2]
            tiles[(0, p, s, a)] = block.matrix_rot90(-rotations[tile_id])
            tile_id += 1

    # Write the cbf file
    destpath = os.path.splitext(filename)[0] + ".cbf"

    write_cspad_cbf(tiles, metro, metro_style, img['TIMESTAMP'], destpath,
                    img['WAVELENGTH'], img['DISTANCE'])
