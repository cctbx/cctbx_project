from __future__ import absolute_import, division, print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cctbx.integration_pickle_viewer
from cctbx.array_family import flex # implicit dependency
from matplotlib import pyplot as plt
from six.moves import cPickle as pickle

def get_CSPAD_active_areas(image, version_phil):
  from libtbx.phil import parse
  from iotbx.detectors.npy import NpyImage
  from spotfinder.applications.xfel.cxi_phil import cxi_basic_start
  data = pickle.load(open(image, "rb"))
  scope = parse(file_name=version_phil)
  basic_scope = cxi_basic_start()
  new_scope = basic_scope.phil_scope.fetch(source=scope)
  phil = new_scope.extract()
  img = NpyImage("dummy", source_data=data)
  img.readHeader(phil)
  tm = img.get_tile_manager(phil)
  return list(tm.effective_tiling_as_flex_int())

LG36_active_areas = [724, 448, 918, 633, 527, 448, 721, 633, 724, 655, 918, 840, 527, 655, 721, 840, 524, 23, 709, 217, 524, 220, 709, 414, 731, 24, 916, 218, 731, 221, 916, 415, 100, 230, 294, 415, 297, 230, 491, 415, 99, 20, 293, 205, 296, 20, 490, 205, 119, 447, 304, 641, 119, 644, 304, 838, 328, 450, 513, 644, 328, 647, 513, 841, 446, 856, 631, 1050, 446, 1053, 631, 1247, 654, 856, 839, 1050, 654, 1053, 839, 1247, 26, 1065, 220, 1250, 223, 1065, 417, 1250, 25, 858, 219, 1043, 222, 858, 416, 1043, 232, 1475, 417, 1669, 232, 1278, 417, 1472, 24, 1476, 209, 1670, 24, 1279, 209, 1473, 448, 1463, 642, 1648, 645, 1463, 839, 1648, 448, 1258, 642, 1443, 645, 1258, 839, 1443, 855, 1143, 1049, 1328, 1052, 1143, 1246, 1328, 855, 937, 1049, 1122, 1052, 937, 1246, 1122, 1060, 1558, 1245, 1752, 1060, 1361, 1245, 1555, 856, 1558, 1041, 1752, 856, 1361, 1041, 1555, 1473, 1359, 1667, 1544, 1276, 1359, 1470, 1544, 1482, 1566, 1676, 1751, 1285, 1566, 1479, 1751, 1463, 1131, 1648, 1325, 1463, 934, 1648, 1128, 1257, 1132, 1442, 1326, 1257, 935, 1442, 1129, 1142, 729, 1327, 923, 1142, 532, 1327, 726, 937, 729, 1122, 923, 937, 532, 1122, 726, 1558, 531, 1752, 716, 1361, 531, 1555, 716, 1557, 734, 1751, 919, 1360, 734, 1554, 919, 1359, 106, 1544, 300, 1359, 303, 1544, 497, 1573, 101, 1758, 295, 1573, 298, 1758, 492, 1133, 127, 1327, 312, 936, 127, 1130, 312, 1133, 335, 1327, 520, 936, 335, 1130, 520]

def plot_preds(pdata, active_areas=LG36_active_areas):
  try:
    preds = pdata['mapped_predictions'][0]
    preds_flat = preds.as_double()
    preds_fast = preds_flat[0::2]
    preds_slow = preds_flat[1::2]
  except KeyError:
    print("pickle may not be an integration pickle! skipping...")
    return
  plt.scatter(preds_slow, preds_fast, c='blue', marker='.')
  for i in range(64):
    aa = active_areas[4*i:4*i+4]
    plt.plot(
      [aa[1], aa[1], aa[3], aa[3], aa[1]],
      [aa[0], aa[2], aa[2], aa[0], aa[0]]
            )
  ax = plt.gca()
  ax.set_xlim(0, 1800)
  ax.set_ylim(0, 1800)
  ax.invert_yaxis()
  plt.axes().set_aspect('equal')
  plt.show()

if __name__ == "__main__":
  import sys
  from libtbx.utils import Sorry
  import libtbx.option_parser
  cmd_line = (libtbx.option_parser.option_parser(
    usage="%s [--detector_version_phil phil] [--image image] integration_pickle(s)" % libtbx.env.dispatcher_name)
    .option(None, "--detector_version_phil", "-d",
            type="string",
            default=None,
            dest="det_phil",
            help="detector version phil for the CSPAD")
    .option(None, "--image", "-i",
            type="string",
            default=None,
            dest="det_image",
            help="image matching the detector version phil")
    ).process(args=sys.argv[1:])
  if cmd_line.options.det_phil is not None and cmd_line.options.det_image is not None:
    print("extracting active areas...")
    active_areas = get_CSPAD_active_areas(
      cmd_line.options.det_image,
      cmd_line.options.det_phil)
  elif cmd_line.options.det_phil is None and cmd_line.options.det_image is None:
    print("using active areas from LG36 CSPAD metrology")
    active_areas = LG36_active_areas
  else:
    raise Sorry("Specify both a detector version phil and an example image to extract active areas.")
  for arg in cmd_line.args:
    file = open(arg, "rb")
    data = pickle.load(file)
    file.close()
    plot_preds(data, active_areas=active_areas)
