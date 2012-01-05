# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# The average_spots jiffy calculates average pickled images of either
# raw, pickled CsPad shots, or spotfinder results on raw, pickled
# CsPad shots.
#
# $Id$

import numpy, os, sys

from libtbx import easy_pickle
from iotbx               import detectors
from iotbx.detectors.npy import NpyImage


def ImageFactory(path):
  if os.path.isfile(path):
    I = NpyImage(path)
    I.readHeader()
    return (I)

detectors.ImageFactory = ImageFactory


from spotfinder.diffraction.imagefiles import FileName
FileName.exts.append("pickle")


# The img_add() function reads the image, distance, and beam energy
# from the pickled file whose name is the string pointed to by path,
# and adds them to img_sum, dist_sum, and nrg_sum.
def img_add(path, img_sum, dist_sum, nrg_sum, nmemb):
  img_cspad = easy_pickle.load(path)

  if (img_sum is None):
    img_sum    = img_cspad["image"].astype(numpy.uint32)
    dist_sum   = img_cspad["distance"]
    nrg_sum    = img_cspad["beamEnrg"]
    nmemb      = 1
  else:
    img_sum   += img_cspad["image"].astype(numpy.uint32)
    dist_sum  += img_cspad["distance"]
    nrg_sum   += img_cspad["beamEnrg"]
    nmemb     += 1

  return (img_sum, dist_sum, nrg_sum, nmemb)


# The spot_add() function reads the image, distance, and beam energy
# from the pickled file whose name is the string pointed to by path.
# The spots found in the image are added to spot_sum, while the
# distance and energy are added to dist_sum and nrg_sum respectively.
def spot_add(path, spot_sum, dist_sum, nrg_sum, nmemb):
  limits="""1479, 1515, 1672, 1699
            1281, 1515, 1474, 1699
            1092, 1506, 1249, 1699
            1065, 1506, 1090, 1699
             853, 1505, 1037, 1698
            1650, 1081, 1672, 1487
            1479, 1303, 1604, 1487
            1281, 1303, 1474, 1487
            1466, 1082, 1650, 1274
            1065, 1308, 1249, 1501
            1253, 1080, 1437, 1273
             853, 1307, 1037, 1500
             622, 1465,  815, 1649
             424, 1465,  617, 1649
             213, 1473,  397, 1666
            1048, 1098, 1241, 1282
             850, 1098, 1043, 1282
             623, 1252,  816, 1436
             425, 1252,  618, 1436
             213, 1275,  397, 1468
            1466,  883, 1650, 1076
            1506,  664, 1699,  848
            1253,  882, 1437, 1075
            1048,  885, 1241, 1069
            1308,  664, 1501,  848
            1099,  656, 1283,  849
             850,  885, 1043, 1069
             634, 1048,  818, 1241
             421, 1048,  605, 1241
             198, 1064,  391, 1248
            1506,  451, 1699,  635
            1308,  451, 1501,  635
            1099,  458, 1283,  651
            1514,  231, 1698,  424
            1085,  262, 1278,  446
            1514,   33, 1698,  226
             885,  656, 1069,  849
             885,  458, 1069,  651
             887,  262, 1080,  446
             634,  850,  818, 1043
             421,  850,  605, 1043
             656,  626,  849,  810
             656,  414,  849,  598
             664,  198,  848,  391
             664,    0,  848,  193
             458,  626,  651,  810
             198,  851,  391, 1035
             458,  414,  651,  598
             451,  198,  635,  391
             451,    0,  635,  193
               1, 1473,  185, 1666
               1, 1275,  185, 1468
               0, 1064,  193, 1248
               0,  851,  193, 1035
             263,  617,  447,  810
              50,  616,  234,  809
             263,  419,  447,  612
              50,  418,  234,  611
             231,  213,  424,  397
              33,  213,  226,  397"""

  args = [path,
          "distl.bins.verbose=True",
          "distl.minimum_spot_area=3",
          "distl.detector_tiling=%s" % limits.replace(" ", "").replace("\n", ","), # XXX should not need to reproduce the tiling on every single use
          "distl.peripheral_margin=1",
          "distl.res.outer=2.1"
          ]
  from spotfinder.command_line.signal_strength import run
  try:
    info = run(args)
  except e:
    return (spot_sum, dist_sum, nrg_sum, nmemb)

  img_cspad = easy_pickle.load(path)

  if (spot_sum is None):
    spot_sum   = numpy.zeros((img_cspad["image"].shape[0],
                              img_cspad["image"].shape[1]), dtype="float")
    dist_sum   = img_cspad["distance"]
    nrg_sum    = img_cspad["beamEnrg"]
    nmemb      = 0
  else:
    dist_sum  += img_cspad["distance"]
    nrg_sum   += img_cspad["beamEnrg"]
    nmemb     += 1

  # XXX Note that x and y are flipped in the summation.
  for spot in info.S.images[info.frames[0]]["spots_inlier"]:
    for i in xrange(len(spot.bodypixels)):
      spot_sum[spot.bodypixels[i].x, spot.bodypixels[i].y] += spot.wts[i]

  return (spot_sum, dist_sum, nrg_sum, nmemb)


# http://www.artima.com/weblogs/viewpost.jsp?thread=4829
def main(argv = None):
  if (argv is None):
    argv = sys.argv

  outpath = "average_prutt_00001.pickle" # XXX Should be argument!

  img_sum  = None
  dist_sum = 0
  nrg_sum  = 0
  nmemb    = 0
  for arg in argv[1:]: # XXX ugly hack!
    if (False): # XXX Should be argument?
      img_sum, dist_sum, nrg_sum, nmemb = img_add(
        arg, img_sum, dist_sum, nrg_sum, nmemb)
    else:
      img_sum, dist_sum, nrg_sum, nmemb = spot_add(
        arg, img_sum, dist_sum, nrg_sum, nmemb)
  if (nmemb == 0):
    return (0)

  # XXX Post-mortem--avoid overflows!  But breaks distance and energy!
  #nmemb = 1.0 * img_sum.max() / (2**14 - 16)

  easy_pickle.dump(outpath,
                   dict(beamEnrg = 1.0 / nmemb * nrg_sum,
                        distance = 1.0 / nmemb * dist_sum,
                        image    = 1.0 / nmemb * img_sum), # XXX implicit cast?
                   )
  print "Wrote average of %d images to '%s'" % (nmemb, outpath)
  return (0)


# Run with "phenix.python average_spots.py `awk '{print $2;}'
# section01.out`".
if (__name__ == "__main__"):
  sys.exit(main())
