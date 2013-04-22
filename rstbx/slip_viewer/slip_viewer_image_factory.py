from __future__ import division
import os
import iotbx.detectors

# store default ImageFactory function
defaultImageFactory = iotbx.detectors.ImageFactory

def SlipViewerImageFactory(filename):
  try:
    return NpyImageFactory(filename)
  except Exception, e:
    return defaultImageFactory(filename)

# The NpyImage requires some special treatment for creating an instance
def NpyImageFactory(filename):

  from iotbx.detectors.npy import NpyImage
  from spotfinder.applications.xfel import cxi_phil
  args = [filename,
          "distl.detector_format_version=CXI 7.1",
          "viewer.powder_arcs.show=False",
          "viewer.powder_arcs.code=3n9c",
         ]
  params = cxi_phil.cxi_versioned_extract(args)
  horizons_phil = params.persist.commands
  if isinstance(filename, basestring) and os.path.isfile(filename):
    I = NpyImage(filename)
  else:
    print "This is not a file; assume the data are in the defined dictionary format"
    I = NpyImage(filename, source_data=params.indexing.data)
  I.readHeader(horizons_phil)
  I.translate_tiles(horizons_phil)
  return I

