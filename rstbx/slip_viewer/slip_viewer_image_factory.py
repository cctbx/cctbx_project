from __future__ import absolute_import, division, print_function
import iotbx.detectors

# store default ImageFactory function
defaultImageFactory = iotbx.detectors.ImageFactory

def SlipViewerImageFactory(filename):
  try:
    return NpyImageFactory(filename)
  except Exception as e:
    return defaultImageFactory(filename)

# Use the dxtbx class as it handles all possible variance of NPY images
def NpyImageFactory(filename):
  from dxtbx.format.FormatPYunspecified import FormatPYunspecified
  img = FormatPYunspecified(filename)
  return img.get_detectorbase()
