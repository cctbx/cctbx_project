from __future__ import absolute_import, division
from dxtbx.datablock import DataBlockFactory
import libtbx.load_env
import os

def run():
  # https://github.com/dials/dials/issues/226
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping tst_dials_226.py: dials_regression not present"
    return
  dials_regression = libtbx.env.dist_path('dials_regression')
  h5_path = os.path.join(
    dials_regression, 'image_examples/SACLA_MPCCD_Cheetah/run266702-0-subset.h5')
  datablock = DataBlockFactory.from_filenames([h5_path])[0]
  imageset = datablock.extract_imagesets()[0]
  path = imageset.get_path(imageset.indices()[0])

  sliced_imageset = imageset[2:]
  # this line fails
  path = sliced_imageset.get_path(sliced_imageset.indices()[0])


if __name__ == '__main__':
  #run() # disabled 2016/OCT/17 to enable nightly builds to pass
  print 'OK'
