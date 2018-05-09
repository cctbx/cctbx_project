from __future__ import absolute_import, division, print_function

import os
import pytest

from dxtbx.datablock import DataBlockFactory

@pytest.mark.skip("disabled by GW 2016/OCT/17 to enable nightly builds to pass")
def test(dials_regression):
  # https://github.com/dials/dials/issues/226
  h5_path = os.path.join(
    dials_regression, 'image_examples/SACLA_MPCCD_Cheetah/run266702-0-subset.h5')
  datablock = DataBlockFactory.from_filenames([h5_path])[0]
  imageset = datablock.extract_imagesets()[0]
  path = imageset.get_path(imageset.indices()[0])

  sliced_imageset = imageset[2:]
  # this line fails
  path = sliced_imageset.get_path(sliced_imageset.indices()[0])
