from __future__ import absolute_import, division, print_function
import os
from dxtbx.datablock import DataBlockFactory, DataBlockDumper

"""
Test deserializing a datablock that has single file indices while using check_format = True or False
"""

def get_indices(datablock):
  imageset = datablock.extract_imagesets()[0]
  return list(imageset.indices())

def test_split_single_image_datablock(dials_regression, tmpdir):
  tmpdir.chdir()
  sacla_file = os.path.join(dials_regression, "image_examples", "SACLA_MPCCD_Cheetah", "run266702-0-subset.h5")
  db = DataBlockFactory.from_filenames([sacla_file])[0]
  assert db.num_images() == 4
  imageset = db.extract_imagesets()[0]
  subset = imageset[2:3]
  subblock = DataBlockFactory.from_imageset(subset)[0]
  assert subblock.num_images() == 1
  assert get_indices(subblock) == [2]

  dumped_filename = "split_datablock.json"
  dump = DataBlockDumper(subblock)
  dump.as_json(dumped_filename)

  db = DataBlockFactory.from_json_file(dumped_filename, check_format = True)[0]
  assert db.num_images() == 1
  assert get_indices(db) == [2]

  db = DataBlockFactory.from_json_file(dumped_filename, check_format = False)[0]
  assert db.num_images() == 1
  assert get_indices(db) == [2]
