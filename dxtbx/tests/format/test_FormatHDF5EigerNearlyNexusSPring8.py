from __future__ import absolute_import, division, print_function

import os
import pytest

from dxtbx.format.FormatHDF5EigerNearlyNexusSPring8 import FormatHDF5EigerNearlyNexusSPring8
from dxtbx.datablock import DataBlockFactory
from dxtbx.model.goniometer import Goniometer

pytestmark = pytest.mark.skipif(
  not os.access('/dls/mx-scratch/rjgildea/zenodo/spring8-ccp4-2018/1443110/ccp4school2018_bl41xu', os.R_OK),
  reason='Test images not available')

def test_spring8_ccp4_2018_zenodo_1443110_data03():
  # https://zenodo.org/record/1443110#.XD8bD5ynzmE
  master_h5 = '/dls/mx-scratch/rjgildea/zenodo/spring8-ccp4-2018/1443110/ccp4school2018_bl41xu/05/data03/data03_master.h5'
  assert FormatHDF5EigerNearlyNexusSPring8.understand(master_h5)

  datablocks = DataBlockFactory.from_filenames([master_h5])
  imageset = datablocks[0].extract_imagesets()[0]
  assert imageset.get_format_class() == FormatHDF5EigerNearlyNexusSPring8

  detector = imageset.get_detector()
  gonio = imageset.get_goniometer()
  scan = imageset.get_scan()
  beam = imageset.get_beam()

  panel = detector[0]
  assert panel.get_pixel_size() == pytest.approx((0.075,0.075))
  assert panel.get_image_size() == (4150,4371)
  assert panel.get_trusted_range() == (-1,2.094707e+06) # XXX really?!
  assert panel.get_fast_axis() == (1,0,0)
  assert panel.get_slow_axis() == (0,-1,0)
  assert panel.get_origin() == pytest.approx(
    (-151.939,169.629,-180), abs=1e-3)
  assert panel.get_distance() == pytest.approx(180)

  assert isinstance(gonio, Goniometer)
  assert gonio.get_rotation_axis() == (-1,0,0) # XXX the rotation as understood by dxtbx is inverted
  assert gonio.get_fixed_rotation() == (1,0,0,0,1,0,0,0,1)
  assert gonio.get_setting_rotation() == (1,0,0,0,1,0,0,0,1)

  assert scan.get_oscillation() == pytest.approx((-10, 1))
  assert scan.get_image_range() == (1, 180)

  assert beam.get_wavelength() == pytest.approx(1.28241, abs=1e-5)
  assert beam.get_s0() == pytest.approx((0,0,-1/beam.get_wavelength()))
