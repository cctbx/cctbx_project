from __future__ import absolute_import, division, print_function

import os
import pytest

from dxtbx.format.FormatNexus import FormatNexus
from dxtbx.datablock import DataBlockFactory
from dxtbx.model.goniometer import Goniometer

pytestmark = pytest.mark.skipif(
  not os.access('/dls/mx/data/mx21314/mx21314-27/VMXi-AB0816/well_7/images', os.R_OK),
  reason='Test images not available')

def test_VMXi_rotation_scan():
  master_h5 = '/dls/mx/data/mx21314/mx21314-27/VMXi-AB0816/well_7/images/image_14364_master.h5'
  assert FormatNexus.understand(master_h5)

  datablocks = DataBlockFactory.from_filenames([master_h5])
  imageset = datablocks[0].extract_imagesets()[0]
  assert imageset.get_format_class() == FormatNexus

  detector = imageset.get_detector()
  gonio = imageset.get_goniometer()
  scan = imageset.get_scan()
  beam = imageset.get_beam()

  panel = detector[0]
  assert panel.get_pixel_size() == (0.075,0.075)
  assert panel.get_image_size() == (2068, 2162)
  assert panel.get_trusted_range() == (-1,4096) # XXX really?!
  assert panel.get_fast_axis() == (1,0,0)
  assert panel.get_slow_axis() == (0,-1,0)
  assert panel.get_origin() == pytest.approx(
    (-78.05999999999999, 87.03, -194.5039999999999))
  assert panel.get_distance() == pytest.approx(194.504)

  assert isinstance(gonio, Goniometer)
  assert gonio.get_rotation_axis() == (0,1,0)
  assert gonio.get_fixed_rotation() == (1,0,0,0,1,0,0,0,1)
  assert gonio.get_setting_rotation() == (1,0,0,0,1,0,0,0,1)

  assert scan.get_oscillation() == pytest.approx((-30,0.1))
  assert scan.get_image_range() == (1, 600)

  assert beam.get_wavelength() == pytest.approx(0.979492)
  assert beam.get_s0() == pytest.approx((0,0,-1), abs=5e-2)

