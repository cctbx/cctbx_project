from __future__ import absolute_import, division, print_function

import os
import pytest

from dxtbx.format.FormatHDF5EigerNearlyNexus import FormatHDF5EigerNearlyNexus
from dxtbx.datablock import DataBlockFactory
from dxtbx.model.goniometer import Goniometer

pytestmark = pytest.mark.skipif(
    not os.access("/dls/mx-scratch/rjgildea/zenodo", os.R_OK),
    reason="Test images not available",
)


def test_soleil_Proxima2A_zenodo_1443110_data03():
    # https://zenodo.org/record/1221344#.XEHr_5ynx2Q
    master_h5 = "/dls/mx-scratch/rjgildea/zenodo/1221344/200Hz/3_5_200Hz_1_master.h5"
    assert FormatHDF5EigerNearlyNexus.understand(master_h5)

    datablocks = DataBlockFactory.from_filenames([master_h5])
    imageset = datablocks[0].extract_imagesets()[0]
    assert imageset.get_format_class() == FormatHDF5EigerNearlyNexus

    detector = imageset.get_detector()
    gonio = imageset.get_goniometer()
    scan = imageset.get_scan()
    beam = imageset.get_beam()

    panel = detector[0]
    assert panel.get_pixel_size() == pytest.approx((0.075, 0.075))
    assert panel.get_image_size() == (3110, 3269)
    assert panel.get_trusted_range() == (-1, 12440)
    assert panel.get_fast_axis() == (1, 0, 0)
    assert panel.get_slow_axis() == (0, -1, 0)
    assert panel.get_thickness() == pytest.approx(0.45)
    assert panel.get_mu() == pytest.approx(3.96763)
    assert panel.get_material() == "Si"
    assert panel.get_origin() == pytest.approx((-120.556, 118.982, -134.255), abs=1e-3)
    assert panel.get_distance() == pytest.approx(134.255)

    assert isinstance(gonio, Goniometer)
    assert gonio.get_rotation_axis() == (1, 0, 0)
    assert gonio.get_fixed_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)
    assert gonio.get_setting_rotation() == (1, 0, 0, 0, 1, 0, 0, 0, 1)

    assert scan.get_oscillation() == pytest.approx((0, 0.5))
    assert scan.get_image_range() == (1, 800)

    assert beam.get_wavelength() == pytest.approx(0.980112, abs=1e-5)
    assert beam.get_s0() == pytest.approx((0, 0, -1 / beam.get_wavelength()))
