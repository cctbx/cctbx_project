from __future__ import absolute_import, division, print_function

import pytest
import os
from dxtbx.datablock import DataBlockFactory
from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.datablock import BeamComparison, DetectorComparison, GoniometerComparison
from dxtbx.datablock import SweepDiff
from dials.util.options import tolerance_phil_scope


@pytest.mark.parametrize(
    "image_file",
    [
        "image_examples/ALS_831/q315r_lyso_001.img",
        #'image_examples/RigakuA200/XV10001.img',
        "image_examples/DLS_I02/X4_wide_M1S4_1_0001.cbf",
    ],
)
def test_cbf_writer(image_file, dials_regression, run_in_tmpdir):
    filename = os.path.join(dials_regression, image_file)
    datablock = DataBlockFactory.from_filenames([filename])[0]
    imageset = datablock.extract_imagesets()[0]

    FormatCBFMini.as_file(
        imageset.get_detector(),
        imageset.get_beam(),
        imageset.get_goniometer(),
        imageset.get_scan(),
        imageset.get_raw_data(0)[0],
        "image_0001.cbf",
    )

    datablock2 = DataBlockFactory.from_filenames(["image_0001.cbf"])[0]
    imageset2 = datablock2.extract_imagesets()[0]

    tolerance = tolerance_phil_scope.extract().tolerance

    diff = SweepDiff(tolerance)
    print("\n".join(diff(imageset, imageset2)))

    assert BeamComparison()(imageset.get_beam(), imageset2.get_beam())
    assert DetectorComparison(origin_tolerance=tolerance.detector.origin)(
        imageset.get_detector(), imageset2.get_detector()
    )
    assert GoniometerComparison()(imageset.get_goniometer(), imageset2.get_goniometer())
    s1 = imageset.get_scan()
    s2 = imageset.get_scan()
    assert s1.get_exposure_times() == s2.get_exposure_times()
    assert s1.get_oscillation() == s2.get_oscillation()
    assert s1.get_image_range() == s2.get_image_range()
    assert imageset.get_raw_data(0) == imageset2.get_raw_data(0)
