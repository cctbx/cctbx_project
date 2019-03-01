# Tests for the scan class, and its helper classes.

from __future__ import absolute_import, division, print_function

import os

from dxtbx.model.scan import ScanFactory
from dxtbx.model.scan_helpers import scan_helper_image_files
from dxtbx.model.scan_helpers import scan_helper_image_formats
import pytest


def test_helper_image_files():
    """Test the static methods in scan_helper_image_files."""

    helper = scan_helper_image_files()

    import dxtbx

    directory = os.path.join(os.path.split(dxtbx.__file__)[0], "tests")

    template = "image_###.dat"

    assert (
        len(scan_helper_image_files.template_directory_to_indices(template, directory))
        == 20
    )

    assert scan_helper_image_files.template_directory_index_to_image(
        template, directory, 1
    ) == os.path.join(directory, "image_001.dat")

    assert (
        scan_helper_image_files.template_index_to_image(template, 1) == "image_001.dat"
    )

    assert scan_helper_image_files.image_to_template_directory(
        os.path.join(directory, "image_001.dat")
    ) == (template, directory)

    assert scan_helper_image_files.image_to_index("image_001.dat") == 1

    assert scan_helper_image_files.image_to_template("image_001.dat") == "image_###.dat"

    assert scan_helper_image_files.image_to_index("image_6.8kev_1_001.cbf") == 1


def test_helper_image_formats():
    """Test the static methods and properties in scan_helper_image_formats."""

    assert scan_helper_image_formats.check_format(scan_helper_image_formats.FORMAT_CBF)
    # suspend this test pending further discussion -- Nick Sauter
    # assert(not(scan_helper_image_formats.check_format('CBF')))


def test_xScanFactory():
    """Test out the ScanFactory."""

    import dxtbx

    directory = os.path.join(os.path.split(dxtbx.__file__)[0], "tests")

    template = "image_###.dat"

    xscans = [
        ScanFactory.single(
            scan_helper_image_files.template_directory_index_to_image(
                template, directory, j + 1
            ),
            scan_helper_image_formats.FORMAT_CBF,
            1.0,
            18 + 0.5 * j,
            0.5,
            j,
        )
        for j in range(20)
    ]

    xscans.reverse()

    with pytest.raises(RuntimeError):
        print(sum(xscans[1:], xscans[0]))

    xscans.sort()

    sum(xscans[1:], xscans[0])

    a = ScanFactory.add(xscans[:10])
    b = ScanFactory.add(xscans[10:])

    a + b

    filename = scan_helper_image_files.template_directory_index_to_image(
        template, directory, 1
    )

    assert len(ScanFactory.search(filename)) == 20

    (a + b)[1:5]
    (a + b)[:10]
