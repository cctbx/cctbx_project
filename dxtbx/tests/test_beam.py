from __future__ import absolute_import, division

#!/usr/bin/env python
# test_beam.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Tests for the beam class.

from dxtbx.model.beam import BeamFactory


def test_beam():
    """A test class for the beam class."""

    import libtbx.load_env
    import os

    dxtbx_dir = libtbx.env.dist_path("dxtbx")

    image = os.path.join(dxtbx_dir, "tests", "phi_scan_001.cbf")
    cbf = BeamFactory.imgCIF(image)


if __name__ == "__main__":

    test_beam()
