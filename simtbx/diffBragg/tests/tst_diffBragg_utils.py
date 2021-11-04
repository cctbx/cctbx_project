from __future__ import division, print_function
import numpy as np
from simtbx.diffBragg import stage_two_utils
from simtbx.diffBragg.phil import phil_scope

PARAMS = phil_scope.extract()

import logging
import sys
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
format = "%(filename)s:%(funcName)s >> %(message)s"
F = logging.Formatter(format)
handler.setFormatter(fmt=F)
logger.setLevel(logging.INFO)
handler.setLevel(logging.INFO)
logger.addHandler(handler)


def test_regions():
    """divides the detector into square regions for doing gain correction refinement"""
    PARAMS.refiner.region_size = 60, 60
    det_shape = 1, 514, 1030
    regions = stage_two_utils.regionize_detector(det_shape, PARAMS.refiner.region_size)
    nregions = len(np.unique(regions))
    assert nregions == 162

    det_shape = 2, 514, 1030
    regions = stage_two_utils.regionize_detector(det_shape, PARAMS.refiner.region_size)
    nregions = len(np.unique(regions))
    assert nregions == 324
    logger.info("OK")


if __name__ == "__main__":
    test_regions()
    print("OK")
