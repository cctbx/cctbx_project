from __future__ import division
#!/usr/bin/env python
# FormatPilatusHelpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Helper methods for class for working with Pilatus images, for instance for
# identifying the regions to be masked.

def pilatus_6M_mask():
    '''Hard coded mask regions for a Pilatus 6M instrument.'''
    # FIX me, the paramters are listed here as f0, f1, s0, s1 but the prototype specifies f0, s0, f1, s1
    return [[488, 494, 1, 2527],
            [982, 988, 1, 2527],
            [1476, 1482, 1, 2527],
            [1970, 1976, 1, 2527],
            [1, 2463, 196, 212],
            [1, 2463, 408, 424],
            [1, 2463, 620, 636],
            [1, 2463, 832, 848],
            [1, 2463, 1044, 1060],
            [1, 2463, 1256, 1272],
            [1, 2463, 1468, 1484],
            [1, 2463, 1680, 1696],
            [1, 2463, 1892, 1908],
            [1, 2463, 2104, 2120],
            [1, 2463, 2316, 2332]]

def pilatus_2M_mask():
    '''Hard coded mask regions for a Pilatus 2M detector.'''

    return [[488, 494, 1, 1679],
            [982, 988, 1, 1679],
            [1, 1475, 196, 212],
            [1, 1475, 408, 424],
            [1, 1475, 620, 636],
            [1, 1475, 832, 848],
            [1, 1475, 1044, 1060],
            [1, 1475, 1256, 1272],
            [1, 1475, 1468, 1484]]

def pilatus_300K_mask():
    '''Hard coded mask regions for a Pilatus 300K instrument.'''

    return [[1, 487, 196, 212],
            [1, 487, 408, 424]]

def determine_pilatus_mask(xdetector):
    '''Return an appropriate pixel mask for a Pilatus detector.'''

    size = xdetector[0].get_image_size()

    if size == (2463, 2527):
        return pilatus_6M_mask()

    if size == (1475, 1679):
        return pilatus_2M_mask()

    if size == (487, 619):
        return pilatus_300K_mask()

    raise RuntimeError, 'unknown image size: %d %d' % size
