from __future__ import division
#!/usr/bin/env python
#
# dxtbx.serialize.goniometer.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def to_dict(gonio):
    ''' Convert the goniometer model to a dictionary

    Params:
        gonio The goniometer model

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    return OrderedDict([
        ('rotation_axis', gonio.get_rotation_axis()),
        ('fixed_rotation', gonio.get_fixed_rotation())])

def from_dict(d):
    ''' Convert the dictionary to a goniometer model

    Params:
        d The dictionary of parameters

    Returns:
        The goniometer model

    '''
    from dxtbx.model import Goniometer

    # If None, return None
    if d == None:
        return None

    # Create the model from the dictionary
    return Goniometer(tuple(d['rotation_axis']),
                      tuple(d['fixed_rotation']))
