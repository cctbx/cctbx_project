#!/usr/bin/env python
#
# dxtbx.serialize.scan.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def to_dict(scan):
    ''' Convert the scan model to a dictionary

    Params:
        scan The scan model

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    return OrderedDict([
        ('image_range', scan.get_image_range()),
        ('oscillation', scan.get_oscillation()),
        ('exposure_time', scan.get_exposure_time()),
        ('epochs', list(scan.get_epochs()))])

def from_dict(d):
    ''' Convert the dictionary to a scan model

    Params:
        d The dictionary of parameters

    Returns:
        The scan model

    '''
    from dxtbx.model import Scan
    from scitbx.array_family import flex

    # If None, return None
    if d == None:
        return None

    # Create the model from the dictionary
    return Scan(tuple(d['image_range']),
                tuple(d['oscillation']),
                float(d['exposure_time']),
                flex.double(d['epochs']))
