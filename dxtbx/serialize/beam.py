#!/usr/bin/env python
#
# dxtbx.serialize.beam.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def to_dict(beam):
    ''' Convert the beam model to a dictionary

    Params:
        beam The beam model

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    return OrderedDict([
        ('direction', beam.get_direction()),
        ('wavelength', beam.get_wavelength()),
        ('divergence', beam.get_divergence()),
        ('sigma_divergence', beam.get_sigma_divergence())])

def from_dict(d):
    ''' Convert the dictionary to a beam model

    Params:
        d The dictionary of parameters

    Returns:
        The beam model

    '''
    from dxtbx.model import Beam

    # If None, return None
    if d == None:
        return None

    # Create the model from the dictionary
    return Beam(tuple(d['direction']),
                float(d['wavelength']),
                float(d['divergence']),
                float(d['sigma_divergence']))
