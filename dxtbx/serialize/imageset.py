#!/usr/bin/env python
#
# dxtbx.serialize.imageset.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class NullSweep(object):
    ''' A class to facilitate reading in without access to images. '''
    template = None
    beam = None
    detector = None
    goniometer = None
    scan = None
    def __init__(self, template):
        self.template = template
    def get_template(self):
        return self.template
    def get_beam(self):
        return self.beam
    def set_beam(self, beam):
        self.beam = beam
    def get_detector(self):
        return self.detector
    def set_detector(self, detector):
        self.detector = detector
    def get_goniometer(self):
        return self.goniometer
    def set_goniometer(self, goniometer):
        self.goniometer = goniometer
    def get_scan(self):
        return self.scan
    def set_scan(self, scan):
        self.scan = scan


def filename_to_absolute(filename):
    ''' Convert filenames to absolute form. '''
    import os
    if isinstance(filename, list):
        return [os.path.abspath(f) for f in filename]

    return os.path.abspath(filename)

def basic_imageset_to_dict(imageset):
    ''' Convert an imageset to a dictionary

    Params:
        imageset The imageset

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    from dxtbx.serialize import beam, detector

    # Return the dictionary representation
    return OrderedDict([
        ("__id__", "imageset"),
        ("filenames", filename_to_absolute(imageset.paths())),
        ("beam", beam.to_dict(imageset.get_beam())),
        ("detector", detector.to_dict(imageset.get_detector()))])

def imagesweep_to_dict(sweep):
    ''' Convert a sweep to a dictionary

    Params:
        sweep The sweep

    Returns:
        A dictionary of the parameters

    '''
    from collections import OrderedDict
    from dxtbx.serialize import beam, detector, goniometer, scan

    # Return the dictionary representation
    return OrderedDict([
        ("__id__", "imageset"),
        ("template", filename_to_absolute(sweep.get_template())),
        ("beam", beam.to_dict(sweep.get_beam())),
        ("detector", detector.to_dict(sweep.get_detector())),
        ("goniometer", goniometer.to_dict(sweep.get_goniometer())),
        ("scan", scan.to_dict(sweep.get_scan()))])

def imageset_to_dict(imageset):
    ''' Convert the imageset to a dictionary

    Params:
        imageset The imageset

    Returns:
        A dictionary of the parameters

    '''
    from dxtbx.imageset import ImageSet, ImageSweep

    # If this is an imageset then return a list of filenames
    if isinstance(imageset, ImageSweep):
        return imagesweep_to_dict(imageset)
    elif isinstance(imageset, ImageSet):
        return basic_imageset_to_dict(imageset)
    elif isinstance(imageset, NullSweep):
        return imagesweep_to_dict(imageset)
    else:
        raise TypeError("Unknown ImageSet Type")

def basic_imageset_from_dict(d):
    ''' Construct an ImageSet class from the dictionary.'''
    from dxtbx.imageset import ImageSetFactory
    from dxtbx.serialize import beam, detector

    # Get the filename list and create the imageset
    filenames = map(str, d['filenames'])
    imageset = ImageSetFactory.new(filenames)[0]

    # Get the existing models as dictionaries
    beam_dict = beam.to_dict(sweep.get_beam())
    detector_dict = detector.to_dict(sweep.get_detector())

    # Set models
    imageset.set_beam(beam.from_dict(d.get('beam'), beam_dict))
    imageset.set_detector(detector.from_dict(d.get('detector'), detector_dict))

    # Return the imageset
    return imageset

def imagesweep_from_dict(d):
    '''Construct and image sweep from the dictionary.'''
    from dxtbx.imageset import ImageSetFactory
    from dxtbx.serialize import beam, detector, goniometer, scan

    # Get the template (required)
    template = str(d['template'])

    # If the scan isn't set, find all available files
    scan_dict = d.get('scan')
    if scan_dict is None:
        image_range = None
    else:
        image_range = scan_dict.get('image_range')

    # Construct the sweep
    try:
        sweep = ImageSetFactory.from_template(template, image_range)[0]
    except Exception:
        sweep = NullSweep(template)
        sweep.set_scan(scan.from_dict(d.get('scan')))

    # Get the existing models as dictionaries
    beam_dict = beam.to_dict(sweep.get_beam())
    gonio_dict = goniometer.to_dict(sweep.get_goniometer())
    detector_dict = detector.to_dict(sweep.get_detector())
    scan_dict = scan.to_dict(sweep.get_scan())

    # Set the models with the exisiting models as templates
    sweep.set_beam(beam.from_dict(d.get('beam'), beam_dict))
    sweep.set_goniometer(goniometer.from_dict(d.get('goniometer'), gonio_dict))
    sweep.set_detector(detector.from_dict(d.get('detector'), detector_dict))
#    sweep.set_scan(scan.from_dict(d.get('scan'), scan_dict))

    # Return the sweep
    return sweep

def imageset_from_dict(d):
    ''' Convert the dictionary to a sweep

    Params:
        d The dictionary of parameters

    Returns:
        The sweep

    '''
    # Check the input
    if d == None:
        return None

    # Check the version and id
    if str(d['__id__']) != "imageset":
        raise ValueError("\"__id__\" does not equal \"imageset\"")

    if "filenames" in d:
        return imageset_from_dict(d)
    elif "template" in d:
        return imagesweep_from_dict(d)
    else:
        raise TypeError("Unable to deserialize given imageset")
