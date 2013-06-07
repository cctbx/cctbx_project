#!/usr/bin/env python
#
# dxtbx.serialize.xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def to_imageset(input_filename, extra_filename=None):
    '''Get an image set from the xds input filename plus an extra filename

    Params:
        input_filename The XDS.INP file
        extra_filename A (G)XPARM.XDS, INTGRATE.HKL or XDS_ASCII.HKL file

    Returns:
        The imageset

    '''
    from iotbx.xds import xds_inp
    from dxtbx.imageset import ImageSetFactory
    import dxtbx

    # Read the input filename
    handle = xds_inp.reader()
    handle.read_file(input_filename)

    # Get the template
    template = handle.name_template_of_data_frames.replace('?', '#')

    # Create the imageset
    imageset = ImageSetFactory.from_template(template)[0]

    # If an extra filename has been specified, try to load models
    if extra_filename:
        models = dxtbx.load(extra_filename)
        imageset.set_beam(models.get_beam())
        imageset.set_detector(models.get_detector())
        imageset.set_goniometer(models.get_goniometer())

    # Return the imageset
    return imageset
