from __future__ import division
#!/usr/bin/env python
#
# dxtbx.serialize.load.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def imageset_from_string(string):
    ''' Load the string and return the models.

    Params:
        string The JSON string to load

    Returns:
        The models

    '''
    import json
    from dxtbx.serialize.imageset import imageset_from_dict
    return imageset_from_dict(json.loads(string))

def imageset(infile):
    ''' Load the given JSON file.

    Params:
        infile The input filename or file object

    Returns:
        The models

    '''
    # If the input is a string then open and read from that file
    if isinstance(infile, str):
        with open(infile, 'r') as infile:
            return imageset_from_string(infile.read())

    # Otherwise assume the input is a file and read from it
    else:
        return imageset_from_string(infile.read())
