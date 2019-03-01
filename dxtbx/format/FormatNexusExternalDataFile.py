#!/usr/bin/env python
#
# FormatNexusExternalDataFile.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division

from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.model import Beam  # import dependency
from dxtbx.model import Detector  # import dependency
from dxtbx.model import Goniometer  # import dependency
from dxtbx.model import Scan  # import dependency


def find_entries(nx_file):
    """
  Find NXmx entries

  """
    if "entry" in nx_file:
        entry = nx_file["entry"]
        if "NX_class" in entry.attrs.keys():
            if entry.attrs["NX_class"] == "NXentry":
                if "definition" not in entry.keys():
                    return entry
    return None


def is_nexus_external_data_file(filename):
    """
  A hacky function to check if this is a nexus file

  """
    import h5py

    # Get the file handle
    handle = h5py.File(filename, "r")

    # Find the NXmx entries
    entry = find_entries(handle)
    if entry is not None:
        if "instrument" not in entry:
            return True
    return False


class FormatNexusExternalDataFile(FormatHDF5):
    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatHDF5.__init__(self, image_file, **kwargs)

    @staticmethod
    def understand(image_file):
        try:
            is_nexus = is_nexus_external_data_file(image_file)
        except IOError:
            return False
        return is_nexus

    @classmethod
    def ignore(cls):
        return True


if __name__ == "__main__":
    import sys

    f = FormatNexusExternalDataFile(sys.argv[1])
