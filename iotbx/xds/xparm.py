#!/usr/bin/env libtbx.python
#
# iotbx.xds.xparm.py
#
#   James Parkhurst, Diamond Light Source, 2012/OCT/16
#
#   Class to read all the data from a (G)XPARM.XDS file
#
class reader:
    """A class to read the XPARM.XDS/GXPARM.XDS file used in XDS"""

    def __init__(self):
        pass

    @staticmethod
    def is_xparm_file(filename):
        """Check if the given file is a (G)XPARM.XDS file.

        Ensure it is named correctly and contains exactly 11 lines and 42
        tokens, otherwise return False.

        Params:
            filename The (G)XPARM.XDS filename

        Returns:
            True/False the file is a (G)XPARM.XDS file

        """
        import os

        # Check filename is (G)XPARM.XDS
        basename = os.path.basename(filename)
        if basename != 'XPARM.XDS' and basename != 'GXPARM.XDS':
            return False

        # Check file contains 11 lines and 42 tokens
        with open(filename, 'r') as file_handle:
            tokens = []
            for count, line in enumerate(file_handle):
                if count+1 > 11:
                    return False
                tokens.extend(line.split())

            if count+1 != 11 or len(tokens) != 42:
                print 3
                return False

        # Is a (G)XPARM.XDS file
        return True

    def _read_lines(self, filename):
        """Read the (G)XPARM.XDS file.

        If the file is not valid then return an IOError.

        Params:
            filename The path to the file

        Returns:
            The lines in the file as a list

        """
        if reader.is_xparm_file(filename):
            return open(filename, 'r').readlines()
        else:
            raise IOError("{0} is not a (G)XPARM.XDS file".format(filename))

    def read_file(self, filename):
        """Read the XPARM.XDS/GXPARAM.XDS file.

        See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
        information about the file format.

        Param:
            filename The path to the file

        """
        # Read the text from the file and split into an array of tokens
        tokens = [l.split() for l in self._read_lines(filename)]

        # Read the parameters from the list of tokens
        self.starting_frame    = int(tokens[0][0])
        self.starting_angle    = float(tokens[0][1])
        self.oscillation_range = float(tokens[0][2])
        self.rotation_axis     = tuple(map(float, tokens[0][3:6]))
        self.wavelength        = float(tokens[1][0])
        self.beam_vector       = tuple(map(float, tokens[1][1:4]))
        self.detector_size     = tuple(map(int, tokens[2][0:2]))
        self.pixel_size        = tuple(map(float, tokens[2][2:4]))
        self.detector_distance = float(tokens[3][0])
        self.detector_origin   = tuple(map(float, tokens[3][1:3]))
        self.detector_x_axis   = tuple(map(float, tokens[4]))
        self.detector_y_axis   = tuple(map(float, tokens[5]))
        self.detector_normal   = tuple(map(float, tokens[6]))
        self.space_group       = int(tokens[7][0])
        self.unit_cell         = tuple(map(float, tokens[7][1:7]))
        self.unit_cell_a_axis  = tuple(map(float, tokens[8]))
        self.unit_cell_b_axis  = tuple(map(float, tokens[9]))
        self.unit_cell_c_axis  = tuple(map(float, tokens[10]))
