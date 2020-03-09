#!/usr/bin/env libtbx.python
#
# iotbx.xds.xparm.py
#
#   Copyright (C) 2013 Diamond Light Source, James Parkhurst & Richard Gildea
#
#   Class to read all the data from a (G)XPARM.XDS file
#
from __future__ import absolute_import, division, print_function

import io
import sys
import warnings

from libtbx import adopt_init_args
from six.moves import range
from six.moves import map

class reader(object):
  """A class to read the XPARM.XDS/GXPARM.XDS file used in XDS"""

  @staticmethod
  def find_version(filename):
    """Check the version if the given file is a (G)XPARM.XDS file.

    If the file contains exactly 11 lines and 42 tokens, it is the old style
    version 1 file. If the file starts with XPARM.XDS it is the new style
    version 2 file. If the file contains segment definitions then it is a
    version 3 file.

    Params:
      filename The (G)XPARM.XDS filename

    Returns:
      The version or None if the file is not recognised

    """

    # Check file contains 11 lines and 42 tokens
    with io.open(filename, 'r', encoding="ascii") as file_handle:
      tokens = []
      version = 1
      count = 0
      for count, line in enumerate(file_handle):
        line_tokens = line.split()
        if count == 0:
          if len(line_tokens) > 0 and line_tokens[0] == 'XPARM.XDS':
              version = 2
        if version == 1:
          if count+1 > 11:
            return None
        elif version == 2:
          if count+1 > 14:
            if len(line_tokens) == 5:
              version = 3
            else:
              return None
        elif version == 3:
          if len(line_tokens) not in (5, 9):
            return None
        tokens.extend(line_tokens)

      if version == 1:
        if count+1 != 11 or len(tokens) != 42:
          return None

    # Is a (G)XPARM.XDS file
    return version

  @staticmethod
  def is_xparm_file(filename, check_filename = True):
    """Check if the given file is a (G)XPARM.XDS file.

    Ensure it is named correctly and contains exactly 11 lines and 42
    tokens, otherwise return False.

    Params:
      filename The (G)XPARM.XDS filename

    Returns:
      True/False the file is a (G)XPARM.XDS file

    """
    try:
      return reader.find_version(filename) is not None
    except UnicodeDecodeError:
      return False

  def read_file(self, filename, check_filename = True):
    """Read the XPARM.XDS/GXPARAM.XDS file.

    See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
    information about the file format.

    Param:
      filename The path to the file

    """

    # Check version and read file
    version = reader.find_version(filename)
    if version is None:
      raise IOError("{} is not a (G)XPARM.XDS file".format(filename))

    with io.open(filename, 'r', encoding="ascii") as fh:
      tokens = [l.split() for l in fh.readlines()]

    # Parse the tokens
    if version == 1:
      self.parse_version_1_tokens(tokens)
    else:
      self.parse_version_2_tokens(tokens)

  def parse_version_1_tokens(self, tokens):
    """Parse the version 1 tokens

    Param:
      tokens The list of tokens

    """
    # Scan and goniometer stuff
    self.starting_frame    = int(tokens[0][0])
    self.starting_angle    = float(tokens[0][1])
    self.oscillation_range = float(tokens[0][2])
    self.rotation_axis     = tuple(map(float, tokens[0][3:6]))

    # Beam stuff
    self.wavelength        = float(tokens[1][0])
    self.beam_vector       = tuple(map(float, tokens[1][1:4]))

    # Detector stuff
    self.num_segments      = 0
    self.detector_size     = tuple(map(int, tokens[2][0:2]))
    self.pixel_size        = tuple(map(float, tokens[2][2:4]))
    self.detector_distance = float(tokens[3][0])
    self.detector_origin   = tuple(map(float, tokens[3][1:3]))
    self.detector_x_axis   = tuple(map(float, tokens[4]))
    self.detector_y_axis   = tuple(map(float, tokens[5]))
    self.detector_normal   = tuple(map(float, tokens[6]))

    # Crystal stuff
    self.space_group       = int(tokens[7][0])
    self.unit_cell         = tuple(map(float, tokens[7][1:7]))
    self.unit_cell_a_axis  = tuple(map(float, tokens[8]))
    self.unit_cell_b_axis  = tuple(map(float, tokens[9]))
    self.unit_cell_c_axis  = tuple(map(float, tokens[10]))

  def parse_version_2_tokens(self, tokens):
    """Parse the version 2 tokens

    Param:
      tokens The list of tokens

    """
    # Scan and goniometer stuff
    self.starting_frame    = int(tokens[1][0])
    self.starting_angle    = float(tokens[1][1])
    self.oscillation_range = float(tokens[1][2])
    self.rotation_axis     = tuple(map(float, tokens[1][3:6]))

    # Beam stuff
    self.wavelength        = float(tokens[2][0])
    self.beam_vector       = tuple(map(float, tokens[2][1:4]))

    # Crystal stuff
    self.space_group       = int(tokens[3][0])
    self.unit_cell         = tuple(map(float, tokens[3][1:7]))
    self.unit_cell_a_axis  = tuple(map(float, tokens[4]))
    self.unit_cell_b_axis  = tuple(map(float, tokens[5]))
    self.unit_cell_c_axis  = tuple(map(float, tokens[6]))

    # Detector stuff
    self.num_segments      = int(tokens[7][0])
    self.detector_size     = tuple(map(int, tokens[7][1:3]))
    self.pixel_size        = tuple(map(float, tokens[7][3:5]))
    self.detector_origin   = tuple(map(float, tokens[8][0:2]))
    self.detector_distance = float(tokens[8][2])
    self.detector_x_axis   = tuple(map(float, tokens[9]))
    self.detector_y_axis   = tuple(map(float, tokens[10]))
    self.detector_normal   = tuple(map(float, tokens[11]))

    # Loop through all the segments
    self.segments = []
    self.orientation = []
    for i in range(self.num_segments):
        self.segments.append(tuple(map(int, tokens[12+i*2])))
        self.orientation.append(tuple(map(float, tokens[13+i*2])))


class writer(object):

  def __init__(self,
               starting_frame,
               starting_angle,
               oscillation_range,
               rotation_axis,
               wavelength,
               beam_vector,
               space_group,
               unit_cell,
               unit_cell_a_axis,
               unit_cell_b_axis,
               unit_cell_c_axis,
               num_segments,
               detector_size,
               pixel_size,
               detector_origin,
               detector_distance,
               detector_x_axis,
               detector_y_axis,
               detector_normal,
               segments=None,
               orientation=None):
    adopt_init_args(self, locals())
    if [num_segments, segments, orientation].count(None) == 3:
      self.num_segments = 1
      self.segments = []
      self.orientation = []
      for i in range(self.num_segments):
        self.segments.append(
          (i+1, 1, self.detector_size[0], 1, self.detector_size[1]))
        self.orientation.append((0, 0, 0, 1, 0, 0, 0, 1, 0))
    warnings.warn("xparm.writer() is deprecated, use xparm.write() instead", DeprecationWarning, stacklevel=2)

  def show(self, out=None):
    """
    http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#XPARM.XDS
    """
    if out is None:
      out = sys.stdout
    print("XPARM.XDS", file=out)
    print("%6i %13.4f %9.4f" %(
      self.starting_frame, self.starting_angle, self.oscillation_range), end=' ', file=out)
    print("%9.6f %9.6f %9.6f" %(self.rotation_axis), file=out)
    print(" %14.6f" %self.wavelength, end=' ', file=out)
    print("%14.6f %14.6f %14.6f" %(self.beam_vector), file=out)
    print("   %3i" %(self.space_group), end=' ', file=out)
    print("%11.4f %11.4f %11.4f %7.3f %7.3f %7.3f" %self.unit_cell, file=out)
    print(" %14.6f  %14.6f  %14.6f" %self.unit_cell_a_axis, file=out)
    print(" %14.6f  %14.6f  %14.6f" %self.unit_cell_b_axis, file=out)
    print(" %14.6f  %14.6f  %14.6f" %self.unit_cell_c_axis, file=out)
    print(" %8i %9i %9i %11.6f %11.6f" %(
      self.num_segments, self.detector_size[0], self.detector_size[1],
      self.pixel_size[0], self.pixel_size[1]), file=out)
    print(" %14.6f %14.6f" %self.detector_origin, end=' ', file=out)
    print(" %14.6f" %self.detector_distance, file=out)
    print(" %14.6f %14.6f %14.6f" %self.detector_x_axis, file=out)
    print(" %14.6f %14.6f %14.6f" %self.detector_y_axis, file=out)
    print(" %14.6f %14.6f %14.6f" %self.detector_normal, file=out)
    for i in range(self.num_segments):
      print(" %9i %9i %9i %9i %9i" %tuple(self.segments[i]), file=out)
      print("".join([" %7.2f"*3] + [" %8.5f"]*6) %tuple(self.orientation[i]), file=out)

  def write_file(self, filename):
    with open(filename, 'w') as f:
      self.show(out=f)

# http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#XPARM.XDS
xparm_xds_template = """XPARM.XDS
{starting_frame:6d} {starting_angle:13.4f} {oscillation_range:9.4f} {rotation_axis[0]:9.6f} {rotation_axis[1]:9.6f} {rotation_axis[2]:9.6f}
 {wavelength:14.6f} {beam_vector[0]:14.6f} {beam_vector[1]:14.6f} {beam_vector[2]:14.6f}
{space_group:6d} {unit_cell[0]:11.4f} {unit_cell[1]:11.4f} {unit_cell[2]:11.4f} {unit_cell[3]:7.3f} {unit_cell[4]:7.3f} {unit_cell[5]:7.3f}
 {unit_cell_a_axis[0]:14.6f}  {unit_cell_a_axis[1]:14.6f}  {unit_cell_a_axis[2]:14.6f}
 {unit_cell_b_axis[0]:14.6f}  {unit_cell_b_axis[1]:14.6f}  {unit_cell_b_axis[2]:14.6f}
 {unit_cell_c_axis[0]:14.6f}  {unit_cell_c_axis[1]:14.6f}  {unit_cell_c_axis[2]:14.6f}
 {num_segments:8d} {detector_size[0]:9d} {detector_size[1]:9d} {pixel_size[0]:11.6f} {pixel_size[1]:11.6f}
 {detector_origin[0]:14.6f} {detector_origin[1]:14.6f}  {detector_distance:14.6f}
 {detector_x_axis[0]:14.6f} {detector_x_axis[1]:14.6f} {detector_x_axis[2]:14.6f}
 {detector_y_axis[0]:14.6f} {detector_y_axis[1]:14.6f} {detector_y_axis[2]:14.6f}
 {detector_normal[0]:14.6f} {detector_normal[1]:14.6f} {detector_normal[2]:14.6f}
"""
xparm_xds_segment_template = """
 (segments[{i}][0]:9d) (segments[{i}][1]:9d) (segments[{i}][2]:9d) (segments[{i}][3]:9d) (segments[{i}][4]:9d)
 (orientation[{i}][0]:7.2f) (orientation[{i}][1]:7.2f) (orientation[{i}][2]:7.2f) (orientation[{i}][3]:8.5f) (orientation[{i}][4]:8.5f) (orientation[{i}][5]:8.5f) (orientation[{i}][6]:8.5f) (orientation[{i}][7]:8.5f) (orientation[{i}][8]:8.5f)
""".lstrip("\n")


def write(
    starting_frame,
    starting_angle,
    oscillation_range,
    rotation_axis,
    wavelength,
    beam_vector,
    space_group,
    unit_cell,
    unit_cell_a_axis,
    unit_cell_b_axis,
    unit_cell_c_axis,
    num_segments,
    detector_size,
    pixel_size,
    detector_origin,
    detector_distance,
    detector_x_axis,
    detector_y_axis,
    detector_normal,
    segments=None,
    orientation=None,
):
    if num_segments is None and segments is None and orientation is None:
        num_segments = 1
        orientation = [(0, 0, 0, 1, 0, 0, 0, 1, 0)]
        segments = [(1, 1, detector_size[0], 1, detector_size[1])]

    template = xparm_xds_template
    for i in range(num_segments):
        template += (
            xparm_xds_segment_template.format(i=i).replace("(", "{").replace(")", "}")
        )

    return template.format(**locals())
