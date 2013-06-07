#!/usr/bin/env libtbx.python
#
# iotbx.xds.xds_inp.py
#
#   James Parkhurst, Diamond Light Source, 2012/OCT/16
#
#   Class to read all the data from a XDS.INP file
#
from __future__ import division

class reader:
  """A class to read the XDS.INP file used in XDS"""

  def __init__(self):
    pass

  @staticmethod
  def is_xds_inp_file(filename):
    """Check if the given file is an XDS.INP file.

    Params:
      filename The XDS.INP filename

    Returns:
      True/False the file is a XDS.INP file

    """
    import os
    return os.path.basename(filename) == 'XDS.INP'

  def read_file(self, filename, check_filename = True):
    """Read the XDS.INP file.

    See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
    information about the file format.

    Param:
      filename The path to the file

    """

    # Check and read file
    if reader.is_xds_inp_file(filename):
      lines = open(filename, 'r').readlines()
    else:
      raise IOError("{0} is not a XDS.INP file".format(filename))

    # Parse the tokens
    self.parse_lines(lines)

  def parse_lines(self, lines):
    """Parse the lines

    Param:
      tokens The list of lines

    """
    import re

    self.untrusted_rectangle = []
    for record in lines:

      record = record.strip()
      comp = [c for c in re.split(r' |=|\t', record) if c != '']

      if record.startswith('DETECTOR='):
        self.detector = comp[1]
        self.minimum_valid_pixel_value = int(comp[3])
        self.overload = int(comp[5])
        continue
      if record.startswith('CORRECTIONS'):
        self.corrections = comp[-1]
        continue
      if record.startswith('DIRECTION_OF_DETECTOR_X-AXIS='):
        self.direction_of_detector_x_axis = map(float, comp[-3:])
        continue
      if record.startswith('DIRECTION_OF_DETECTOR_Y-AXIS='):
        self.direction_of_detector_y_axis = map(float, comp[-3:])
        continue
      if record.startswith('TRUSTED_REGION='):
        self.trusted_region = map(float, comp[-2:])
        continue
      if record.startswith('SENSOR_THICKNESS='):
        self.sensor_thickness = float(comp[-1])
        continue

      if record.startswith('UNTRUSTED_RECTANGLE='):
        self.untrusted_rectangle.append(map(int, comp[-4:]))
        continue

      if record.startswith('MAXIMUM_NUMBER_OF_PROCESSORS='):
        self.maximum_number_of_processor = int(comp[-1])
        continue
      if record.startswith('NY='):
        self.nx = int(comp[1])
        self.ny = int(comp[3])
        self.px = float(comp[5])
        self.py = float(comp[7])
        continue
      if record.startswith('ORGX='):
        self.orgx = float(comp[1])
        self.orgy = float(comp[3])
        continue
      if record.startswith('ROTATION_AXIS='):
        self.rotation_axis = map(float, comp[-3:])
        continue
      if record.startswith('DETECTOR_DISTANCE='):
        self.detector_distance = float(comp[-1])
        continue
      if record.startswith('X-RAY_WAVELENGTH='):
        self.xray_wavelength = float(comp[-1])
        continue
      if record.startswith('INCIDENT_BEAM_DIRECTION='):
        self.incident_beam_direction = map(float, comp[-3:])
        continue
      if record.startswith('FRACTION_OF_POLARIZATION='):
        self.fraction_of_polarization = float(comp[-1])
        continue
      if record.startswith('POLARIZATION_PLANE_NORMAL='):
        self.polarization_plane_normal = map(float, comp[-3:])
        continue
      if record.startswith('FRIEDEL\'S_LAW='):
        self.friedels_law = bool(comp[-1])
        continue
      if record.startswith('NAME_TEMPLATE_OF_DATA_FRAMES='):
        self.name_template_of_data_frames = comp[-1]
        continue
      if record.startswith('STARTING_ANGLE='):
        self.starting_angle = float(comp[1])
        self.starting_frane = float(comp[3])
        continue
      if record.startswith('INCLUDE_RESOLUTION_RANGE='):
        self.include_resolution_range = map(float, comp[-2:])
        continue
      if record.startswith('UNIT_CELL_CONSTANTS='):
        self.unit_cell_constants = map(float, comp[-6:])
        continue
      if record.startswith('SPACE_GROUP_NUMBER='):
        self.space_group_number = int(comp[-1])
        continue
      if record.startswith('MAX_FAC_Rmeas='):
        self.max_fac_rmeas = float(comp[-1])
        continue
      if record.startswith('DATA_RANGE='):
        self.data_range = map(int, comp[-2])
        continue
