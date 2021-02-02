from __future__ import absolute_import, division, print_function
from six.moves import range

class reader:
  """A class to read the INTEGRATE.HKL file used in XDS"""

  def __init__(self):
    """Initialise the file contents."""
    self._header = {}
    self.hkl = []
    self.iobs = []
    self.sigma = []
    self.xyzcal = []
    self.rlp = []
    self.peak = []
    self.corr = []
    self.maxc = []
    self.xyzobs = []
    self.alfbet0 = []
    self.alfbet1 = []
    self.psi = []
    self.iseg = []

  @staticmethod
  def is_integrate_hkl_file(filename):
    '''Check that the file is identified as an INTEGRATE.HKL

    Params:
      filename The path to the file

    Returns:
      True/False Is the file an INTEGRATE.HKL file

    '''
    with open(filename, 'r') as fh:
      return fh.read(26) == '!OUTPUT_FILE=INTEGRATE.HKL'

  def read_file(self, filename):
    """Read the INTEGRATE.HKL file.

    See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
    information about the file format.

    Params:
      filename The path to the file

    """

    # Check the file is an INTEGRATE.HKL file
    if not reader.is_integrate_hkl_file(filename):
      raise IOError("{0} is not an INTEGRATE.HKL file".format(filename))

    # Read the lines from the file
    with open(filename, 'r') as fh:
      lines = fh.readlines()

    # Loop through the lines in the file. First off, parse the header
    # lines until we reach !END_OF_HEADER. Then parse the data lines
    # until we read !END_OF_DATA
    in_header = True
    for l in lines:
      if in_header:
        if l.strip().startswith('!END_OF_HEADER'):
          in_header = False
          continue
        else:
          if not l.strip().startswith('!'):
            continue
          self._parse_header_line(l.strip()[1:])
      else:
        if l.strip().startswith('!END_OF_DATA'):
          break
        else:
          self._parse_data_line(l)

    # Set the header parameters
    self._set_header_parameters()

  def _parse_str(self, s):
    """Parse a string to either an int, float or string

    Params:
      s The input string

    Returns:
      The parsed value

    """
    try:
      return int(s)
    except ValueError:
      try:
        return float(s)
      except ValueError:
        return str(s)

  def _parse_value(self, value):
    """Parse the value or array of values contained in the string

    Params:
      value The value to parse

    Returns:
      The parsed value

    """
    values = value.split()
    if len(values) == 1:
      return self._parse_str(values[0])
    else:
      return tuple([self._parse_str(s) for s in values])

  def _set_header_parameters(self):
    """Get the parameters from the header dict

    """
    self.space_group       = self._header['SPACE_GROUP_NUMBER']
    self.unit_cell         = self._header['UNIT_CELL_CONSTANTS']
    self.detector_size     = (self._header['NX'], self._header['NY'])
    self.pixel_size        = (self._header['QX'], self._header['QY'])
    self.starting_frame    = self._header['STARTING_FRAME']
    self.starting_angle    = self._header['STARTING_ANGLE']
    self.oscillation_range = self._header['OSCILLATION_RANGE']
    self.rotation_axis     = self._header['ROTATION_AXIS']
    self.wavelength        = self._header['X-RAY_WAVELENGTH']
    self.beam_vector       = self._header['INCIDENT_BEAM_DIRECTION']
    self.detector_x_axis   = self._header['DIRECTION_OF_DETECTOR_X-AXIS']
    self.detector_y_axis   = self._header['DIRECTION_OF_DETECTOR_Y-AXIS']
    self.detector_origin   = (self._header['ORGX'], self._header['ORGY'])
    self.detector_distance = self._header['DETECTOR_DISTANCE']
    self.unit_cell_a_axis  = self._header['UNIT_CELL_A-AXIS']
    self.unit_cell_b_axis  = self._header['UNIT_CELL_B-AXIS']
    self.unit_cell_c_axis  = self._header['UNIT_CELL_C-AXIS']
    self.sigma_divergence  = self._header['BEAM_DIVERGENCE_E.S.D.']
    self.sigma_mosaicity   = self._header['REFLECTING_RANGE_E.S.D.']
    self.template          = self._header['NAME_TEMPLATE_OF_DATA_FRAMES']
    self.detector_type     = self._header['DETECTOR']
    self.minpk             = self._header['MINPK']
    self.cut               = self._header['CUT']
    self.variance_model    = self._header['VARIANCE_MODEL']
    del(self._header)

  def _parse_header_line(self, line):
    """Parse a line that has been identified as a header line

    Params:
      line The line to parse

    """
    name_value = line.split('=')
    if (len(name_value) < 2):
      return

    name = name_value[0]
    if (len(name_value) > 2):
      for i in range(1, len(name_value)-1):
        value_name = name_value[i].split()
        value = ''.join(value_name[:-1])
        self._header[name] = self._parse_value(value)
        name = value_name[-1]

    value = name_value[-1]
    self._header[name] = self._parse_value(value)

  def _parse_data_line(self, line):
    """Parse a data line from the Integrate.hkl file

    Params:
      line The line to parse

    """
    # Split the tokens
    tokens = line.split()
    tokens = [int(t) for t in tokens[0:3]] + [float(t) for t in tokens[3:]]

    # Get the reflection information and append to the lists
    self.hkl.append(tuple(tokens[0:3]))
    self.iobs.append(tokens[3])
    self.sigma.append(tokens[4])
    self.xyzcal.append(tuple(tokens[5:8]))
    self.rlp.append(tokens[8])
    self.peak.append(tokens[9])
    self.corr.append(tokens[10])
    self.maxc.append(tokens[11])
    self.xyzobs.append(tuple(tokens[12:15]))
    self.alfbet0.append(tuple(tokens[15:17]))
    self.alfbet1.append(tuple(tokens[17:19]))
    self.psi.append(tokens[19])
    if len(tokens) > 20:
      self.iseg.append(int(tokens[20]))

  def as_miller_arrays(self,
                       crystal_symmetry=None,
                       force_symmetry=False,
                       merge_equivalents=True,
                       base_array_info=None,
                       anomalous=None):
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="xds_integrate_hkl")
    from cctbx.array_family import flex
    from cctbx import crystal, miller, sgtbx
    crystal_symmetry = crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group_info=sgtbx.space_group_info(number=self.space_group))
    indices = flex.miller_index(self.hkl)
    miller_set = miller.set(crystal_symmetry, indices, anomalous_flag=anomalous)
    return (miller.array(
      miller_set, data=flex.double(self.iobs), sigmas=flex.double(self.sigma))
            .set_info(base_array_info.customized_copy(
              labels=["iobs", "sigma_iobs"])).set_observation_type_xray_intensity(),
            miller.array(miller_set, data=flex.vec3_double(self.xyzcal))
            .set_info(base_array_info.customized_copy(
              labels=["xyzcal"])),
            miller.array(miller_set, data=flex.vec3_double(self.xyzobs))
            .set_info(base_array_info.customized_copy(
              labels=["xyzobs"])),
            miller.array(miller_set, data=flex.double(self.rlp))
            .set_info(base_array_info.customized_copy(
              labels=["rlp"])),
            miller.array(miller_set, data=flex.double(self.peak))
            .set_info(base_array_info.customized_copy(
              labels=["peak"])),
            miller.array(miller_set, data=flex.double(self.corr))
            .set_info(base_array_info.customized_copy(
              labels=["corr"])),
            miller.array(miller_set, data=flex.double(self.maxc))
            .set_info(base_array_info.customized_copy(
              labels=["maxc"])),
            miller.array(miller_set, data=flex.vec2_double(self.alfbet0))
            .set_info(base_array_info.customized_copy(
              labels=["alfbet0"])),
            miller.array(miller_set, data=flex.vec2_double(self.alfbet1))
            .set_info(base_array_info.customized_copy(
              labels=["alfbet1"])),
            miller.array(miller_set, data=flex.double(self.psi))
            .set_info(base_array_info.customized_copy(
              labels=["psi"])))
