class error(RuntimeError):
  """ Base class of all errors
      self.args is a tuple starting with the line number and the line
  """

class illegal_command_or_atom_name_error(error):
  """ """

class illegal_continuation_line_error(error):
  """ """

class illegal_argument_error(error):
  """ self.args has an extra trailing item: the illegal argument """

class illegal_scatterer_error(error):
  """ """

class missing_sfac_error(error):
  """ """
