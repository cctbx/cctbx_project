from iotbx import crystal_symmetry_from_any
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from libtbx.optparse_wrapper import OptionParser, OptionError, make_option
import sys

class iotbx_option_parser(OptionParser):

  def __init__(self, usage=None, description=None):
    OptionParser.__init__(self, usage=usage, description=description)
    self.symmetry_callback = symmetry_callback()
    self.chunk_callback = chunk_callback()

  def format_help (self, formatter=None):
    if formatter is None:
      formatter = self.formatter
    result = []
    if self.usage:
      result.append(self.get_usage() + "\n")
    result.append(self.format_option_help(formatter))
    if self.description:
      result.append("\n")
      result.append(self.format_description(formatter) + "\n")
      result.append("\n")
    return "".join(result)

  def enable_unit_cell(self):
    self.add_option(make_option(None, "--unit_cell",
      action="callback",
      type="string",
      callback=self.symmetry_callback,
      help="External unit cell parameters",
      metavar="10,10,20,90,90,120|FILENAME"))
    return self

  def enable_space_group(self):
    self.add_option(make_option(None, "--space_group",
      action="callback",
      type="string",
      callback=self.symmetry_callback,
      help="External space group symbol",
      metavar="P212121|FILENAME"))
    return self

  def enable_symmetry(self):
    self.add_option(make_option(None, "--symmetry",
      action="callback",
      type="string",
      callback=self.symmetry_callback,
      help="External file with symmetry information",
      metavar="FILENAME"))
    return self

  def enable_symmetry_comprehensive(self):
    self.enable_unit_cell()
    self.enable_space_group()
    self.enable_symmetry()
    return self

  def enable_chunk(self):
    self.add_option(make_option(None, "--chunk",
      action="callback",
      type="string",
      callback=self.chunk_callback,
      help="Number of chunks for parallel execution and index for one process",
      metavar="n,i"))
    return self

  def process(self, args=None, nargs=None, min_nargs=None, max_nargs=None):
    assert nargs is None or (min_nargs is None and max_nargs is None)
    try:
      (options, args) = self.parse_args(args)
    except OptionError, e:
      print >> sys.stderr, e
      sys.exit(1)
    if (min_nargs is None): min_nargs = nargs
    if (min_nargs is not None):
      if (len(args) < min_nargs):
        self.error("Not enough arguments (at least %d required, %d given)." % (
          min_nargs, len(args)))
    if (max_nargs is None): max_nargs = nargs
    if (max_nargs is not None):
      if (len(args) > max_nargs):
        self.error("Too many arguments (at most %d allowed, %d given)." % (
          max_nargs, len(args)))
    return processed_options(self, options, args,
      self.symmetry_callback.get(),
      self.chunk_callback)

class processed_options:

  def __init__(self, parser, options, args, symmetry, chunk_callback):
    self.parser = parser
    self.options = options
    self.args = args
    self.symmetry = symmetry
    self.chunk_n = chunk_callback.n
    self.chunk_i = chunk_callback.i

class symmetry_callback:

  def __init__(self):
    self.unit_cell = None
    self.space_group_info = None

  def __call__(self, option, opt, value, parser):
    if (opt == "--unit_cell"):
      unit_cell = None
      try: unit_cell = uctbx.unit_cell(value)
      except: pass
      if (unit_cell is not None):
        self.unit_cell = unit_cell
      else:
        crystal_symmetry = crystal_symmetry_from_any.extract_from(value)
        if (   crystal_symmetry is None
            or crystal_symmetry.unit_cell() is None):
          raise OptionError("cannot read parameters: " + value, opt)
        self.unit_cell = crystal_symmetry.unit_cell()
    elif (opt == "--space_group"):
      space_group_info = None
      space_group_info = sgtbx.space_group_info(symbol=value)
      try: space_group_info = sgtbx.space_group_info(symbol=value)
      except: pass
      if (space_group_info is not None):
        self.space_group_info = space_group_info
      else:
        crystal_symmetry = crystal_symmetry_from_any.extract_from(value)
        if (   crystal_symmetry is None
            or crystal_symmetry.space_group_info() is None):
          raise OptionError("unknown space group: " + value, opt)
        self.space_group_info = crystal_symmetry.space_group_info()
    elif (opt == "--symmetry"):
      crystal_symmetry = crystal_symmetry_from_any.extract_from(value)
      if (   crystal_symmetry is None
          or crystal_symmetry.space_group_info() is None):
        raise OptionError("cannot read symmetry: " + value, opt)
      if (crystal_symmetry.unit_cell() is not None):
        self.unit_cell = crystal_symmetry.unit_cell()
      if (crystal_symmetry.space_group_info() is not None):
        self.space_group_info = crystal_symmetry.space_group_info()
    else:
      raise RuntimeError, "Programming error."

  def get(self):
    return crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group_info=self.space_group_info)

class chunk_callback:

  def __init__(self):
    self.n = 1
    self.i = 0

  def __call__(self, option, opt, value, parser):
    assert opt == "--chunk"
    try:
      self.n, self.i = [int(i) for i in value.split(",")]
    except:
      raise OptionError(
        "Two comma-separated positive integers required.",
        opt)
    if (self.n < 1):
      raise OptionError(
        "First integer (number of chunks) must be greater than 0 (%d given)."
        % self.n, opt)
    if (self.i < 0):
      raise OptionError(
        "Second integer (index of chunks) must be positive (%d given)."
        % self.i, opt)
    if (self.n < self.i):
      raise OptionError(
        ("First integer (number of chunks, %d given) must be greater"
        + " than second integer (index of chunks, %d given).")%(self.n,self.i),
        opt)
