from iotbx import crystal_symmetry_from_any
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from libtbx.option_parser import libtbx_option_parser, OptionError, make_option

class option_parser(libtbx_option_parser):

  def __init__(self, usage=None, description=None, more_help=None):
    libtbx_option_parser.__init__(self,
      usage=usage, description=description, more_help=more_help)
    self.symmetry_callback = symmetry_callback()

  def enable_unit_cell(self):
    self.add_option(make_option(None, "--unit_cell",
      action="callback",
      type="string",
      callback=self.symmetry_callback,
      help="External unit cell parameters",
      metavar="10,10,20,90,90,120|FILENAME"))
    self.symmetry_callback.is_enabled = True
    return self

  def enable_space_group(self):
    self.add_option(make_option(None, "--space_group",
      action="callback",
      type="string",
      callback=self.symmetry_callback,
      help="External space group symbol",
      metavar="P212121|FILENAME"))
    self.symmetry_callback.is_enabled = True
    return self

  def enable_symmetry(self):
    self.add_option(make_option(None, "--symmetry",
      action="callback",
      type="string",
      callback=self.symmetry_callback,
      help="External file with symmetry information",
      metavar="FILENAME"))
    self.symmetry_callback.is_enabled = True
    return self

  def enable_symmetry_comprehensive(self):
    self.enable_unit_cell()
    self.enable_space_group()
    self.enable_symmetry()
    return self

  def enable_resolution(self, default=None):
    self.add_option(make_option(None, "--resolution",
      action="store",
      default=default,
      type="float",
      help="High resolution limit (minimum d-spacing, d_min)",
      metavar="FLOAT"))
    return self

  def enable_low_resolution(self, default=None):
    self.add_option(make_option(None, "--low_resolution",
      action="store",
      default=default,
      type="float",
      help="Low resolution limit (maximum d-spacing, d_max)",
      metavar="FLOAT"))
    return self

  def enable_resolutions(self, default_low=None, default_high=None):
    self.enable_resolution(default=default_high)
    self.enable_low_resolution(default=default_low)
    return self

  def process(self, args=None, nargs=None, min_nargs=None, max_nargs=None):
    result = libtbx_option_parser.process(self,
      args=args, nargs=nargs, min_nargs=min_nargs, max_nargs=max_nargs)
    result.symmetry = self.symmetry_callback.get()
    return result

iotbx_option_parser = option_parser

class symmetry_callback(object):

  def __init__(self):
    self.is_enabled = False
    self.unit_cell = None
    self.space_group_info = None

  def __call__(self, option, opt, value, parser):
    if (opt == "--unit_cell"):
      unit_cell = None
      try: unit_cell = uctbx.unit_cell(value)
      except Exception: pass
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
      except Exception: pass
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
