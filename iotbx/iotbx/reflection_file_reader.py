from iotbx import mtz
from iotbx.scalepack import reader as scalepack_reader
from iotbx.cns import reflection_reader as cns_reflection_reader
from iotbx import crystal_symmetry_from_any
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
import sys

class any_reflection_file:

  def __init__(self, file_name):
    self._file_name = file_name
    open(file_name) # test read access
    self._file_type = None
    if (self._file_type == None):
      try: self._file_content = mtz.Mtz(file_name)
      except RuntimeError: pass
      else: self._file_type = "ccp4_mtz"
    if (self._file_type == None):
      try: self._file_content = cns_reflection_reader.cns_reflection_file(
        open(file_name))
      except cns_reflection_reader.CNS_input_Error: pass
      else: self._file_type = "cns_reflection_file"
    if (self._file_type == None):
      try: self._file_content = scalepack_reader.scalepack_file(
        open(file_name))
      except scalepack_reader.ScalepackFormatError: pass
      else: self._file_type = "scalepack_merged"

  def file_name(self):
    return self._file_name

  def file_type(self):
    return self._file_type

  def file_content(self):
    return self._file_content

  def as_miller_arrays(self, crystal_symmetry=None, force_symmetry=00000):
    if (self.file_type() == None):
      return []
    info_prefix = self.file_name() + ":"
    if (info_prefix.startswith("./") or info_prefix.startswith(".\\")):
      info_prefix = info_prefix[2:]
    return self._file_content.as_miller_arrays(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      info_prefix=info_prefix)

def usage():
  return (  "usage: iotbx.any_reflection_file_reader.py"
          + " [--unit_cell=1,1,1,90,90,90]"
          + " [--space_group=P212121]"
          + " [--extract_symmetry=any_file_format]"
          + " [--force_symmetry]"
          + " any_reflection_file_format ...")

def run(args):
  unit_cell = None
  space_group_info = None
  force_symmetry = 00000
  remaining_args = []
  for arg in args:
    if (arg.startswith("--unit_cell=")):
      params = arg.split("=", 1)[1]
      unit_cell = uctbx.unit_cell(params)
    elif (arg.startswith("--space_group=")):
      symbol = arg.split("=", 1)[1]
      space_group_info = sgtbx.space_group_info(symbol=symbol)
    elif (arg.startswith("--extract_symmetry=")):
      file_name = arg.split("=", 1)[1]
      crystal_symmetry = crystal_symmetry_from_any.extract_from(file_name)
      unit_cell = crystal_symmetry.unit_cell()
      space_group_info = crystal_symmetry.space_group_info()
    elif (arg == "--force_symmetry"):
      force_symmetry = 0001
    elif (arg.startswith("--")):
      print usage()
      raise RuntimeError, "Unknown option: " + arg
    else:
      remaining_args.append(arg)
  args = remaining_args
  for file_name in args:
    print "file_name:", file_name
    sys.stdout.flush()
    reflection_file = any_reflection_file(file_name)
    print "file_type:", reflection_file.file_type()
    miller_arrays = reflection_file.as_miller_arrays(
      crystal_symmetry=crystal.symmetry(
        unit_cell=unit_cell,
        space_group_info=space_group_info),
      force_symmetry=force_symmetry)
    for miller_array in miller_arrays:
      miller_array.show_comprehensive_summary()
      print
    print
