from iotbx import mtz
from iotbx.scalepack import reader as scalepack_reader
from iotbx.cns import reflection_reader as cns_reflection_reader
from iotbx import crystal_symmetry_from_any
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from scitbx.python_utils import dicts

class any_reflection_file:

  def __init__(self, file_name):
    self.file_name = file_name
    open(file_name, "r") # test read access
    self.file_type = None
    if (self.file_type == None):
      try: self.mtz_file = mtz.Mtz(file_name)
      except RuntimeError: pass
      else: self.file_type = "ccp4_mtz"
    if (self.file_type == None):
      try: self.cns_file = cns_reflection_reader.cns_reflection_file(
        open(file_name, "r"))
      except cns_reflection_reader.CNS_input_Error: pass
      else: self.file_type = "cns_reflection_file"
    if (self.file_type == None):
      try: self.scalepack_file = scalepack_reader.scalepack_file(
        open(file_name, "r"))
      except scalepack_reader.ScalepackFormatError: pass
      else: self.file_type = "scalepack_merged"

  def all_miller_arrays(self, crystal_symmetry=None,
                              unit_cell=None,
                              space_group_info=None):
    assert crystal_symmetry == None or (unit_cell == None and space_group_info == None)
    if (crystal_symmetry != None):
      unit_cell = crystal_symmetry.unit_cell()
      space_group_info = crystal_symmetry.space_group_info()
    result = dicts.easy()
    if (self.file_type == "scalepack_merged"):
      if (unit_cell == None):
        unit_cell = self.scalepack_file.unit_cell
      if (space_group_info == None):
        space_group_info = self.scalepack_file.space_group_info
      crystal_symmetry = crystal.symmetry(
        unit_cell=unit_cell,
        space_group_info=space_group_info)
      result[self.file_name+":f_obs,sigma"] \
      = self.scalepack_file.as_miller_array(
          info="Scalepack file: "+self.file_name,
          crystal_symmetry=crystal_symmetry).f_sq_as_f()
    elif (self.file_type == "cns_reflection_file"):
      crystal_symmetry = crystal.symmetry(
        unit_cell=unit_cell,
        space_group_info=space_group_info)
      miller_arrays = self.cns_file.as_miller_arrays(
        crystal_symmetry=crystal_symmetry)
      for key,miller_array in miller_arrays.items():
        result[self.file_name+":"+key] = miller_array
    elif (self.file_type == "ccp4_mtz"):
      pass
    else:
      raise RuntimeError, "Internal error."
    return result

def usage():
  return (  "usage: iotbx.any_reflection_file_reader.py"
          + " [--unit_cell=1,1,1,90,90,90]"
          + " [--space_group=P212121]"
          + " [--extract_symmetry=any_file_format]"
          + " any_reflection_file_format ...")

def run(args):
  unit_cell = None
  space_group_info = None
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
    elif (arg.startswith("--")):
      print usage()
      raise RuntimeError, "Unknown option: " + arg
    else:
      remaining_args.append(arg)
  args = remaining_args
  for file_name in args:
    print "file_name:", file_name
    reflection_file = any_reflection_file(file_name)
    print "file_type:", reflection_file.file_type
    miller_arrays = reflection_file.all_miller_arrays(
      unit_cell=unit_cell,
      space_group_info=space_group_info)
    for key,miller_array in miller_arrays.items():
      print "Key:", key
      miller_array.show_comprehensive_summary()
      print
