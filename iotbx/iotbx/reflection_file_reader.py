from iotbx import mtz
from iotbx.scalepack import merge as scalepack_merge
from iotbx.scalepack import no_merge_original_index as scalepack_no_merge
from iotbx.cns import reflection_reader as cns_reflection_reader
from iotbx.cns import index_fobs_sigma_reader as cns_index_fobs_sigma_reader
from iotbx.dtrek import reflnlist_reader as dtrek_reflnlist_reader
from iotbx.shelx import hklf as shelx_hklf
from iotbx.xds.read_ascii import reader as xds_ascii_reader
from iotbx.solve.fpfm_reader import reader as solve_fpfm_reader
from iotbx import crystal_symmetry_from_any
from iotbx.option_parser import iotbx_option_parser
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from scitbx.python_utils import easy_pickle
from libtbx.utils import UserError
import sys, os

class any_reflection_file:

  def __init__(self, file_name):
    self._file_name = file_name
    open(file_name) # test read access
    self._file_type = None
    if (self._file_type is None):
      try: self._file_content = mtz.Mtz(file_name)
      except RuntimeError: pass
      else: self._file_type = "ccp4_mtz"
    if (self._file_type is None):
      try: self._file_content = cns_reflection_reader.cns_reflection_file(
        open(file_name))
      except cns_reflection_reader.CNS_input_Error: pass
      else: self._file_type = "cns_reflection_file"
    if (self._file_type is None):
      try: self._file_content = cns_index_fobs_sigma_reader.reader(
        file_name=file_name)
      except RuntimeError: pass
      else: self._file_type = "cns_index_fobs_sigma"
    if (self._file_type is None):
      try: self._file_content = scalepack_merge.reader(
        open(file_name))
      except scalepack_merge.FormatError: pass
      else: self._file_type = "scalepack_merge"
    if (self._file_type is None):
      try: self._file_content = scalepack_no_merge.reader(file_name)
      except: pass
      else: self._file_type = "scalepack_no_merge_original_index"
    if (self._file_type is None):
      try: self._file_content = dtrek_reflnlist_reader.reflnlist(
        open(file_name))
      except: pass
      else: self._file_type = "dtrek_reflnlist"
    if (self._file_type is None):
      try: self._file_content = shelx_hklf.reader(
        open(file_name))
      except: pass
      else: self._file_type = "shelx_hklf"
    if (self._file_type is None):
      try: self._file_content = xds_ascii_reader(
        open(file_name))
      except: pass
      else: self._file_type = "xds_ascii"
    if (self._file_type is None):
      try: self._file_content = solve_fpfm_reader(file_name=file_name)
      except: pass
      else: self._file_type = "solve_fpfm"
    if (self._file_type is None):
      try: self._file_content = easy_pickle.load(file_name)
      except: pass
      else:
        if (isinstance(self._file_content, miller.array)):
          self._file_content = [self._file_content]
        else:
          miller_arrays = []
          try:
            for miller_array in self._file_content:
              if (isinstance(miller_array, miller.array)):
                miller_arrays.append(miller_array)
          except:
            pass
          else:
            if (len(miller_arrays) == 0):
              self._file_content = None
            else:
              self._file_content = miller_arrays
        if (self._file_content is not None):
          self._file_type = "cctbx.miller.array"

  def file_name(self):
    return self._file_name

  def file_type(self):
    return self._file_type

  def file_content(self):
    return self._file_content

  def as_miller_arrays(self, crystal_symmetry=None, force_symmetry=00000):
    if (self.file_type() is None):
      return []
    if (self.file_type() == "cctbx.miller.array"):
      return self.file_content()
    info_prefix = self.file_name() + ":"
    if (info_prefix.startswith("./") or info_prefix.startswith(".\\")):
      info_prefix = info_prefix[2:]
    return self._file_content.as_miller_arrays(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      info_prefix=info_prefix)

def run(args):
  command_line = (iotbx_option_parser(
    usage="iotbx.reflection_file_reader [options] reflection_file ...",
    description="Example: iotbx.reflection_file_reader w1.sca w2.mtz w3.cns")
    .enable_symmetry_comprehensive()
    .option(None, "--weak_symmetry",
      action="store_true",
      default=00000,
      dest="weak_symmetry",
      help="symmetry on command line is weaker than symmetry found in files")
    .option(None, "--pickle",
      action="store",
      type="string",
      dest="pickle",
      help="write all data to FILE ('--pickle .' copies name of input file)",
      metavar="FILE")
  ).process()
  if (command_line.options.pickle is None):
    all_miller_arrays = None
  else:
    all_miller_arrays = []
  for file_name in command_line.args:
    print
    print "file_name:", file_name
    sys.stdout.flush()
    reflection_file = any_reflection_file(file_name)
    print "file_type:", reflection_file.file_type()
    miller_arrays = reflection_file.as_miller_arrays(
      crystal_symmetry=command_line.symmetry,
      force_symmetry=not command_line.options.weak_symmetry)
    for miller_array in miller_arrays:
      print
      miller_array.show_comprehensive_summary()
    if (all_miller_arrays is not None):
      all_miller_arrays.extend(miller_arrays)
  if (all_miller_arrays is not None and len(all_miller_arrays) > 0):
    if (len(all_miller_arrays) == 1):
      all_miller_arrays = all_miller_arrays[0]
    pickle_file_name = command_line.options.pickle
    if (pickle_file_name == "."):
      if (len(command_line.args) > 1):
        raise UserError(
          "Ambiguous name for pickle file (more than one input file).")
      pickle_file_name = os.path.basename(command_line.args[0])
      if (pickle_file_name.lower().endswith(".pickle")):
        raise UserError("Input file is already a pickle file.")
    if (not pickle_file_name.lower().endswith(".pickle")):
      pickle_file_name += ".pickle"
    print
    print "Writing all Miller arrays to file:", pickle_file_name
    easy_pickle.dump(pickle_file_name, all_miller_arrays)
    print
