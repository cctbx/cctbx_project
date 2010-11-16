from iotbx import mtz
from iotbx.scalepack import merge as scalepack_merge
from iotbx.scalepack import no_merge_original_index as scalepack_no_merge
from iotbx.cns import reflection_reader as cns_reflection_reader
from iotbx.cns import index_fobs_sigma_reader as cns_index_fobs_sigma_reader
from iotbx.dtrek import reflnlist_reader as dtrek_reflnlist_reader
from iotbx.shelx import hklf as shelx_hklf
from iotbx.shelx import crystal_symmetry_from_ins
from iotbx.xds.read_ascii import reader as xds_ascii_reader
from iotbx.solve.fpfm_reader import reader as solve_fpfm_reader
from iotbx.option_parser import option_parser
from cctbx import miller
from cctbx import crystal
from libtbx import easy_pickle, smart_open
from libtbx.utils import Sorry, detect_binary_file
import sys, os, os.path
import re

def unpickle_miller_arrays(file_name):
  result = easy_pickle.load(file_name)
  if (isinstance(result, miller.array)):
    return [result]
  result = list(result)
  for miller_array in result:
    if (not isinstance(miller_array, miller.array)):
      return None
  return result

def try_all_readers(file_name):
  try: content = mtz.object(file_name=file_name)
  except RuntimeError: pass
  else: return ("ccp4_mtz", content)
  if (detect_binary_file.from_initial_block(file_name=file_name)):
    try: content = unpickle_miller_arrays(file_name=file_name)
    except KeyboardInterrupt: raise
    except: pass
    else: return ("cctbx.miller.array", content)
    return (None, None)
  try: content = cns_reflection_reader.cns_reflection_file(
    open(file_name))
  except cns_reflection_reader.CNS_input_Error: pass
  else: return ("cns_reflection_file", content)
  try: content = cns_index_fobs_sigma_reader.reader(
    file_name=file_name)
  except RuntimeError: pass
  else: return ("cns_index_fobs_sigma", content)
  try: content = scalepack_merge.reader(
    open(file_name))
  except scalepack_merge.FormatError: pass
  else: return ("scalepack_merge", content)
  try: content = scalepack_no_merge.reader(file_name)
  except KeyboardInterrupt: raise
  except: pass
  else: return ("scalepack_no_merge_original_index", content)
  try: content = dtrek_reflnlist_reader.reflnlist(
    open(file_name))
  except KeyboardInterrupt: raise
  except: pass
  else: return ("dtrek_reflnlist", content)
  try: content = shelx_hklf.reader(
    open(file_name))
  except KeyboardInterrupt: raise
  except: pass
  else: return ("shelx_hklf", content)
  try: content = xds_ascii_reader(
    open(file_name))
  except KeyboardInterrupt: raise
  except: pass
  else: return ("xds_ascii", content)
  try: content = solve_fpfm_reader(file_name=file_name)
  except KeyboardInterrupt: raise
  except: pass
  else: return ("solve_fpfm", content)
  return (None, None)

class any_reflection_file(object):

  def __init__(self, file_name, ensure_read_access=True, strict=True):
    if (   file_name.startswith("amplitudes=")
        or file_name.startswith("hklf3=")
        or file_name.startswith("intensities=")
        or file_name.startswith("hklf4=")
        or file_name.startswith("hklf+ins/res=") ):
      self._observation_type, self._file_name = file_name.split("=", 1)
    elif (   file_name.endswith("=amplitudes")
          or file_name.endswith("=hklf3")
          or file_name.endswith("=intensities")
          or file_name.endswith("=hklf4")
          or file_name.endswith("=hklf+ins/res")):
      self._file_name, self._observation_type = file_name.split("=", 1)
    else:
      self._file_name = file_name
      self._observation_type = None
    if (self._observation_type == "hklf3"):
      self._observation_type = "amplitudes"
    elif (self._observation_type == "hklf4"):
      self._observation_type = "intensities"
    self._file_type = None
    file_name = self._file_name
    try:
      open(file_name) # test read access
    except IOError, e:
      if (ensure_read_access):
        raise Sorry(str(e))
      return
    if (self._observation_type is not None):
      try: self._file_content = shelx_hklf.reader(
        smart_open.for_reading(file_name),
        strict=strict)
      except KeyboardInterrupt: raise
      except:
        raise Sorry("Not a SHELX reflection file: %s\n"
          "  =%s can only be used for SHELX reflection files."
          % (file_name, self._observation_type))
      else: self._file_type = "shelx_hklf"
    else:
      self._file_type, self._file_content = try_all_readers(
        file_name=file_name)

  def file_name(self):
    return self._file_name

  def file_type(self):
    return self._file_type

  def file_content(self):
    return self._file_content

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    if (self._file_type is None):
      return []
    info_source = self._file_name
    if (info_source.startswith("./") or info_source.startswith(".\\")):
      info_source = info_source[2:]
    if (base_array_info is None):
      base_array_info = miller.array_info(
        source=info_source,
        source_type=self._file_type)
    if (self._file_type == "cctbx.miller.array"):
      result = []
      for miller_array in self._file_content:
        info = miller_array.info()
        if (info is None or not isinstance(info, miller.array_info)):
          info = base_array_info
        info.source = info_source
        info.crystal_symmetry_from_file = crystal.symmetry(
          unit_cell=miller_array.unit_cell(),
          space_group_info=miller_array.space_group_info())
        result.append(miller_array.customized_copy(
          crystal_symmetry=miller_array.join_symmetry(
            other_symmetry=crystal_symmetry,
            force=force_symmetry))
              .set_info(info)
              .set_observation_type(miller_array.observation_type()))
      return result
    if ((   crystal_symmetry is None
         or crystal_symmetry.unit_cell() is None)
        and self._observation_type == 'hklf+ins/res'
        ):
        name, ext = os.path.splitext(self._file_name)
        if ext != '.hkl': # it may be compressed: name.hkl.gz
          name, ext = os.path.splitext(name)
        for shelx_file_name in ('%s.ins' % name, '%s.res' % name):
          try:
            shelx_file = open(shelx_file_name)
            break
          except IOError:
            continue
        else:
          raise Sorry("Can't open files %s.ins or %s.res"
                      "required by the option hklf+ins/res" % ((name,)*2))
        crystal_symmetry = crystal_symmetry_from_ins.extract_from(
          file=shelx_file)
        shelx_file.seek(0)
        remaining = shelx_file.read()
        shelx_file.close()
        m = re.search("^HKLF\s*(\d)", remaining, re.X|re.M|re.S)
        if m is None:
          raise Sorry("%s does not contain the mandatory HKLF instruction"
                      % shelx_file.name)
        if m.group(1) == "4":
          self._observation_type = "intensities"
        elif m.group(1) == "3":
          self._observation_type = "amplitudes"
        else:
          raise Sorry("HKLF %s not supported" % m.group(1))
    result = self._file_content.as_miller_arrays(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      merge_equivalents=merge_equivalents,
      base_array_info=base_array_info)
    if (self.file_type() == "shelx_hklf"):
      if (self._observation_type == "intensities"):
        result[0].set_info(result[0].info().customized_copy(
          labels=["Iobs", "SigIobs"]))
        result[0].set_observation_type_xray_intensity()
      elif (self._observation_type == "amplitudes"):
        result[0].set_info(result[0].info().customized_copy(
          labels=["Fobs", "SigFobs"]))
        result[0].set_observation_type_xray_amplitude()
      else:
        raise Sorry("Unresolved amplitude/intensity ambiguity: %s\n"
          "  SHELX reflection files may contain amplitudes or intensities.\n"
          "  Please append   =amplitudes\n"
          "             or   =hklf3\n"
          "             or   =intensities\n"
          "             or   =hklf4\n"
          "  to the file name argument or parameter to resolve the"
          " ambiguity.\n"
          "  If a corresponding .ins file is available, look for the"
          " HKLF codeword.\n"
          "  Alternatively, run the phenix.reflection_statistics"
          " command twice,\n"
          "  once with =amplitudes and once with =intensities. Inspect"
          " the <I^2>/(<I>)^2\n"
          "  statistics. For acentric structures the values should"
          " fluctuate around\n"
          "  2.0, for centric structures around 3.0. If the statistics"
          " are not conclusive\n"
          "  it will be best to recover the original reflection data, such"
          " as SCALEPACK,\n"
          "  SCALA MTZ, XDS, or d*TREK files." % self._file_name)
    return result

def collect_arrays(
      file_names,
      crystal_symmetry,
      force_symmetry,
      merge_equivalents=True,
      discard_arrays=False,
      verbose=2,
      report_out=None):
  if (report_out is None):
    report_out = sys.stdout
  if (discard_arrays):
    result = None
  else:
    result = []
  for file_name in file_names:
    if (verbose > 0):
      print >> report_out
      print >> report_out, "file_name:", file_name
      report_out.flush()
    reflection_file = any_reflection_file(file_name)
    if (verbose > 0):
      print >> report_out, "file_type:", reflection_file.file_type()
      print >> report_out
    miller_arrays = reflection_file.as_miller_arrays(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      merge_equivalents=merge_equivalents)
    for miller_array in miller_arrays:
      if (verbose > 0):
        if (verbose > 1):
          miller_array.show_comprehensive_summary(f=report_out)
        else:
          miller_array.show_summary(f=report_out)
        if (verbose > 2):
          miller_array.show_array(f=report_out)
        print >> report_out
    if (not discard_arrays):
      result.extend(miller_arrays)
  return result

def run(args):
  command_line = (option_parser(
    usage="iotbx.reflection_file_reader [options] reflection_file ...",
    description="Example: iotbx.reflection_file_reader w1.sca w2.mtz w3.cns")
    .enable_symmetry_comprehensive()
    .option(None, "--weak_symmetry",
      action="store_true",
      default=False,
      help="symmetry on command line is weaker than symmetry found in files")
    .option(None, "--show_data",
      action="store_true",
      default=False,
      help="show Miller indices and data of all arrays")
    .option(None, "--pickle",
      action="store",
      type="string",
      help="write all data to FILE ('--pickle .' copies name of input file)",
      metavar="FILE")
  ).process(args=args)
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return
  if (command_line.options.show_data):
    verbose = 3
  else:
    verbose = 2
  all_miller_arrays = collect_arrays(
    file_names=command_line.args,
    crystal_symmetry=command_line.symmetry,
    force_symmetry=not command_line.options.weak_symmetry,
    discard_arrays=command_line.options.pickle is None,
    verbose=verbose,
    report_out=sys.stdout)
  if (all_miller_arrays is not None and len(all_miller_arrays) > 0):
    if (len(all_miller_arrays) == 1):
      all_miller_arrays = all_miller_arrays[0]
    pickle_file_name = command_line.options.pickle
    if (pickle_file_name == "."):
      if (len(command_line.args) > 1):
        raise Sorry(
          "Ambiguous name for pickle file (more than one input file).")
      pickle_file_name = os.path.basename(command_line.args[0])
      if (pickle_file_name.lower().endswith(".pickle")):
        raise Sorry("Input file is already a pickle file.")
    if (not pickle_file_name.lower().endswith(".pickle")):
      pickle_file_name += ".pickle"
    print
    print "Writing all Miller arrays to file:", pickle_file_name
    easy_pickle.dump(pickle_file_name, all_miller_arrays)
    print
