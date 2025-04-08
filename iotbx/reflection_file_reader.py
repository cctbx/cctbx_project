"""
This module provides a generic frontend to all of the reflection file formats
supported in ``iotbx``.  Note that this module can also be used indirectly via
the even more generic :py:mod:`iotbx.file_reader` module, which provides
a unified API for reading in any file (but calls
:py:class:`iotbx.reflection_file_reader.any_reflection_file` internally).
Currently, the supported formats include:

- **CIF**: Crystallographic Information Format, the common syntax for
  specifying most structured data encountered in crystallography (but more
  widely used in small-molecule versus macromolecular crystallography), usually
  as ASCII text.  May encapsulate a variety of other data types, but only
  reflection data (of any type) is handled by this particular module.  Uses
  :py:mod:`iotbx.cif` internally.
- **MTZ**: Binary file format established by `CCP4 <http://www.ccp4.ac.uk>`_
  capable of storing any numerical data, and used by most major macromolecular
  crystallography software packages.
  Because of its speed and broad compatibility, this is the primary interchange
  format for reflection data in Phenix.  Uses :py:mod:`iotbx.mtz` internally.
- **Scalepack**: Fixed-format ASCII text produced by the program of the same
  name and the HKL2000 graphical interface.  This is actually two separate
  formats: one for merged intensities (with or without Friedel mates), another
  for unmerged intensities and associated processing metadata.  Uses either
  :py:mod:`iotbx.scalepack.merge` or :py:mod:`iotbx.scalepack.no_merge_original_index` internally.
- **CNS**: ASCII format, not as flexible as MTZ or CIF but able to store
  either amplitudes or intensities, R-free flags, and Hendrickson-Lattman
  coefficients.  Uses :py:mod:`iotbx.cns.reflection_reader` internally.
- **SHELX**: Fixed-format ASCII used by the eponymous software suite.  This
  format has significan disadvantages, discussed below.  Uses
  :py:mod:`iotbx.shelx.hklf` internally.
- **XDS**: ASCII format for processed intensities (both merged and unmerged).
  Uses :py:mod:`iotbx.xds.read_ascii` internally.
- **D*Trek**: ASCII format produced by software sold by Rigaku.

Independently, the :py:class:`cctbx.miller.array` class defines output methods
for CIF, MTZ, CNS, SHELX, and unmerged Scalepack files (although only the first
two are recommended for routine use).

Note that the underlying formats do not always contain complete information
about the crystal or even the data type.  SHELX format is especially
problematic as it not only omits crystal symmetry, but the same format may be
used to store either amplitudes or intensities, without any distinguishing
features.  As a crude workaround for the latter problem, the data type may be
specified as part of the file name::

  hkl_file = any_reflection_file("data.hkl=hklf4")
  hkl_file = any_reflection_file("data.hkl=intensities")

Other formats (CNS, unmerged Scalepack) may have incomplete or missing crystal
symmetry.  MTZ, XDS, and (usually) CIF files will be more complete.
"""

from __future__ import absolute_import, division, print_function
from iotbx import mtz
from iotbx.scalepack import merge as scalepack_merge
from iotbx.scalepack import no_merge_original_index as scalepack_no_merge
from iotbx.cif import reader as cif_reader
from iotbx.cns import reflection_reader as cns_reflection_reader
from iotbx.cns import index_fobs_sigma_reader as cns_index_fobs_sigma_reader
from iotbx.dtrek import reflnlist_reader as dtrek_reflnlist_reader
from iotbx.shelx import hklf as shelx_hklf
from iotbx.shelx import crystal_symmetry_from_ins
from iotbx.xds.read_ascii import reader as xds_ascii_reader
from iotbx.xds.integrate_hkl import reader as xds_integrate_hkl_reader
from iotbx.solve.fpfm_reader import reader as solve_fpfm_reader
from iotbx.option_parser import option_parser
from cctbx import miller
from cctbx import crystal
from libtbx import easy_pickle
from libtbx.utils import Sorry, detect_binary_file
import sys, os, os.path
import re

def unpickle_miller_arrays(file_name):
  result = easy_pickle.load(file_name)
  # Python 3 pickle fix
  # =========================================================================
  if sys.version_info.major == 3:
    result = easy_pickle.fix_py2_pickle(result)
  # =========================================================================
  if (isinstance(result, miller.array)):
    return [result]
  result = list(result)
  for miller_array in result:
    if (not isinstance(miller_array, miller.array)):
      return None
  return result

def _cif_prefilter(file_name):
  f_root, f_ext = os.path.splitext(file_name)
  if f_ext.lower() == '.gz': f_root, f_ext = os.path.splitext(f_root)
  if f_ext.lower() in ['.cif', '.mmcif', '.dic']: return True
  with open(file_name) as f:
    for l in f:
      if l.strip().startswith('#'): continue
      if not l.strip(): continue
      return l.strip().lower().startswith('data_')
  return False

def try_all_readers(file_name):
  try: content = mtz.object(file_name=file_name)
  except RuntimeError: pass
  else: return ("ccp4_mtz", content)
  if (detect_binary_file.from_initial_block(file_name=file_name)):
    try: content = unpickle_miller_arrays(file_name=file_name)
    except KeyboardInterrupt: raise
    except Exception: pass
    else: return ("cctbx.miller.array", content)
  try:
    with open(file_name) as fh:
      content = cns_reflection_reader.cns_reflection_file(fh)
  except cns_reflection_reader.CNS_input_Error: pass
  else: return ("cns_reflection_file", content)
  try: content = cns_index_fobs_sigma_reader.reader(
    file_name=file_name)
  except RuntimeError: pass
  else: return ("cns_index_fobs_sigma", content)
  try:
    with open(file_name) as fh:
      content = scalepack_merge.reader(fh)
  except scalepack_merge.FormatError: pass
  else: return ("scalepack_merge", content)
  try: content = scalepack_no_merge.reader(file_name)
  except KeyboardInterrupt: raise
  except Exception: pass
  else: return ("scalepack_no_merge_original_index", content)
  try:
    with open(file_name) as fh:
      content = dtrek_reflnlist_reader.reflnlist(fh)
  except KeyboardInterrupt: raise
  except Exception: pass
  else: return ("dtrek_reflnlist", content)
  try: content = shelx_hklf.reader(
    file_name=file_name)
  except KeyboardInterrupt: raise
  except Exception: pass
  else: return ("shelx_hklf", content)
  try:
    # The cif parser uses a lot of memory when reading a file with millions
    # of words (like an xds_ascii file). Thus we filter out obvious non-cif
    # files.
    assert _cif_prefilter(file_name)
    content = cif_reader(file_path=file_name)
    looks_like_a_reflection_file = False
    for block in content.model().values():
      if '_refln_index_h' in block or '_refln.index_h' in block:
        looks_like_a_reflection_file = True
        break
    if not looks_like_a_reflection_file:
      raise RuntimeError
  except KeyboardInterrupt: raise
  except Exception: pass
  else: return ("cif", content)
  try:
    with open(file_name) as fh:
      content = xds_ascii_reader(fh)
  except KeyboardInterrupt: raise
  except Exception: pass
  else: return ("xds_ascii", content)
  try:
    content = xds_integrate_hkl_reader()
    content.read_file(file_name)
  except KeyboardInterrupt: raise
  except Exception: pass
  else: return ("xds_integrate_hkl", content)
  try: content = solve_fpfm_reader(file_name=file_name)
  except KeyboardInterrupt: raise
  except Exception: pass
  else: return ("solve_fpfm", content)
  return (None, None)


class any_reflection_file(object):
  """
  Proxy object for reading a reflection file of unspecified format, and
  extracting the Miller arrays contained therein.

  Examples
  --------
  >>> from iotbx.reflection_file_reader import any_reflection_file
  >>> hkl_file = any_reflection_file("data.mtz")
  >>> print hkl_file.file_type()
  'ccp4_mtz'
  >>> print type(hkl_file.file_content())
  <class 'iotbx_mtz_ext.object'>
  >>> miller_arrays = hkl_file.as_miller_arrays()
  """

  def __init__(self, file_name, ensure_read_access=True, strict=True):
    # strict is no longer used
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
      with open(file_name): # test read access
        pass
    except IOError as e:
      if (ensure_read_access):
        raise Sorry(str(e))
      return
    if (self._observation_type is not None):
      try: self._file_content = shelx_hklf.reader(
        file_name=file_name)
      except KeyboardInterrupt: raise
      except Exception:
        raise Sorry("Not a SHELX reflection file: %s\n"
          "  =%s can only be used for SHELX reflection files."
          % (file_name, self._observation_type))
      else: self._file_type = "shelx_hklf"
    else:
      self._file_type, self._file_content = try_all_readers(
        file_name=file_name)

  def file_name(self):
    """Returns the file name."""
    return self._file_name

  def file_type(self):
    """Return a string specifying the format type (e.g. 'ccp4_mtz')."""
    return self._file_type

  def file_content(self):
    """Return the underlying format-specific object."""
    return self._file_content

  def as_miller_arrays(self,
                       crystal_symmetry=None,
                       force_symmetry=False,
                       merge_equivalents=True,
                       base_array_info=None,
                       assume_shelx_observation_type_is=None,
                       enforce_positive_sigmas=False,
                       anomalous=None,
                       reconstruct_amplitudes=True
     ):
    """
    Convert the contents of the reflection file into a list of
    :py:class:`cctbx.miller.array` objects, each of which may contain multiple
    columns of data from the underlying file.  By default this will
    automatically merge redundant observations to obtain a unique set under
    symmetry.

    :param crystal_symmetry: :py:class:`cctbx.crystal.symmetry` object
      (defaults to using internally specified symmetry, if any)
    :param force_symmetry: TODO
    :param merge_equivalents: merge redundant obervations (default=True)
    :param base_array_info: :py:class:`cctbx.miller.array_info` object
      containing basic information to be propagated to the arrays
    :param assume_shelx_observation_type_is: if specified, instead of raising
      an exception if the SHELX file type is not known from the file name plus
      data type tag, the function will force the specified data type.
    :param reconstruct_amplitudes: ignored by all other readers than mtz reader.
      If set to True mean amplitudes and adjacent anomalous diffference columns will be
      fused into anomalous miller_array object.
      If False, tells the reader not to fuse mean amplitude and adjacent anomalous
      difference columns into anomalous miller_array objects.
    """
    assert (assume_shelx_observation_type_is in
            [None, "amplitudes", "intensities"])
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
          space_group_info=miller_array.space_group_info(),
          raise_sorry_if_incompatible_unit_cell=True)
        result.append(miller_array.customized_copy(
          crystal_symmetry=miller_array.join_symmetry(
            other_symmetry=crystal_symmetry,
            force=force_symmetry,
            raise_sorry_if_incompatible_unit_cell=True))
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
          file=shelx_file, close_file=False)
        shelx_file.seek(0)
        remaining = shelx_file.read()
        shelx_file.close()
        m = re.search(r"^HKLF\s*(\d)", remaining, re.X|re.M|re.S)
        if m is None:
          raise Sorry("%s does not contain the mandatory HKLF instruction"
                      % shelx_file.name)
        if m.group(1) == "4":
          self._observation_type = "intensities"
        elif m.group(1) == "3":
          self._observation_type = "amplitudes"
        else:
          raise Sorry("HKLF %s not supported" % m.group(1))
    if (self._file_type == "ccp4_mtz"):
      result = self._file_content.as_miller_arrays(
        crystal_symmetry=crystal_symmetry,
        force_symmetry=force_symmetry,
        merge_equivalents=merge_equivalents,
        base_array_info=base_array_info,
        anomalous=anomalous,
        reconstruct_amplitudes=reconstruct_amplitudes
      )
    else:
      result = self._file_content.as_miller_arrays(
        crystal_symmetry=crystal_symmetry,
        force_symmetry=force_symmetry,
        merge_equivalents=merge_equivalents,
        base_array_info=base_array_info,
        anomalous=anomalous,
      )
    if (self.file_type() == "shelx_hklf"):
      if ((self._observation_type == "intensities") or
          (assume_shelx_observation_type_is == "intensities")):
        result[0].set_info(result[0].info().customized_copy(
          labels=["Iobs", "SigIobs"]))
        result[0].set_observation_type_xray_intensity()
      elif ((self._observation_type == "amplitudes") or
            (assume_shelx_observation_type_is == "amplitudes")):
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
    # discard reflections where sigma <= 0
    # XXX note that this will happen after data merging, so for unmerged data
    # it is better to specify merge_equivalents=False!
    if (enforce_positive_sigmas):
      result_ = []
      for array in result :
        result_.append(array.enforce_positive_sigmas())
      result = result_
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
      print(file=report_out)
      print("file_name:", file_name, file=report_out)
      report_out.flush()
    reflection_file = any_reflection_file(file_name)
    if (verbose > 0):
      print("file_type:", reflection_file.file_type(), file=report_out)
      print(file=report_out)
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
        print(file=report_out)
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
    print()
    print("Writing all Miller arrays to file:", pickle_file_name)
    easy_pickle.dump(pickle_file_name, all_miller_arrays)
    print()
