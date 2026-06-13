'''
iotbx.file_io.reader: read_file, the single place a data file is parsed for
content (parse once, type known), plus FileIOResult. Delegates parsing to the
existing format readers; does not reimplement any parser.
'''
from libtbx.utils import Sorry
from iotbx.file_io._maps import any_file_type
from iotbx.file_io.detection import get_file_type
import os
from iotbx.file_reader import any_file, splitext


class FileIOResult(object):
  '''
  Uniform container for a parsed data file.

  Attributes
  ----------
  data_type : str
      DataManager datatype str (e.g. 'model', 'restraint', 'json')
  file_object : object
      Parsed content. For any_file-backed types this is exactly
      any_file(...).file_object (e.g. model -> the pdb hierarchy input,
      whose .input feeds mmtbx.model.manager). For 'restraint' it is the
      iotbx.cif reader (.model() yields the cif model). For json/yaml it is
      the parsed Python object.
  reader : object or None
      The underlying any_file_input, or None.
  '''

  def __init__(self, data_type, file_object, reader=None):
    self.data_type = data_type
    self.file_object = file_object
    self.reader = reader

  @property
  def file_content(self):
    '''
    Alias for file_object, mirroring iotbx.file_reader.any_file.

    Returns
    -------
    object
        The parsed content (the same value as file_object)
    '''
    # alias, mirrors iotbx.file_reader.any_file
    return self.file_object


def _read_with_any_file(filename, datatype, force=False):
  '''
  Parse an any_file-backed datatype and wrap it in a FileIOResult.

  Parameters
  ----------
  filename : str
      The filepath to parse
  datatype : str
      The expected DataManager datatype
  force : bool, optional
      If True, force any_file to parse filename as datatype (used to extract a
      specific type from a combined CIF that auto-detect would resolve to the
      precedence winner)

  Returns
  -------
  FileIOResult
      The parsed file_object together with the underlying any_file reader

  Raises
  ------
  Sorry
      If the file is missing or does not parse as the expected datatype
  '''
  if force:
    a = any_file(filename, force_type=any_file_type[datatype],
                 raise_sorry_if_errors=True)
    if a.file_object is None:
      raise Sorry('%s is not a valid %s file' % (filename, datatype))
    return FileIOResult(datatype, a.file_object, reader=a)
  # auto-detect + parse in one pass (mirrors the historical process_<type>_file),
  # then validate so a wrong-type file yields the canonical clean message rather
  # than a parser-internal error. any_file raises "Couldn't find the file" when
  # the path is missing.
  a = any_file(filename)
  if (a.file_object is None) or (a.file_type != any_file_type[datatype]):
    raise Sorry('%s is not a recognized %s file' % (filename, datatype))
  return FileIOResult(datatype, a.file_object, reader=a)


def _read_restraint(filename, cif_engine, force=False):
  '''
  Read a restraint CIF and wrap the iotbx.cif reader in a FileIOResult.

  Parameters
  ----------
  filename : str
      The filepath of the restraint CIF
  cif_engine : str
      The CIF engine passed to iotbx.cif.reader (e.g. 'xcif' or 'ucif')
  force : bool, optional
      If True, skip the restraint type check (a combined CIF will not detect as
      'restraint' but still carries restraint data)

  Returns
  -------
  FileIOResult
      A result whose file_object is the iotbx.cif reader

  Raises
  ------
  Sorry
      If the file is missing or is not a recognized restraints file
  '''
  if not os.path.isfile(filename):
    raise Sorry("Couldn't find the file %s" % filename)
  # cheaply confirm this really is a restraint CIF (reject model/reflection
  # CIFs) before the parse, preserving the historical message
  if not force and get_file_type(filename) != 'restraint':
    raise Sorry('%s is not a recognized restraints file' % filename)
  import iotbx.cif
  reader = iotbx.cif.reader(file_path=filename, strict=False, engine=cif_engine)
  return FileIOResult('restraint', reader, reader=None)


def _read_serialized(filename, datatype, load, type_label):
  '''
  Parse a text-serialized file (JSON or YAML) and wrap the result.

  Parameters
  ----------
  filename : str
      The filepath of the serialized file
  datatype : str
      The DataManager datatype ('json' or 'yaml')
  load : callable
      A loader taking an open text handle and returning a Python object
      (e.g. json.load or yaml.safe_load)
  type_label : str
      Human-readable format name used in the error message (e.g. 'JSON', 'YAML')

  Returns
  -------
  FileIOResult
      A result whose file_object is the parsed Python object

  Raises
  ------
  Sorry
      If the file is missing or does not parse as the given format
  '''
  if not os.path.isfile(filename):
    raise Sorry("Couldn't find the file %s" % filename)
  try:
    with open(filename, 'r', encoding='utf-8') as f:
      obj = load(f)
  except Exception:
    raise Sorry('%s is not a valid %s file' % (filename, type_label))
  return FileIOResult(datatype, obj, reader=None)


def _read_json(filename):
  '''
  Parse a JSON file and wrap the resulting Python object in a FileIOResult.

  Parameters
  ----------
  filename : str
      The filepath of the JSON file

  Returns
  -------
  FileIOResult
      A result whose file_object is the parsed Python object

  Raises
  ------
  Sorry
      If the file is missing or is not valid JSON
  '''
  import json
  return _read_serialized(filename, 'json', json.load, 'JSON')


def _read_yaml(filename):
  '''
  Parse a YAML file and wrap the resulting Python object in a FileIOResult.

  Parameters
  ----------
  filename : str
      The filepath of the YAML file

  Returns
  -------
  FileIOResult
      A result whose file_object is the parsed Python object

  Raises
  ------
  Sorry
      If the file is missing, PyYAML is not installed, or the file is not
      valid YAML
  '''
  try:
    import yaml
  except ImportError:
    raise Sorry('PyYAML is required to read YAML files. Install it with '
                '"libtbx.conda install pyyaml" or "libtbx.pip install PyYAML".')
  return _read_serialized(filename, 'yaml', yaml.safe_load, 'YAML')


def _read_phil(filename, datatype):
  '''
  Parse a phil file directly via iotbx.phil (no any_file).

  Parameters
  ----------
  filename : str
      The filepath to parse.
  datatype : str
      The DataManager datatype ('phil').

  Returns
  -------
  FileIOResult
      A result whose file_object is the parsed phil scope.

  Raises
  ------
  Sorry
      If the file has no phil objects. A phil syntax error raises the parser's
      RuntimeError, which _dispatch_read routes to the any_file fallback.
  '''
  from iotbx.phil import parse as parse_phil
  phil_object = parse_phil(file_name=filename, process_includes=True)
  if len(phil_object.objects) == 0:
    raise Sorry('%s is an empty parameter file' % filename)
  return FileIOResult(datatype, phil_object, reader=None)


def _read_model(filename, datatype):
  '''
  Parse a PDB or mmCIF model directly via iotbx.pdb.hierarchy.input.

  Parameters
  ----------
  filename : str
      The filepath to parse.
  datatype : str
      The DataManager datatype ('model').

  Returns
  -------
  FileIOResult
      A result whose file_object is the pdb hierarchy input, whose `.input`
      feeds mmtbx.model.manager.

  Raises
  ------
  Sorry
      If the file does not parse or has no ATOM/HETATM records.
  '''
  import iotbx.pdb.hierarchy
  try:
    pdb_inp = iotbx.pdb.hierarchy.input(filename)
  except ValueError as e:
    raise Sorry(str(e))
  if pdb_inp.hierarchy.models_size() == 0:
    raise Sorry('%s has no ATOM or HETATM records' % filename)
  return FileIOResult(datatype, pdb_inp, reader=None)


def _read_real_map(filename, datatype):
  '''
  Read a CCP4/MRC map directly via iotbx.map_manager.

  Parameters
  ----------
  filename : str
      The filepath to read.
  datatype : str
      The DataManager datatype ('real_map').

  Returns
  -------
  FileIOResult
      A result whose file_object is the map_manager.

  Raises
  ------
  Sorry
      If iotbx.map_manager cannot read the file as a map. _dispatch_read catches
      any reader-level failure and falls back to any_file.
  '''
  from iotbx.map_manager import map_manager
  from libtbx.utils import null_out
  map_object = map_manager(file_name=str(filename), log=null_out())
  return FileIOResult(datatype, map_object, reader=None)


def _read_reflections(filename, datatype):
  '''
  Read reflections (mtz/sca/cns/cif/...) directly via any_reflection_file.

  Serves both the miller_array and map_coefficients datatypes. str() guards
  Boost.Python from a unicode path.

  Parameters
  ----------
  filename : str
      The filepath to read.
  datatype : str
      The DataManager datatype ('miller_array' or 'map_coefficients').

  Returns
  -------
  FileIOResult
      A result whose file_object is the any_reflection_file object.

  Raises
  ------
  Sorry
      If the file is not a recognized reflection file.
  '''
  from iotbx.reflection_file_reader import any_reflection_file
  hkl_file = any_reflection_file(str(filename))
  if hkl_file.file_type() is None:
    raise Sorry('%s is not a valid reflection file' % filename)
  return FileIOResult(datatype, hkl_file, reader=None)


def _read_ncs_spec(filename, datatype):
  '''
  Read an ncs object directly via mmtbx.ncs.

  Tries BIOMT records from a PDB input first; if that does not parse, reads the
  ncs_spec/phil text format.

  Parameters
  ----------
  filename : str
      The filepath to read.
  datatype : str
      The DataManager datatype ('ncs_spec').

  Returns
  -------
  FileIOResult
      A result whose file_object is the mmtbx.ncs.ncs object.

  Raises
  ------
  Sorry
      If the file yields no NCS operators.
  '''
  from mmtbx.ncs.ncs import ncs
  from libtbx.utils import null_out
  import iotbx.pdb
  ncs_object = ncs()
  try:
    pdb_inp = iotbx.pdb.input(file_name=filename)
    ncs_object.ncs_from_pdb_input_BIOMT(pdb_inp=pdb_inp, log=null_out(), quiet=True)
  except Exception:
    ncs_object.read_ncs(file_name=filename, log=null_out(), quiet=True)
  if ncs_object.max_operators() == 0:
    raise Sorry('%s has no NCS operators' % filename)
  return FileIOResult(datatype, ncs_object, reader=None)


def _read_sequence(filename, datatype):
  '''
  Read sequences directly via iotbx.bioinformatics.any_sequence_format.

  Mirrors any_file._try_as_seq: rejects misformatted or gapped input and drops
  empty sequences. (any_file's .ccp4/.ncs endswith guards are unnecessary here
  because the datatype is already known.)

  Parameters
  ----------
  filename : str
      The filepath to read.
  datatype : str
      The DataManager datatype ('sequence').

  Returns
  -------
  FileIOResult
      A result whose file_object is the list of sequence objects.

  Raises
  ------
  Sorry
      If the file has no sequence data, is misformatted, or contains gaps (an
      alignment rather than a sequence).
  '''
  from iotbx.bioinformatics import any_sequence_format
  objects, non_compliant = any_sequence_format(filename)
  if objects is None:
    raise Sorry('%s has no sequence data' % filename)
  if len(non_compliant) != 0:
    raise Sorry('%s has misformatted sequence data' % filename)
  for seq_obj in objects:
    if '-' in seq_obj.sequence:
      raise Sorry('%s contains gaps (an alignment, not a sequence)' % filename)
  objects = [o for o in objects if len(o) != 0]
  return FileIOResult(datatype, objects, reader=None)


def _read_direct(filename, datatype, force=False):
  '''
  Read an any_file-backed datatype by calling the underlying format reader
  directly (no any_file).

  The object returned is identical to any_file(force_type=...).file_object for
  that datatype, so consumer contracts are unchanged.

  Parameters
  ----------
  filename : str
      The filepath to read.
  datatype : str
      The DataManager datatype to read it as.
  force : bool, optional
      Accepted for signature parity with the any_file fallback; the direct
      readers parse the file as `datatype` regardless (a combined CIF reads as
      the requested block via the underlying reader).

  Returns
  -------
  FileIOResult
      The parsed result.

  Raises
  ------
  Sorry
      On a reader-level failure (the exception type is reader-specific);
      _dispatch_read catches any such failure and falls back to any_file during
      the deprecation window.
  '''
  if datatype == 'sequence':
    return _read_sequence(filename, datatype)
  if datatype == 'ncs_spec':
    return _read_ncs_spec(filename, datatype)
  if datatype in ('miller_array', 'map_coefficients'):
    return _read_reflections(filename, datatype)
  if datatype == 'real_map':
    return _read_real_map(filename, datatype)
  if datatype == 'model':
    return _read_model(filename, datatype)
  if datatype == 'phil':
    return _read_phil(filename, datatype)
  raise Sorry('no direct reader for datatype %r' % datatype)


def _dispatch_read(filename, datatype, cif_engine, force=False):
  '''
  Route filename to the reader for its datatype.

  Parameters
  ----------
  filename : str
      The filepath to read
  datatype : str
      The DataManager datatype to read it as
  cif_engine : str
      The CIF engine forwarded to the restraint reader
  force : bool, optional
      If True, force the datatype-specific reader to extract datatype from a
      combined CIF (json/yaml ignore it)

  Returns
  -------
  FileIOResult
      The parsed result from the datatype-specific reader

  Raises
  ------
  Sorry
      If datatype cannot be read by any reader, or the underlying reader fails
  '''
  if datatype == 'restraint':
    return _read_restraint(filename, cif_engine, force=force)
  if datatype == 'json':
    return _read_json(filename)
  if datatype == 'yaml':
    return _read_yaml(filename)
  if datatype in any_file_type:
    if not os.path.isfile(filename):
      raise Sorry("Couldn't find the file %s" % filename)
    try:
      return _read_direct(filename, datatype, force=force)
    except Exception:
      # any_file fallback for the deprecation window. The direct readers call the
      # same parsers any_file does, so any failure here is a direct-reader gap (a
      # datatype not yet migrated, or a compression/edge case) -- route it to
      # any_file, which raises a clean Sorry (and read_file's outer handler then
      # retries via decompression). Broad by design: a genuine bug on a VALID file
      # is still caught by each datatype's happy-path `reader is None` test, and
      # this mirrors _read_json/_read_yaml's except Exception. TODO(any_file-
      # deprecation): drop this try/except once direct readers + tests cover every
      # any_file-backed datatype.
      return _read_with_any_file(filename, datatype, force=force)
  raise Sorry('Cannot read %s as datatype %r' % (filename, datatype))


def _read_decompressed(filename, datatype, cif_engine, force=False):
  '''
  Decompress filename to a temp file carrying the inner extension and read that.

  Fallback for compressed files whose underlying format reader cannot read the
  compressor directly (the MTZ/MRC readers, and .bz2 for the text readers). The
  format readers materialize their data in memory, so the temp file is removed
  before returning.

  Parameters
  ----------
  filename : str
      The filepath of the compressed file
  datatype : str
      The DataManager datatype to read it as
  cif_engine : str
      The CIF engine forwarded to the restraint reader
  force : bool, optional
      If True, force the reader to extract datatype from a combined CIF

  Returns
  -------
  FileIOResult
      The parsed result read from the decompressed temp file

  Raises
  ------
  Sorry
      If the file cannot be decompressed or the decompressed content does not
      read as datatype (with the original filename restored in the message)
  '''
  import tempfile
  from libtbx import smart_open
  if filename.endswith('.zip'):
    # smart_open has no .zip reader; report it clearly instead of letting the
    # raw bytes fall through to a misleading "could not be decompressed".
    raise Sorry('%s: ZIP archives are not supported by the I/O layer' % filename)
  _, inner_ext, _ = splitext(filename)
  fd, tmp = tempfile.mkstemp(suffix=inner_ext)
  os.close(fd)
  try:
    if filename.endswith('.Z'):
      # Python has no LZW (.Z) decompressor, and smart_open pipes .Z through
      # gunzip as text (which corrupts a binary payload). Decompress to the
      # temp file as raw bytes via gunzip instead.
      import subprocess
      try:
        with open(tmp, 'wb') as out:
          result = subprocess.run(['gunzip', '-c', filename], stdout=out,
                                  stderr=subprocess.PIPE)
      except FileNotFoundError:
        raise Sorry('gunzip is required to read %s (.Z) but was not found on PATH'
                    % filename)
      if result.returncode != 0:
        detail = result.stderr.decode('utf-8', 'replace').strip()
        raise Sorry('%s could not be decompressed%s'
                    % (filename, (': ' + detail) if detail else ''))
    else:
      try:
        with smart_open.for_reading(filename, mode='r') as handle:
          data = handle.read()
      except RuntimeError as e:
        raise Sorry(str(e))             # e.g. the zstandard package is not installed
      except Exception:
        raise Sorry('%s could not be decompressed' % filename)
      # smart_open's compressor readers yield bytes for mode='r'
      with open(tmp, 'wb') as out:
        out.write(data)
    try:
      return _dispatch_read(tmp, datatype, cif_engine, force)
    except Sorry as e:
      # preserve the specific message but report the original filename, not tmp
      raise Sorry(str(e).replace(tmp, filename))
  finally:
    try:
      os.remove(tmp)
    except OSError:
      pass


def read_file(filename, file_type=None, cif_engine='xcif', force=False):
  '''
  Parse filename once and return a FileIOResult.

  file_type, when given, is the DataManager datatype to read it as (skips
  detection); otherwise it is detected. A compressed file whose format reader
  cannot read the compressor directly is transparently decompressed to a temp
  file and retried.

  Parameters
  ----------
  filename : str
      The filepath to parse
  file_type : str, optional
      The DataManager datatype to read it as. If None, the datatype is
      detected.
  cif_engine : str, optional
      The CIF engine forwarded to the restraint reader
  force : bool, optional
      If True, force the reader to extract file_type from a combined CIF that
      auto-detect would otherwise resolve to the precedence winner

  Returns
  -------
  FileIOResult
      The parsed result

  Raises
  ------
  Sorry
      If the file is missing, unrecognized, or cannot be read as that type
  '''
  datatype = file_type
  if datatype is None:
    datatype = get_file_type(filename)
  if datatype is None:
    if not os.path.isfile(filename):
      raise Sorry("Couldn't find the file %s" % filename)
    raise Sorry('%s is not a recognized file type' % filename)
  try:
    return _dispatch_read(filename, datatype, cif_engine, force)
  except Sorry:
    # the file may use a compressor the underlying reader cannot read directly
    # (text readers handle .gz/.Z; mrcfile handles .gz; the MTZ reader handles
    # none) -> decompress to a temp file and retry once.
    _, _, compress_ext = splitext(filename)
    if (compress_ext is not None) and os.path.isfile(filename):
      return _read_decompressed(filename, datatype, cif_engine, force)
    raise
