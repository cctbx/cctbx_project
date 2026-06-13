'''
iotbx.file_io.detection: fast, robust determination of a file's DataManager
datatype without a full parse. Never raises (returns None on failure).

Tiers: extension map -> cheap content verification (magic + text/binary, plus an
authoritative xcif parse for CIF) -> any_file fallback. See the design spec.
'''
import os

from libtbx.utils import detect_binary_file, Sorry
from iotbx.file_io._maps import (
  extension_to_datatype, data_manager_type, CIF_SENTINEL)
from iotbx.file_reader import any_file, splitext, strip_shelx_format_extension
from iotbx.file_io.text_detection import sniff_text_datatype

# datatypes whose payload is text (not binary)
_TEXT_DATATYPES = frozenset(
  ['model', 'sequence', 'phil', 'ncs_spec', 'restraint', 'json', 'yaml'])

# text datatypes with a positive content sniffer. json/yaml are excluded: any_file
# cannot classify them, and they have dedicated readers, so they keep
# extension-trust. restraint never reaches the text branch (.cif/.mmcif resolve as
# CIF_SENTINEL before it).
_SNIFFABLE_TEXT_DATATYPES = frozenset(['model', 'sequence', 'phil', 'ncs_spec'])

# decompression errors from a bounded prefix read of a corrupt OR truncated
# compressed file. A bad header raises a format-specific error (LZMAError,
# ZstdError, zlib.error); a valid-header but truncated stream (interrupted
# download, partial copy) raises EOFError. None of these subclass OSError (unlike
# gzip's BadGzipFile, already caught), so without listing them they miss the
# caught set in get_file_type and trip the unexpected-error warning, mislabeling
# a corrupt input as an internal bug. lzma/zlib are optional on minimal builds
# and zstandard is an optional dependency, so collect them behind guarded imports.
_DECOMPRESS_ERRORS = (EOFError,)
try:
  import zlib as _zlib
  _DECOMPRESS_ERRORS = _DECOMPRESS_ERRORS + (_zlib.error,)
except ImportError:
  pass
try:
  import lzma as _lzma
  _DECOMPRESS_ERRORS = _DECOMPRESS_ERRORS + (_lzma.LZMAError,)
except ImportError:
  pass
try:
  import zstandard as _zstandard
  _DECOMPRESS_ERRORS = _DECOMPRESS_ERRORS + (_zstandard.ZstdError,)
except ImportError:
  pass

# compression suffixes (lowercase) we can cheaply stream a prefix from
_PEEKABLE = frozenset(['.gz', '.bz2', '.xz', '.lzma', '.zst', '.zstd'])

_PREFIX_BYTES = 65536  # 64 KB cap for content sniffing


def _read_prefix(file_name, n=None):
  '''
  Read up to n decompressed bytes from the start of file_name (binary).

  n defaults to _PREFIX_BYTES, read at call time (not bound as a default
  argument) so the cap stays tunable.

  Parameters
  ----------
  file_name : str
      The filepath to read a prefix from (may be compressed)
  n : int, optional
      The maximum number of bytes to read. If None, _PREFIX_BYTES is used.

  Returns
  -------
  bytes
      Up to n decompressed bytes from the start of the file
  '''
  if n is None:
    n = _PREFIX_BYTES
  from libtbx import smart_open
  # bz2_open is binary-only and rejects 'rb' (asserts mode in ('r','w')); its
  # 'r' mode already yields bytes, so use that for .bz2 while every other path
  # (compressed or plain) takes 'rb'.
  mode = 'r' if file_name.lower().endswith('.bz2') else 'rb'
  handle = smart_open.for_reading(file_name, mode=mode)
  try:
    return handle.read(n)
  finally:
    try:
      handle.close()
    except Exception:
      pass


def _looks_like_mtz(prefix):
  '''
  Test whether a byte prefix carries the MTZ magic stamp.

  Parameters
  ----------
  prefix : bytes
      A prefix of the file content

  Returns
  -------
  bool
      True if the prefix starts with the MTZ magic 'MTZ '
  '''
  return prefix[:4] == b'MTZ '


def _looks_like_ccp4_map(prefix):
  '''
  Test whether a byte prefix carries the CCP4/MRC map magic stamp.

  Parameters
  ----------
  prefix : bytes
      A prefix of the file content

  Returns
  -------
  bool
      True if the 'MAP ' stamp is present at byte offset 208
  '''
  # CCP4/MRC: the 'MAP ' stamp sits at byte offset 208
  return prefix[208:212] == b'MAP '


def _is_binary(prefix):
  '''
  Test whether a byte prefix looks like binary (non-text) content.

  Parameters
  ----------
  prefix : bytes
      A prefix of the file content

  Returns
  -------
  bool
      True if the prefix is classified as binary
  '''
  if not prefix:
    return False
  # detect_binary_file renders a verdict only after monitor_initial (default
  # 1000) characters; cap the window at the prefix length so short binary files
  # (e.g. truncated data misnamed .pdb) are still classified as binary instead
  # of slipping through as None -> not-binary.
  detector = detect_binary_file(monitor_initial=min(len(prefix), 1000))
  return bool(detector.is_binary_file(prefix))


def _is_cif(filename):
  '''
  Whether filename has a .cif/.mmcif extension (compression-aware).

  Parameters
  ----------
  filename : str

  Returns
  -------
  bool
  '''
  _, ext, _ = splitext(strip_shelx_format_extension(filename))
  return ext.lower() in ('.cif', '.mmcif')


def _cif_datatypes(filename, cif_engine='xcif'):
  '''
  Parse a CIF (fast via xcif at any size) and return the set of DataManager
  datatypes it contains.

  Parameters
  ----------
  filename : str
      The CIF filepath
  cif_engine : str, optional
      The iotbx.cif engine ("xcif" or "ucif")

  Returns
  -------
  set of str
      A subset of {'model', 'miller_array', 'restraint'}; empty if the file
      cannot be parsed as a CIF
  '''
  import iotbx.cif
  types = set()
  try:
    cif_model = iotbx.cif.reader(file_path=filename, engine=cif_engine).model()
  except Exception:
    return types
  for block_name, block in cif_model.items():
    keys = list(block.keys())
    if any(k.startswith('_atom_site.') or k.startswith('_atom_site_')
           for k in keys):
      types.add('model')
    if any(k.startswith('_refln.') or k.startswith('_refln_') for k in keys):
      types.add('miller_array')
    if (block_name.startswith('comp_')
        or any(k.startswith('_chem_comp_atom') or k.startswith('_chem_comp_bond')
               or k.startswith('_chem_link') for k in keys)):
      types.add('restraint')
  return types


def _primary_cif_type(types):
  '''
  Choose a single representative datatype from a set of CIF datatypes.

  Parameters
  ----------
  types : set of str
      Datatypes present in a CIF (from _cif_datatypes)

  Returns
  -------
  str or None
      The representative datatype by precedence (miller_array > model >
      restraint, mirroring any_file._try_as_cif), or None if the set is empty
  '''
  for datatype in ('miller_array', 'model', 'restraint'):
    if datatype in types:
      return datatype
  return None


def _filtered(datatype, valid_types):
  '''
  Restrict a detected datatype to an allowed set.

  Parameters
  ----------
  datatype : str or None
      The detected datatype
  valid_types : iterable of str or None
      The allowed datatypes. If None, any datatype is allowed.

  Returns
  -------
  str or None
      datatype if it is allowed (and not None), otherwise None
  '''
  if datatype is None:
    return None
  if (valid_types is not None) and (datatype not in valid_types):
    return None
  return datatype


def _fallback(file_name, valid_types):
  '''
  Authoritative fallback: a full parse via any_file (current behavior).

  Parameters
  ----------
  file_name : str
      The filepath to parse
  valid_types : iterable of str or None
      The allowed datatypes used to filter the result

  Returns
  -------
  str or None
      The DataManager datatype from a full any_file parse, filtered by
      valid_types, or None
  '''
  return _filtered(data_manager_type.get(any_file(file_name).file_type),
                   valid_types)


def get_file_type(filename, valid_types=None, verify=True, logger=None,
                  cif_engine='xcif'):
  '''
  Return the DataManager datatype for filename (e.g. 'model', 'miller_array',
  'json'), or None if unrecognized/unsupported. Never raises.

  Detection is tiered: extension map -> cheap content verification (binary-format
  magic, text/binary gate, authoritative xcif parse for CIF) -> any_file
  fallback. The body is wrapped so the public contract (return None on failure)
  holds; expected I/O/parse errors map to None, and any unexpected error is
  surfaced (via logger or warnings.warn, or re-raised when IOTBX_FILE_IO_DEBUG
  is set) before degrading to None.

  Parameters
  ----------
  filename : str
      The filepath to detect the datatype of
  valid_types : iterable of str, optional
      The allowed datatypes. If None, any datatype is allowed.
  verify : bool, optional
      If True, confirm the extension-based guess against the file content;
      if False, trust the extension.
  logger : file-like, optional
      Destination for an unexpected-error message. If None, the message is
      emitted via warnings.warn instead.
  cif_engine : str, optional
      The iotbx.cif engine ("xcif" or "ucif") used to resolve a CIF datatype.

  Returns
  -------
  str or None
      The DataManager datatype, or None if unrecognized or not in valid_types
  '''
  try:
    clean = strip_shelx_format_extension(filename)
    if not os.path.isfile(clean):
      return None
    _, file_ext, compress_ext = splitext(clean)
    ext = file_ext[1:].lower() if file_ext.startswith('.') else file_ext.lower()
    candidate = extension_to_datatype.get(ext)

    # compression gate (decides whether a cheap prefix peek is possible)
    can_peek = True
    if compress_ext is not None:
      ce = compress_ext.lower()
      if ce == '.zip':
        return None                  # not readable by the I/O layer
      # .Z reads at population (gunzip) but offers no cheap prefix peek; any
      # non-streaming compressor likewise degrades to extension-trust
      can_peek = ce in _PEEKABLE

    # a .cif/.mmcif is resolved by an authoritative xcif parse (fast at any
    # size), which also reveals combined CIFs (model + reflections + restraints)
    if candidate == CIF_SENTINEL:
      primary = _primary_cif_type(_cif_datatypes(filename, cif_engine))
      if primary is None:
        return _fallback(filename, valid_types)
      return _filtered(primary, valid_types)

    if candidate is None:
      return _fallback(filename, valid_types)

    # extension-trust (verify disabled, .Z, or a compressor we cannot peek)
    if (not verify) or (not can_peek):
      return _filtered(candidate, valid_types)

    # bounded prefix; a missing optional decompressor (zstandard) degrades to
    # extension-trust rather than a hard failure
    try:
      prefix = _read_prefix(clean)
    except RuntimeError:
      return _filtered(candidate, valid_types)

    # binary-format magic
    if (candidate == 'miller_array') and (ext == 'mtz'):
      if _looks_like_mtz(prefix):
        return _filtered('miller_array', valid_types)
      return _fallback(filename, valid_types)
    if candidate == 'real_map':
      if _looks_like_ccp4_map(prefix):
        return _filtered('real_map', valid_types)
      return _fallback(filename, valid_types)

    # text gate
    if candidate in _TEXT_DATATYPES:
      if _is_binary(prefix):
        return _fallback(filename, valid_types)
      # positive content verification for the sniffable text types: a confident
      # content match wins over the extension (reclassify); fall back to any_file
      # only when no signature matches. json/yaml have no sniffer (any_file cannot
      # classify them) and keep extension-trust.
      if candidate in _SNIFFABLE_TEXT_DATATYPES:
        sniffed = sniff_text_datatype(prefix.decode('utf-8', errors='replace'))
        if sniffed is not None:
          return _filtered(sniffed, valid_types)
        return _fallback(filename, valid_types)

    return _filtered(candidate, valid_types)
  except (IOError, OSError, UnicodeDecodeError, RuntimeError, ValueError,
          Sorry) + _DECOMPRESS_ERRORS:
    return None
  except Exception as e:  # safety net: surface bugs, keep the contract
    if os.environ.get('IOTBX_FILE_IO_DEBUG'):
      raise
    message = 'iotbx.file_io.get_file_type: unexpected error for %r: %s' % (
      filename, e)
    if logger is not None:
      print(message, file=logger)
    else:
      import warnings
      warnings.warn(message)
    return None
