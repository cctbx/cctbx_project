from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/__init__.py
#
# xcif — fast CIF parser with iotbx.cif.reader-compatible API.
#
# The reader class parses CIF files via the C++ xcif_ext engine and
# returns a model() that provides the same dict-like interface as
# iotbx.cif.model (case-insensitive block/tag lookup, flex.std_string
# column access, loop iteration, show() serialization).
#
# The adapter layer (_cif_adapter, _block_adapter, _loop_adapter) wraps
# the C++ Document/Block/Loop objects lazily — Python strings and flex
# arrays are materialized only on access, not during parsing.

import sys
import xcif_ext
import copy
from cctbx.array_family import flex

try:
  from collections.abc import MutableMapping
except ImportError:
  from collections import MutableMapping

try:
  from six.moves import StringIO
except ImportError:
  from io import StringIO


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------

class reader(object):
  """Drop-in replacement for iotbx.cif.reader (read path only).

  Accepts the same constructor signature.  The model() method returns
  a _cif_adapter that supports the dict-like iotbx.cif.model API.
  """

  def __init__(self,
               file_path=None,
               file_object=None,
               input_string=None,
               cif_object=None,
               builder=None,
               raise_if_errors=True,
               strict=True):
    assert [file_path, file_object, input_string].count(None) == 2
    self.file_path = file_path
    if file_path is not None:
      self._doc = xcif_ext.parse_file(file_path)
    else:
      if file_object is not None:
        input_string = file_object.read()
        if hasattr(file_object, 'close'):
          file_object.close()
        if file_path is None:
          file_path = "memory"
      self._doc = xcif_ext.parse(input_string)
    self._model = _cif_adapter(self._doc)

  def model(self):
    return self._model

  def error_count(self):
    return 0

  def show_errors(self, max_errors=50, out=None):
    pass


# ---------------------------------------------------------------------------
# CIF adapter (wraps Document — dict of block_name -> block)
# ---------------------------------------------------------------------------

class _cif_adapter(MutableMapping):

  def __init__(self, doc):
    self._doc = doc
    self._names = [doc[i].name for i in range(len(doc))]
    self._names_lower = {}
    for n in self._names:
      self._names_lower[n.lower()] = n

  def __len__(self):
    return len(self._names)

  def __iter__(self):
    for n in self._names:
      yield n

  def __getitem__(self, key):
    real_key = self._names_lower.get(key.lower())
    if real_key is None:
      raise KeyError('Unknown CIF data block name: "%s"' % key)
    blk = self._doc.find_block(real_key)
    if blk is None:
      raise KeyError('Unknown CIF data block name: "%s"' % key)
    return _block_adapter(blk)

  def __setitem__(self, key, value):
    raise NotImplementedError("xcif model is read-only")

  def __delitem__(self, key):
    raise NotImplementedError("xcif model is read-only")

  def keys(self):
    return list(self._names)

  def get(self, key, default=None):
    try:
      return self[key]
    except KeyError:
      return default

  def items(self):
    for n in self._names:
      yield n, self[n]

  def show(self, out=None, indent="  ", indent_row=None,
           data_name_field_width=34,
           loop_format_strings=None,
           align_columns=True):
    if out is None:
      out = sys.stdout
    for name in self._names:
      print("data_%s" % name, file=out)
      self[name].show(
        out=out, indent=indent,
        data_name_field_width=data_name_field_width,
        align_columns=align_columns)

  def __str__(self):
    s = StringIO()
    self.show(out=s)
    return s.getvalue()

  def __deepcopy__(self, memo):
    # Materialize to iotbx.cif.model only on copy
    return copy.deepcopy(self._to_cif_model(), memo)

  def __getstate__(self):
    return self._to_cif_model()

  def __setstate__(self, state):
    # After unpickle, we're just a plain model.cif
    self.__class__ = state.__class__
    self.__dict__ = state.__dict__

  def _to_cif_model(self):
    """Materialize to iotbx.cif.model.cif (only called on copy/pickle)."""
    from iotbx.cif import model
    cif_obj = model.cif()
    for name in self._names:
      blk = self._doc.find_block(name)
      cif_obj[name] = _block_adapter(blk)._to_block_model()
    return cif_obj


# ---------------------------------------------------------------------------
# Block adapter (wraps Block — dict of tag -> value or flex.std_string)
# ---------------------------------------------------------------------------

class _block_adapter(MutableMapping):

  def __init__(self, xcif_block):
    self._b = xcif_block
    self._pair_tags = None
    self._loop_map = None
    self._tag_to_loop = None

  def _ensure_pair_tags(self):
    if self._pair_tags is None:
      self._pair_tags = list(self._b.pair_tags)

  def _ensure_loop_map(self):
    if self._loop_map is not None:
      return
    self._loop_map = {}
    self._tag_to_loop = {}
    for lp in self._b.loops:
      wrapped = _loop_adapter(lp)
      loop_tags = list(lp.tags)
      if loop_tags:
        cat = _loop_category(loop_tags)
        self._loop_map[cat] = wrapped
      for t in loop_tags:
        self._tag_to_loop[t.lower()] = wrapped

  @property
  def loops(self):
    self._ensure_loop_map()
    return dict(self._loop_map)

  def iterloops(self):
    self._ensure_loop_map()
    for lp in self._loop_map.values():
      yield lp

  def __len__(self):
    return len(self.keys())

  def __iter__(self):
    for k in self.keys():
      yield k

  def keys(self):
    self._ensure_pair_tags()
    self._ensure_loop_map()
    result = list(self._pair_tags)
    for lp in self._b.loops:
      result.extend(list(lp.tags))
    return result

  def __contains__(self, key):
    if self._b.has_tag(key):
      return True
    self._ensure_loop_map()
    return key.lower() in self._tag_to_loop

  def __getitem__(self, key):
    # Single tag-value pair?
    val = self._b.find_value(key)
    if val is not None:
      return val
    # Looped tag?
    self._ensure_loop_map()
    lp = self._tag_to_loop.get(key.lower())
    if lp is not None:
      return lp[key]
    # Loop category name?
    kl = key.lower()
    for cat, loop_adapt in self._loop_map.items():
      if cat.lower() == kl:
        return loop_adapt
    raise KeyError(key)

  def __setitem__(self, key, value):
    raise NotImplementedError("xcif model is read-only")

  def __delitem__(self, key):
    raise NotImplementedError("xcif model is read-only")

  def get_loop(self, loop_name, default=None):
    self._ensure_loop_map()
    kl = loop_name.lower()
    for cat, lp in self._loop_map.items():
      if cat.lower() == kl:
        return lp
    return default

  def get_loop_or_row(self, loop_name, default=None):
    lp = self.get_loop(loop_name)
    if lp is not None:
      return lp
    # Synthesize a single-row loop from scalar tag-value pairs
    # that share this category prefix (matches iotbx.cif.model behavior).
    prefix = loop_name.lower()
    if not prefix.endswith('.'):
      prefix += '.'
    self._ensure_pair_tags()
    found = {}
    for tag in self._pair_tags:
      if tag.lower().startswith(prefix):
        val = self._b.find_value(tag)
        found[tag] = flex.std_string([val])
    if found:
      return _synthetic_loop_adapter(found)
    return default

  def show(self, out=None, indent="  ", indent_row=None,
           data_name_field_width=34,
           loop_format_strings=None,
           align_columns=True):
    if out is None:
      out = sys.stdout
    fmt = "%%-%is" % (data_name_field_width - 1)
    self._ensure_pair_tags()
    self._ensure_loop_map()
    # Tag-value pairs
    for tag in self._pair_tags:
      val = self._b.find_value(tag)
      print(fmt % tag, _format_value(val), file=out)
    # Loops
    for lp in self._b.loops:
      _loop_adapter(lp).show(
        out=out, indent=indent,
        align_columns=align_columns)
      print(file=out)

  def __str__(self):
    s = StringIO()
    self.show(out=s)
    return s.getvalue()

  def __deepcopy__(self, memo):
    return copy.deepcopy(self._to_block_model(), memo)

  def _to_block_model(self):
    from iotbx.cif import model
    blk = model.block()
    self._ensure_pair_tags()
    for tag in self._pair_tags:
      blk[tag] = self._b.find_value(tag)
    for lp in self._b.loops:
      blk.add_loop(_loop_adapter(lp)._to_loop_model())
    for sf_block in self._b.save_frames:  # needs the new binding
      sf = model.save()
      for tag in sf_block.pair_tags:
        sf[tag] = sf_block.find_value(tag)
      for xloop in sf_block.loops:
        sf.add_loop(_loop_adapter(xloop)._to_loop_model())
      blk[sf_block.name] = sf
    return blk


# ---------------------------------------------------------------------------
# Loop adapter (wraps Loop — dict of tag -> flex.std_string column)
# ---------------------------------------------------------------------------

class _loop_adapter(MutableMapping):

  def __init__(self, xcif_loop):
    self._lp = xcif_loop
    self._tags = None
    self._keys_lower = None

  def _ensure_tags(self):
    if self._tags is None:
      self._tags = list(self._lp.tags)
      self._keys_lower = {}
      for t in self._tags:
        self._keys_lower[t.lower()] = t

  def __len__(self):
    self._ensure_tags()
    return len(self._tags)

  def __iter__(self):
    self._ensure_tags()
    for t in self._tags:
      yield t

  def __getitem__(self, key):
    self._ensure_tags()
    real_key = self._keys_lower.get(key.lower())
    if real_key is None:
      raise KeyError(key)
    return self._lp.column_as_flex_string(real_key)

  def __setitem__(self, key, value):
    raise NotImplementedError("xcif model is read-only")

  def __delitem__(self, key):
    raise NotImplementedError("xcif model is read-only")

  def keys(self):
    self._ensure_tags()
    return list(self._tags)

  def name(self):
    self._ensure_tags()
    return _loop_category(self._tags)

  def size(self):
    return self._lp.length

  def n_rows(self):
    return self._lp.length

  def n_columns(self):
    return self._lp.width

  def show(self, out=None, indent="  ", indent_row=None,
           fmt_str=None, align_columns=True):
    if out is None:
      out = sys.stdout
    if indent_row is None:
      indent_row = indent
    self._ensure_tags()
    n_rows = self._lp.length
    n_cols = self._lp.width
    if n_rows == 0 or n_cols == 0:
      return
    print("loop_", file=out)
    for tag in self._tags:
      print(indent + tag, file=out)
    # Collect columns for alignment
    columns = []
    for tag in self._tags:
      col = list(self._lp.column_as_flex_string(tag))
      columns.append([_format_value(v) for v in col])
    if align_columns and fmt_str is None:
      fmt_parts = []
      for col in columns:
        non_multi = [v for v in col if "\n" not in v]
        w = max((len(v) for v in non_multi), default=1)
        # Right-align numeric columns, left-align string columns
        # (matching iotbx.cif.model.loop.show behavior)
        non_null = [v for v in non_multi if v != "." and v != "?"]
        is_numeric = True
        for v in non_null:
          try:
            float(v)
          except ValueError:
            is_numeric = False
            break
        if not is_numeric:
          w = -w
        fmt_parts.append("%%%is" % w)
      row_fmt = indent_row + "  ".join(fmt_parts)
      for r in range(n_rows):
        vals = tuple(columns[c][r] for c in range(n_cols))
        print((row_fmt % vals).rstrip(), file=out)
    else:
      for r in range(n_rows):
        parts = [columns[c][r] for c in range(n_cols)]
        print(indent_row + " ".join(parts), file=out)

  def __str__(self):
    s = StringIO()
    self.show(out=s)
    return s.getvalue()

  def __deepcopy__(self, memo):
    return copy.deepcopy(self._to_loop_model(), memo)

  def _to_loop_model(self):
      from iotbx.cif import model
      tags = list(self._lp.tags)

      # Initialize an entirely empty loop. DO NOT pass header=tags.
      lp = model.loop()

      for tag in tags:
        # On the first iteration, len(lp) is 0, so the assertion is bypassed.
        # On the second iteration, len(lp) is 1, and the assertion safely checks 28 == 28.
        lp[tag] = self._lp.column_as_flex_string(tag)

      return lp

# ---------------------------------------------------------------------------
# Synthetic loop adapter (wraps grouped scalar pairs as a single-row loop)
# ---------------------------------------------------------------------------

class _synthetic_loop_adapter(MutableMapping):
  """Single-row loop synthesized from scalar tag-value pairs that share
  a common category prefix.  Returned by _block_adapter.get_loop_or_row()
  when no real CIF loop_ exists for the requested category, matching
  iotbx.cif.model.block.get_loop_or_row() behavior."""

  def __init__(self, tag_value_dict):
    self._data = tag_value_dict
    self._tags = list(tag_value_dict.keys())
    self._keys_lower = {}
    for t in self._tags:
      self._keys_lower[t.lower()] = t

  def __len__(self):
    return len(self._tags)

  def __iter__(self):
    for t in self._tags:
      yield t

  def __getitem__(self, key):
    real_key = self._keys_lower.get(key.lower())
    if real_key is None:
      raise KeyError(key)
    return self._data[real_key]

  def __setitem__(self, key, value):
    raise NotImplementedError("xcif model is read-only")

  def __delitem__(self, key):
    raise NotImplementedError("xcif model is read-only")

  def keys(self):
    return list(self._tags)

  def name(self):
    return _loop_category(self._tags)

  def size(self):
    return 1

  def n_rows(self):
    return 1

  def n_columns(self):
    return len(self._tags)

  def show(self, out=None, indent="  ", indent_row=None,
           fmt_str=None, align_columns=True):
    if out is None:
      out = sys.stdout
    if indent_row is None:
      indent_row = indent
    n_cols = len(self._tags)
    if n_cols == 0:
      return
    print("loop_", file=out)
    for tag in self._tags:
      print(indent + tag, file=out)
    columns = []
    for tag in self._tags:
      col = list(self._data[tag])
      columns.append([_format_value(v) for v in col])
    if align_columns:
      fmt_parts = []
      for col in columns:
        w = max((len(v) for v in col), default=1)
        non_null = [v for v in col if v != "." and v != "?"]
        is_numeric = True
        for v in non_null:
          try:
            float(v)
          except ValueError:
            is_numeric = False
            break
        if not is_numeric:
          w = -w
        fmt_parts.append("%%%is" % w)
      row_fmt = indent_row + "  ".join(fmt_parts)
      vals = tuple(columns[c][0] for c in range(n_cols))
      print((row_fmt % vals).rstrip(), file=out)
    else:
      parts = [columns[c][0] for c in range(n_cols)]
      print(indent_row + " ".join(parts), file=out)

  def __str__(self):
    s = StringIO()
    self.show(out=s)
    return s.getvalue()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _loop_category(tags):
  """Extract the common category prefix from a list of loop tags."""
  if not tags:
    return ""
  if len(tags) == 1:
    return tags[0]
  # Find common prefix
  prefix = tags[0]
  for t in tags[1:]:
    i = 0
    while i < len(prefix) and i < len(t) and prefix[i] == t[i]:
      i += 1
    prefix = prefix[:i]
  # Trim to last '.' or '_' boundary
  if not (prefix.endswith('.') or prefix.endswith('_')):
    dot = prefix.rfind('.')
    if dot >= 0:
      prefix = prefix[:dot]
    else:
      under = prefix.rfind('_')
      if under > 0:
        prefix = prefix[:under]
  return prefix.rstrip('.')


def _format_value(value):
  """Format a CIF value for output (add quotes if needed)."""
  if not value:
    return "''"
  if '\n' in value:
    if value[0] != '\n':
      value = '\n' + value
    if value[-1] != '\n':
      value = value + '\n'
    return "\n;%s;\n" % value
  has_space = ' ' in value or '\t' in value
  if has_space or value[0] in ('#', '$', '[', ']', '_', ';'):
    if "'" not in value:
      return "'%s'" % value
    if '"' not in value:
      return '"%s"' % value
    return "\n;\n%s\n;\n" % value
  return value
