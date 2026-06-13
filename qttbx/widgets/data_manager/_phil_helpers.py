"""Internal helpers for the DataManagerWidget.

Path normalization, .style file_type parsing, data-type detection, and
walking a PhilModel for compatible path parameters. These live in their
own module so the table model and widget code can both call them, and
so they can be unit-tested without spinning up Qt.
"""

import os


def normalize_path(filename):
  """Canonical absolute path form used for all internal comparisons.

  Applies ``expanduser``, ``normpath``, and ``abspath`` in that order.
  The result is the form used everywhere a filename leaves the user
  gesture / PHIL boundary: comparisons against PHIL values, keys in
  the table model, stale-row records, and arguments to the
  DataManager's ``process_<type>_file`` methods.

  ``realpath`` is deliberately *not* applied: two symlinks to the same
  physical file are treated as distinct files. The DataManager and
  PHIL string values use the path the user provided; the widget
  follows the same convention.

  Parameters
  ----------
  filename : str
    Any path string (absolute, relative, with ``~``, or with ``.``/``..``).

  Returns
  -------
  str
    The canonical absolute path.
  """
  return os.path.abspath(os.path.normpath(os.path.expanduser(filename)))


def parse_file_type_style(style_str):
  """Extract the ``<suffix>`` from a PHIL ``.style`` token of the form
  ``file_type:<suffix>``.

  PHIL ``.style`` values are space-separated tokens; this helper finds
  the first ``file_type:<suffix>`` token and returns the suffix.
  Returns ``None`` if no such token is present (or the suffix is empty).

  Parameters
  ----------
  style_str : str or None
    The contents of a PHIL definition's ``.style`` attribute.

  Returns
  -------
  str or None
    The suffix (e.g. ``"pdb"``), or ``None``.
  """
  if not style_str:
    return None
  for token in style_str.split():
    if token.startswith("file_type:"):
      suffix = token[len("file_type:"):]
      if suffix:
        return suffix
      return None
  return None


def detect_data_type(filename):
  """Detect the DataManager data type for ``filename``.

  Delegates to :func:`iotbx.file_io.get_file_type` (fast detection that
  supersedes the old any_file pass). Returns ``None`` if the file is
  missing, unreadable, or unrecognized. Never raises.

  Parameters
  ----------
  filename : str
    A filesystem path. Normalized internally before detection.

  Returns
  -------
  str or None
    A DataManager data type (e.g. ``"model"``, ``"miller_array"``), or
    ``None``.
  """
  from iotbx.file_io import get_file_type
  filename = normalize_path(filename)
  if not os.path.exists(filename):
    return None
  return get_file_type(filename)


def compatible_phil_params(phil_model, file_data_type):
  """Walk a PhilModel for path definitions matching a data type.

  Parameters
  ----------
  phil_model : qttbx.phil.PhilModel
    The active PHIL model.
  file_data_type : str
    A DataManager data type name (e.g. ``"model"``, ``"miller_array"``).

  Returns
  -------
  list of (str, libtbx.phil.definition)
    Pairs ``(phil_path, definition)`` for every definition with
    ``.type == path`` and a ``.style file_type:<suffix>`` whose suffix
    maps (via :data:`iotbx.data_manager.data_manager_type`) to
    ``file_data_type``. Empty list when nothing matches or when
    ``phil_model`` is None.
  """
  if phil_model is None:
    return []
  from iotbx.data_manager import data_manager_type
  out = []
  for phil_path, defn in phil_model.iter_definitions():
    phil_type = getattr(getattr(defn, "type", None), "phil_type", None)
    if phil_type != "path":
      continue
    suffix = parse_file_type_style(getattr(defn, "style", None) or "")
    if suffix is None:
      continue
    if data_manager_type.get(suffix) == file_data_type:
      out.append((phil_path, defn))
  return out
