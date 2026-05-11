"""DataManagerTableModel -- read-only table model aggregating a DataManager
(file list) and a PhilModel (path-parameter bindings).

Row identity is (normalized_filename, data_type) for DataManager-derived
rows; stale rows (PHIL bindings whose target file is missing or type-
mismatched) are appended after normal rows. Columns are Filename, Type,
Used for, Delete.

Cache invariants (see spec section 4.7) are maintained on every PhilModel
signal (dataChanged, rowsInserted, rowsRemoved, modelReset). The
:meth:`DataManagerTableModel.reconcile_stale` hook removes matching stale
rows whose underlying PHIL value still points at the newly-added file,
keeping the table consistent with successful ``add_file`` gestures.
"""

import os

from PySide2.QtCore import QAbstractTableModel, Qt, QModelIndex


_PRETTY_DATA_TYPE = {
  "model":            "Model",
  "miller_array":     "Miller array",
  "real_map":         "Real map",
  "map_coefficients": "Map coefficients",
  "restraint":        "Restraint",
  "sequence":         "Sequence",
  "ncs_spec":         "NCS spec",
  "phil":             "PHIL",
}


def pretty_data_type(data_type):
  """Return the user-facing label for a DataManager data type.

  Maps the canonical snake_case names from ``iotbx.data_manager`` to a
  human-readable form (Title Case with spaces; acronyms preserved).
  Unknown data types fall back to the snake_case string with
  underscores turned into spaces and the first letter capitalized,
  so future ``DataManager`` additions still render reasonably without
  this map being updated.
  """
  if data_type is None:
    return ""
  pretty = _PRETTY_DATA_TYPE.get(data_type)
  if pretty is not None:
    return pretty
  s = data_type.replace("_", " ")
  return s[:1].upper() + s[1:] if s else s


def _trash_icon():
  """Return the system's standard trash icon, or ``None`` if unavailable.

  Shared by the Delete column's header decoration (see
  :py:meth:`DataManagerTableModel.headerData`) and the delegate's
  per-row paint, so header and cell icons stay in sync. Returns
  ``None`` when the running Qt build does not expose
  ``QStyle.SP_TrashIcon`` (older than Qt 5.12) or when no
  :class:`QApplication` has been constructed yet.
  """
  from PySide2.QtWidgets import QApplication, QStyle
  sp = getattr(QStyle, "SP_TrashIcon", None)
  if sp is None:
    return None
  app = QApplication.instance()
  if app is None:
    return None
  icon = app.style().standardIcon(sp)
  if icon.isNull():
    return None
  return icon

from qttbx.widgets.data_manager._phil_helpers import normalize_path


class StaleRow(object):
  """A row representing a PHIL ``path`` value whose file cannot be added
  to the DataManager (missing or type mismatch).

  One stale row per PHIL binding -- if three PHIL parameters all
  reference the same missing file, the table holds three stale rows.

  Attributes
  ----------
  filename : str
    Normalized path of the missing/mismatched file.
  expected_type : str
    DataManager data type the slot expects.
  phil_path : str
    The PHIL path whose value pointed at this file.
  message : str
    Tooltip text, e.g. "file not found" or "type mismatch: expected
    miller_array".
  """

  __slots__ = ("filename", "expected_type", "phil_path", "message")

  def __init__(self, filename, expected_type, phil_path, message):
    self.filename = filename
    self.expected_type = expected_type
    self.phil_path = phil_path
    self.message = message

  def __repr__(self):
    return ("StaleRow(filename=%r, expected_type=%r, phil_path=%r, "
            "message=%r)" % (self.filename, self.expected_type,
                              self.phil_path, self.message))


class DataManagerTableModel(QAbstractTableModel):
  """Read-only Qt model aggregating DataManager files and PhilModel bindings.

  Each DataManager-derived row corresponds to a single
  ``(normalized_filename, data_type)`` pair enumerated from the attached
  :class:`iotbx.data_manager.DataManager`. When a
  :class:`qttbx.phil.PhilModel` is also attached, a fourth column shows
  which PHIL path parameters reference the file. Stale rows (PHIL
  bindings whose target file is missing or type-mismatched) are appended
  *after* the normal DM-derived rows so DM-row indexing remains stable.
  Live updates flow from PhilModel signals (see :meth:`attach`) and from
  explicit :meth:`refresh`; :meth:`reconcile_stale` removes matching
  stale rows after a successful ``add_file`` so the table stays
  consistent without callers having to manage stale-row lifetimes by
  hand.

  Attributes
  ----------
  COL_FILENAME, COL_TYPE, COL_USED_FOR, COL_DELETE : int
    Column indices used by callers and the item delegate. All four
    columns are always present; without a PhilModel attached, the
    ``Used for`` column stays empty (no chips, no ``+`` button), and
    the Delete column remains visible so files can still be removed.
  """

  COL_FILENAME = 0
  COL_TYPE = 1
  COL_USED_FOR = 2
  COL_DELETE = 3

  def __init__(self, parent=None):
    QAbstractTableModel.__init__(self, parent)
    self._data_manager = None
    self._phil_model = None
    self._rows = []                  # list[(normalized_filename, data_type)]
    # Binding caches built by _rebuild_caches; kept O(1) for data() lookups.
    self._bindings_by_file_and_type = {}  # (norm_fn, dtype) -> list[phil_path]
    self._compatible_by_type = {}         # dtype -> list[(phil_path, defn)]
    self._target_by_phil_path = {}        # phil_path -> (norm_fn|None, dtype)
    self._stale = []                      # list[StaleRow]
    # Optional root directory: when set, filenames inside the root render
    # as paths relative to it; filenames outside the root keep their
    # full normalized form. Internal state (cache keys, _rows, stale
    # rows) always uses the full normalized path.
    self._root = None

  def attach(self, data_manager, phil_model=None):
    """Attach a DataManager and optional PhilModel; rebuild rows and caches.

    Wrapped in ``beginResetModel`` / ``endResetModel`` so any attached
    view is correctly invalidated even on re-attach. After the reset,
    connects to the PhilModel's ``dataChanged``, ``rowsInserted``,
    ``rowsRemoved``, and ``modelReset`` signals so external PHIL
    mutations refresh the binding caches incrementally and the
    ``Used for`` column stays in sync without a full ``refresh()``.

    Parameters
    ----------
    data_manager : iotbx.data_manager.DataManager
      The DataManager whose files this table reflects.
    phil_model : qttbx.phil.PhilModel, optional
      The PhilModel whose path parameters drive the ``Used for`` column.
      When ``None`` the binding caches stay empty; the ``Used for``
      column is still present but renders blank, and the widget runs
      as a pure file manager.
    """
    # Disconnect from any previously attached PhilModel before swapping.
    self.detach()
    self.beginResetModel()
    self._data_manager = data_manager
    self._phil_model = phil_model
    self._rebuild_rows()
    self._rebuild_caches()
    self.endResetModel()
    if phil_model is not None:
      phil_model.dataChanged.connect(self._on_phil_data_changed)
      phil_model.rowsInserted.connect(self._on_phil_rows_inserted)
      phil_model.rowsRemoved.connect(self._on_phil_rows_removed)
      phil_model.modelReset.connect(self._on_phil_model_reset)

  def detach(self):
    """Drop references to attached models and clear table + caches.

    Disconnects from the previously attached PhilModel's signals (if
    any), wrapped in ``try`` / ``except`` so a partially-attached or
    already-disposed PhilModel does not prevent cleanup.
    """
    if self._phil_model is not None:
      try:
        self._phil_model.dataChanged.disconnect(self._on_phil_data_changed)
        self._phil_model.rowsInserted.disconnect(self._on_phil_rows_inserted)
        self._phil_model.rowsRemoved.disconnect(self._on_phil_rows_removed)
        self._phil_model.modelReset.disconnect(self._on_phil_model_reset)
      except (TypeError, RuntimeError):
        # Not connected, or the underlying QObject has been destroyed.
        pass
    self.beginResetModel()
    self._data_manager = None
    self._phil_model = None
    self._rows = []
    self._bindings_by_file_and_type = {}
    self._compatible_by_type = {}
    self._target_by_phil_path = {}
    self._stale = []
    self.endResetModel()

  def refresh(self):
    """Re-enumerate the DataManager and rebuild binding caches.

    Call after external code mutates the attached DataManager outside the
    widget's own gestures. The model emits ``modelReset`` so any
    attached view repopulates from scratch.
    """
    self.beginResetModel()
    self._rebuild_rows()
    self._rebuild_caches()
    self.endResetModel()

  def _rebuild_rows(self):
    """Populate ``self._rows`` from the attached DataManager.

    Iterates ``self._data_manager.datatypes`` and calls the
    ``get_<datatype>_names`` method exposed by each per-datatype
    DataManager mixin (see ``iotbx/data_manager/model.py``,
    ``miller_array.py``, etc.). Datatypes whose mixin does not provide a
    names accessor, and accessors that raise, are silently skipped --
    enumeration must not fail because of a single misbehaving mixin.
    """
    self._rows = []
    if self._data_manager is None:
      return
    for data_type in self._data_manager.datatypes:
      getter = getattr(self._data_manager,
                       "get_%s_names" % data_type, None)
      if getter is None:
        continue
      try:
        names = getter()
      except Exception:
        names = []
      for n in names:
        self._rows.append((normalize_path(n), data_type))

  def _rebuild_caches(self):
    """Walk the attached PhilModel for path definitions and rebuild caches.

    Populates three dicts used by :meth:`used_for` and (in later tasks)
    by the binding popup:

    * ``_compatible_by_type`` -- ``dtype -> list[(phil_path, defn)]``,
      the universe of PHIL paths that *could* point at a file of each
      DataManager data type.
    * ``_target_by_phil_path`` -- ``phil_path -> (norm_fn|None, dtype)``,
      where each non-``.multiple`` path-typed definition currently points
      (or ``None`` if unset). For ``.multiple = True`` path definitions
      the slot stores ``(None, dtype)`` as a placeholder; live targets
      live entirely in ``_bindings_by_file_and_type``. Incremental updates
      for ``.multiple`` paths are handled in Task 7 by re-querying
      ``value_at_path`` and diffing.
    * ``_bindings_by_file_and_type`` -- ``(norm_fn, dtype) -> list[phil_path]``,
      the inverse mapping that gives a row-keyed O(1) lookup for the
      ``Used for`` column.

    Notes on ``.multiple = True`` path definitions
    ----------------------------------------------
    :meth:`PhilModel.iter_definitions` yields each ``.multiple`` template
    once per existing instance (e.g. three times when the extract has
    two instances plus the master template), so this method deduplicates
    by ``phil_path`` and queries :meth:`PhilModel.value_at_path` exactly
    once per definition. For ``.multiple`` paths ``value_at_path`` returns
    a list (``scope_extract_list``); each non-empty entry is normalized
    and added to ``_bindings_by_file_and_type`` so the same template path
    can legitimately appear in multiple file rows' ``Used for`` columns.

    Empty strings are treated identically to ``None`` (unbound) to keep
    this method symmetric with :meth:`_on_phil_data_changed` in Task 7,
    which uses ``if new_value else None`` when reacting to edits that
    clear a path via the empty string.

    Idempotent: each call clears and refills the caches. Safe to call
    without a PhilModel attached (just leaves the caches empty).
    """
    self._bindings_by_file_and_type = {}
    self._compatible_by_type = {}
    self._target_by_phil_path = {}
    if self._phil_model is None:
      return
    from iotbx.data_manager import data_manager_type
    from qttbx.widgets.data_manager._phil_helpers import parse_file_type_style
    # Dedupe: ``.multiple = True`` definitions are yielded once per
    # instance, but we want each template path processed exactly
    # once. Matches the early-add pattern in widget._auto_import.
    seen_paths = set()
    for phil_path, defn in self._phil_model.iter_definitions():
      if phil_path in seen_paths:
        continue
      seen_paths.add(phil_path)
      phil_type = getattr(getattr(defn, "type", None), "phil_type", None)
      if phil_type != "path":
        continue
      suffix = parse_file_type_style(getattr(defn, "style", None) or "")
      if suffix is None:
        continue
      dtype = data_manager_type.get(suffix)
      if dtype is None:
        continue
      self._compatible_by_type.setdefault(dtype, []).append((phil_path, defn))
      value = self._phil_model.value_at_path(phil_path)
      if getattr(defn, "multiple", False):
        # ``.multiple = True`` path: value is a list of strings. The
        # template path can be bound to many files simultaneously, so the
        # _target_by_phil_path entry is just a placeholder; the real
        # mapping lives in _bindings_by_file_and_type per-instance.
        # Dedupe repeated values so external code that put the same
        # path in the list twice does not produce duplicate chips.
        self._target_by_phil_path[phil_path] = (None, dtype)
        if not value:
          continue
        seen_values = set()
        for entry in value:
          if not entry:
            continue
          nv = normalize_path(entry)
          if nv in seen_values:
            continue
          seen_values.add(nv)
          self._bindings_by_file_and_type.setdefault(
            (nv, dtype), []).append(phil_path)
        continue
      # Non-``.multiple`` path: value is a scalar (str or None). Treat
      # empty string as unbound to match Task 7's signal handler, which
      # collapses empty edits to the unbound state.
      if not value:
        self._target_by_phil_path[phil_path] = (None, dtype)
        continue
      nv = normalize_path(value)
      self._target_by_phil_path[phil_path] = (nv, dtype)
      self._bindings_by_file_and_type.setdefault((nv, dtype), []).append(
        phil_path)

  def _affected_rows_for_keys(self, keys):
    """Return the set of table-row indices touched by ``keys``.

    ``keys`` is an iterable of ``(normalized_filename, data_type)``
    tuples. Walks ``self._rows`` once and returns each matching row
    index. Used by the dataChanged handler to know which Used-for
    cells to repaint without emitting on every row.
    """
    affected = set()
    if not keys:
      return affected
    keys = set(keys)
    for i, row_key in enumerate(self._rows):
      if row_key in keys:
        affected.add(i)
    return affected

  def _emit_used_for_changed(self, rows):
    """Emit ``dataChanged`` for the Used-for cell of each row in ``rows``.

    Issued one row at a time (rather than a single bracketing range)
    so we never spuriously invalidate the Used-for cells of rows whose
    bindings did not actually change.
    """
    for r in sorted(rows):
      idx = self.index(r, self.COL_USED_FOR)
      self.dataChanged.emit(idx, idx, [Qt.DisplayRole])

  def _definition_for_phil_path(self, phil_path):
    """Return the libtbx.phil.definition for ``phil_path`` from caches.

    Looks in ``self._compatible_by_type``; returns ``None`` if no
    matching definition has been catalogued (which means it is not a
    file-typed path definition we track).
    """
    for entries in self._compatible_by_type.values():
      for p, d in entries:
        if p == phil_path:
          return d
    return None

  def _on_phil_data_changed(self, topLeft, bottomRight, roles=None):
    """Incremental cache update when PHIL value(s) change externally.

    For each cell in the changed range whose ``phil_path`` is one of
    our tracked file-typed definitions: drop the old binding entry,
    install the new one (using :func:`normalize_path` on the new
    value, treating ``None`` / empty string as unbound), and remember
    which ``(file, data_type)`` keys may have flipped so we can emit
    ``dataChanged`` for the corresponding Used-for cells.

    ``.multiple = True`` path definitions cannot be diffed with the
    scalar old/new comparison this method uses -- their value is a
    list, and adding ``[A]`` -> ``[A, B]`` doesn't change a single
    target. When at least one ``.multiple`` definition is in the
    range, fall back to :py:meth:`_rebuild_caches` and broadcast a
    Used-for change for every row. This matches the conservative
    rebuild used by the rows-inserted / rows-removed handlers.
    """
    if self._phil_model is None:
      return
    affected_keys = set()
    multi_seen = False
    for phil_row in range(topLeft.row(), bottomRight.row() + 1):
      idx = self._phil_model.index(
        phil_row, topLeft.column(), topLeft.parent())
      phil_path = self._phil_model.path_for_index(idx)
      if phil_path is None or phil_path not in self._target_by_phil_path:
        continue
      defn = self._definition_for_phil_path(phil_path)
      if defn is not None and getattr(defn, "multiple", False):
        multi_seen = True
        continue
      old_fn, dtype = self._target_by_phil_path[phil_path]
      if old_fn is not None:
        bucket = self._bindings_by_file_and_type.get((old_fn, dtype))
        if bucket is not None and phil_path in bucket:
          bucket.remove(phil_path)
          if not bucket:
            self._bindings_by_file_and_type.pop((old_fn, dtype), None)
        affected_keys.add((old_fn, dtype))
      new_value = self._phil_model.value_at_path(phil_path)
      new_fn = normalize_path(new_value) if new_value else None
      self._target_by_phil_path[phil_path] = (new_fn, dtype)
      if new_fn is not None:
        self._bindings_by_file_and_type.setdefault(
          (new_fn, dtype), []).append(phil_path)
        affected_keys.add((new_fn, dtype))
    if multi_seen:
      # Rebuild caches conservatively; .multiple path lists can change
      # in ways the scalar diff above can't represent.
      self._rebuild_caches()
      if self._rows:
        self.dataChanged.emit(
          self.index(0, self.COL_USED_FOR),
          self.index(len(self._rows) - 1, self.COL_USED_FOR),
          [Qt.DisplayRole])
      return
    self._emit_used_for_changed(self._affected_rows_for_keys(affected_keys))

  def _on_phil_rows_inserted(self, parent, first, last):
    """Conservatively rebuild caches when PhilModel grows.

    A new ``.multiple`` instance (or any structural addition) can
    surface phil_paths that weren't in the universe before. Rebuilding
    is cheap relative to the cost of getting wrong bindings, and is
    rare in practice.
    """
    self._rebuild_caches()
    if self._rows:
      self.dataChanged.emit(
        self.index(0, self.COL_USED_FOR),
        self.index(len(self._rows) - 1, self.COL_USED_FOR),
        [Qt.DisplayRole])

  def _on_phil_rows_removed(self, parent, first, last):
    """Conservatively rebuild caches when PhilModel shrinks."""
    self._rebuild_caches()
    if self._rows:
      self.dataChanged.emit(
        self.index(0, self.COL_USED_FOR),
        self.index(len(self._rows) - 1, self.COL_USED_FOR),
        [Qt.DisplayRole])

  def _on_phil_model_reset(self):
    """Rebuild caches and reset on a full PhilModel reset."""
    self.beginResetModel()
    self._rebuild_caches()
    self.endResetModel()

  def used_for(self, row):
    """Return the list of phil_paths bound to the file at ``row``.

    O(1) lookup via :attr:`_bindings_by_file_and_type`. Returns an empty
    list when ``row`` is out of range, when no PhilModel is attached, or
    when no PHIL path currently points at the file.

    Parameters
    ----------
    row : int
      Row index, 0 .. ``rowCount() - 1``.

    Returns
    -------
    list of str
      Dotted PHIL paths whose current value matches the row's file. A
      fresh list -- callers may mutate it without aliasing the cache.
    """
    if row < 0 or row >= len(self._rows):
      return []
    filename, data_type = self._rows[row]
    return list(self._bindings_by_file_and_type.get((filename, data_type), ()))

  def _label_for_phil_path(self, phil_path):
    """Return the human label for ``phil_path``, or ``None``.

    Looks up the cached definition in :attr:`_compatible_by_type` via
    :meth:`_definition_for_phil_path` and applies
    :func:`qttbx.phil.label_for_definition`. Returns ``None`` when the
    path is not one of our tracked file-typed definitions, so callers
    can fall back to the raw phil_path string.
    """
    from qttbx.phil import label_for_definition
    defn = self._definition_for_phil_path(phil_path)
    if defn is None:
      return None
    return label_for_definition(defn)

  def used_for_with_labels(self, row):
    """Return ``[(phil_path, label), ...]`` for the bindings at ``row``.

    ``label`` is :func:`qttbx.phil.label_for_definition` applied to the
    cached definition (i.e. ``.short_caption`` when present, otherwise
    the prettified parameter name). When the definition is not in the
    cache (an out-of-universe phil_path), ``label`` falls back to the
    raw ``phil_path``. Stale rows return a single tuple for their bound
    phil_path so the delegate can render them identically.

    Spec section 5 / section 6.1: chips show
    ``label_for_definition(definition)``, not the raw dotted PHIL path.

    Parameters
    ----------
    row : int
      Row index, 0 .. ``rowCount() - 1``.

    Returns
    -------
    list of (str, str)
      Pairs of ``(phil_path, label)``. A fresh list; callers may mutate
      it freely.
    """
    if 0 <= row < len(self._rows):
      out = []
      for phil_path in self.used_for(row):
        label = self._label_for_phil_path(phil_path) or phil_path
        out.append((phil_path, label))
      return out
    sr = self._stale_at(row)
    if sr is not None:
      label = self._label_for_phil_path(sr.phil_path) or sr.phil_path
      return [(sr.phil_path, label)]
    return []

  def has_compatible_params(self, row):
    """Return True iff the row's data type has at least one compatible
    PHIL path definition (i.e., a definition with ``.type = path`` and
    a ``.style = file_type:<x>`` token that maps to this row's data
    type via :data:`iotbx.data_manager.data_manager_type`).

    Stale rows always report ``True`` so the user can still see the
    stale chip and click the unbind glyph.
    """
    if 0 <= row < len(self._rows):
      _fn, data_type = self._rows[row]
      return bool(self._compatible_by_type.get(data_type))
    return self._stale_at(row) is not None

  # ----- root directory -----

  def root(self):
    """Return the current root directory (normalized), or None."""
    return self._root

  def set_root(self, root):
    """Set the root directory used for filename display.

    Filenames inside ``root`` (including ``root`` itself) render as
    paths relative to it; filenames outside continue to render in their
    full normalized form. Internal state (cache keys, ``_rows``, stale
    rows) is unaffected -- only the displayed string changes.

    Parameters
    ----------
    root : str or None
      A directory path, or None to clear the root and render all
      filenames in their full form.
    """
    new_root = normalize_path(root) if root else None
    if new_root == self._root:
      return
    self._root = new_root
    # Repaint every Filename cell (sorting may also change for the
    # Filename column; the view will re-poll on dataChanged).
    n = self.rowCount()
    if n > 0:
      self.dataChanged.emit(
        self.index(0, self.COL_FILENAME),
        self.index(n - 1, self.COL_FILENAME),
        [Qt.DisplayRole])

  def filename_for_row(self, row):
    """Return the full normalized filename for ``row``.

    Always returns the internal canonical (absolute, normalized) path,
    independent of any root setting. Use this -- not ``data(row,
    COL_FILENAME)`` -- whenever you need a filename to pass back into
    the DataManager, PhilModel, or any of the widget's public API
    methods. Stale rows return their stored normalized filename.

    Returns
    -------
    str or None
      The full normalized filename, or ``None`` for an out-of-range row.
    """
    if 0 <= row < len(self._rows):
      return self._rows[row][0]
    sr = self._stale_at(row)
    return sr.filename if sr is not None else None

  def data_type_for_row(self, row):
    """Return the canonical DataManager data type for ``row``.

    Always returns the snake_case key (e.g. ``"miller_array"``,
    ``"model"``), independent of the user-facing pretty label rendered
    by ``data(row, COL_TYPE)``. Use this -- not ``data(row, COL_TYPE)``
    -- whenever you need a data type to feed back into the DataManager
    or :func:`compatible_phil_params`. Stale rows return their
    ``expected_type``.

    Returns
    -------
    str or None
      The canonical data type, or ``None`` for an out-of-range row.
    """
    if 0 <= row < len(self._rows):
      return self._rows[row][1]
    sr = self._stale_at(row)
    return sr.expected_type if sr is not None else None

  def display_filename(self, filename):
    """Return the display string for ``filename`` under the current root.

    When the root is set and ``filename`` is inside it (or equal to it),
    return the relative form. Otherwise return the full normalized form.
    """
    if not self._root:
      return filename
    try:
      common = os.path.commonpath([self._root, filename])
    except ValueError:
      # Different drives on Windows, or one of the paths is empty.
      return filename
    if common != self._root:
      return filename
    rel = os.path.relpath(filename, self._root)
    # commonpath has already excluded the upward-traversal case, so
    # rel never starts with "..". Treat the root itself as ".".
    return rel

  # ----- Stale-row API -----

  def _stale_at(self, row):
    """Return the StaleRow at the given row index, or None for a non-stale row."""
    n_real = len(self._rows)
    if row < n_real:
      return None
    sr_idx = row - n_real
    if 0 <= sr_idx < len(self._stale):
      return self._stale[sr_idx]
    return None

  def is_stale(self, row):
    """Return True if ``row`` is a stale row (vs a normal DataManager row)."""
    return self._stale_at(row) is not None

  def stale_message(self, row):
    """Return the tooltip message for a stale row, or None for non-stale."""
    sr = self._stale_at(row)
    return sr.message if sr is not None else None

  def stale_phil_path(self, row):
    """Return the phil_path the stale row pointed at, or None."""
    sr = self._stale_at(row)
    return sr.phil_path if sr is not None else None

  def append_stale_row(self, stale_row):
    """Append a StaleRow to the table. Emits the appropriate row-insertion signal.

    Guards against ``None`` so accidental calls (e.g. from a reconciliation
    path that already cleared the slot) do not append a ghost entry that
    would silently corrupt :meth:`_stale_at` row arithmetic.
    """
    if stale_row is None:
      return
    new_row_index = self.rowCount()
    self.beginInsertRows(QModelIndex(), new_row_index, new_row_index)
    self._stale.append(stale_row)
    self.endInsertRows()

  def remove_stale_row(self, row):
    """Remove the stale row at the given index. No-op if ``row`` is not stale.

    Uses indexed deletion (``del self._stale[sr_idx]``) rather than
    ``list.remove`` so a future ``StaleRow.__eq__`` override cannot
    cause the wrong sibling stale row to be removed when two stale
    rows happen to compare equal.
    """
    sr = self._stale_at(row)
    if sr is None:
      return
    sr_idx = row - len(self._rows)
    self.beginRemoveRows(QModelIndex(), row, row)
    del self._stale[sr_idx]
    self.endRemoveRows()

  def clear_stale_rows(self):
    """Drop every stale row.

    Spec section 6.1 lists three ways stale rows are removed: the user
    clicks the chip ✕, the user-corrected file is added (handled via
    :meth:`reconcile_stale`), or :meth:`DataManagerWidget.refresh`
    rebuilds everything. This method implements the third path.

    Kept separate from :meth:`refresh` so internal callers
    (e.g. :meth:`DataManagerWidget.add_file`) that re-enumerate the
    DataManager but want any unresolved stale rows preserved for a
    subsequent :meth:`reconcile_stale` pass can do so. No-op when no
    stale rows are present.
    """
    if not self._stale:
      return
    self.beginResetModel()
    self._stale = []
    self.endResetModel()

  def reconcile_stale(self, filename, data_type):
    """Remove stale rows that match (filename, data_type) and whose
    PHIL value still points there. See spec section 4.5.1.

    Parameters
    ----------
    filename : str
      Filename as added (will be normalized internally).
    data_type : str
      DataManager data type the file was added as.
    """
    if not self._stale:
      return
    norm = normalize_path(filename)
    n_real = len(self._rows)
    # Walk stale rows in descending index order so begin/endRemoveRows
    # invariants hold (rowCount must drop monotonically across the
    # remove sequence). _stale.pop() inside the begin/end window is
    # what makes rowCount() consistent with what Qt expects.
    for i in range(len(self._stale) - 1, -1, -1):
      sr = self._stale[i]
      if (sr.filename == norm
          and sr.expected_type == data_type
          and self._phil_path_still_points_to(sr.phil_path, norm)):
        row = n_real + i
        self.beginRemoveRows(QModelIndex(), row, row)
        self._stale.pop(i)
        self.endRemoveRows()

  def _phil_path_still_points_to(self, phil_path, normalized_filename):
    """Return True if the PhilModel value at ``phil_path`` still
    references ``normalized_filename``.

    Handles both scalar (single-valued) and ``.multiple = True`` path
    definitions: for the latter, :py:meth:`qttbx.phil.PhilModel.value_at_path`
    returns a list of strings, any of which may match. Strings and lists
    are distinguished by duck-typing (strings have ``.encode``) so that
    ``isinstance`` is avoided per cctbx convention.

    Parameters
    ----------
    phil_path : str
      Dotted PHIL path of the binding being reconciled.
    normalized_filename : str
      Already-normalized filename to compare against; see
      :func:`qttbx.widgets.data_manager._phil_helpers.normalize_path`.

    Returns
    -------
    bool
      ``True`` when the PhilModel's current value at ``phil_path``
      (after normalization) matches ``normalized_filename``. ``False``
      when no PhilModel is attached, the value is ``None`` / empty, or
      no entry of a ``.multiple`` list matches.
    """
    if self._phil_model is None:
      return False
    val = self._phil_model.value_at_path(phil_path)
    if not val:
      return False
    # ``.multiple = True`` path definitions return a list of values;
    # everything else is a string. Distinguish by duck-typing: strings
    # expose ``.encode``, lists do not.
    if not hasattr(val, "encode"):
      for v in val:
        if v and normalize_path(v) == normalized_filename:
          return True
      return False
    return normalize_path(val) == normalized_filename

  # ----- Qt model API -----

  def rowCount(self, parent=QModelIndex()):
    """Number of rows (DataManager-derived + stale); ``0`` for any
    non-root parent index."""
    if parent.isValid():
      return 0
    return len(self._rows) + len(self._stale)

  def columnCount(self, parent=QModelIndex()):
    """Column count -- always ``4`` (Filename, Type, Used for, Delete).

    When no :class:`qttbx.phil.PhilModel` is attached the widget runs as
    a pure file manager: the Used-for column is still present but stays
    empty (``used_for`` returns ``[]`` and ``has_compatible_params``
    returns ``False``), so no chips or ``+`` button render. The Delete
    column remains visible so files can still be removed.
    """
    if parent.isValid():
      return 0
    return 4

  def data(self, index, role=Qt.DisplayRole):
    """Return the value for ``index`` under ``role``.

    Only ``Qt.DisplayRole`` is handled in the skeleton; later tasks add
    tooltip, decoration, and edit roles. The ``Used for`` column returns
    a list of bound phil_paths (empty when the file is unused); the
    ``Delete`` column returns the empty string and is rendered as a
    button by the item delegate added in Task 11. Stale rows (PHIL
    bindings whose target file is missing or type-mismatched) render
    the missing path, the expected data type, and the bound phil_path
    as a single-element list.
    """
    if not index.isValid():
      return None
    sr = self._stale_at(index.row())
    if sr is not None:
      col = index.column()
      if role == Qt.DisplayRole:
        if col == self.COL_FILENAME:
          return self.display_filename(sr.filename)
        if col == self.COL_TYPE:
          return pretty_data_type(sr.expected_type)
        if col == self.COL_USED_FOR:
          return [sr.phil_path]
        if col == self.COL_DELETE:
          return ""
      return None
    if index.row() >= len(self._rows):
      return None
    filename, data_type = self._rows[index.row()]
    col = index.column()
    if role == Qt.DisplayRole:
      if col == self.COL_FILENAME:
        return self.display_filename(filename)
      if col == self.COL_TYPE:
        return pretty_data_type(data_type)
      if col == self.COL_USED_FOR:
        return self.used_for(index.row())
      if col == self.COL_DELETE:
        return ""
    return None

  def headerData(self, section, orientation, role=Qt.DisplayRole):
    """Column header labels and decorations for the horizontal header.

    The Delete column shows the system's standard trash icon as its
    decoration so the header matches the per-row glyph rendered by the
    delegate. Its display text is empty.
    """
    if orientation != Qt.Horizontal:
      return None
    if role == Qt.DisplayRole:
      if section == self.COL_FILENAME:
        return "Filename"
      if section == self.COL_TYPE:
        return "Type"
      if section == self.COL_USED_FOR:
        return "Used for"
      if section == self.COL_DELETE:
        return ""
      return None
    if role == Qt.DecorationRole and section == self.COL_DELETE:
      return _trash_icon()
    return None

  def sort(self, column, order=Qt.AscendingOrder):
    """Sort normal rows by ``column``.

    Sortable columns:
      - ``COL_FILENAME``: alphabetical by normalized filename.
      - ``COL_TYPE``: by data type, secondary by filename for stability.
      - ``COL_USED_FOR``: by binding count (number of chips). Ties are
        broken by filename so the order is deterministic.

    ``COL_DELETE`` is not sortable; clicks are silently ignored.

    Stale rows are not reordered -- they continue to render after the
    sorted normal rows in their own insertion order, which keeps the
    per-binding stale-row identity intact (spec section 6.1).

    Parameters
    ----------
    column : int
      One of the sortable columns above. Other values are no-ops.
    order : Qt.SortOrder
      ``Qt.AscendingOrder`` or ``Qt.DescendingOrder``.
    """
    if column not in (self.COL_FILENAME, self.COL_TYPE, self.COL_USED_FOR):
      return
    reverse = (order == Qt.DescendingOrder)
    self.layoutAboutToBeChanged.emit()
    if column == self.COL_TYPE:
      # Sort by the user-facing pretty label so the order matches what
      # the user sees; break ties by displayed filename.
      self._rows.sort(
        key=lambda r: (pretty_data_type(r[1]), self.display_filename(r[0])),
        reverse=reverse)
    elif column == self.COL_USED_FOR:
      # Sort by chip count; break ties by displayed filename.
      self._rows.sort(
        key=lambda r: (
          len(self._bindings_by_file_and_type.get(r, ())),
          self.display_filename(r[0])),
        reverse=reverse)
    else:
      # COL_FILENAME: sort by the string the user sees.
      self._rows.sort(key=lambda r: self.display_filename(r[0]),
                      reverse=reverse)
    self.layoutChanged.emit()
