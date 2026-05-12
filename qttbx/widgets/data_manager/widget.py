"""DataManagerWidget -- widget for adding/deleting files via drag-drop
and binding them to PHIL ``path`` parameters with ``.style file_type:<x>``.
"""

import os

from qttbx.qt.QtCore import Qt, QEvent, Signal
from qttbx.qt.QtWidgets import (
  QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QTableView,
  QHeaderView, QMessageBox, QFileDialog)

from libtbx.utils import Sorry

from qttbx.widgets.data_manager._table_model import (
  DataManagerTableModel, StaleRow)
from qttbx.widgets.data_manager._delegate import DataManagerItemDelegate
from qttbx.widgets.data_manager._binding_popup import (
  DataManagerBindingPopup)
from qttbx.widgets.data_manager._phil_helpers import (
  normalize_path, detect_data_type, compatible_phil_params)


class DataManagerWidget(QWidget):
  """File-manager widget bound to a DataManager and (optionally) a PhilModel.

  Combines a ``QTableView`` of the attached
  :class:`iotbx.data_manager.DataManager`'s files with a ``+ add`` button
  and a chip-based ``Used for`` column. When a
  :class:`qttbx.phil.PhilModel` is attached, the widget keeps PHIL
  ``path`` definitions (any with a ``.style = file_type:<x>`` token)
  in sync with user gestures: adding a file emits
  :py:attr:`fileAdded`; binding to a PHIL path emits
  :py:attr:`bindingChanged`; removing a file cascades to clear any
  bindings then emits :py:attr:`fileRemoved`.

  See spec section 3 for the full public contract.

  Signals
  -------
  fileAdded : (str, str)
    ``(normalized_filename, data_type)`` on every successful add.
  fileRemoved : (str, str)
    ``(normalized_filename, data_type)`` on every successful remove.
  bindingChanged : (str, str, bool)
    ``(normalized_filename, phil_path, is_bound)`` on every successful
    bind or unbind.
  """

  fileAdded = Signal(str, str)
  fileRemoved = Signal(str, str)
  bindingChanged = Signal(str, str, bool)

  def __init__(self, parent=None, phil_model=None, data_manager=None,
               root=None, root_label="Root"):
    """Construct the widget.

    Parameters
    ----------
      parent : QWidget, optional
      phil_model : qttbx.phil.PhilModel, optional
        When ``None``, the widget runs in pure file-manager mode and
        the ``Used for`` column is hidden. When provided, the widget
        walks the model for ``file_type:*`` path definitions and keeps
        PHIL <-> DataManager in sync.
      data_manager : iotbx.data_manager.DataManager, optional
        When ``None`` a fresh :class:`iotbx.data_manager.DataManager`
        is constructed and exposed via :py:attr:`data_manager`.
      root : str, optional
        Optional root directory. When set, filenames inside the root
        render as paths relative to it; filenames outside the root
        keep their full normalized form. The root is shown in a small
        bar at the top of the widget with a button to change it.
      root_label : str, optional
        The text shown next to the root path. Defaults to ``"Root"``;
        can be changed at construction time or runtime via
        :py:meth:`set_root_label`.
    """
    QWidget.__init__(self, parent)
    if data_manager is None:
      from iotbx.data_manager import DataManager
      data_manager = DataManager()
    self._data_manager = data_manager
    self._phil_model = phil_model

    self._table_model = DataManagerTableModel(self)
    self._table_model.attach(data_manager, phil_model=phil_model)
    if root is not None:
      self._table_model.set_root(root)
    self._root_label_text = root_label
    self._auto_import()

    self._delegate = DataManagerItemDelegate(self)
    self._delegate.unbindClicked.connect(self._on_unbind_clicked)
    self._delegate.addChipClicked.connect(self._on_add_clicked)
    self._delegate.deleteClicked.connect(self._on_delete_clicked)

    self._table = QTableView(self)
    self._table.setModel(self._table_model)
    self._table.setItemDelegate(self._delegate)
    # Per-column resize modes:
    #   Filename / Type -- Interactive (user can drag header dividers),
    #   Used for        -- Stretch (absorbs remaining width),
    #   Delete          -- Fixed at a narrow width that fits the icon.
    h = self._table.horizontalHeader()
    h.setSectionResizeMode(DataManagerTableModel.COL_FILENAME,
                           QHeaderView.Interactive)
    h.setSectionResizeMode(DataManagerTableModel.COL_TYPE,
                           QHeaderView.Interactive)
    h.setSectionResizeMode(DataManagerTableModel.COL_USED_FOR,
                           QHeaderView.Stretch)
    h.setSectionResizeMode(DataManagerTableModel.COL_DELETE,
                           QHeaderView.Fixed)
    self._table.setColumnWidth(DataManagerTableModel.COL_DELETE, 28)
    # Enable click-to-sort. Filename, Type, and Used for are all
    # sortable; Delete is not. We catch attempts to put the sort
    # indicator on Delete and revert to the previous valid column.
    self._table.setSortingEnabled(True)
    h.setSortIndicatorShown(True)
    self._last_sort_section = -1
    self._last_sort_order = Qt.AscendingOrder
    self._reverting_sort = False
    h.sortIndicatorChanged.connect(self._on_sort_indicator_changed)
    # Row height tracks the delegate's sizeHint so Used-for cells with
    # multiple chips grow vertically. Setting it on the vertical header
    # works even though the header itself is hidden.
    self._table.verticalHeader().setSectionResizeMode(
      QHeaderView.ResizeToContents)
    self._table.verticalHeader().setVisible(False)
    self._table.setSelectionMode(QTableView.SingleSelection)
    self._table.setSelectionBehavior(QTableView.SelectRows)
    # Clear the selection when the user clicks empty space below the
    # last row (Qt's default keeps the selection indefinitely).
    self._table.viewport().installEventFilter(self)
    # Belt-and-suspenders: explicitly recompute affected row heights when
    # binding state changes, so chip-stack rows grow immediately.
    self._table_model.dataChanged.connect(self._on_table_data_changed)
    self._table_model.rowsInserted.connect(
      lambda *_a: self._table.resizeRowsToContents())
    self._table_model.modelReset.connect(
      self._table.resizeRowsToContents)

    self._popup = DataManagerBindingPopup(self)
    self._popup.bindToggled.connect(self._on_popup_bind_toggled)
    self._popup_row = -1

    self._add_btn = QPushButton(u"Add files…", parent=self)
    self._add_btn.clicked.connect(self._on_add_button)

    # Top bar: root label, Set... button, stretch, Add files... button.
    self._root_label_widget = QLabel(self._render_root_label(), parent=self)
    self._root_set_btn = QPushButton(u"Set…", parent=self)
    self._root_set_btn.clicked.connect(self._on_root_button)
    top = QHBoxLayout()
    top.addWidget(self._root_label_widget)
    top.addWidget(self._root_set_btn)
    top.addStretch(1)
    top.addWidget(self._add_btn)
    layout = QVBoxLayout(self)
    layout.addLayout(top)
    layout.addWidget(self._table)
    self.setAcceptDrops(True)

  # ----- properties -----

  @property
  def data_manager(self):
    """The :class:`iotbx.data_manager.DataManager` this widget reflects."""
    return self._data_manager

  @property
  def phil_model(self):
    """The attached :class:`qttbx.phil.PhilModel`, or ``None``."""
    return self._phil_model

  @property
  def root(self):
    """The current root directory (normalized), or ``None``."""
    return self._table_model.root()

  @property
  def root_label(self):
    """The label text shown next to the root path in the root bar."""
    return self._root_label_text

  # ----- public API -----

  def set_root(self, root):
    """Set the root directory.

    Filenames inside ``root`` render as relative paths; filenames
    outside continue to render in their full normalized form. Pass
    ``None`` to clear.
    """
    self._table_model.set_root(root)
    self._root_label_widget.setText(self._render_root_label())

  def set_root_label(self, label):
    """Change the text shown next to the root path."""
    self._root_label_text = label
    self._root_label_widget.setText(self._render_root_label())

  def _render_root_label(self):
    """Return the displayed text for the root bar."""
    r = self._table_model.root()
    return "%s: %s" % (self._root_label_text, r if r else "(none)")

  def _on_root_button(self):
    """Slot: prompt the user for a new root directory."""
    start = self._table_model.root() or os.path.expanduser("~")
    chosen = QFileDialog.getExistingDirectory(
      self, "Choose %s directory" % self._root_label_text, start)
    if chosen:
      self.set_root(chosen)

  def add_file(self, filename, data_type=None):
    """Add a file to the DataManager and emit :py:attr:`fileAdded`.

    Detects the data type (when not supplied) via
    :func:`qttbx.widgets.data_manager._phil_helpers.detect_data_type`,
    dispatches to the appropriate ``process_<type>_file`` method on
    the DataManager, refreshes the table model, and reconciles any
    matching stale rows.

    Parameters
    ----------
      filename : str
        Path to add. Will be normalized internally.
      data_type : str, optional
        DataManager data type. When omitted, ``detect_data_type`` is
        called on the normalized path.

    Returns
    -------
      str
        The resolved data type.

    Raises
    ------
      libtbx.utils.Sorry
        When the data type cannot be detected, or when the DataManager
        does not provide a ``process_<data_type>_file`` method.
    """
    norm = normalize_path(filename)
    if data_type is None:
      data_type = detect_data_type(norm)
      if data_type is None:
        raise Sorry(
          "Could not detect data type for %s" % os.path.basename(norm))
    processor = getattr(
      self._data_manager, "process_%s_file" % data_type, None)
    if processor is None:
      raise Sorry(
        "DataManager does not support data type %r" % data_type)
    processor(norm)
    self._table_model.refresh()
    self._table_model.reconcile_stale(norm, data_type)
    self.fileAdded.emit(norm, data_type)
    return data_type

  def remove_file(self, filename):
    """Remove a file from the DataManager and clear its bindings.

    Walks the table model for rows matching the normalized filename,
    cascades :py:meth:`unbind` for every PHIL path that currently
    references the file, calls the DataManager's
    ``remove_<data_type>`` method (when present), and emits
    :py:attr:`fileRemoved`. No-op when no row matches.

    Parameters
    ----------
      filename : str
        Path to remove. Will be normalized internally.
    """
    norm = normalize_path(filename)
    rows = [r for r in range(self._table_model.rowCount())
            if not self._table_model.is_stale(r)
            and self._table_model.filename_for_row(r) == norm]
    for r in rows:
      data_type = self._table_model.data_type_for_row(r)
      for phil_path in list(self._table_model.used_for(r)):
        self.unbind(norm, phil_path)
      remover = getattr(self._data_manager, "remove_%s" % data_type, None)
      if remover is not None:
        try:
          remover(norm)
        except Exception:
          # DataManager remove_<type> implementations vary across
          # mixins (e.g. by-name vs. by-handle); swallow individual
          # failures so the table row still disappears from view.
          pass
      self.fileRemoved.emit(norm, data_type)
    self._table_model.refresh()

  def bind(self, filename, phil_path):
    """Bind ``filename`` to the PHIL parameter at ``phil_path``.

    For a non-``.multiple`` path definition: sets the scalar to the
    normalized filename. Raises :class:`libtbx.utils.Sorry` when the
    target is already bound to a different filename. Idempotent on
    the same filename.

    For a ``.multiple = True`` path definition: appends the normalized
    filename to the definition's list of instances. Idempotent on the
    same filename (does not append duplicates).

    **File existence is not verified.** ``bind`` writes the PHIL value
    even when ``filename`` is missing on disk or has not been added to
    the DataManager. The widget will *not* show a row for the binding
    in that case; on the next construction, :py:meth:`_auto_import`
    will surface the missing file as a stale row. If you want
    immediate visual feedback, call :py:meth:`add_file` first and let
    its :class:`Sorry` propagate.

    Parameters
    ----------
      filename : str
        Path to bind. Will be normalized internally.
      phil_path : str
        Dotted PHIL path of a path-typed definition.

    Raises
    ------
      libtbx.utils.Sorry
        When no PhilModel is attached, or the non-``.multiple`` target
        is already bound to a different filename.
    """
    if self._phil_model is None:
      raise Sorry("No PhilModel attached; cannot bind.")
    norm = normalize_path(filename)
    defn = self._phil_model.definition_for_path(phil_path)
    is_multiple = bool(getattr(defn, "multiple", False))
    if is_multiple:
      instances = self._phil_model.instances_for_path(phil_path)
      if any(v is not None and normalize_path(v) == norm
             for v in instances):
        return
      self._phil_model.append_instance(phil_path, norm)
    else:
      current = self._phil_model.value_at_path(phil_path)
      if current == "":
        current = None
      if current is not None and normalize_path(current) == norm:
        return
      if current is not None:
        raise Sorry(
          "%s is already bound to %s" % (phil_path, current))
      self._phil_model.set_value_at_path(phil_path, norm)
    self.bindingChanged.emit(norm, phil_path, True)

  def unbind(self, filename, phil_path):
    """Remove a binding. Idempotent.

    For a non-``.multiple`` path: clears the scalar when it currently
    equals ``filename`` (after normalization). Silent no-op otherwise.

    For a ``.multiple = True`` path: removes the matching instance.
    Silent no-op when no instance matches.

    Parameters
    ----------
      filename : str
        Path to unbind. Will be normalized internally.
      phil_path : str
        Dotted PHIL path of a path-typed definition.
    """
    if self._phil_model is None:
      return
    norm = normalize_path(filename)
    defn = self._phil_model.definition_for_path(phil_path)
    is_multiple = bool(getattr(defn, "multiple", False))
    if is_multiple:
      removed = self._phil_model.remove_instance_with_value(phil_path, norm)
      if not removed:
        return
    else:
      current = self._phil_model.value_at_path(phil_path)
      if current is None or normalize_path(current) != norm:
        return
      self._phil_model.set_value_at_path(phil_path, None)
    self.bindingChanged.emit(norm, phil_path, False)

  def refresh(self):
    """Re-enumerate the attached DataManager and rebuild everything.

    Per spec section 6.1, an external :py:meth:`refresh` removes stale
    rows in addition to rebuilding from the DataManager. Internal
    callers that need the file list refreshed but want stale rows
    preserved for a follow-up :py:meth:`DataManagerTableModel.reconcile_stale`
    (e.g. :py:meth:`add_file`) call ``_table_model.refresh()`` directly.
    """
    self._table_model.refresh()
    self._table_model.clear_stale_rows()

  # ----- internal handlers -----

  def eventFilter(self, obj, event):
    """Clear the table selection on clicks that land outside any row.

    Installed on ``self._table.viewport()`` so we see mouse presses
    before the view consumes them. Qt's default selection model
    leaves the previously-selected row highlighted after empty-space
    clicks; we explicitly clear it to match modern UI conventions.
    """
    if (obj is self._table.viewport()
        and event.type() == QEvent.MouseButtonPress):
      # QMouseEvent.pos() is deprecated in Qt6; use position().toPoint() when
      # available, fall back to pos() on PySide2.
      if hasattr(event, "position"):
        p = event.position().toPoint()
      else:
        p = event.pos()
      index = self._table.indexAt(p)
      if not index.isValid():
        self._table.clearSelection()
    return False  # don't consume; normal event flow continues

  def _on_sort_indicator_changed(self, section, order):
    """Revert the indicator to the previous valid section when the user
    clicks the Delete column header (Delete is not sortable). For valid
    sections, remember the (section, order) pair so we can revert to it
    on the next bad click."""
    if self._reverting_sort:
      return
    if section == DataManagerTableModel.COL_DELETE:
      self._reverting_sort = True
      self._table.horizontalHeader().setSortIndicator(
        self._last_sort_section, self._last_sort_order)
      self._reverting_sort = False
      return
    self._last_sort_section = section
    self._last_sort_order = order

  def _on_table_data_changed(self, topLeft, bottomRight, roles=None):
    """Recompute row heights for cells whose Used-for chip count may
    have changed. Without this, ``QHeaderView.ResizeToContents`` does
    not always re-poll on per-cell ``dataChanged`` and rows that gain
    chips stay at their initial single-line height.
    """
    for row in range(topLeft.row(), bottomRight.row() + 1):
      self._table.resizeRowToContents(row)

  def _on_unbind_clicked(self, row, phil_path):
    """Delegate-signal slot: ``unbindClicked`` -> :py:meth:`unbind`.

    Stale rows clear the PHIL value directly and remove the stale row
    from the table model (per spec section 6.1); live rows route through
    the public :py:meth:`unbind`.
    """
    filename = self._table_model.filename_for_row(row)
    if self._table_model.is_stale(row):
      if self._phil_model is None:
        return
      self._phil_model.set_value_at_path(phil_path, None)
      self._table_model.remove_stale_row(row)
      return
    self.unbind(filename, phil_path)

  def _on_delete_clicked(self, row):
    """Delegate-signal slot: ``deleteClicked`` -> :py:meth:`remove_file`.

    For stale rows, clears the bound PHIL value and removes the stale
    row from the table model; live rows route through
    :py:meth:`remove_file`.
    """
    filename = self._table_model.filename_for_row(row)
    if self._table_model.is_stale(row):
      if self._phil_model is None:
        return
      phil_path = self._table_model.stale_phil_path(row)
      if phil_path:
        self._phil_model.set_value_at_path(phil_path, None)
      self._table_model.remove_stale_row(row)
      return
    self.remove_file(filename)

  def _on_add_clicked(self, row, anchor_rect):
    """Delegate-signal slot: ``addChipClicked`` -> open binding popup.

    Builds the list of compatible PHIL paths (via
    :func:`compatible_phil_params`), populates the popup, and shows it
    anchored at the ``+ add`` chip's bottom-left in viewport coords.
    No-op when no PhilModel is attached or no compatible parameters
    exist.
    """
    if self._phil_model is None:
      return
    filename = self._table_model.filename_for_row(row)
    data_type = self._table_model.data_type_for_row(row)
    candidates = self._build_popup_candidates(filename, data_type)
    if not candidates:
      return
    self._popup.populate(candidates)
    self._popup_row = row
    viewport = self._table.viewport()
    anchor_global = viewport.mapToGlobal(anchor_rect.bottomLeft())
    self._popup.move(anchor_global)
    self._popup.show()

  def _build_popup_candidates(self, filename, data_type):
    """Return ``(phil_path, label, checked, disabled, tooltip)`` candidates.

    Combines :func:`compatible_phil_params` (the universe of paths for
    the row's data type) with the row's currently-bound paths
    (:py:meth:`DataManagerTableModel.used_for`) to produce a
    five-tuple per candidate. Non-``.multiple`` paths already bound
    to a different filename are reported as ``disabled=True`` with a
    tooltip naming the other filename.
    """
    from qttbx.phil import label_for_definition
    norm = normalize_path(filename)
    candidates = []
    row = self._row_for(norm, data_type)
    bound_set = set(self._table_model.used_for(row)) if row >= 0 else set()
    for phil_path, defn in compatible_phil_params(
        self._phil_model, data_type):
      checked = phil_path in bound_set
      is_multiple = bool(getattr(defn, "multiple", False))
      disabled = False
      tooltip = ""
      if not is_multiple and not checked:
        cur = self._phil_model.value_at_path(phil_path)
        if cur is not None and cur != "" and normalize_path(cur) != norm:
          disabled = True
          tooltip = "already bound to %s" % cur
      label = label_for_definition(defn) or phil_path
      candidates.append((phil_path, label, checked, disabled, tooltip))
    return candidates

  def _row_for(self, normalized_filename, data_type):
    """Return the live row index matching ``(filename, data_type)``, or ``-1``."""
    for r in range(self._table_model.rowCount()):
      if self._table_model.is_stale(r):
        continue
      fn = self._table_model.filename_for_row(r)
      dt = self._table_model.data_type_for_row(r)
      if fn == normalized_filename and dt == data_type:
        return r
    return -1

  def _on_popup_bind_toggled(self, phil_path, want_bound):
    """Popup-signal slot: ``bindToggled`` -> :py:meth:`bind` /
    :py:meth:`unbind`.

    Catches :class:`libtbx.utils.Sorry` from :py:meth:`bind` and
    surfaces it as a non-modal :class:`QMessageBox` warning so the
    user can correct the conflict.
    """
    row = self._popup_row
    if row < 0:
      return
    filename = self._table_model.filename_for_row(row)
    try:
      if want_bound:
        self.bind(filename, phil_path)
      else:
        self.unbind(filename, phil_path)
    except Sorry as e:
      QMessageBox.warning(self, "Cannot bind", str(e))

  def _on_add_button(self):
    """``Add files…`` button slot: open a :class:`QFileDialog`.

    Each selected path is fed through :py:meth:`add_file`; failures
    are collected and reported in a single :class:`QMessageBox`
    warning after the loop.
    """
    paths, _filter = QFileDialog.getOpenFileNames(
      self, "Add files", "", "All files (*)")
    failed = []
    for p in paths:
      try:
        self.add_file(p)
      except Sorry as e:
        failed.append((p, str(e)))
    if failed:
      msg = "\n".join("%s: %s" % (os.path.basename(p), m)
                      for p, m in failed)
      QMessageBox.warning(self, "Could not add files", msg)

  # ----- drag/drop -----

  def dragEnterEvent(self, event):
    """Accept drags that carry local file URLs.

    Highlights the table view with a dashed red border (per qttbx
    convention; see ``qttbx.widgets.phil.text_base``) while the drag
    is hovering over the widget.
    """
    md = event.mimeData()
    if md.hasUrls() and all(u.isLocalFile() for u in md.urls()):
      event.acceptProposedAction()
      self._set_drop_highlight(True)

  def dragMoveEvent(self, event):
    """Accept the proposed action while the drag is inside the widget."""
    event.acceptProposedAction()

  def dragLeaveEvent(self, event):
    """Remove the drop highlight when the drag exits the widget."""
    self._set_drop_highlight(False)

  def dropEvent(self, event):
    """Accept dropped local files and add them to the DataManager.

    Each URL is fed through :py:meth:`add_file`; failures (both
    :class:`libtbx.utils.Sorry` for unrecognized types and any other
    :class:`Exception` raised by the DataManager processors) are
    collected and reported in a single :class:`QMessageBox.warning`
    at the end so a multi-file drop with one bad file does not spawn
    a popup per failure.
    """
    md = event.mimeData()
    paths = [u.toLocalFile() for u in md.urls() if u.isLocalFile()]
    failed = []
    for p in paths:
      try:
        self.add_file(p)
      except Sorry as e:
        failed.append((p, str(e)))
      except Exception as e:
        failed.append((p, str(e)))
    self._set_drop_highlight(False)
    if failed:
      msg = "\n".join("%s: %s" % (os.path.basename(p), m)
                      for p, m in failed)
      QMessageBox.warning(self, "Could not add files", msg)
    event.acceptProposedAction()

  def _set_drop_highlight(self, on):
    """Toggle the dashed drop-target border on the table view.

    Palette-aware: bright red in light mode, softer pink-red in dark
    mode so the border remains visible without screaming."""
    if on:
      from qttbx.qt.QtGui import QPalette
      palette = self._table.palette()
      is_dark = palette.color(QPalette.Base).lightness() < 128
      color = "rgb(255, 100, 100)" if is_dark else "rgb(220, 50, 50)"
      self._table.setStyleSheet(
        "QTableView { border: 2px dashed %s; }" % color)
    else:
      self._table.setStyleSheet("")

  # ----- auto-import -----

  def _auto_import(self):
    """Walk the PhilModel and seed the DataManager from non-None paths.

    For every path-typed PHIL definition tagged with a
    ``.style = file_type:<suffix>`` whose ``<suffix>`` maps to a known
    DataManager data type (per :data:`iotbx.data_manager.data_manager_type`),
    resolve the current value via
    :py:meth:`qttbx.phil.PhilModel.value_at_path` and either add the file
    to the DataManager or append a stale row when the file is missing or
    has the wrong type. ``.multiple = True`` paths return a list and are
    iterated. Implements spec section 4.5: PHIL state takes precedence
    over an empty DataManager, but a partially-populated DataManager
    survives intact because :py:meth:`_auto_import_one` calls the
    DataManager's ``process_<type>_file`` method, which is idempotent.

    Never raises; per-value failures during DataManager processing are
    captured by :py:meth:`_auto_import_one` and surfaced as stale rows.
    Final step is a table-model :py:meth:`refresh` so any newly-added
    files appear as normal rows and the binding cache picks them up.
    Stale rows survive this refresh by design (per spec section 6.1;
    only :py:meth:`refresh` on the widget clears them).
    """
    if self._phil_model is None:
      return
    from iotbx.data_manager import data_manager_type
    from qttbx.widgets.data_manager._phil_helpers import (
      parse_file_type_style)
    # iter_definitions() yields .multiple = True paths once per instance
    # (template + N instances); dedupe by phil_path so a .multiple slot
    # with one missing file produces a single stale row, not N rows.
    seen = set()
    for phil_path, defn in self._phil_model.iter_definitions():
      if phil_path in seen:
        continue
      seen.add(phil_path)
      phil_type = getattr(getattr(defn, "type", None), "phil_type", None)
      if phil_type != "path":
        continue
      suffix = parse_file_type_style(getattr(defn, "style", None) or "")
      if suffix is None:
        continue
      expected_dt = data_manager_type.get(suffix)
      if expected_dt is None:
        continue
      value = self._phil_model.value_at_path(phil_path)
      if not value:
        continue
      # For .multiple paths, value is a list; iterate. Distinguish by
      # duck-typing per cctbx convention (no isinstance).
      if hasattr(value, "encode"):
        values = [value]
      else:
        values = list(value)
      for v in values:
        if not v:
          continue
        self._auto_import_one(phil_path, expected_dt, v)
    # Re-enumerate so newly-added DataManager files appear as normal
    # rows. The cache rebuild inside refresh() repopulates bindings
    # against the (now-present) files.
    self._table_model.refresh()

  def _auto_import_one(self, phil_path, expected_dt, value):
    """Resolve a single PHIL path value, adding to DM or appending stale.

    Three outcomes:

    1. The file does not exist on disk -- append a "file not found"
       stale row.
    2. The file exists but :func:`detect_data_type` reports a type
       different from the PHIL slot's declared ``expected_dt`` -- append
       a "type mismatch" stale row.
    3. The file exists and the type matches -- dispatch to the
       DataManager's ``process_<expected_dt>_file`` method. If the
       DataManager lacks the processor or the call raises, fall back to
       a stale row carrying the exception text.

    Never raises; all error paths surface as stale rows. Caller is
    expected to call :py:meth:`DataManagerTableModel.refresh` afterward
    to repaint the regular rows.

    Parameters
    ----------
    phil_path : str
      The PHIL path being processed (passed through into any stale row
      so the user knows which slot triggered the failure).
    expected_dt : str
      The DataManager data type the slot expects.
    value : str
      The raw PHIL value (will be normalized internally).
    """
    norm = normalize_path(value)
    if not os.path.exists(norm):
      self._table_model.append_stale_row(StaleRow(
        filename=norm, expected_type=expected_dt,
        phil_path=phil_path, message="file not found"))
      return
    detected = detect_data_type(norm)
    if detected != expected_dt:
      self._table_model.append_stale_row(StaleRow(
        filename=norm, expected_type=expected_dt,
        phil_path=phil_path,
        message="type mismatch: expected %s" % expected_dt))
      return
    try:
      processor = getattr(
        self._data_manager, "process_%s_file" % expected_dt, None)
      if processor is None:
        return
      processor(norm)
    except Exception as e:
      self._table_model.append_stale_row(StaleRow(
        filename=norm, expected_type=expected_dt,
        phil_path=phil_path, message=str(e)))
