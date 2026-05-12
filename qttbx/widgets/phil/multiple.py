"""MultipleWidget, RepeatableScopeWidget -- handle PHIL .multiple = True.

``MultipleWidget`` wraps a definition with ``.multiple = True`` (and a
non-list-native type) -- value is a Python list, edited via N inner widgets.

``RepeatableScopeWidget`` is the form-side handler for a SCOPE with
``.multiple = True``. See its class docstring for details.
"""

from qttbx.qt.QtCore import QPersistentModelIndex, QSignalBlocker
from qttbx.qt.QtWidgets import (
  QFormLayout, QHBoxLayout, QListWidget, QListWidgetItem, QPushButton,
  QTabWidget, QVBoxLayout, QWidget,
)

from qttbx.widgets.phil import PhilWidget


class MultipleWidget(PhilWidget):
  """PHIL definition-multiple editor.

  Wraps N inner widgets (one per element) in a ``QListWidget`` with a
  ``[+] [-] [up] [down]`` toolbar. ``inner_cls`` is a single-value
  ``PhilWidget`` subclass to instantiate per row; if None, looked up from
  the registry by ``definition.type.phil_type``.

  The ``value()`` is a Python list of per-row values. ``setValue(list)``
  rebuilds the row list. ``isValid()`` is True iff every row is valid AND
  the row count is in ``[size_min, size_max]``.
  """

  def __init__(self, definition, inner_cls=None, parent=None,
               size_min=None, size_max=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with ``.multiple == True`` and a type not in
        ``{ints, floats, strings, words}``.
      inner_cls: type, optional
        ``PhilWidget`` subclass to use per row. Looked up from the registry
        if None.
      parent: QWidget, optional
      size_min, size_max: int, optional
        Inclusive bounds on row count; not enforced if None.
    """
    PhilWidget.__init__(self, definition, parent)
    if inner_cls is None:
      from qttbx.widgets.phil import _widget_registry
      inner_cls = _widget_registry[definition.type.phil_type]
    self._inner_cls = inner_cls
    self._size_min = size_min
    self._size_max = size_max
    # Toolbar buttons.
    self._add_btn = QPushButton("+", parent=self)
    self._del_btn = QPushButton("−", parent=self)
    self._up_btn = QPushButton("↑", parent=self)
    self._down_btn = QPushButton("↓", parent=self)
    for b in (self._add_btn, self._del_btn, self._up_btn, self._down_btn):
      b.setMaximumWidth(32)
    # List of rows (each row is an inner_cls widget, hosted in a QListWidget item).
    self._list = QListWidget(self)
    # Layout.
    bar = QHBoxLayout()
    bar.setContentsMargins(0, 0, 0, 0)
    bar.addWidget(self._add_btn)
    bar.addWidget(self._del_btn)
    bar.addWidget(self._up_btn)
    bar.addWidget(self._down_btn)
    bar.addStretch(1)
    layout = QVBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._list)
    layout.addLayout(bar)
    # Wire toolbar.
    self._add_btn.clicked.connect(self._on_add)
    self._del_btn.clicked.connect(self._on_remove)
    self._up_btn.clicked.connect(self._on_move_up)
    self._down_btn.clicked.connect(self._on_move_down)
    # Initial value from the definition.
    # For a multi definition, ``extract()`` may return:
    #   - a list-shaped object (e.g. when the master phil has multiple instances)
    #   - a scalar default (e.g. ``val = 1`` with .multiple = True)
    #   - None (no default specified)
    # Normalize to a Python list.
    ex = self.definition.extract()
    if ex is None:
      initial = []
    else:
      try:
        initial = list(ex)
      except TypeError:
        initial = [ex]
    self.setValue(initial)

  def value(self):
    """Return the list of per-row values."""
    out = []
    for i in range(self._list.count()):
      row_widget = self._list.itemWidget(self._list.item(i))
      out.append(row_widget.value())
    return out

  def setValue(self, values):
    """Rebuild the row list to match ``values`` (a Python list)."""
    if values is None:
      values = []
    blocker = QSignalBlocker(self._list)
    for i in range(self._list.count()):
      item = self._list.item(i)
      inner = self._list.itemWidget(item)
      self._list.removeItemWidget(item)
      if inner is not None:
        inner.deleteLater()
    self._list.clear()
    for v in values:
      self._append_row(v)
    del blocker

  def isValid(self):
    """True iff every row is valid AND the row count is in [size_min, size_max]."""
    n = self._list.count()
    if self._size_min is not None and n < self._size_min:
      return False
    if self._size_max is not None and n > self._size_max:
      return False
    for i in range(n):
      row_widget = self._list.itemWidget(self._list.item(i))
      if not row_widget.isValid():
        return False
    return True

  def errorString(self):
    """Return a human-readable error referencing the offending row(s)."""
    n = self._list.count()
    if self._size_min is not None and n < self._size_min:
      return "{n} rows below minimum {m}".format(n=n, m=self._size_min)
    if self._size_max is not None and n > self._size_max:
      return "{n} rows exceed maximum {m}".format(n=n, m=self._size_max)
    errs = []
    for i in range(n):
      row_widget = self._list.itemWidget(self._list.item(i))
      if not row_widget.isValid():
        errs.append("row {i}: {e}".format(i=i, e=row_widget.errorString()))
    return "; ".join(errs)

  def setReadOnly(self, ro=True):
    """Disable editing of all rows and the add/remove/move buttons."""
    for i in range(self._list.count()):
      row_widget = self._list.itemWidget(self._list.item(i))
      row_widget.setReadOnly(bool(ro))
    for b in (self._add_btn, self._del_btn, self._up_btn, self._down_btn):
      b.setEnabled(not bool(ro))

  def _append_row(self, value):
    """Append one row hosting an ``inner_cls`` widget initialized to ``value``."""
    row_widget = self._inner_cls(self.definition, parent=self._list)
    blocker = QSignalBlocker(row_widget)
    row_widget.setValue(value)
    del blocker
    row_widget.valueChanged.connect(self._on_row_changed)
    item = QListWidgetItem(self._list)
    item.setSizeHint(row_widget.sizeHint())
    self._list.addItem(item)
    self._list.setItemWidget(item, row_widget)

  def _on_row_changed(self, _v):
    """Forward a row commit as the wrapper's valueChanged."""
    self.valueChanged.emit(self.value())

  def _on_add(self):
    """[+] toolbar action: append a default-valued row."""
    blocker = QSignalBlocker(self._list)
    # Append with whatever the inner class produces by default
    # (calls setValue(None) which most widgets treat as the empty/Auto state).
    self._append_row(None)
    del blocker
    self.valueChanged.emit(self.value())

  def _on_remove(self):
    """[-] toolbar action: remove the currently-selected row."""
    row = self._list.currentRow()
    if row < 0:
      return
    blocker = QSignalBlocker(self._list)
    item = self._list.item(row)
    inner = self._list.itemWidget(item)
    self._list.removeItemWidget(item)
    self._list.takeItem(row)
    if inner is not None:
      inner.deleteLater()
    del blocker
    self.valueChanged.emit(self.value())

  def _on_move_up(self):
    """[up] toolbar action: swap the selected row with the one above."""
    row = self._list.currentRow()
    if row <= 0:
      return
    self._swap_rows(row, row - 1)
    self._list.setCurrentRow(row - 1)
    self.valueChanged.emit(self.value())

  def _on_move_down(self):
    """[down] toolbar action: swap the selected row with the one below."""
    row = self._list.currentRow()
    if row < 0 or row >= self._list.count() - 1:
      return
    self._swap_rows(row, row + 1)
    self._list.setCurrentRow(row + 1)
    self.valueChanged.emit(self.value())

  def _swap_rows(self, a, b):
    """Swap row widgets at positions ``a`` and ``b``."""
    wa = self._list.itemWidget(self._list.item(a))
    wb = self._list.itemWidget(self._list.item(b))
    va = wa.value()
    vb = wb.value()
    blocker = QSignalBlocker(self._list)
    wa.setValue(vb)
    wb.setValue(va)
    del blocker


class RepeatableScopeWidget(QWidget):
  """Form-side editor for a PHIL scope with ``.multiple = True``.

  Each tab is a ``QFormLayout`` of ``PhilField``s constructed with
  ``persistent_index=`` so they stay bound to the right instance even
  after sibling tabs are removed.

  Corner-widget buttons:
    [+] add a new instance via ``model.add_scope_instance(full_path)``
    [x] remove the current tab's instance via
        ``model.remove_scope_instance(full_path, tab_index)``

  External sync: listens to ``model.rowsInserted`` / ``model.rowsRemoved``
  on the parent index so externally-driven adds/removes (from another view,
  programmatic mutations) keep the tab strip in sync.
  """

  def __init__(self, model, full_path, parent=None):
    """
    Parameters
    ----------
      model: qttbx.phil.PhilModel
      full_path: str
        Dotted path of the multi-scope (e.g. "ncs_group").
      parent: QWidget, optional
    """
    QWidget.__init__(self, parent)
    self._model = model
    self._full_path = full_path
    # Locate the parent QModelIndex (the master scope or an outer scope).
    self._parent_qmi = self._compute_parent_qmi()
    # Tabs.
    self._tabs = QTabWidget(self)
    # Corner widget: [+] and [x] buttons.
    self._add_btn = QPushButton("+", parent=self)
    self._del_btn = QPushButton("×", parent=self)   # MULTIPLICATION SIGN
    for b in (self._add_btn, self._del_btn):
      b.setMaximumWidth(32)
    corner = QWidget(self)
    corner_layout = QHBoxLayout(corner)
    corner_layout.setContentsMargins(0, 0, 0, 0)
    corner_layout.addWidget(self._add_btn)
    corner_layout.addWidget(self._del_btn)
    self._tabs.setCornerWidget(corner)
    # Layout.
    layout = QVBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._tabs)
    # Wire buttons.
    self._add_btn.clicked.connect(self._on_add)
    self._del_btn.clicked.connect(self._on_remove)
    # Listen for external adds/removes on the parent index.
    self._model.rowsInserted.connect(self._on_rows_inserted)
    self._model.rowsRemoved.connect(self._on_rows_removed)
    # Build initial tabs.
    self._rebuild_tabs()

  def _compute_parent_qmi(self):
    """Compute the QModelIndex of the multi-scope's parent in the model tree.

    Delegates to ``PhilModel._parent_item_for_path`` which handles both
    top-level multi-scopes (parent is the master-scope item) and
    intermediate scopes via ``scope_indices`` for nested multi paths.
    """
    parent_item, _ = self._model._parent_item_for_path(self._full_path)
    return self._model.createIndex(parent_item.row(), 0, parent_item)

  def _scope_name(self):
    """Return the last segment of full_path (the multi-scope's name)."""
    return self._full_path.rsplit(".", 1)[-1]

  def _rebuild_tabs(self):
    """Clear and rebuild all tabs from the current model state.

    Each removed tab body is explicitly ``deleteLater()``-ed so the
    ``PhilField`` instances it hosts (which connected to
    ``model.dataChanged`` in their constructor) are destroyed and
    Qt's auto-disconnect cleans up the slot chain. Without this, prior
    generations of ``PhilField`` continue firing on every ``dataChanged``
    emission, costing O(n^2) per ``setData`` over time.
    """
    while self._tabs.count():
      body = self._tabs.widget(0)
      self._tabs.removeTab(0)
      if body is not None:
        body.deleteLater()
    from qttbx.phil import label_for_definition
    parent_item = self._parent_qmi.internalPointer()
    name = self._scope_name()
    for i in range(parent_item.childCount()):
      child = parent_item.child(i)
      if child.definition.name == name and child.instance_index is not None:
        label = label_for_definition(child.definition)
        self._tabs.addTab(self._build_tab_widget(child),
                          "{n} {i}".format(n=label, i=child.instance_index + 1))

  def _build_tab_widget(self, scope_item):
    """Build a tab body holding a QFormLayout of PhilField for each child def."""
    from qttbx.widgets.phil import PhilField
    body = QWidget(self._tabs)
    form = QFormLayout(body)
    for i in range(scope_item.childCount()):
      child = scope_item.child(i)
      if child.definition.is_definition:
        qpi = QPersistentModelIndex(
          self._model.createIndex(child.row(), 1, child))
        field = PhilField(self._model, persistent_index=qpi)
        form.addRow(field)
    return body

  def _on_add(self):
    """Handler for the [+] button: add a scope instance via the model."""
    self._model.add_scope_instance(self._full_path)
    # rowsInserted will trigger _on_rows_inserted -> _rebuild_tabs.

  def _on_remove(self):
    """Handler for the [x] button: remove the active tab's instance."""
    idx = self._tabs.currentIndex()
    if idx < 0:
      return
    self._model.remove_scope_instance(self._full_path, idx)

  def _on_rows_inserted(self, parent, _first, _last):
    """Slot for model.rowsInserted; rebuild if the change touches our parent."""
    if parent == self._parent_qmi:
      self._rebuild_tabs()

  def _on_rows_removed(self, parent, _first, _last):
    """Slot for model.rowsRemoved; rebuild if the change touches our parent."""
    if parent == self._parent_qmi:
      self._rebuild_tabs()
