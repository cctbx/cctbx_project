"""ChoiceWidget -- a PHIL choice editor backed by QComboBox."""

from qttbx.qt.QtCore import Qt, QSignalBlocker
from qttbx.qt.QtWidgets import QComboBox, QHBoxLayout, QListWidget, QListWidgetItem, QVBoxLayout

from qttbx.widgets.phil import PhilWidget


_NONE_LABEL = "(none)"

# Marker stored at Qt.UserRole on the synthetic "(none)" combo entry.
# Distinguishes it from a real choice whose text happens to equal _NONE_LABEL.
_NONE_SENTINEL = object()


class ChoiceWidget(PhilWidget):
  """PHIL single-choice editor backed by a QComboBox.

  The initial choice list is read from ``definition.words`` (with the
  leading ``*`` of the default-marked entry stripped). When ``_allow_none``
  is True a leading ``"(none)"`` entry is prepended so users can clear the
  selection. ``setChoices(list)`` may be called at runtime to swap the
  choice list -- e.g. populate from an MTZ file's array names.
  """

  def __init__(self, definition, parent=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with .type.phil_type == 'choice'. The available
        choices are taken from ``definition.words`` at construction; use
        ``setChoices`` to update them later.
      parent: QWidget, optional
    """
    PhilWidget.__init__(self, definition, parent)
    self._combo = QComboBox(self)
    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._combo)
    # Build initial choice list from definition.words.
    choices = [w.value.lstrip("*") for w in self.definition.words]
    self._populate(choices)
    self._combo.currentIndexChanged.connect(self._on_index_changed)
    # Initial value from definition extract.
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)

  def value(self):
    """Return the current selection.

    Returns
    -------
      value: str or None
        The selected choice text, or None when the synthetic ``(none)`` entry
        (identified via Qt.UserRole sentinel) is selected.
    """
    idx = self._combo.currentIndex()
    if idx < 0:
      return None
    if self._combo.itemData(idx, Qt.UserRole) is _NONE_SENTINEL:
      return None
    return self._combo.currentText()

  def setValue(self, value):
    """Select the entry matching ``value``. Idempotent and signal-blocked.

    Parameters
    ----------
      value: str or None
        ``None`` selects the synthetic ``(none)`` entry (identified by sentinel,
        not by text). Any other string must match an existing concrete choice
        exactly.

    Raises
    ------
      ValueError
        If ``value`` is not None and does not match any concrete choice.
    """
    target_idx = -1
    if value is None:
      for i in range(self._combo.count()):
        if self._combo.itemData(i, Qt.UserRole) is _NONE_SENTINEL:
          target_idx = i
          break
    else:
      text = str(value)
      for i in range(self._combo.count()):
        if (self._combo.itemData(i, Qt.UserRole) is not _NONE_SENTINEL
            and self._combo.itemText(i) == text):
          target_idx = i
          break
    if target_idx < 0:
      raise ValueError(
        "{v!r} is not one of {c!r}".format(v=value, c=self._current_choices()))
    if target_idx == self._combo.currentIndex():
      return
    blocker = QSignalBlocker(self._combo)
    self._combo.setCurrentIndex(target_idx)
    del blocker

  def setChoices(self, choices):
    """Replace the available options.

    Preserves the current selection if it is still present in the new
    list; otherwise selects the first concrete entry; otherwise selects
    ``"(none)"`` when allow_none is set.

    Parameters
    ----------
      choices: list of str
        The new option strings, in display order. Leading ``"(none)"``
        is added automatically when allow_none is set.
    """
    saved = self.value()
    self._populate(choices)
    if saved is not None and saved in choices:
      self.setValue(saved)
    elif choices:
      self.setValue(choices[0])
    elif self._allow_none:
      self.setValue(None)

  def _populate(self, choices):
    """Refill the combo box with the given choices, signal-blocked.

    Adds the leading ``"(none)"`` entry first when allow_none is set; that
    synthetic entry is tagged at ``Qt.UserRole`` with ``_NONE_SENTINEL`` so
    ``value()`` / ``setValue(None)`` can identify it without text comparison
    (a real choice whose text equals ``"(none)"`` is then unambiguous).

    Parameters
    ----------
      choices: list of str
    """
    blocker = QSignalBlocker(self._combo)
    self._combo.clear()
    if self._allow_none:
      self._combo.addItem(_NONE_LABEL)
      self._combo.setItemData(0, _NONE_SENTINEL, Qt.UserRole)
    for c in choices:
      self._combo.addItem(c)
    del blocker

  def _current_choices(self):
    """Return the current concrete choice texts (excluding the sentinel entry).

    Returns
    -------
      list of str
    """
    out = []
    for i in range(self._combo.count()):
      if self._combo.itemData(i, Qt.UserRole) is _NONE_SENTINEL:
        continue
      out.append(self._combo.itemText(i))
    return out

  def isValid(self):
    """Always True; combo selection is always one of the registered entries.

    Returns
    -------
      bool
    """
    return True

  def _on_index_changed(self, _i):
    """Slot for QComboBox.currentIndexChanged; emits valueChanged."""
    self.valueChanged.emit(self.value())


class ChoiceMultiWidget(PhilWidget):
  """PHIL choice(multi=True) editor backed by a QListWidget with checkboxes.

  ``value()`` returns a list of the selected choice strings (an empty list
  when nothing is checked, except when the choice's ``.optional == False``,
  in which case the widget reports invalid until at least one item is
  checked).
  """

  def __init__(self, definition, parent=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with .type.phil_type == 'choice' and
        .type.multi == True.
      parent: QWidget, optional
    """
    PhilWidget.__init__(self, definition, parent)
    self._list = QListWidget(self)
    layout = QVBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._list)
    choices = [w.value.lstrip("*") for w in self.definition.words]
    self._populate(choices)
    self._list.itemChanged.connect(self._on_item_changed)
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    if initial:
      self.setValue(initial)

  def value(self):
    """Return the list of currently checked choice strings."""
    out = []
    for i in range(self._list.count()):
      item = self._list.item(i)
      if item.checkState() == Qt.Checked:
        out.append(item.text())
    return out

  def setValue(self, value):
    """Check the items matching ``value`` (a list of choice strings)."""
    desired = set() if value is None else set(value)
    blocker = QSignalBlocker(self._list)
    for i in range(self._list.count()):
      item = self._list.item(i)
      want = (item.text() in desired)
      have = (item.checkState() == Qt.Checked)
      if want != have:
        item.setCheckState(Qt.Checked if want else Qt.Unchecked)
    del blocker

  def setChoices(self, choices):
    """Replace the available options; keep matching selections checked."""
    saved = self.value()
    self._populate(choices)
    self.setValue([s for s in saved if s in choices])

  def isValid(self):
    """True unless choice is non-optional and selection is empty."""
    optional = getattr(self.definition, "optional", True)
    if optional is False and not self.value():
      return False
    return True

  def errorString(self):
    """Return the validation error message, or '' when valid.

    Returns
    -------
      str
    """
    if self.isValid():
      return ""
    return "at least one selection required"

  def _populate(self, choices):
    """Refill the list widget with the given choices, signal-blocked.

    Each choice becomes a checkable ``QListWidgetItem`` with initial
    ``Qt.Unchecked`` state.

    Parameters
    ----------
      choices: list of str
    """
    blocker = QSignalBlocker(self._list)
    self._list.clear()
    for c in choices:
      item = QListWidgetItem(c)
      item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
      item.setCheckState(Qt.Unchecked)
      self._list.addItem(item)
    del blocker

  def _on_item_changed(self, _item):
    """Slot for list itemChanged; emits valueChanged."""
    self.valueChanged.emit(self.value())
