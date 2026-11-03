"""AtomSelectionWidget, AtomSelectionTextWidget -- PHIL atom_selection editors.

The PHIL ``atom_selection`` type is a quoted string at the parser level
(`iotbx.phil.atom_selection_converters` extends `libtbx.phil.qstr_converters`).
v4 widgets do shape-only validation: non-empty and no literal ``$`` (matches
qstr semantics). Real syntactic / semantic validation against a model is
deferred to a future revision that integrates with the DataManager widget.

Both widgets register themselves for ``phil_type == 'atom_selection'`` at
module import time. ``AtomSelectionWidget`` is the registered default; the
multi-line ``AtomSelectionTextWidget`` is opt-in via
``PhilField(model, path, widget=AtomSelectionTextWidget)``.
"""

from qttbx.qt.QtCore import QSignalBlocker
from qttbx.qt.QtWidgets import QHBoxLayout, QVBoxLayout

from libtbx import Auto
from qttbx.widgets.phil import PhilWidget, register_widget
from qttbx.widgets.phil.text_base import ValidatedLineEdit, ValidatedTextEdit


def _parse_atom_selection(text, allow_none, use_auto):
  """Parse the displayed text into a Python value.

  Parameters
  ----------
    text: str
      The raw line-edit / text-edit content (NOT pre-stripped).
    allow_none: bool
    use_auto: bool

  Returns
  -------
    str or None or libtbx.Auto

  Raises
  ------
    ValueError
      If text is empty when neither allow_none nor use_auto is set, or
      if the text contains a literal '$' (matches the qstr_converters
      parent semantics).
  """
  if text == "":
    if use_auto:
      return Auto
    if allow_none:
      return None
    raise ValueError("value required")
  if "$" in text:
    raise ValueError("'$' is not allowed")
  return text


def _format_atom_selection(value):
  """Format an atom-selection value for display.

  Parameters
  ----------
    value: str or None or libtbx.Auto

  Returns
  -------
    str
      Empty for None / Auto; ``str(value)`` otherwise.
  """
  if value is None or value is Auto:
    return ""
  return str(value)


class AtomSelectionWidget(PhilWidget):
  """PHIL atom_selection editor (single line).

  Validates only the string shape (non-empty, no ``$``). Future revision will
  add syntactic validation against a model when the DataManager widget exists.
  """

  def __init__(self, definition, parent=None, width=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      width: int, optional
        Pixel hint passed to ``setMinimumWidth`` on the line edit.
    """
    PhilWidget.__init__(self, definition, parent)
    self._line_edit = ValidatedLineEdit(parse=self._parse, parent=self)
    if width is not None:
      self._line_edit.setMinimumWidth(int(width))
    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._line_edit)
    self._line_edit.valueChanged.connect(self._emit_value_changed)
    self._line_edit.validityChanged.connect(self.validityChanged)
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)
    self._line_edit.validate()

  def value(self):
    """Return the current value (str / None / Auto)."""
    return self._parse(self._line_edit.text())

  def setValue(self, value):
    """Display ``value`` in the line edit. Idempotent and signal-blocked."""
    new_text = _format_atom_selection(value)
    if new_text == self._line_edit.text():
      return
    blocker = QSignalBlocker(self._line_edit)
    self._line_edit.setText(new_text)
    del blocker
    self._line_edit.validate()

  def isValid(self):
    """Forward to the inner ValidatedLineEdit."""
    return self._line_edit.isValid()

  def errorString(self):
    """Forward to the inner ValidatedLineEdit."""
    return self._line_edit.errorString()

  def setReadOnly(self, ro=True):
    """Forward read-only state to the inner ValidatedLineEdit."""
    self._line_edit.setReadOnly(bool(ro))

  def _parse(self, text):
    """Parse callable used by the inner ValidatedLineEdit."""
    return _parse_atom_selection(text, self._allow_none, self._use_auto)

  def _emit_value_changed(self, _v):
    """Re-derive the current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


class AtomSelectionTextWidget(PhilWidget):
  """PHIL atom_selection editor (multi-line)."""

  def __init__(self, definition, parent=None, width=None, lines=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      width: int, optional
      lines: int, optional
        Visible-row hint; converted to a fixed height via fontMetrics.
    """
    PhilWidget.__init__(self, definition, parent)
    self._text_edit = ValidatedTextEdit(parse=self._parse, parent=self)
    if width is not None:
      self._text_edit.setMinimumWidth(int(width))
    if lines is not None:
      h = self._text_edit.fontMetrics().lineSpacing() * int(lines) + 8
      self._text_edit.setFixedHeight(h)
    layout = QVBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._text_edit)
    self._text_edit.valueChanged.connect(self._emit_value_changed)
    self._text_edit.validityChanged.connect(self.validityChanged)
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)
    self._text_edit.validate()

  def value(self):
    """Return the current value (str / None / Auto)."""
    return self._parse(self._text_edit.toPlainText())

  def setValue(self, value):
    """Display ``value`` in the text edit. Idempotent and signal-blocked."""
    new_text = _format_atom_selection(value)
    if new_text == self._text_edit.toPlainText():
      return
    blocker = QSignalBlocker(self._text_edit)
    self._text_edit.setPlainText(new_text)
    del blocker
    self._text_edit.validate()

  def isValid(self):
    """Forward to the inner ValidatedTextEdit."""
    return self._text_edit.isValid()

  def errorString(self):
    """Forward to the inner ValidatedTextEdit."""
    return self._text_edit.errorString()

  def setReadOnly(self, ro=True):
    """Forward read-only state to the inner ValidatedTextEdit."""
    self._text_edit.setReadOnly(bool(ro))

  def _parse(self, text):
    """Parse callable used by the inner ValidatedTextEdit."""
    return _parse_atom_selection(text, self._allow_none, self._use_auto)

  def _emit_value_changed(self, _v):
    """Re-derive the current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


register_widget("atom_selection", AtomSelectionWidget)
