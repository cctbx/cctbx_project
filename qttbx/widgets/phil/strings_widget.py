"""StringsWidget, StringsTextWidget -- PHIL strings list editors."""

from qttbx.qt.QtCore import QSignalBlocker
from qttbx.qt.QtWidgets import QHBoxLayout, QVBoxLayout

from libtbx import Auto
from qttbx.widgets.phil import PhilWidget
from qttbx.widgets.phil.text_base import ValidatedLineEdit, ValidatedTextEdit


def _check_size(items, size_min, size_max):
  """Raise ValueError if ``items`` has fewer than ``size_min`` or more than
  ``size_max`` entries.

  Parameters
  ----------
    items: sequence
    size_min: int or None
    size_max: int or None
  """
  n = len(items)
  if size_min is not None and n < size_min:
    raise ValueError("{n} items below minimum {m}".format(n=n, m=size_min))
  if size_max is not None and n > size_max:
    raise ValueError("{n} items exceed maximum {m}".format(n=n, m=size_max))


class StringsWidget(PhilWidget):
  """PHIL strings editor (whitespace-separated single line).

  The PHIL ``strings`` type is a ``list[str]``. This short-form variant
  accepts whitespace-separated tokens via ``str.split()``; for free-form
  strings with internal spaces use ``StringsTextWidget`` (one item per
  line).
  """

  def __init__(self, definition, parent=None,
               size_min=None, size_max=None, width=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      size_min, size_max: int, optional
        Inclusive bounds on the number of items.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._size_min = size_min
    self._size_max = size_max
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
    """Return the current value (list of strings / None / Auto)."""
    return self._parse(self._line_edit.text())

  def setValue(self, value):
    """Display ``value`` in the line edit. Idempotent and signal-blocked."""
    if value is None or value is Auto:
      new_text = ""
    else:
      new_text = " ".join(str(s) for s in value)
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
    tokens = text.split()
    if not tokens:
      if self._use_auto:
        return Auto
      if self._allow_none:
        return None
      raise ValueError("at least one item required")
    _check_size(tokens, self._size_min, self._size_max)
    return tokens

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


class StringsTextWidget(PhilWidget):
  """PHIL strings editor (one item per line).

  Long-form variant of ``StringsWidget``: each non-blank line of the text
  edit is one list item. Blank / whitespace-only lines are dropped.
  Suitable for free-form strings whose items contain internal whitespace.
  """

  def __init__(self, definition, parent=None,
               size_min=None, size_max=None, width=None, lines=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      size_min, size_max: int, optional
        Inclusive bounds on the number of items.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
      lines: int, optional
        Visible-row hint; converted to a fixed height via ``fontMetrics``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._size_min = size_min
    self._size_max = size_max
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
    """Return the current value (list of strings / None / Auto)."""
    return self._parse(self._text_edit.toPlainText())

  def setValue(self, value):
    """Display ``value`` in the text edit. Idempotent and signal-blocked."""
    if value is None or value is Auto:
      new_text = ""
    else:
      new_text = "\n".join(str(s) for s in value)
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
    items = [ln for ln in text.splitlines() if ln.strip() != ""]
    if not items:
      if self._use_auto:
        return Auto
      if self._allow_none:
        return None
      raise ValueError("at least one item required")
    _check_size(items, self._size_min, self._size_max)
    return items

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)
