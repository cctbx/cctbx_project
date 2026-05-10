"""QstrWidget, QstrTextWidget -- PHIL qstr editors (short and long variants)."""

from PySide2.QtCore import QSignalBlocker
from PySide2.QtWidgets import QHBoxLayout, QVBoxLayout

from libtbx import Auto
from qttbx.widgets.phil import PhilWidget
from qttbx.widgets.phil.text_base import ValidatedLineEdit, ValidatedTextEdit


def _make_qstr_parse(allow_none, use_auto, min_length, max_length):
  """Return a (text -> str/None/Auto) parse callable closing over the limits.

  Rejects literal ``$`` to match the wxtbx convention (the character has
  shell-substitution semantics elsewhere in PHENIX).
  """
  def _parse(text):
    """Parse the displayed text into a Python str/None/Auto."""
    if text == "":
      if use_auto:
        return Auto
      if allow_none:
        return None
      raise ValueError("value required")
    if "$" in text:
      raise ValueError("'$' is not allowed")
    if min_length is not None and len(text) < min_length:
      raise ValueError(
        "length {n} below minimum {m}".format(n=len(text), m=min_length))
    if max_length is not None and len(text) > max_length:
      raise ValueError(
        "length {n} exceeds maximum {m}".format(n=len(text), m=max_length))
    return text
  return _parse


def _format_qstr(value):
  """Format a str/None/Auto value for display in a text input."""
  if value is None or value is Auto:
    return ""
  return str(value)


class QstrWidget(PhilWidget):
  """PHIL qstr (quoted string) editor (single line).

  PHIL's ``qstr`` type is the same as ``str`` except the value is round-tripped
  through ``libtbx.phil.tokenizer.word(quote_token='"')`` on serialization.
  In-memory the value is a plain Python ``str``. Literal ``$`` is rejected
  to match the wxtbx convention.
  """

  def __init__(self, definition, parent=None,
               min_length=None, max_length=None, width=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with .type.phil_type == 'qstr'.
      parent: QWidget, optional
      min_length, max_length: int, optional
        Inclusive bounds on string length when non-empty. PHIL has no
        equivalent attributes; these are form-author choices.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._min_length = min_length
    self._max_length = max_length
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
    new_text = _format_qstr(value)
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
    fn = _make_qstr_parse(self._allow_none, self._use_auto,
                          self._min_length, self._max_length)
    return fn(text)

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


class QstrTextWidget(PhilWidget):
  """PHIL qstr editor (multi-line).

  Long-text variant of ``QstrWidget`` wrapping a ``ValidatedTextEdit``.
  Newlines are preserved verbatim in the value. Literal ``$`` is rejected.
  """

  def __init__(self, definition, parent=None,
               min_length=None, max_length=None, width=None, lines=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      min_length, max_length: int, optional
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
      lines: int, optional
        Visible-row hint; converted to a fixed height via ``fontMetrics``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._min_length = min_length
    self._max_length = max_length
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
    new_text = _format_qstr(value)
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
    fn = _make_qstr_parse(self._allow_none, self._use_auto,
                          self._min_length, self._max_length)
    return fn(text)

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)
