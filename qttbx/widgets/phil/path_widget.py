"""PathWidget, PathsTextWidget -- PHIL path editors."""

import os

from qttbx.qt.QtCore import QSignalBlocker
from qttbx.qt.QtWidgets import (
  QFileDialog, QHBoxLayout, QPushButton, QVBoxLayout,
)

from libtbx import Auto
from qttbx.widgets.phil import PhilWidget
from qttbx.widgets.phil.text_base import ValidatedLineEdit, ValidatedTextEdit


def _make_path_parse(allow_none, use_auto, must_exist):
  """Return a (text -> str/None/Auto) parse callable for a single path.

  Parameters
  ----------
    allow_none: bool
    use_auto: bool
    must_exist: bool
      When True, the parse fails if ``os.path.exists(text)`` is False.

  Returns
  -------
    callable
      A closure suitable for ``ValidatedLineEdit``'s ``parse=`` argument.
  """
  def _parse(text):
    """Parse a single path string; raises ValueError on bad input."""
    text = text.strip()
    if text == "":
      if use_auto:
        return Auto
      if allow_none:
        return None
      raise ValueError("value required")
    if must_exist and not os.path.exists(text):
      raise ValueError("path does not exist: {p}".format(p=text))
    return text
  return _parse


def _make_paths_parse(allow_none, use_auto, must_exist, size_min, size_max):
  """Return a (text -> list[str]/None/Auto) parse callable for newline-separated paths.

  Parameters
  ----------
    allow_none: bool
    use_auto: bool
    must_exist: bool
      When True, the parse fails if any path doesn't exist on disk.
    size_min, size_max: int, optional
      Inclusive bounds on the number of paths.

  Returns
  -------
    callable
      A closure suitable for ``ValidatedTextEdit``'s ``parse=`` argument.
  """
  def _parse(text):
    """Parse newline-separated paths; raises ValueError on bad input."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
      if use_auto:
        return Auto
      if allow_none:
        return None
      raise ValueError("at least one path required")
    if size_min is not None and len(lines) < size_min:
      raise ValueError(
        "{n} paths below minimum {m}".format(n=len(lines), m=size_min))
    if size_max is not None and len(lines) > size_max:
      raise ValueError(
        "{n} paths exceed maximum {m}".format(n=len(lines), m=size_max))
    if must_exist:
      for p in lines:
        if not os.path.exists(p):
          raise ValueError("path does not exist: {p}".format(p=p))
    return lines
  return _parse


class PathWidget(PhilWidget):
  """PHIL path editor: line edit + browse button.

  Optional ``must_exist=True`` validates the path with ``os.path.exists``.
  """

  def __init__(self, definition, parent=None,
               must_exist=False, file_filter="All files (*)", width=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      must_exist: bool
        If True, validation fails when the path does not exist.
      file_filter: str
        Filter string passed to ``QFileDialog.getOpenFileName``.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth`` on the line edit.
    """
    PhilWidget.__init__(self, definition, parent)
    self._must_exist = bool(must_exist)
    self._file_filter = file_filter
    self._line_edit = ValidatedLineEdit(parse=self._parse, parent=self)
    if width is not None:
      self._line_edit.setMinimumWidth(int(width))
    self._browse_btn = QPushButton("...", parent=self)
    self._browse_btn.setMaximumWidth(32)
    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._line_edit, stretch=1)
    layout.addWidget(self._browse_btn)
    self._line_edit.valueChanged.connect(self._emit_value_changed)
    self._line_edit.validityChanged.connect(self.validityChanged)
    self._browse_btn.clicked.connect(self._on_browse)
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
    if value is None or value is Auto:
      new_text = ""
    else:
      new_text = str(value)
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
    """Forward read-only state to the inner ValidatedLineEdit and disable browse."""
    self._line_edit.setReadOnly(bool(ro))
    self._browse_btn.setEnabled(not bool(ro))

  def _parse(self, text):
    """Parse callable used by the inner ValidatedLineEdit."""
    fn = _make_path_parse(self._allow_none, self._use_auto, self._must_exist)
    return fn(text)

  def _on_browse(self):
    """Slot for the [...] button; opens QFileDialog and updates the line edit."""
    start = self._line_edit.text() or ""
    path, _ = QFileDialog.getOpenFileName(self, "Select file", start,
                                          self._file_filter)
    if path:
      self._apply_browse_result(path)

  def _apply_browse_result(self, path):
    """Set the text from a browse result. Factored for testability."""
    blocker = QSignalBlocker(self._line_edit)
    self._line_edit.setText(path)
    del blocker
    self._line_edit.validate()
    if self._line_edit.isValid():
      self.valueChanged.emit(self.value())

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


class PathsTextWidget(PhilWidget):
  """PHIL multi-path editor: one path per line, no browse button."""

  def __init__(self, definition, parent=None,
               must_exist=False, size_min=None, size_max=None,
               width=None, lines=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      must_exist: bool
        If True, validation fails when any path does not exist.
      size_min, size_max: int, optional
        Inclusive bounds on the number of non-empty lines.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
      lines: int, optional
        Visible-row hint; converted to a fixed height via ``fontMetrics``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._must_exist = bool(must_exist)
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
    """Return the current value (list of paths / None / Auto)."""
    return self._parse(self._text_edit.toPlainText())

  def setValue(self, value):
    """Display ``value`` in the text edit. Idempotent and signal-blocked.

    Accepts ``None``, ``Auto``, a single string, or a list-like of strings.
    """
    if value is None or value is Auto:
      new_text = ""
    elif hasattr(value, "__iter__") and not hasattr(value, "lower"):
      new_text = "\n".join(str(p) for p in value)
    else:
      new_text = str(value)
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
    fn = _make_paths_parse(self._allow_none, self._use_auto,
                           self._must_exist, self._size_min, self._size_max)
    return fn(text)

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)
