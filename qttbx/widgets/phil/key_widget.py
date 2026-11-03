"""KeyWidget -- a PHIL key editor for identifier-shaped strings."""

from qttbx.qt.QtCore import QSignalBlocker
from qttbx.qt.QtWidgets import QHBoxLayout

import libtbx.phil
from qttbx.widgets.phil import PhilWidget
from qttbx.widgets.phil.text_base import ValidatedLineEdit


class KeyWidget(PhilWidget):
  """PHIL key editor.

  Wraps a ValidatedLineEdit whose parse callable validates that the text
  matches ``libtbx.phil.is_standard_identifier`` (begins with a letter or
  underscore, then letters / digits / underscores / dots). Empty input
  resolves to None when allow_none is True (the default).
  """

  def __init__(self, definition, parent=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with .type.phil_type == 'key'.
      parent: QWidget, optional
    """
    PhilWidget.__init__(self, definition, parent)
    self._line_edit = ValidatedLineEdit(parse=self._parse, parent=self)
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
    """Return the current Python value.

    Returns
    -------
      value: str or None
        The displayed identifier text, or None when empty and allow_none
        is set.

    Raises
    ------
      ValueError
        If the text is empty and allow_none is not set.
    """
    text = self._line_edit.text().strip()
    if text == "":
      if self._allow_none:
        return None
      raise ValueError("value required")
    return text

  def setValue(self, value):
    """Display ``value`` in the line edit.

    Idempotent and signal-blocked.

    Parameters
    ----------
      value: str or None
    """
    new_text = "" if value is None else str(value)
    if new_text == self._line_edit.text():
      return
    blocker = QSignalBlocker(self._line_edit)
    self._line_edit.setText(new_text)
    del blocker
    self._line_edit.validate()

  def isValid(self):
    """Forward to the inner ValidatedLineEdit.

    Returns
    -------
      bool
    """
    return self._line_edit.isValid()

  def errorString(self):
    """Forward to the inner ValidatedLineEdit.

    Returns
    -------
      str
    """
    return self._line_edit.errorString()

  def setReadOnly(self, ro=True):
    """Forward read-only state to the inner ValidatedLineEdit.

    Parameters
    ----------
      ro: bool
    """
    self._line_edit.setReadOnly(bool(ro))

  def _parse(self, text):
    """Parse callable used by the inner ValidatedLineEdit.

    Parameters
    ----------
      text: str

    Returns
    -------
      str or None

    Raises
    ------
      ValueError
        If text is not a standard identifier, or is empty when allow_none
        is not set.
    """
    text = text.strip()
    if text == "":
      if self._allow_none:
        return None
      raise ValueError("value required")
    if not libtbx.phil.is_standard_identifier(text):
      raise ValueError("{t!r} is not a valid identifier".format(t=text))
    return text

  def _emit_value_changed(self, _v):
    """Re-derive the current value (resolving None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)
