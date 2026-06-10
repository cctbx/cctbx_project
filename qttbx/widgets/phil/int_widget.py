"""IntWidget -- a PHIL int editor wrapping ValidatedLineEdit."""

from qttbx.qt.QtCore import QSignalBlocker
from qttbx.qt.QtWidgets import QHBoxLayout

from libtbx import Auto
from qttbx.widgets.phil import PhilWidget
from qttbx.widgets.phil.text_base import ValidatedLineEdit


class IntWidget(PhilWidget):
  """PHIL int editor.

  Wraps a ValidatedLineEdit whose parse callable does ``int(text)`` and
  enforces ``definition.type.value_min`` and ``definition.type.value_max``.
  Empty input is interpreted as Auto when ``setUseAuto(True)``, as None
  when ``setAllowNone(True)`` (the default), and as a validation error
  otherwise.
  """

  def __init__(self, definition, parent=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with .type.phil_type == 'int'. Optional attributes
        ``value_min``, ``value_max`` are honored.
      parent: QWidget, optional
    """
    PhilWidget.__init__(self, definition, parent)
    self._line_edit = ValidatedLineEdit(parse=self._parse, parent=self)
    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._line_edit)
    # Forward signals up.
    self._line_edit.valueChanged.connect(self._emit_value_changed)
    self._line_edit.validityChanged.connect(self.validityChanged)
    # Initial value from the PHIL definition.
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)
    # Prime validation (empty string is invalid for non-optional, etc.).
    self._line_edit.validate()

  # -- PhilWidget overrides ---------------------------------------------------

  def value(self):
    """Return the current Python value.

    Returns
    -------
      value: int or None or libtbx.Auto
        The parsed integer; None when text is empty and allow_none is True;
        Auto when text is empty and use_auto is True.

    Raises
    ------
      ValueError
        If the text is empty and neither allow_none nor use_auto is set.
    """
    text = self._line_edit.text().strip()
    if text == "":
      if self._use_auto:
        return Auto
      if self._allow_none:
        return None
      raise ValueError("value required")
    return int(text)

  def setValue(self, value):
    """Display ``value`` in the line edit.

    Idempotent: returns early without firing signals if the formatted value
    matches the current text. Programmatic update is wrapped in a
    QSignalBlocker so ``valueChanged`` is not re-emitted.

    Parameters
    ----------
      value: int or None or libtbx.Auto
    """
    new_text = self._format(value)
    if new_text == self._line_edit.text():
      return                                       # idempotent
    blocker = QSignalBlocker(self._line_edit)
    self._line_edit.setText(new_text)
    del blocker                                    # release block
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

  # -- helpers ----------------------------------------------------------------

  def _parse(self, text):
    """Parse callable used by the inner ValidatedLineEdit.

    Parameters
    ----------
      text: str
        Raw text from the line edit.

    Returns
    -------
      int or None

    Raises
    ------
      ValueError
        If text cannot be parsed as int, falls outside [value_min, value_max],
        or is empty when neither allow_none nor use_auto is set.
    """
    text = text.strip()
    if text == "":
      if self._use_auto or self._allow_none:
        return None                  # validation passes; value() resolves sentinel
      raise ValueError("value required")
    n = int(text)                    # raises ValueError on bad text
    vmin = getattr(self.definition.type, 'value_min', None)
    vmax = getattr(self.definition.type, 'value_max', None)
    if vmin is not None and n < vmin:
      raise ValueError(
        "value {n} is below minimum {m}".format(n=n, m=vmin))
    if vmax is not None and n > vmax:
      raise ValueError(
        "value {n} exceeds maximum {m}".format(n=n, m=vmax))
    return n

  def _format(self, value):
    """Format a Python value as the displayed text.

    Parameters
    ----------
      value: int or None or libtbx.Auto

    Returns
    -------
      str
        Empty for None / Auto; ``str(value)`` otherwise.
    """
    if value is None or value is Auto:
      return ""
    return str(value)

  def _emit_value_changed(self, _v):
    """Re-derive the current value (resolving Auto/None) and emit.

    Used as the slot for the inner ValidatedLineEdit's valueChanged so that
    callers see a value that reflects the widget's Auto/None semantics
    rather than whatever the line-edit parse callable returned.
    """
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)
