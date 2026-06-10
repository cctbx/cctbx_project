"""FloatWidget -- a PHIL float editor wrapping ValidatedLineEdit."""

from qttbx.qt.QtCore import QSignalBlocker
from qttbx.qt.QtWidgets import QHBoxLayout

from libtbx import Auto
from qttbx.widgets.phil import PhilWidget
from qttbx.widgets.phil.text_base import ValidatedLineEdit


class FloatWidget(PhilWidget):
  """PHIL float editor.

  Wraps a ValidatedLineEdit whose parse callable does ``float(text)`` and
  enforces ``definition.type.value_min`` / ``value_max``. Auto/None
  handling mirrors IntWidget. Display formatting uses ``%g`` by default;
  pass ``decimals=N`` to use ``%.Nf`` instead.
  """

  def __init__(self, definition, parent=None, decimals=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with .type.phil_type == 'float'. Optional attributes
        ``value_min``, ``value_max`` are honored.
      parent: QWidget, optional
      decimals: int, optional
        If given, format the displayed value as ``"%.{decimals}f"``.
        Default ``None`` uses ``"%g"``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._decimals = decimals
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
      value: float or None or libtbx.Auto
        The parsed float; None when text is empty and allow_none is True;
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
    return float(text)

  def setValue(self, value):
    """Display ``value`` in the line edit.

    Idempotent and signal-blocked, mirroring IntWidget.setValue.

    Parameters
    ----------
      value: float or None or libtbx.Auto
    """
    new_text = self._format(value)
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
      float or None

    Raises
    ------
      ValueError
        If text cannot be parsed as float, falls outside [value_min, value_max],
        or is empty when neither allow_none nor use_auto is set.
    """
    text = text.strip()
    if text == "":
      if self._use_auto or self._allow_none:
        return None
      raise ValueError("value required")
    x = float(text)
    vmin = getattr(self.definition.type, 'value_min', None)
    vmax = getattr(self.definition.type, 'value_max', None)
    if vmin is not None and x < vmin:
      raise ValueError("value {x} is below minimum {m}".format(x=x, m=vmin))
    if vmax is not None and x > vmax:
      raise ValueError("value {x} exceeds maximum {m}".format(x=x, m=vmax))
    return x

  def _format(self, value):
    """Format a Python value as the displayed text.

    Parameters
    ----------
      value: float or None or libtbx.Auto

    Returns
    -------
      str
        Empty for None / Auto; ``"%.<decimals>f" % value`` if ``decimals``
        was set on the constructor; otherwise ``"%g" % value``.
    """
    if value is None or value is Auto:
      return ""
    if self._decimals is not None:
      return ("%." + str(self._decimals) + "f") % value
    return "%g" % value

  def _emit_value_changed(self, _v):
    """Re-derive the current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)
