"""FloatsWidget, FloatsTextWidget -- PHIL floats list editors."""

from PySide2.QtCore import QSignalBlocker
from PySide2.QtWidgets import QHBoxLayout, QVBoxLayout

from libtbx import Auto
from qttbx.widgets.phil import PhilWidget
from qttbx.widgets.phil.text_base import ValidatedLineEdit, ValidatedTextEdit


def _parse_float_token(tok, value_min, value_max, allow_none_elements):
  """Convert a single token to a float, ``None``, or raise.

  Parameters
  ----------
    tok: str
      Whitespace-stripped token to parse.
    value_min: float or None
      Inclusive lower bound for non-None values.
    value_max: float or None
      Inclusive upper bound for non-None values.
    allow_none_elements: bool
      If True, the literal token ``"None"`` or ``"*"`` parses to ``None``.

  Returns
  -------
    f: float or None

  Raises
  ------
    ValueError
      If ``tok`` does not parse as a number or violates the bounds.
  """
  if tok in ("None", "*") and allow_none_elements:
    return None
  try:
    f = float(tok)
  except ValueError:
    raise ValueError("{t!r} is not a number".format(t=tok))
  if value_min is not None and f < value_min:
    raise ValueError(
      "value {f} below minimum {m}".format(f=f, m=value_min))
  if value_max is not None and f > value_max:
    raise ValueError(
      "value {f} exceeds maximum {m}".format(f=f, m=value_max))
  return f


def _format_float_one(v, decimals):
  """Format a single float (or ``None``) for display.

  Parameters
  ----------
    v: float or None
    decimals: int or None
      If ``None``, format with ``"%g"``; otherwise format as fixed-point
      with the given number of decimal places.

  Returns
  -------
    text: str
  """
  if v is None:
    return "None"
  if decimals is None:
    return "%g" % v
  return ("{:." + str(int(decimals)) + "f}").format(v)


def _format_float_list(value, decimals, sep=" "):
  """Format a list of floats (with optional ``None`` elements) for display.

  Parameters
  ----------
    value: list, None, or libtbx.Auto
    decimals: int or None
      Per-element formatting precision (see ``_format_float_one``).
    sep: str, optional
      Token separator. ``" "`` for the single-line widget, ``"\\n"`` for
      the multi-line variant.

  Returns
  -------
    text: str
  """
  if value is None or value is Auto:
    return ""
  return sep.join(_format_float_one(v, decimals) for v in value)


class _FloatsBase(PhilWidget):
  """Shared parse/validate logic for ``FloatsWidget`` and ``FloatsTextWidget``.

  This base subclasses ``PhilWidget`` directly so that the two concrete
  classes (single-line and multi-line) can share ``_parse_tokens`` without
  duplicating the per-element bounds, allow-None-element, and size-bounds
  logic.
  """

  def __init__(self, definition, parent=None,
               size_min=None, size_max=None, allow_none_elements=False,
               decimals=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      size_min, size_max: int, optional
        Inclusive bounds on the number of items.
      allow_none_elements: bool, optional
        If True, the literal token ``"None"`` (or ``"*"``) in the input
        produces a ``None`` element in the output list.
      decimals: int, optional
        If given, format displayed values as fixed-point with that many
        decimal places. Default ``None`` uses ``"%g"``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._size_min = size_min
    self._size_max = size_max
    self._allow_none_elements = bool(allow_none_elements)
    self._decimals = decimals

  def _parse_tokens(self, tokens):
    """Convert an iterable of tokens into a validated list of floats / Nones.

    Empty input resolves through the standard empty-input ladder
    (Auto if ``_use_auto``, else ``None`` if ``_allow_none``, else raise).

    Parameters
    ----------
      tokens: list of str

    Returns
    -------
      out: list, None, or libtbx.Auto

    Raises
    ------
      ValueError
        On bad token, per-element bound violation, or size violation.
    """
    if not tokens:
      if self._use_auto:
        return Auto
      if self._allow_none:
        return None
      raise ValueError("at least one item required")
    vmin = getattr(self.definition.type, "value_min", None)
    vmax = getattr(self.definition.type, "value_max", None)
    out = []
    for tok in tokens:
      out.append(_parse_float_token(
        tok, vmin, vmax, self._allow_none_elements))
    n = len(out)
    if self._size_min is not None and n < self._size_min:
      raise ValueError(
        "{n} items below minimum {m}".format(n=n, m=self._size_min))
    if self._size_max is not None and n > self._size_max:
      raise ValueError(
        "{n} items exceed maximum {m}".format(n=n, m=self._size_max))
    return out


class FloatsWidget(_FloatsBase):
  """PHIL floats editor (whitespace-separated single line).

  The PHIL ``floats`` type is a ``list[float]``. Per-element bounds come
  from ``definition.type.value_min`` / ``value_max`` (the same attributes
  scalar ``float`` uses). For multi-line entry use ``FloatsTextWidget``.
  """

  def __init__(self, definition, parent=None,
               size_min=None, size_max=None, allow_none_elements=False,
               decimals=None, width=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      size_min, size_max: int, optional
        Inclusive bounds on the number of items.
      allow_none_elements: bool, optional
        If True, the literal token ``"None"`` (or ``"*"``) parses to a
        ``None`` element in the output list.
      decimals: int, optional
        If given, format displayed values as fixed-point with that many
        decimal places. Default ``None`` uses ``"%g"``.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
    """
    _FloatsBase.__init__(self, definition, parent,
                         size_min=size_min, size_max=size_max,
                         allow_none_elements=allow_none_elements,
                         decimals=decimals)
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
    """Return the current value (list of float / None / Auto)."""
    return self._parse(self._line_edit.text())

  def setValue(self, value):
    """Display ``value`` in the line edit. Idempotent and signal-blocked."""
    new_text = _format_float_list(value, self._decimals)
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
    return self._parse_tokens(text.split())

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


class FloatsTextWidget(_FloatsBase):
  """PHIL floats editor (one item per line).

  Long-form variant of ``FloatsWidget``. Each non-blank line contributes
  one token; per-element bounds and size bounds match the single-line
  variant.
  """

  def __init__(self, definition, parent=None,
               size_min=None, size_max=None, allow_none_elements=False,
               decimals=None, width=None, lines=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      size_min, size_max: int, optional
        Inclusive bounds on the number of items.
      allow_none_elements: bool, optional
        If True, the literal token ``"None"`` (or ``"*"``) on a line parses
        to a ``None`` element in the output list.
      decimals: int, optional
        If given, format displayed values as fixed-point with that many
        decimal places. Default ``None`` uses ``"%g"``.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
      lines: int, optional
        Visible-row hint; converted to a fixed height via ``fontMetrics``.
    """
    _FloatsBase.__init__(self, definition, parent,
                         size_min=size_min, size_max=size_max,
                         allow_none_elements=allow_none_elements,
                         decimals=decimals)
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
    """Return the current value (list of float / None / Auto)."""
    return self._parse(self._text_edit.toPlainText())

  def setValue(self, value):
    """Display ``value`` in the text edit. Idempotent and signal-blocked."""
    new_text = _format_float_list(value, self._decimals, sep="\n")
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
    tokens = [ln.strip() for ln in text.splitlines() if ln.strip() != ""]
    return self._parse_tokens(tokens)

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)
