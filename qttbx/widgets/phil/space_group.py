"""SpaceGroupWidget -- PHIL space_group editor backed by cctbx.sgtbx.

Accepts space group numbers (e.g. ``"19"``), Hermann-Mauguin symbols
(``"P 21 21 21"``, ``"P21"``), and Hall symbols (``"P 2yb"``). The widget
``value()`` is a ``cctbx.sgtbx.space_group_info`` object (or None / libtbx.Auto).
A small label to the right of the line edit shows a canonical form when the
input is valid -- by default the full Hermann-Mauguin symbol, or the Hall
symbol when ``side_label=SIDE_LABEL_HALL`` is passed at construction.

Registers itself for ``phil_type == 'space_group'`` at module import time.
"""

from PySide2.QtCore import QSignalBlocker
from PySide2.QtWidgets import QHBoxLayout, QLabel

from cctbx import sgtbx
from libtbx import Auto
from qttbx.widgets.phil import PhilWidget, register_widget
from qttbx.widgets.phil.text_base import ValidatedLineEdit


SIDE_LABEL_FULL_HERMANN_MAUGUIN = "full_hermann_mauguin"
SIDE_LABEL_HALL = "hall"
_VALID_SIDE_LABELS = (SIDE_LABEL_FULL_HERMANN_MAUGUIN, SIDE_LABEL_HALL)


class SpaceGroupWidget(PhilWidget):
  """PHIL space_group editor backed by ``cctbx.sgtbx.space_group_info``.

  Accepts three input forms:
    - space group number (e.g. ``"19"``)
    - Hermann-Mauguin (full or short, e.g. ``"P 21 21 21"`` / ``"P21"``)
    - Hall symbol (e.g. ``"P 2yb"``) via fallback parsing

  Empty input resolves to ``Auto`` when ``setUseAuto(True)``, ``None`` when
  ``setAllowNone(True)`` (the default), or is rejected.
  """

  def __init__(self, definition, parent=None, width=None,
               side_label=SIDE_LABEL_FULL_HERMANN_MAUGUIN):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with ``.type.phil_type == 'space_group'``.
      parent: QWidget, optional
      width: int, optional
        Pixel hint passed to ``setMinimumWidth`` on the line edit.
      side_label: str, optional
        Selects what the side label shows when input is valid. One of
        ``SIDE_LABEL_FULL_HERMANN_MAUGUIN`` (the canonical full HM symbol;
        the default) or ``SIDE_LABEL_HALL`` (the canonical Hall symbol).
    """
    PhilWidget.__init__(self, definition, parent)
    if side_label not in _VALID_SIDE_LABELS:
      raise ValueError(
        "side_label must be one of {opts!r}; got {got!r}".format(
          opts=_VALID_SIDE_LABELS, got=side_label))
    self._side_label_kind = side_label
    self._line_edit = ValidatedLineEdit(parse=self._parse, parent=self)
    if width is not None:
      self._line_edit.setMinimumWidth(int(width))
    self._side_label = QLabel("", parent=self)
    self._side_label.setStyleSheet(
      "QLabel { color: #555; font-style: italic; }")
    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._line_edit, stretch=1)
    layout.addWidget(self._side_label)
    self._line_edit.valueChanged.connect(self._emit_value_changed)
    self._line_edit.validityChanged.connect(self.validityChanged)
    self._line_edit.validityChanged.connect(self._on_validity_changed)
    # Refresh the side label whenever the text changes (covers both user
    # edits via textEdited and programmatic setText). The line edit's own
    # validityChanged only fires on transitions, which is too coarse here.
    self._line_edit.textChanged.connect(self._refresh_side_label)
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)
    self._line_edit.validate()
    self._on_validity_changed(self._line_edit.isValid())

  def value(self):
    """Return the current value (space_group_info / None / Auto)."""
    return self._parse(self._line_edit.text())

  def setValue(self, value):
    """Display ``value`` in the line edit. Idempotent and signal-blocked."""
    new_text = self._format(value)
    if new_text == self._line_edit.text():
      return
    blocker = QSignalBlocker(self._line_edit)
    self._line_edit.setText(new_text)
    del blocker
    self._line_edit.validate()
    self._on_validity_changed(self._line_edit.isValid())

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
    """Parse callable used by the inner ValidatedLineEdit.

    Tries ``sgtbx.space_group_info(symbol=text)`` first (handles SG number
    and Hermann-Mauguin forms). On RuntimeError, falls back to the Hall
    path: ``sgtbx.space_group_info(group=sgtbx.space_group(hall_symbol=text))``.
    When both parses fail, the raised ValueError combines both diagnostic
    messages so a user typing Hall sees the Hall-form error too.

    Returns
    -------
      sg: cctbx.sgtbx.space_group_info or None or libtbx.Auto

    Raises
    ------
      ValueError
        On empty input (when neither allow_none nor use_auto is set), or
        when both the symbol-based and Hall-based parses raise.
    """
    text = text.strip()
    if text == "":
      if self._use_auto:
        return Auto
      if self._allow_none:
        return None
      raise ValueError("value required")
    # First try the standard symbol form (number or Hermann-Mauguin).
    try:
      return sgtbx.space_group_info(symbol=text)
    except RuntimeError as symbol_err:
      symbol_msg = str(symbol_err)
    # Fall back to Hall input.
    try:
      group = sgtbx.space_group(hall_symbol=text)
    except RuntimeError as hall_err:
      raise ValueError(
        "not a space group symbol ({s}) and not a Hall symbol ({h})".format(
          s=symbol_msg, h=str(hall_err)))
    return sgtbx.space_group_info(group=group)

  def _format(self, value):
    """Format a space_group_info value for display in the line edit."""
    if value is None or value is Auto:
      return ""
    return str(value)

  def _on_validity_changed(self, is_valid):
    """Update the side label on a validity transition."""
    self._refresh_side_label()

  def _refresh_side_label(self, *_args):
    """Update the side label per the configured side_label kind.

    Reads the current value (resolving Auto/None and parsing failures) and
    sets the side label to the canonical full Hermann-Mauguin (default) or
    Hall symbol when valid; clears it otherwise.
    """
    if not self._line_edit.isValid():
      self._side_label.setText("")
      return
    try:
      sg = self.value()
    except Exception:
      self._side_label.setText("")
      return
    if sg is None or sg is Auto:
      self._side_label.setText("")
      return
    if self._side_label_kind == SIDE_LABEL_HALL:
      self._side_label.setText(sg.type().hall_symbol())
    else:
      # Default: full Hermann-Mauguin.
      self._side_label.setText(sg.type().lookup_symbol())

  def _emit_value_changed(self, _v):
    """Re-derive the current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


register_widget("space_group", SpaceGroupWidget)
