"""UnitCellWidget -- PHIL unit_cell editor backed by cctbx.uctbx.

Accepts six whitespace-separated floats ``a b c alpha beta gamma`` and
validates via ``cctbx.uctbx.unit_cell((...))``. The widget ``value()`` is a
``cctbx.uctbx.unit_cell`` object (or None / libtbx.Auto).

Registers itself for ``phil_type == 'unit_cell'`` at module import time.
"""

from qttbx.qt.QtCore import QSignalBlocker
from qttbx.qt.QtWidgets import QHBoxLayout

from cctbx import uctbx
from libtbx import Auto
from qttbx.widgets.phil import PhilWidget, register_widget
from qttbx.widgets.phil.text_base import ValidatedLineEdit


class UnitCellWidget(PhilWidget):
  """PHIL unit_cell editor backed by ``cctbx.uctbx.unit_cell``.

  Accepts six whitespace-separated floats. Empty input resolves to ``Auto``
  when ``setUseAuto(True)``, ``None`` when ``setAllowNone(True)`` (default),
  or is rejected.
  """

  def __init__(self, definition, parent=None, width=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with ``.type.phil_type == 'unit_cell'``.
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
    self._line_edit.validityChanged.connect(self._on_validity_changed)
    self._line_edit.textChanged.connect(self._refresh_tooltip)
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)
    self._line_edit.validate()
    self._on_validity_changed(self._line_edit.isValid())

  def value(self):
    """Return the current value (unit_cell / None / Auto)."""
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

    Returns
    -------
      uc: cctbx.uctbx.unit_cell or None or libtbx.Auto

    Raises
    ------
      ValueError
        On empty input (when neither allow_none nor use_auto), bad token
        count, non-numeric tokens, or rejected geometry.
    """
    stripped = text.strip()
    if stripped == "":
      if self._use_auto:
        return Auto
      if self._allow_none:
        return None
      raise ValueError("value required")
    tokens = stripped.split()
    if len(tokens) != 6:
      raise ValueError(
        "expected 6 numbers (a b c alpha beta gamma); got {n}".format(
          n=len(tokens)))
    try:
      params = tuple(float(t) for t in tokens)
    except ValueError:
      raise ValueError("could not parse all six values as numbers")
    try:
      return uctbx.unit_cell(params)
    except RuntimeError as e:
      raise ValueError(str(e))

  def _format(self, value):
    """Format a unit_cell value as the displayed text."""
    if value is None or value is Auto:
      return ""
    params = value.parameters()
    return " ".join("%g" % v for v in params)

  def _on_validity_changed(self, _is_valid):
    """Slot for ``validityChanged``; refreshes the tooltip on transitions."""
    self._refresh_tooltip()

  def _refresh_tooltip(self, *_args):
    """Refresh the line-edit tooltip based on current validity and value.

    Wired to ``validityChanged`` (fires on transitions) AND ``textChanged``
    (fires on every keystroke and on programmatic setText) so the tooltip
    stays accurate even when the user types one valid input after another
    without crossing an invalid state.
    """
    if not self._line_edit.isValid():
      self._line_edit.setToolTip("")
      return
    try:
      uc = self.value()
    except Exception:
      self._line_edit.setToolTip("")
      return
    if uc is None or uc is Auto:
      self._line_edit.setToolTip("")
    else:
      params = uc.parameters()
      self._line_edit.setToolTip(
        "a = {0:g}\nb = {1:g}\nc = {2:g}\n"
        "alpha = {3:g}\nbeta = {4:g}\ngamma = {5:g}".format(*params))

  def _emit_value_changed(self, _v):
    """Re-derive the current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


register_widget("unit_cell", UnitCellWidget)
