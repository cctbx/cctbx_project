"""BoolWidget -- a PHIL bool editor wrapping QCheckBox."""

from PySide2.QtCore import Qt, QSignalBlocker
from PySide2.QtWidgets import QCheckBox, QHBoxLayout

from libtbx import Auto
from qttbx.widgets.phil import PhilWidget


class BoolWidget(PhilWidget):
  """PHIL bool editor backed by a QCheckBox.

  The checkbox is tristate when allow_none or use_auto is set on the
  widget; the partially-checked state then represents None or Auto
  respectively. Tristate enablement flows only through ``setAllowNone`` /
  ``setUseAuto`` -- ``setValue(None)`` does NOT enable tristate as a side
  effect (programmatic ``setCheckState(Qt.PartiallyChecked)`` on a
  non-tristate checkbox is a Qt-side no-op, leaving the widget in its
  prior state).
  """

  def __init__(self, definition, parent=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        A PHIL definition with .type.phil_type == 'bool'.
      parent: QWidget, optional
    """
    PhilWidget.__init__(self, definition, parent)
    self._check = QCheckBox(self)
    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._check)
    # Tristate flows from allow_none / use_auto, not from setValue.
    self._check.setTristate(self._allow_none or self._use_auto)
    self._check.stateChanged.connect(self._on_state_changed)
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)

  def value(self):
    """Return the current Python value.

    Returns
    -------
      value: bool or None or libtbx.Auto
        True / False from the checkbox's Checked / Unchecked state;
        Auto from PartiallyChecked when use_auto is set; otherwise None
        for PartiallyChecked.
    """
    s = self._check.checkState()
    if s == Qt.Checked:
      return True
    if s == Qt.Unchecked:
      return False
    if self._use_auto:
      return Auto
    return None

  def setValue(self, value):
    """Set the checkbox state from a Python value.

    Idempotent and signal-blocked. None and Auto map to PartiallyChecked
    (which is only visible when tristate has been enabled via
    ``setAllowNone`` / ``setUseAuto``).

    Parameters
    ----------
      value: bool or None or libtbx.Auto
    """
    new_state = self._state_for(value)
    if new_state == self._check.checkState():
      return
    blocker = QSignalBlocker(self._check)
    self._check.setCheckState(new_state)
    del blocker

  def _state_for(self, value):
    """Map a Python value to the corresponding Qt.CheckState.

    Parameters
    ----------
      value: bool or None or libtbx.Auto

    Returns
    -------
      state: Qt.CheckState
    """
    if value is None or value is Auto:
      return Qt.PartiallyChecked
    return Qt.Checked if value else Qt.Unchecked

  def setAllowNone(self, enable=True):
    """Toggle whether None is a permitted value.

    Also updates the underlying checkbox's tristate flag, since None is
    represented by PartiallyChecked.

    Parameters
    ----------
      enable: bool
    """
    PhilWidget.setAllowNone(self, enable)
    self._check.setTristate(enable or self._use_auto)

  def setUseAuto(self, enable=True):
    """Toggle whether Auto is a permitted value.

    Also updates the underlying checkbox's tristate flag, since Auto is
    represented by PartiallyChecked.

    Parameters
    ----------
      enable: bool
    """
    PhilWidget.setUseAuto(self, enable)
    self._check.setTristate(enable or self._allow_none)

  def isValid(self):
    """Return True iff the current state is a permitted value.

    Returns
    -------
      bool
        True for Checked / Unchecked unconditionally. PartiallyChecked is
        only valid when allow_none or use_auto is set.
    """
    if self._check.checkState() == Qt.PartiallyChecked:
      return self._allow_none or self._use_auto
    return True

  def _on_state_changed(self, _state):
    """Slot for QCheckBox.stateChanged; emits the widget's valueChanged."""
    self.valueChanged.emit(self.value())
