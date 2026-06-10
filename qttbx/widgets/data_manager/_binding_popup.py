"""DataManagerBindingPopup -- a frameless checkbox-list popup for
choosing which PHIL parameters a file is bound to."""

from qttbx.qt.QtCore import Qt, Signal
from qttbx.qt.QtWidgets import (
  QFrame, QVBoxLayout, QScrollArea, QCheckBox, QWidget)


class DataManagerBindingPopup(QFrame):
  """Popup-style frame with a vertical list of checkboxes, one per
  compatible PHIL parameter for the file the popup was opened for.

  Use ``populate(candidates)`` then position the popup with ``move()``
  and call ``show()`` -- the consumer chooses the anchor (the
  DataManagerWidget anchors at the ``+`` chip's bottom-left in
  viewport-global coordinates). Each candidate is a tuple
  ``(phil_path, label, checked, disabled, tooltip)``. The popup emits
  ``bindToggled(phil_path, want_bound)`` immediately on each user
  click; no OK button.
  """

  bindToggled = Signal(str, bool)

  def __init__(self, parent=None):
    QFrame.__init__(self, parent, Qt.Popup)
    self.setFrameShape(QFrame.StyledPanel)
    layout = QVBoxLayout(self)
    layout.setContentsMargins(4, 4, 4, 4)
    self._scroll = QScrollArea(self)
    self._scroll.setWidgetResizable(True)
    self._content = QWidget(self._scroll)
    self._content_layout = QVBoxLayout(self._content)
    self._content_layout.setContentsMargins(0, 0, 0, 0)
    self._scroll.setWidget(self._content)
    layout.addWidget(self._scroll)
    self._checkboxes = []   # list[QCheckBox]
    self._paths = []        # list[str]

  def populate(self, candidates):
    """Replace the checkbox list with the given candidates.

    Parameters
    ----------
    candidates : list of tuple
      ``(phil_path, label, checked, disabled, tooltip)``.
    """
    # Tear down old children.
    while self._content_layout.count():
      item = self._content_layout.takeAt(0)
      w = item.widget()
      if w is not None:
        w.blockSignals(True)       # prevent stale signal emissions during deleteLater window
        w.deleteLater()
    self._checkboxes = []
    self._paths = []
    for phil_path, label, checked, disabled, tooltip in candidates:
      cb = QCheckBox(label, parent=self._content)
      cb.setChecked(bool(checked))
      cb.setEnabled(not disabled)
      if tooltip:
        cb.setToolTip(tooltip)
      idx = len(self._paths)
      cb.toggled.connect(
        lambda state, i=idx: self._emit_toggle(i, state))
      self._content_layout.addWidget(cb)
      self._checkboxes.append(cb)
      self._paths.append(phil_path)
    self._content_layout.addStretch(1)

  def candidate_count(self):
    """Return the number of currently populated candidates."""
    return len(self._paths)

  def _toggle_at(self, idx):
    """Test hook: toggle the checkbox at ``idx`` programmatically.

    No-op if the checkbox is disabled (matches user behavior --
    disabled checkboxes ignore clicks).
    """
    cb = self._checkboxes[idx]
    if cb.isEnabled():
      cb.toggle()

  def _emit_toggle(self, idx, state):
    self.bindToggled.emit(self._paths[idx], bool(state))
