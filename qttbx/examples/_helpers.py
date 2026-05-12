"""Internal helpers shared by qttbx example apps."""

from qttbx.qt.QtCore import QEvent, QObject, QTimer
from qttbx.qt.QtWidgets import QHeaderView


class ProportionalColumns(QObject):
  """Apply fractional column widths to a ``QTableView`` and reapply on resize.

  Used by the qttbx example apps to give the DataManagerWidget's table
  a stable proportional layout (e.g. 60/20/20 across Filename, Type,
  and Used for). Columns not listed in ``proportions`` keep their
  current widths so a Delete column at a fixed pixel width still works.

  All listed columns are switched to ``QHeaderView.Interactive`` so the
  helper is free to rewrite their widths on every resize. The initial
  sizing and every subsequent reapply are deferred via
  ``QTimer.singleShot(0, ...)`` because Qt's ``Resize`` event fires
  before the viewport reflects the new width.
  """

  def __init__(self, table, proportions):
    """
    Parameters
    ----------
    table : QTableView
      The table whose columns should be sized proportionally.
    proportions : list of (int, float)
      ``(column_index, fraction)`` pairs. Fractions should sum to
      approximately 1.0; the last entry receives the rounding
      remainder so the sum exactly matches the available width.
    """
    QObject.__init__(self, table)
    self._table = table
    self._proportions = list(proportions)
    self._proportional_cols = set(c for c, _ in self._proportions)
    h = table.horizontalHeader()
    for col, _frac in self._proportions:
      h.setSectionResizeMode(col, QHeaderView.Interactive)
    table.installEventFilter(self)
    QTimer.singleShot(0, self.apply)

  def eventFilter(self, obj, event):
    if obj is self._table and event.type() == QEvent.Resize:
      QTimer.singleShot(0, self.apply)
    return False

  def apply(self):
    """Recompute and apply column widths so the proportions hold."""
    table = self._table
    model = table.model()
    if model is None:
      return
    n_cols = model.columnCount()
    fixed_w = sum(
      table.columnWidth(c)
      for c in range(n_cols)
      if c not in self._proportional_cols)
    available = table.viewport().width() - fixed_w
    if available <= 0:
      return
    assigned = 0
    for i, (col, frac) in enumerate(self._proportions):
      if i == len(self._proportions) - 1:
        w = available - assigned
      else:
        w = int(available * frac)
        assigned += w
      table.setColumnWidth(col, w)
