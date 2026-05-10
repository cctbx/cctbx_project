"""QStyledItemDelegate that hosts PhilWidget editors in a QTreeView."""

from PySide2.QtCore import QEvent, Qt
from PySide2.QtWidgets import QStyledItemDelegate

from qttbx.widgets.phil import widget_for_definition


class PhilItemDelegate(QStyledItemDelegate):
  """Tree-view controller dispatching the right PhilWidget per cell.

  Used with a QTreeView bound to a PhilModel. Each editable cell asks the
  delegate for an editor; the delegate looks up the right PhilWidget
  subclass via the registry, installs itself as an event filter on the
  editor to keep it open while the user has invalid input, and round-trips
  the value through the model on commit.
  """

  def createEditor(self, parent, option, index):
    """Build the editor widget for the given model index.

    Parameters
    ----------
      parent: QWidget
        Qt-supplied parent for the new editor.
      option: QStyleOptionViewItem
        Qt-supplied style option (unused).
      index: QModelIndex
        The model index whose cell is being edited.

    Returns
    -------
      editor: PhilWidget
        A fresh PhilWidget instance dispatched from the registry by the
        underlying PhilItem's definition type.
    """
    item = index.internalPointer()
    editor = widget_for_definition(item.definition, parent=parent)
    editor.installEventFilter(self)
    return editor

  def setEditorData(self, editor, index):
    """Push the model's current value into the editor.

    Parameters
    ----------
      editor: PhilWidget
        The editor returned by createEditor.
      index: QModelIndex
        The model index being edited.
    """
    editor.setValue(index.data(Qt.EditRole))

  def setModelData(self, editor, model, index):
    """Commit the editor's value back to the model.

    Returns early without writing if the editor is currently invalid; the
    model is left unchanged.

    Parameters
    ----------
      editor: PhilWidget
        The editor returned by createEditor.
      model: PhilModel
      index: QModelIndex
    """
    if not editor.isValid():
      return                                   # leave model unchanged
    model.setData(index, editor.value(), Qt.EditRole)

  def eventFilter(self, obj, event):
    """Suppress focus-out / Enter close while the editor is invalid.

    Esc still cancels normally (Qt convention) because the parent class's
    eventFilter handles QEvent.KeyPress with Qt.Key_Escape itself; this
    override only intercepts the close-on-commit events when invalid.

    Parameters
    ----------
      obj: QObject
        The watched object (the editor returned by createEditor).
      event: QEvent

    Returns
    -------
      bool
        True to consume (suppress) the event; otherwise the parent's
        eventFilter result.
    """
    if event.type() == QEvent.FocusOut:
      if hasattr(obj, "isValid") and not obj.isValid():
        return True                            # eat the event
    elif event.type() == QEvent.KeyPress:
      key = event.key()
      if key in (Qt.Key_Return, Qt.Key_Enter):
        if hasattr(obj, "isValid") and not obj.isValid():
          return True
    return QStyledItemDelegate.eventFilter(self, obj, event)
