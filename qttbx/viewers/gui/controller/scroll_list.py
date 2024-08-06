
from .controller import Controller

class ScrollableListController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.entries = [] # list of ScrollEntryControllers

  def update(self):
    raise NotImplementedError

  @property
  def refs(self):
    return [entry.ref for entry in self.entries]

  @property
  def model_entries(self):
    return [entry for entry in self.entries if "ModelEntryController" in str(type(entry))]

  def add_entry(self, entry):
    assert entry.parent_list == self, 'This entry was not initialized with this parent list'
    self.entries.append(entry)
    self.view.container.layout.addWidget(entry.view)




  def remove_entry(self,entry):
    assert entry in self.entries, f"Cannot remove an entry {entry} not in the list: {self.entries}"
    self.entries.remove(entry)
    layout_index = self._find_widget_index(self.view.container.layout,entry.view)
    self._remove_widget_at(self.view.container.layout,layout_index)
    del entry
    self.update()

  def _find_widget_index(self,layout, target_widget):
    for i in range(layout.count()):
      widget = layout.itemAt(i).widget()
      if widget == target_widget:
        return i
    return None

  def _remove_widget_at(self,layout, index):
    item = layout.itemAt(index)
    if item is not None:
      widget = item.widget()
      if widget is not None:
        layout.removeWidget(widget)
        #widget.is_destroyed = True # Workaround for slots called after delete
        widget.deleteLater()
