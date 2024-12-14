
from .controller import Controller
from PySide2.QtCore import QObject

class ScrollEntryController(Controller, QObject):
  """
  Controls views/scroll_entry.py
  """
  def __init__(self,parent=None,view=None,ref=None,show=True):
    QObject.__init__(self) # QObject
    Controller.__init__(self,parent=parent,view=view)
    assert ref is not None, f"An Entry is a gui analog to a Ref. Must provide ref"
    self.ref = ref
    self.ref.entry = self
    self.show = True
    self._is_destroyed = False

    # Connect signals from view

    # Active
    if self.view.active_toggle is not None:
      self.view.active_toggle.stateChanged.connect(self._toggle_active_func)  
      #self.state.signals.new_active_ref.connect(self._check_active)


    # Close
    self.view.button_close.clicked.connect(self.remove)

    # Display name
    self.view.label_name.setText(self.ref.label)


  @property
  def active(self):
    return self.view.active_toggle.is_checked

  @active.setter
  def active(self,value):
    #if not self.view.is_destroyed and not self.is_destroyed:
    assert isinstance(value,bool), "Active must be boolean"
    self.view.active_toggle.is_checked=value

  @property
  def is_destroyed(self):
    return self._is_destroyed

  @is_destroyed.setter
  def is_destroyed(self,value):
    self._is_destroyed = value
    self.view.is_destroyed = value

  @property
  def parent_list(self):
    return self.parent


  def remove(self):

    self.ref.show = False

    # harsh reset of viewer, the problem is that the viewer will send back old pairings if same model
    self.is_destroyed = True
    self.parent_list.remove_entry(self) # remove from gui
    #self.state.signals.remove_ref.emit(self.ref)
    #self.state.signals.clear.emit("Resetting....")



  def _toggle_active_func(self,is_checked):
    # 1. set active (change visual state)
    # 2. call toggle_active_func (implemented by subclasses)
    if is_checked:
      self.active = True
      self.toggle_active_func(True)
    else:
      self.active = False
      self.toggle_active_func(False)

  def toggle_active_func(self,is_checked):
    # implement for subclasses. Called when toggle is switched
    raise NotImplementedError


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
