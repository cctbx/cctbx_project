from PySide2.QtCore import QObject

from .controller import Controller


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
    return self.ref.active

  @active.setter
  def active(self,value):
    if not self.view.is_destroyed and not self.is_destroyed:
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

    # # remove from active
    # if (self.ref == self.state.active_model_ref):
    #   self.state.active_model_entry = None
    # if (self.ref == self.state.active_map_ref):
    #   self.state.active_map_entry = None
    # if (self.ref == self.state.active_selection_ref):
    #   self.state.active_selection_entry = None

    # # delete from data manager
    # if isinstance(self.ref,ModelRef):
    #   if self.ref.data.filepath in self.state.data_manager.get_model_names():
    #     name = self.ref.data.filepath
    #   elif self.ref.data.filename in self.state.data_manager.get_model_names():
    #     name = self.ref.data.filepath
    #   else:
    #     name = None
    #   if name:
    #     self.state.data_manager.remove_model(name)

    # if isinstance(self.ref,MapRef):
    #   if self.ref.data.filepath in self.state.data_manager.get_real_map_names():
    #     name = self.ref.data.filepath
    #   elif self.ref.data.filename in self.state.data_manager.get_real_map_names():
    #     name = self.ref.data.filepath
    #   if name:
    #     self.state.data_manager.remove_real_map(name)


    # # delete children

    # for ref_id,ref in list(self.state.references.items()):
    #   if not isinstance(ref,ModelRef):
    #     if hasattr(ref,"model_ref") and ref.model_ref == self.ref:
    #       if ref.entry is not None:
    #         if ref.entry in ref.entry.parent_list.entries:
    #           ref.entry.parent_list.remove_entry(ref.entry)
    # if self.ref.id in self.state.references:
    #   del self.state.references[self.ref.id]

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