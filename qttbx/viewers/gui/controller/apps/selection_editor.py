"""
The top level Controller for for the molstar base app, which implements a Molstar viewer running 
in a QtWebView with very rudimentary GUI controls. 
"""
from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QEvent

from qttbx.viewers.gui.controller.selection import SelectionTabController
from qttbx.viewers.gui.controller.molstar import MolstarController
from qttbx.viewers.gui.controller import Controller


class SelectionEditorAppController(Controller):
  """
  This is the top level Controller instance for the base molstar app
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)


    # Main Level Components
    self.molstar = MolstarController(parent=self,view=self.view.viewer_tab_view)
    self.selection = SelectionTabController(parent=self,view=self.view.selection_tab_view)


    # signals
    self.view.signal_close.connect(self.close_event)
    self.state.signals.tab_change.connect(self.change_tab_to)



    # Finish
    self.state.start()
    if hasattr(self,"molstar"):
      self.molstar._update_state_from_remote()


  def update_from_remote(self, json_data):
    # Used by flask connector
    # this is connected up in main
    # Each message should have a clear connection to a function in the viewer
    # message is a json string
    assert isinstance(json_data,dict), (
        f"Expected json data ({json_data}) to be parsed to a dict, got: {type(json_data)}")
    #try:
    message_dict = json_data # alias
    command = message_dict.get("command")
    args = message_dict.get("args", [])
    kwargs = message_dict.get("kwargs", {})

    # Ensure command is allowed to prevent unauthorized access
    if command not in self.molstar.api_function_names:
      self.log(f"Command {command} is not allowed.")
      return

    # Look for the function on the controller
    func = getattr(self.molstar, command, None)
    if func is None:
      # Look for it on the viewer
      func = getattr(self.molstar.viewer, command, None)
      if func is None:
        self.log(f"No such command: {command}")
        return

    # Ensure args and kwargs are lists and dicts respectively
    if not isinstance(args, list) or not isinstance(kwargs, dict):
      self.log("Invalid args or kwargs format.")
      return
    func(**kwargs)



  def close_application(self):
    # manually call this function to close
    closeEvent = QEvent(QEvent.Close)
    QApplication.sendEvent(self.view, closeEvent)

  def close_event(self):
    # slot for built in close event
    if hasattr(self,"molstar"):
      self.molstar.close_viewer()
    if hasattr(self,"chimerax"):
      self.chimerax.close_viewer()

  def change_tab_to(self,name):
    """
    Change the active tab by a name
    """
    # TODO: Clean up tab change behavior for Selection and Viewer as child
    #if self.view._has_child_window:
    #self.log("Change to tab: ",name)
    tab_widget = self.view.tabs
    for index in range(tab_widget.count()):
      if tab_widget.tabText(index) == name:
        tab_widget.setCurrentIndex(index)
        return
