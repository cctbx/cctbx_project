from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QEvent


from ..selection import SelectionTabController
from ..molstar_controller import MolstarController
from ..geometry.top_tab import GeometryTableTopTabController
from ..restraint_edits.top_tab import EditsTableTopTabController
from ..cif import CifTabController
from ..restraint import RestraintTabController
from ..controller import Controller


class ViewerGUIController(Controller):
  def __init__(self,parent=None,view=None,params=None,log=None):
    super().__init__(parent=parent,view=view)
    self.log = log
    if params and params.viewer_choice:
      self.viewer_choice = params.viewer_choice
    else:
      self.viewer_choice = 'molstar'



    # Main Level Components


    self.molstar = MolstarController(parent=self,view=self.view.viewer_tab_view)
    # else:
    #   self.chimerax = ChimeraXController(parent=self,view=self.view.viewer_tab_view)

    self.selection = SelectionTabController(parent=self,view=self.view.selection_tab_view)


    #self.data = DataTabController(parent=self,view=self.view.data_tab_view)

    # if 'all' in show_tab  or 'atoms' in show_tab:
    #   self.sites = SitesTabController(parent=self,view=self.view.sites_tab_view)
    # if 'all' in show_tab  or 'cif' in show_tab:
    self.cif = CifTabController(parent=self,view=self.view.cif_tab_view)

    # if 'all' in show_tab  or 'restraints' in show_tab:
    #   #self.restraints = GeometryTopTabController(parent=self,view=self.view.restraints_tab_view)
    #   self.restraints = RestraintTableTopTabController(parent=self,view=self.view.restraints)

    self.geometry = GeometryTableTopTabController(parent=self,view=self.view.geo)
    self.edits = EditsTableTopTabController(parent=self,view=self.view.edits)
    self.restraints = RestraintTabController(parent=self,view=self.view.restraints_tab_view)
    #self.qscore = QscoreTabController(parent=self,view=self.view.qscore_tab_view)


    # signals
    self.view.signal_close.connect(self.close_event)
    self.state.signals.tab_change.connect(self.change_tab_to)

    




    # Finish
    self.state.start()
    if hasattr(self,"molstar"):
      self.molstar._update_state_from_remote()


  def update_from_remote(self, json_data):
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
      print(f"Command {command} is not allowed.")
      return

    # Look for the function on the controller
    func = getattr(self.molstar, command, None)
    if func is None:
      # Look for it on the viewer
      func = getattr(self.molstar.viewer, command, None)
      if func is None:
        print(f"No such command: {command}")
        return

    # Ensure args and kwargs are lists and dicts respectively
    if not isinstance(args, list) or not isinstance(kwargs, dict):
      print("Invalid args or kwargs format.")
      return
    func(**kwargs)
    # except json.JSONDecodeError as e:
    #   print(f"Error decoding JSON: {e}")
    # except Exception as e:
    #   print(f"An error occurred: {e}")


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
    #print("Change to tab: ",name)
    tab_widget = self.view.tabs
    for index in range(tab_widget.count()):
      if tab_widget.tabText(index) == name:
        tab_widget.setCurrentIndex(index)
        return