
from ..selection import SelectionTabController
from ..data import DataTabController
from ..sites import SitesTabController
from ..chimerax_controller import ChimeraXController
from ..molstar_controller import MolstarController
from ..restraints import RestraintsTopTabController
from ..restraints_table import RestraintsTableTopTabController
from ..qscore import QscoreTabController
from ..cif import CifTabController
from ..controller import Controller


from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QEvent


class ViewerGUIController(Controller):
  def __init__(self,parent=None,view=None,params=None,log=None):
    super().__init__(parent=parent,view=view)
    self.log = log
    if params and params.viewer_choice:
      self.viewer_choice = params.viewer_choice
    else:
      self.viewer_choice = 'molstar'



    # Main Level Components

    if params and params.show_tab:
      show_tab = params.show_tab
      if 'all' in params.show_tab:
        show_tab = 'all'
    else:
      show_tab = []


    if self.viewer_choice == 'molstar':
      self.molstar = MolstarController(parent=self,view=self.view.viewer_tab_view)
    else:
      self.chimerax = ChimeraXController(parent=self,view=self.view.viewer_tab_view)

    self.selection = SelectionTabController(parent=self,view=self.view.selection_tab_view)
    self.data = DataTabController(parent=self,view=self.view.data_tab_view)

    if 'all' in show_tab  or 'atoms' in show_tab:
      self.sites = SitesTabController(parent=self,view=self.view.sites_tab_view)
    if 'all' in show_tab  or 'cif' in show_tab:
      self.cif = CifTabController(parent=self,view=self.view.cif_tab_view)
    if 'all' in show_tab  or 'restraints' in show_tab:
      #self.restraints = RestraintsTopTabController(parent=self,view=self.view.restraints_tab_view)
      self.restraints_table = RestraintsTableTopTabController(parent=self,view=self.view.restraints_table_tab_view)
    if 'all' in show_tab  or 'qscore' in show_tab:
      self.qscore = QscoreTabController(parent=self,view=self.view.qscore_tab_view)


    # signals
    self.view.signal_close.connect(self.close_event)
    self.state.signals.tab_change.connect(self.change_tab_to)

    # Finish
    self.state._sync()

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
