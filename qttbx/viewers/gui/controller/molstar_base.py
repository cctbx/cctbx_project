"""
The controller for the molstar web app
"""
import time

from PySide2.QtCore import QTimer, QUrl
from PySide2.QtWidgets import QMessageBox

from .controller import Controller
from molstar_adaptbx.phenix.api import MolstarState
from molstar_adaptbx.phenix.molstar import MolstarGraphics
from molstar_adaptbx.phenix.utils import get_conda_env_directory
from molstar_adaptbx.phenix.server_utils import NodeHttpServer
from qttbx.viewers.gui.controller.viewer_controls_base import ViewerControlsBaseController
from qttbx.viewers.gui.model.refs import ModelRef


class MolstarBaseController(Controller):
  """
  A controller which manages external graphics for Molstar. However it does not actually
    produce any JavaScript code. That is delegated to its 'graphics' member variable, and 
    instance of the MolstarGraphics class. 

  So an example program flow would be
    1. This class responds to a signal or has a method directly called
    2. This class calls the corresponding method on self.graphics (MolstarGraphics instance)
    3. The MolstarGraphics instance runs Javascript code in the QtWebView to affect the Molstar webapp

  In this way there is a separation between the Qt GUI infrastructure (signals and slots) and 
    the Python-to-Javascript api present in MolstarGraphics. No other class besides this 
    one (MolstarController) should interact directly with the MolstarGraphics instance.

  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Start server and molstar graphics viewer
    env_name = "molstar_env"
    env_dir = get_conda_env_directory(env_name)


    env_bin_dir = f"{env_dir}/bin"
    molstar_install_dir = "/Users/user/software/debug/modules/molstar"
    server = NodeHttpServer([
      f"{env_bin_dir}/node",
      f"{molstar_install_dir}/src/phenix/server.js"
    ],port=8081,allow_port_change=True)

    self.graphics = MolstarGraphics(
      web_view = self.view.web_view,
      dm=self.state.data_manager,
      server = server,
    )
    self.graphics.state = self.state
    self.graphics_controls = ViewerControlsBaseController(parent=self,view=self.view.viewer_controls)


    self._blocking_commands = True
    self._picking_granularity = "atom"

    # Signals
    self.graphics.web_view.loadStarted.connect(self._on_load_started)
    self.graphics.web_view.loadFinished.connect(self._on_load_finished_pre_sync)
    self.state.signals.has_synced.connect(self._on_load_finished_post_sync)

    # Selections
    self.state.signals.picking_level.connect(self.set_picking_level)
    self.state.signals.select_all.connect(self.select_all)
    self.state.signals.deselect_all.connect(self.deselect_all)

    #timer for update
    self.sync_timer = QTimer()
    self.sync_timer.setInterval(2000)
    self.sync_timer.timeout.connect(self._update_state_from_remote)
    self.timer_accumulator = 0
    self.timer_max_retries = 10

    # Start by default
    self.start_viewer()

  # Functions used during the initial sync
  def _on_load_started(self):
    """
    Initial sync step 1: The web view window emits
      a signal 'loadStarted' which calls this function.
    """
    self._blocking_commands = True
    self.graphics._blocking_commands = True
    self.view.parent_explicit.setEnabled(False)

  def _on_load_finished_pre_sync(self, ok):
    """
    Initial sync step 2: The web view window emits
      a signal 'loadFinished' which calls this function.

      In turn, a timer is started which periodically polls to
        see if the Molstar web app has loaded (distinct from the web view)
    """
    if ok:
      self.log("Page loaded successfully. Waiting for sync")
      self.sync_timer.start()
    else:
      self.log("An error occurred while loading web app.")
      self._blocking_commands = True

  def _on_load_finished_post_sync(self, ok):
    """
    Initial sync step 3: The timer periodically calls
    '_update_state_from_remote', which in turn calls 'sync_remote'

    If the results are satisfactory, this function unblocks everything
      and loads data.
    """
    self.view.parent_explicit.setEnabled(True) # enable gui
    self._blocking_commands = False
    self.graphics._blocking_commands = False
    active_model_ref = self.state.get_active_ref(ref_class=ModelRef)
    if active_model_ref:
      self.load_model_from_ref(active_model_ref)

  def _update_state_from_remote(self):
    """
    The function that decides, at each call, if the viewer has synced.
    """
    molstar_state = self.graphics.sync_remote()
    if (molstar_state and 
    isinstance(molstar_state,MolstarState) and 
    molstar_state.connection_id == self.graphics.connection_id):
      self.sync_timer.stop()
      self.state.signals.has_synced.emit(True)
    else:
      self.timer_accumulator+=1
    if self.timer_accumulator>self.timer_max_retries:
      self.sync_timer.stop() # give up

  # API
  def start_viewer(self):
    self.graphics.start_viewer()

  def log_message(self,message):
    self.graphics.log_message(message)


  def load_model_from_ref(self,ref,label=None,format='pdb'):
    if not isinstance(ref,ModelRef):
      return None
    if self._blocking_commands:
      print("Skipping load_model_from_ref due to self._blocking_commands=True")
      return
    if ref is not None:
      # Check if already loaded
      if "molstar" in ref.identifiers and "molstar":
        if ref.identifiers["molstar"] is not None and ref.identifiers["molstar"] in self.state.external_loaded["molstar"]:
          self.log("Not loading model because already loaded into molstar")
          self.log(self.state.external_loaded)
          self.log(ref.identifiers["molstar"])
        return

      # Ref needs to be loaded
      self.graphics.load_model(
        filename=ref.filename, # mmtbx model
        #ref_id=ref.uuid,
      )

  def toggle_selection_mode(self,value):
    self.graphics._toggle_selection_mode(value)

  def select_all(self):
    self.graphics.select_all()

  def deselect_all(self):
    self.graphics.deselect_all()

  def close_viewer(self):
    self.graphics.close_viewer()

  def set_picking_level(self,picking_level):
    if 'atom' in picking_level:
      self.set_granularity("element")
    elif "residue" in picking_level:
      self.set_granularity("residue")
    else:
      pass

  def set_granularity(self,value="residue"):
    assert value in ['element','residue'], 'Provide one of the implemented picking levels'
    self._picking_granularity = value
    self.graphics._set_granularity(value=value)