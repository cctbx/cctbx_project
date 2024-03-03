import json
import sys
import threading
from pathlib import Path
#import argparse

from PySide2.QtGui import QIcon
from PySide2.QtWidgets import QApplication, QWidget, QApplication, QWidget, QPushButton, QVBoxLayout
from PySide2.QtWidgets import QApplication, QDialog, QVBoxLayout, QPushButton, QMessageBox, QMainWindow

from PySide2.QtCore import QObject, QEvent, Qt,  QEvent, QSize, Signal
from PySide2.QtSvg import QSvgRenderer

from iotbx.data_manager import DataManager
from ..view.apps.main import ViewerGUIView
from ..controller.apps.main import ViewerGUIController
from ..state.state import State
from ...last.selection_utils import Selection, SelectionQuery
from . import ViewerChoiceDialog, check_program_access
from ...last.python_utils import DotDict
from flask import Flask, request, jsonify

QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

class FlaskSignalEmitter(QObject):
    flask_signal = Signal(dict)

# Create an instance of the signal emitter
flask_signal_emitter = FlaskSignalEmitter()



flask_app = Flask(__name__)

@flask_app.route('/process_json', methods=['POST'])
def process_json():
    json_data = request.get_json()
    print("type:",type(json_data))
    print(json_data)
    # Continue with actual json message
    if isinstance(json_data,str):
      json_data = json.loads(json_data)
    print("JSON Message recieved:")
    has_command = False
    has_kwargs = False
    for key,value in json_data.items():
      if key != "kwargs":
        print(key,":",value)
    if "kwargs" in json_data:
      has_kwargs = True
      for k,value in json_data["kwargs"].items():
        print(k,":")
        print(value)
        print()
    if "command" in json_data:
      has_command = True
    assert isinstance(json_data,dict), (
        f"Expected json data ({json_data}) to be parsed to a dict, got: {type(json_data)}")

    if has_command and has_kwargs:
      flask_signal_emitter.flask_signal.emit(json_data)

    return jsonify({'status_ok': "ok"})


class ViewerGUIApp:
  def __init__(self,state,view,controller):
    self.controller = controller
    self.view = view
    self.state = state




def main(dm=None,params=None,log=None):

  # first check that the necessary programs are available
  programs_to_check = ['npm', 'http-server']
  inaccessible_programs = check_program_access(programs_to_check)

  if inaccessible_programs:
    print(f"The following required programs are inaccessible or not found: {', '.join(inaccessible_programs)}")
    sys.exit()
  else:
    #print("All programs are accessible.")
    pass

  choice = None
  if params:
    if params.viewer_choice:
      choice = params.viewer_choice



  # start app
  QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
  qapp = QApplication(sys.argv)

  # get icon
  icon_path =  Path(__file__).parent / '../view/assets/icons/phenix/icon.icns'
  icon = QIcon(str(icon_path))
  qapp.setWindowIcon(icon)

  # choose viewer
  if not choice:
    choice_dialog = ViewerChoiceDialog()
    res = choice_dialog.exec_()

    # If a choice was made, show the main window
    if res != QDialog.Accepted:
      QMessageBox.warning(None, "No Choice", "No option was selected. Exiting application.")
      qapp.quit()
      sys.exit()



    choice = choice_dialog.choice
    if not params:
      params = DotDict()
      params.show_tab = "all"
    params.viewer_choice = choice
  # Set up a data manager if no provided
  if not dm:
    dm = DataManager()


  # DEBUG: load some data automatically
  #dm.process_model_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp.pdb")
  #dm.process_real_map_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp_calc.mrc")
  #dm.process_model_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp.cif")

  # Core top level object initialization
  state = State(dm)
  view = ViewerGUIView(params=params)
  controller = ViewerGUIController(parent=state,view=view,params=params)
  app = ViewerGUIApp(state,view,controller)

  # DEBUG: Sync references for test data
  state.signals.references_change.emit()


  # Reach into the Console tab to make variables accessible
  if params and params.show_tab:
    if 'all' in params.show_tab or 'console' in params.show_tab:
      try:
        import qtconsole
        app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'app': app})
        #include Selection dataclasses to build querys in console
        app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'Selection':Selection})
        app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'SelectionQuery':SelectionQuery})

      except:
        print("No qtconsole found")

  # # Create an instance of the event filter
  # globalEventFilter = GlobalEventFilter()

  # # Install the event filter on the QApplication instance
  # qapp.installEventFilter(globalEventFilter)


  # Flask set up
  controller.view.show()
  controller.flask_signal_emitter = flask_signal_emitter
  controller.flask_signal_emitter.flask_signal.connect(controller.update_from_remote)

  def run_flask_app():
    print("Running flask app on port: ",params.rest_server_port)
    flask_app.run(debug=True, port=params.rest_server_port, use_reloader=False)

  threading.Thread(target=run_flask_app, daemon=True).start()




  sys.exit(qapp.exec_())

if __name__ == '__main__':
  main()
