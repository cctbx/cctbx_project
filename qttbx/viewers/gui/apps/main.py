import json
import sys
import os
import threading
from pathlib import Path

from flask import Flask, request, jsonify
from PySide2.QtGui import QIcon
from PySide2.QtCore import QObject, Qt, Signal
from PySide2.QtWidgets import (
    QApplication,
    QDialog,
    QMessageBox
)
from iotbx.data_manager import DataManager

from ..view.apps.main import ViewerGUIView
from ..controller.apps.main import ViewerGUIController
from ..state.state import State
from ...core.selection import Selection
from . import ViewerChoiceDialog, check_program_access
from ...core.python_utils import DotDict

QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)


# # Start flask app to communicate with old Phenix GUI
# class FlaskSignalEmitter(QObject):
#     flask_signal = Signal(dict)

# # Create an instance of the signal emitter
# flask_signal_emitter = FlaskSignalEmitter()

# flask_app = Flask(__name__)

# @flask_app.route('/process_json', methods=['POST'])
# def process_json():
#     json_data = request.get_json()
#     self.log("type:",type(json_data))
#     self.log(json_data)
#     # Continue with actual json message
#     if isinstance(json_data,str):
#       json_data = json.loads(json_data)
#     self.log("JSON Message recieved:")
#     has_command = False
#     has_kwargs = False
#     for key,value in json_data.items():
#       if key != "kwargs":
#         self.log(key,":",value)
#     if "kwargs" in json_data:
#       has_kwargs = True
#       for k,value in json_data["kwargs"].items():
#         self.log(k,":")
#         self.log(value)
#         self.log()
#     if "command" in json_data:
#       has_command = True
#     assert isinstance(json_data,dict), (
#         f"Expected json data ({json_data}) to be parsed to a dict, got: {type(json_data)}")

#     if has_command and has_kwargs:
#       flask_signal_emitter.flask_signal.emit(json_data)

#     return jsonify({'status_ok': "ok"})


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
    #self.log("All programs are accessible.")
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

  # Core top level object initialization
  state = State(dm,params=params)
  view = ViewerGUIView(params=params)
  controller = ViewerGUIController(parent=state,view=view)
  app = ViewerGUIApp(state,view,controller)

  # DEBUG: Sync references for test data
  state.signals.references_change.emit()


  # Reach into the Console tab to make variables accessible
  app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'app': app})
  controller.view.show()


  # # Flask set up
  # controller.flask_signal_emitter = flask_signal_emitter
  # controller.flask_signal_emitter.flask_signal.connect(controller.update_from_remote)

  # def run_flask_app():
  #   state.log("Running flask app on port: ",params.rest_server_port)
  #   flask_app.run(debug=True, port=params.rest_server_port, use_reloader=False)

  # threading.Thread(target=run_flask_app, daemon=True).start()




  sys.exit(qapp.exec_())

if __name__ == '__main__':
  main()
