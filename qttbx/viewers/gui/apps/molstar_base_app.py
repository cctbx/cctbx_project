"""
This is the root file for the molstar base app, which implements the minimum usable
  Molstar viewer running in a QtWebView.
"""
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
from libtbx import group_args

from qttbx.viewers.gui.view.apps.molstar_base_app import MolstarBaseAppView
from qttbx.viewers.gui.controller.apps.molstar_base_app import MolstarBaseAppController
from qttbx.viewers.gui.modelstate import State
#from ...core.selection import Selection
from qttbx.viewers.gui.apps import check_program_access

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


class MolstarBaseApp:
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



  # start app
  QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
  qapp = QApplication(sys.argv)

  # get icon
  icon_path =  Path(__file__).parent / '../view/assets/icons/phenix/icon.icns'
  icon = QIcon(str(icon_path))
  qapp.setWindowIcon(icon)

  # Set up a data manager if no provided
  if not dm:
    dm = DataManager()

  # Core top level object initialization
  state = State(dm,params=params)
  view = MolstarBaseAppView(params=params)
  controller = MolstarBaseAppController(parent=state,view=view)
  app = MolstarBaseApp(state,view,controller)

  # Start
  controller.view.show()
  sys.exit(qapp.exec_())

if __name__ == '__main__':
  main()
