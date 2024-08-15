"""
The Controller base class in the model-view-controller (MVC) architecture as implemented here
"""
import time
from PySide2.QtWidgets import QMessageBox
from ..state.state import State

class Controller:
  """
  The 'Controller' in MVC architecture. 'State' is used as an alternative to 'Model' in 
    documentation to differentiate from a molecular model.

  Principles:

  1. The controller is 'aware' of both the Views and the Model (State)
  2. The Views and the Model (State) have no awareness of each other (ie, Views cannot read/write on state)
  3. The Views and the Model (State) have no awareness of the Controller. State/Views define signals and emit.
     Controllers implement slots to react to changes, and then directly call methods on state/views

  """
  @classmethod
  def from_empty(cls):
    state = State.from_empty()
    return cls(parent=state,view=None)

  def __init__(self,parent=None,view=None):
    assert parent is not None, 'Controller must have parent'
    assert hasattr(parent,'state'), 'Parent must implement "state" property'
    #assert view is not None, 'Controller must have a view'
    self.parent = parent
    self.view = view
    self.debug = True
    self.log = None
    if parent:
      self.log = self.parent.log
    self.notifications = {}
    # self.state.signals.warning.connect(self.warning)
    # self.state.signals.notification.connect(self.notification)
    # self.state.signals.error.connect(self.error)

  @property
  def params(self):
    return self.parent.params

  @property
  def state(self):
    return self.parent.state

  def warning(self,message):
    self._notify(message,title="Warning")
  def notification(self,message):
    self._notify(message,title="Notification")
  def error(self,message):
    self._notify(message,title="Error")

  def _notify(self, notification, title="Notification"):
    message =  notification["msg"]
    uuid =  notification["uuid"]
    if uuid not in self.notifications:
      self.notifications[uuid] =  notification
      # Check if a message box is already open
      if hasattr(self, 'message_box') and self.message_box.isVisible():
        return  # Exit if a message box is already showing

      # Create the message box
      self.message_box = QMessageBox()
      self.message_box.setWindowTitle(title)
      self.message_box.setText(message)
      self.message_box.setIcon(QMessageBox.Information)
      
      # Add an OK button
      self.message_box.addButton(QMessageBox.Ok)
      
      # Set the default button (optional)
      self.message_box.setDefaultButton(QMessageBox.Ok)

      # Show the message box
      self.message_box.exec_()

      # Clean up after the message box is closed
      self.message_box = None

      