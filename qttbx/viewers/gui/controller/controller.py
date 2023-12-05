from ..state.state import State

class Controller:

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

  @property
  def state(self):
    return self.parent.state


class ViewController:
  """
  TODO: Remove this file
  The 'Controller' in MVC architecture. 

  Principles:
  
  1. The controller is 'aware' of both the Views (Molstar viewer/ QT GUI) and the State
  2. The Views/State have no awareness of each other (ie, Views cannot read/write on state)
  3. The Views/State have no awareness of the Controller. State/Views define signals and emit. 
     Controller has slots to react to changes, and then directly calls methods on state/views

  """
  def __init__(self,state,viewer,gui):
    self.state = state
    self.viewer = viewer
    self.gui = gui
    
  #
  # Slots for Selections
  #





  




  
