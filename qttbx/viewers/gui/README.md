## Phenix Viewer
A molecular viewer using Molstar in a QtWebEngineView window.

### Installation
Assuming you are already working on an active Phenix/CCTBX development branch:
1. Get the branch
```sh
git fetch -a
```
2. Checkout the branch
```
git checkout molstar-veiwer-base
```
3. Install dependencies from conda (or mamba)
```
mamba -c conda-forge install pyside2 nodejs qtawesome ipykernel qt-webengine qtconsole-base
```
4. Change to Phenix modules directory
5. Download the Molstar branch
```
git clone git@github.com:phenix-project/phenix-molstar.git
```
6. Change to the `phenix-molstar` directory
7. Install nodejs dependencies
```
npm install
```
8. Compile Typescript to Javascript
```
npm run build
9. Refresh the build to register the command line files.
```
10. Start the viewer with a model as an argument
```
phenix.start_molstar 1yjp.cif
```
Alternatively, run directly if command line files not recognized:
```
phenix.python <MODULES_DIRECTORY>/cctbx_project/qttbx/command_line/start_molstar.py 1yjp.cif
```

### Development guide
The code is structured with the intention to follow a strict model-view-controller (MVC) architecture, to the degree that each component has its own subdirectory in the source code in the ```qttbx/viewers/gui``` subdirectory. 
```
cctbx_project/
└── qttbx/
    ├── command_line/
    │   └── start_molstar.py
    ├── programs/
    │   └── start_molstar.py
    └── viewers/
        ├── molstar/
        │   └── molstar.py
        └── gui/
            ├── README.md
            ├── controller/
            ├── model/
            └── view/
```
The contents of each directory are described below:

1. **`qttbx/command_line/start_molstar.py`**: A stub to register the program as a command line tool.
2. **`qttbx/programs/start_molstar.py`**: The starting point of the program using the `ProgramTemplate` class. All parameter parsing and file IO should be done here before the GUI is instantiated.
3. **`qttbx/viewers/molstar/`**: This directory contains the python code necessary to interface with the Molstar web app running in a QtWebEngineView window. The `molstar.py` file defines a Python api that can be used to affect the Molstar graphics by running corresponding Javascript snippets. Anything in the rest of the GUI that affects the Molstar graphics must come through the `MolstarGraphics` class in this file.
4. **`qttbx/viewers/gui/`**:
   * **`controller`**: Objects in this directory are a subclass of the `Controller` base class. The controllers are responsible for the application logic. All controller instances have access to a shared program state. For modularity, each controller instance tends to be responsible for controlling a single view.
   * **`model`**: Objects in this directory are designed to store data and manage the program state. This includes:
     * **`Data`**: Objects to hold data, often these are just a thin wrapper around external data structures like an `mmtbx.model.manager`
     * **`Ref`**: Objects which manage data for the lifetime of the GUI. The purpose of these objects is to separate GUI specific state (stored in Ref objects) from non-GUI data.
     * **`State`**: A single instance which contains all the 'Ref' objects managed by the GUI. All controllers have access to the same 'State' singleton. The state object facilitates controller communication through a "signal bus", whcihmakes use of the Qt signals/slots paradigm to enable controllers to comunicate in a decoupled and asynchronous manner. All the model data structures including `State` are unable to access any of the controllers or the view objects

   * **`view`**: Objects in this directory are a subclass of `QObject`. This means all the GUI widgets the user "sees" are defined here. The view objects are unable to access any of the controllers or any of the model objects. An example program flow: a user pushes a button. The view emits a signal. For anything to happen, a controller must be connected and listening to that signal. If so, then the controller will determine what to do. If the final step is to update the view object, it is the _controller_ who calls the appropriate methods on the view. A major advantage of this setup is GUI elements can be added without modifying ANY code in the `controller` or `model` subdirectories.

### Programming Example
Here is a simplified example showing how different object might interact using code. We will define behavior to: 
1. Respond to the user (Display the active selection string if a user presses a button)
2. Respond to changes in state (Display the active selection string whenever it changes)
```python
#
# Class definitions
#

class SignalBus(QObject):                        # Define a list of signals for the state
  active_selection = Signal()                    # A signal for a new active selection

class State:                                     # Define an object to store program state
  def __init__(self,dm):                         # Initialize the state with only the data manager
    self.signals = SignalBus()
    self.dm = dm
    self._active_selection_ref = None            # The only state this class manages is the active selection

  @property
  def active_selection_ref(self):
    return self._active_selection_ref

  @active_selection_ref.setter                    # When the active selection is changed, the state emits a signal
  def active_selection_ref(self,value):
    self._active_selection_ref = value
    self.signals.active_selection.emit(self.active_selection_ref)



class SelectionView(QWidget):                     # QWidget is a GUI widget, ultimately a subclass of QObject
  def __init__(self,parent=None):                 # Initialization can be with or without a parent QObject
    super().__init__(parent=parent)

    self.print_sel_button = QPushButton()         # define a button using a built-in Qt button widget
    self.active_sel_label = QLabel()              # Define a label widget to hold text


class SelectionController(Controller):            # The controller base class ensures all controllers behave similarly
 def __init__(self,parent=None,view=None):        #  such as having access to self.state and self.view
    super().__init__(parent=parent,view=view)

                                                  # connect to the view's 'clicked' signal (Defined inside QPushButton)
    self.view.print_sel_button.clicked.connect(self.print_active_selection)

                                                  # connect to the state's signal to get notified of a new active selection
    self.state.signals.active_selection.connect(self.show_active_selection) 

def show_active_selection(self):                  # This function is a "Slot" in Qt terminology, connected to signals.
  selection_ref = self.state.active_selection_ref # Get the active selection from the shared state object
  selection = selection_ref.data                  # Access the raw data, a selection object.
  sel_str = selection.string                      
  self.view.active_sel_label.setText(sel_str)     # Update the view in response

def print_active_selection(self):                 # Another "Slot" connected to the print_sel_button
  selection_ref = self.state.active_selection_ref
  selection = selection_ref.data                  
  sel_str = selection.string                      
  print(sel_str)                                  # Print. A different response


class SelectionRef(Ref):                          # Define a new SelectionRef using the Ref base class
  def __init__(self,data,show=True):              # Storing GUI specific data, such as whether it should be visible
    super().__init__(data=data, show=show)        #   would be undesireable to store on a raw data object


                                                  
class Selection:                                  # Define a data object
  def __init__(self,selection_string):            # For data use something that is useful outside of the GUI,
    self.string = selection_string                #   potentially an existing external datastructure.                          
  

#
# Instantiation
#

# Initialize the GUI
view = SelectionView()                            # No arguments needed for the view
state = State(dm)                                 # Initialize State with the datamanager from the program template
controller = SelectionController(state,view)      # initialize the controller with the state and the view it should manage


# Add some new data to the state
selection = Selection("resname ALA and resseq 10:20")
selection_ref = SelectionRef(selection)
state.active_selection_ref = selection_ref        # causes the label to be updated in the view

# Additionally, whenever the user clicks on the 'print_sel_button', the active selection will be printed.
```

