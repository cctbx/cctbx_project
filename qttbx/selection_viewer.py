import sys
import types
import requests
import json
from collections import defaultdict
from PySide2.QtCore import QAbstractTableModel, Qt, Slot
from PySide2.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QTabWidget, QTableView
import pandas as pd
from urllib.parse import urlencode
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation as R
from PySide2 import QtCore
from PySide2.QtCore import Qt, QTimer
from PySide2.QtWidgets import QApplication, QGraphicsView, QGraphicsScene, QGraphicsRectItem
from PySide2.QtCore import Qt, QRectF
from PySide2.QtGui import QBrush, QColor


from PySide2.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QLabel,
    QMessageBox,
    QLineEdit,
    QMainWindow,
    QProgressBar,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
    QFileDialog,
    QTabWidget
)
from phenix.programs.hotspots import Program as HotSpots
from phenix.program_template import call_program
from iotbx import phil
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager
from qttbx.viewers.chimerax import ChimeraXViewer
from cctbx.array_family import flex
from iotbx.pdb.atom_selection import selection_string_from_selection
from mmtbx.validation.ramalyze import ramalyze

from qttbx.sel_convert_chimera import (
  translate_phenix_selection_string,
  convert_selection_str_to_int,
  convert_selection_int_to_str
)

from cctbx_utils import model_to_df
from cctbx_utils import get_restraint_dfs_from_model

def convert_id_str_to_selection_str(id_str):
  # convert pdb id_str to selection string 
  values = [id_str[:2], id_str[2:8], id_str[8:11], id_str[11:]]
  values = [v.strip() for v in values]
  keys = ["chain","resid","resname","name"]
  
  parts = []
  parts_res = []
  for key,value in zip(keys,values):
      if len(value)>0:
          parts.append(key+" "+value)
          if key != "name":
              parts_res.append(key+" "+value)
              
  selection_string = " and ".join(parts)
  return selection_string


def get_rama_dfs(model):
  ramalyze_result = ramalyze(model.get_hierarchy())
  rama_type_dict = {
  0:"General",
  1:"Glycine",
  2:"CisPro",
  3:"TransPro",
  4:"PrePro",
  5:"IleVal",
  }
  ramalyze_type_dict = {0:"Outlier",
  1:"Allowed",
  2:"Favored",
  3:"Any",
  4:"Not Favored",
  }
    



  data = {
  "Chain":[],
  "Residue Name":[],
  "Residue Seq":[],
  "Residue Type":[],
  "Ramalyze Type":[],
  "Score":[],
  "Phi":[],
  "Psi":[],
  "Integer Selection":[]
  }

  for result in ramalyze_result.results:
    
    selection_string = convert_id_str_to_selection_str(result.id_str())
    selection_int = convert_selection_str_to_int(model,selection_string)
    data["Integer Selection"].append(selection_int)
    data["Residue Type"].append(rama_type_dict[result.rama_type])
    data["Ramalyze Type"].append(result.ramalyze_type())
    data["Chain"].append(result.chain_id)
    data["Residue Seq"].append(result.resseq.strip())
    data["Residue Name"].append(result.resname)
    data["Score"].append(result.score)
    data["Phi"].append(result.phi)
    data["Psi"].append(result.psi)

  #print(data)
  return pd.DataFrame(data)
    


def rotation_from_abc(a, b, c):
    # Compute vectors based on the given points
    x = a - c
    y =  b - ((a + c) / 2)
    z = np.cross(x, y)

    # Normalize each vector
    x_norm = x / np.linalg.norm(x)
    y_norm = y / np.linalg.norm(y)
    z_norm = z / np.linalg.norm(z)

    # Create a rotation matrix using the normalized vectors
    rotation_matrix = np.column_stack((x_norm, y_norm, z_norm))

    # Convert the rotation matrix to a Rotation object
    rotation = R.from_matrix(rotation_matrix)
    
    return rotation

def translate_along_z(R, a, b):
    # Get the z-axis direction from R
    z_direction = R.as_matrix()[:, 2]
    
    # Scale the z-direction by b units
    offset = b * z_direction
    
    # Translate a by the offset
    new_a = a + offset
    
    return new_a

def matrix_string_from_R(R,center=np.array([0,0,0]),offset=10):
    
    center = translate_along_z(R,center,offset)
    matrix_flattened = np.vstack([R.as_matrix().T,center]).T.flatten()
    return ",".join([str(round(v,4)) for v in matrix_flattened])

def chimerax_composition_from_selection_response(response):
    
    def split_into_pairs(lst):
        return [[lst[i], lst[i + 1]] for i in range(0, len(lst), 2)]


    content =  response.json()["log messages"]["note"][0]
    composition = {}
    try:
        for value,key in split_into_pairs(content.split()[:-1]):
            composition[key.strip(",")] = int(value.strip(","))
    except:
        assert False, "Unable to parse composition from :"+"".join(content)
    return composition



def group_by_columns(df, columns):
    # Make a copy of the DataFrame to avoid modifying the original
    df_copy = df.copy()

    # Reset index if any of the groupby columns is the index
    df_copy.reset_index(drop=True, inplace=True)

    # Create a column with the original index
    df_copy['Atom Index'] = df.index
    
    # Group by the specified columns and aggregate
    grouped_df = df_copy.groupby(columns).agg({
        col: 'first' if col in columns else list
        for col in df_copy.columns
    })
    
    # Here the grouped columns are the index, so we convert the index to a normal column
    grouped_df.reset_index(drop=True, inplace=True)

    # Remove columns where rows differ
    for col in grouped_df.columns:
        if col not in columns and col != 'Atom Index':
            grouped_df.drop(col, axis=1, inplace=True)

    # Sort by the columns
    grouped_df.sort_values(by=columns, inplace=True)
    
    return grouped_df




  
class FastTableView(QTableView):
    mouseReleased = QtCore.Signal()

    def __init__(self, parent=None,default_col_width=150):
        super(FastTableView, self).__init__(parent)
        self.horizontalHeader().setDefaultSectionSize(default_col_width)

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        self.mouseReleased.emit()

      

class PandasModel(QAbstractTableModel):
  def __init__(self, df=pd.DataFrame(), parent=None, suppress_columns=[]):
    QAbstractTableModel.__init__(self, parent=parent)
    self._df = df
    self.suppress_columns = set(suppress_columns)
    self.visible_columns = [col for col in self._df.columns if col not in self.suppress_columns]

  def headerData(self, section, orientation, role=Qt.DisplayRole):
    if role != Qt.DisplayRole:
      return None
    if orientation == Qt.Horizontal:
      return self.visible_columns[section]
    else:
      return self._df.index[section]

  def data(self, index, role=Qt.DisplayRole):
    if role != Qt.DisplayRole:
      return None
    return str(self._df.loc[index.row(), self.visible_columns[index.column()]])

  def setData(self, index, value, role):
    if not index.isValid() or role != Qt.EditRole:
      return False
    self._df.loc[index.row(), self.visible_columns[index.column()]] = value
    self.dataChanged.emit(index, index, (Qt.DisplayRole,))
    return True

  def rowCount(self, parent=None):
    return len(self._df.values)

  def columnCount(self, _, parent=None):
    return len(self.visible_columns)

  def flags(self, index):
    return Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable

class ClickableRect(QGraphicsRectItem):
  def __init__(self, x, y, width, height, color):
    super(ClickableRect, self).__init__(x, y, width, height)
    self.setBrush(QBrush(QColor(color)))

  def mousePressEvent(self, event):
    print(f"Clicked on rectangle with color: {self.brush().color().name()}")

class Histogram(QGraphicsView):
  def __init__(self, data):
    super(Histogram, self).__init__()
    self.scene = QGraphicsScene()
    self.setScene(self.scene)

    x = 0
    for height, color in data:
      rect = ClickableRect(x, 0, 50, -height, color)
      self.scene.addItem(rect)
      x += 60



class SelectionViewerWindow(QMainWindow):
    todo = """

    Organization:
    The data manager stores filepaths. Each filepath has a "name", Path(filepath).name. 
    This name is the intermediate identifier between cctbx and chimera. Chimera has a "spec" identifier

    At any moment, this class has:
      1. a current model path
      2. a current model name
      3. a current spec
    

    TODO: phenix logo and larger title
    TODO: focusing
    TODO: clear selection button
    TODO: periodically check if chimera disconnected
    TODO: Remove redundant phenix selection windows, only need textedit, with enter key
    TODO: 'through' bug
    TODO: warn about through changed to :
    TODO: altloc bug
    TODO: Checkbox to only translate even when chimera open
    TODO: check if within includes symmetry, re-box to P1
    TODO: sync model list of chimera with model list data manager
    TODO: use self.notes as intermediate for logging
    TODO: accept integer selections from gui
    TODO: view debug atom pairs

    
    
    """
  
    
  
    def __init__(self):
        super().__init__()
        self.debug = True
        self.strict_selections = False
        self._data_manager = None # Stores models
        self._viewer = None # chimerax viewer
        self._current_model_path = None
        self._current_model_name = None
        self._current_model_spec = None
      
        self._model_name_to_path = {}
        self._model_name_to_spec = {}

        self.dfs = {}
        self.df_residues_q = None
        self.df_residues_cc = None
        self.dfs_restraints = {}
        self.restraint_types = ["bond","nonbonded","angle","dihedral","chirality"]#,"planarity","parallelity"]
        self.dfs_validation = defaultdict(dict)
        self.validation_types = ["ramachandran"]
        self.validation_types_map = ["q-score","ccatoms"]
        self.process = True # make restraints?
        self.notes = []
        self.truncation_limit = 100
        self.setWindowTitle("Phenix Viewer")
        self.setGeometry(100, 100, 400, 200) # Set window position and size
        
        


        layout = QVBoxLayout()
        
        # Header
        self.label_header = QLabel("""
        This utility translates Phenix selection language to ChimeraX selection language.
        It enables viewing Phenix selections graphically.
        """
                                  )
        
        layout.addWidget(self.label_header)
        
        # files
        self.filename = None
        self.upload_button = QPushButton("Load Model File")
        self.upload_button.setFixedWidth(150) 
        self.upload_button.clicked.connect(self.upload_file)
        layout.addWidget(self.upload_button)
        
        self.filename_map = None
        self.upload_button_map = QPushButton("Load Map File")
        self.upload_button_map.setFixedWidth(150) 
        self.upload_button_map.clicked.connect(self.upload_file_map)
        layout.addWidget(self.upload_button_map)
        

        # Label
        self.filename_label = QLabel()
        layout.addWidget(self.filename_label)
        
        # model number
        self.label_model_number = QLabel("Model in ChimeraX")
        layout.addWidget(self.label_model_number)
        
        self.combo_box_models = QComboBox()
        self.update_model_list()
        layout.addWidget(self.combo_box_models)
        self.combo_box_models.currentIndexChanged.connect(self.on_combobox_models_changed)
            
        # # workaround checkbox
        # self.checkbox_workaround = QCheckBox("If selection fails, try with simpler vocabulary.")
        # self.checkbox_workaround.setChecked(True)
        # layout.addWidget(self.checkbox_workaround)
        
        # Label
        self.title = QLabel("\nEnter Phenix selection string")
        layout.addWidget(self.title)

        # Text box
        self.text_input = QLineEdit()
        self.text_input.returnPressed.connect(self.button_clicked)
        layout.addWidget(self.text_input)
        
        # Button
        self.button = QPushButton("Select",self)
        self.button.setFixedWidth(100) 
        self.button.clicked.connect(self.button_clicked)
        layout.addWidget(self.button)
        
        # output
        self.phenix_selection_result_title = QLabel("\nPhenix selection used:")
        layout.addWidget(self.phenix_selection_result_title)
        
        self.phenix_selection_string = QTextEdit()
        layout.addWidget(self.phenix_selection_string)
        
        self.phenix_selection_result = QLabel()
        layout.addWidget(self.phenix_selection_result)
        
        self.text_output_title = QLabel("\nChimeraX selection used:")
        layout.addWidget(self.text_output_title)
        
        self.text_output = QTextEdit()
        layout.addWidget(self.text_output)
        
        self.text_response = QLabel()
        layout.addWidget(self.text_response)
        
        # this is where we notify users of any errors or workarounds used
        self.output_notes = QLabel()
        layout.addWidget(self.output_notes)
        
        self.output_widgets = [self.phenix_selection_result_title,
                          self.phenix_selection_string,
                          self.phenix_selection_result,
                          self.text_output_title,
                          self.text_output,
                          self.text_response,
                          self.output_notes]
                          
        
        # connect to chimera
        self.checkbox_chimerax = QCheckBox("Remote ChimeraX connected")
        # Connect the checkbox to a function that updates the label
        self.checkbox_chimerax.stateChanged.connect(self.checkbox_changed)
        layout.addWidget(self.checkbox_chimerax)
        

        # Create the QTabWidget.
        self.tab_widget = QTabWidget()
        
        
        self.tab1_content = QWidget()
        self.tab1_content.setLayout(layout)
        self.tab_widget.addTab(self.tab1_content, "Selections")

        # Atoms tab
        self.tab2_content = FastTableView()
        #self.tab2_content.clicked.connect(self.on_atom_select)
        self.tab2_content.setSelectionBehavior(QTableView.SelectRows)
        self.tab2_content.mouseReleased.connect(self.on_mouse_released_atoms)
        self.tab_widget.addTab(self.tab2_content, "Atoms")
      
        # Residues tab
        self.tab_comp_content = FastTableView()
        self.tab_comp_content.setSelectionBehavior(QTableView.SelectRows)
        self.tab_comp_content.mouseReleased.connect(self.on_mouse_released_comp)
        self.tab_widget.addTab(self.tab_comp_content, "Residues")

        # restraints tab
        self.restraints_tab = QTabWidget()
        self.tab_widget.addTab(self.restraints_tab, "Restraints")

        for restraint_type in self.restraint_types:
            # Create a FastTableView for each restraint_type
            setattr(self, f"tab_content_{restraint_type}", FastTableView())
            tab_content = getattr(self, f"tab_content_{restraint_type}")
            tab_content.setSelectionBehavior(QTableView.SelectRows)
            tab_content.mouseReleased.connect(self.on_mouse_released_restraints)

            # Add the FastTableView to the restraints_tab
            self.restraints_tab.addTab(tab_content, restraint_type.capitalize())



        # validation tab
        self.validation_tab = QTabWidget()
        self.tab_widget.addTab(self.validation_tab, "Validation")

        for validation_type in self.validation_types:
            # Create a FastTableView for each restraint_type
            setattr(self, f"tab_content_{validation_type}", FastTableView())
            tab_content = getattr(self, f"tab_content_{validation_type}")
            tab_content.setSelectionBehavior(QTableView.SelectRows)
            tab_content.mouseReleased.connect(self.on_mouse_released_validation)

            # Add the FastTableView to the validation tab
            self.validation_tab.addTab(tab_content, validation_type.capitalize())



       # Hotspots
        self.tab_hot_content = FastTableView()
        self.tab_hot_content.setSelectionBehavior(QTableView.SelectRows)
        self.tab_hot_content.mouseReleased.connect(self.on_mouse_released_hotspots)
        self.tab_widget.addTab(self.tab_hot_content, "Hotspots")
    
        # Settings
        self.settings_layout = QVBoxLayout()

        # Create and add a button to the layout
        self.hotspot_color_button = QPushButton("Color Hotspots")
        self.hotspot_color_button.clicked.connect(self.on_hotspot_button_click)  # Connect to a new slot function
        self.settings_layout.addWidget(self.hotspot_color_button)

        # Create a new QWidget to house the settings_layout
        self.settings_tab_content = QWidget()
        self.settings_tab_content.setLayout(self.settings_layout)

        # Add the new QWidget to the tab_widget
        self.tab_widget.addTab(self.settings_tab_content, "Settings")


        # # hotspots tab
        # self.tab_hot_content = FastTableView(default_col_width=75)
        # self.tab_hot_content.setSelectionBehavior(QTableView.SelectRows)
        # self.tab_hot_content.mouseReleased.connect(self.on_mouse_released_hotspots)
        # self.tab_widget.addTab(self.tab_hot_content, "Hotspots")
        # # Create Button
        # self.hotspot_button = QPushButton("Click Me!")
        # self.hotspot_button.clicked.connect(self.on_hotspot_button_click)

        # # Create Vertical Layout for hotspots
        # self.hotspot_layout = QVBoxLayout()

        # # Add Button and Table to the layout
        # self.hotspot_layout.addWidget(self.hotspot_button)
        # self.hotspot_layout.addWidget(self.tab_hot_content)  # Re-use the table you already created, instead of creating a new one

        # # Create a new QWidget to house the hotspot_layout
        # self.hotspot_tab_content = QWidget()
        # self.hotspot_tab_content.setLayout(self.hotspot_layout)

        # # Add the new QWidget to the tab_widget
        # self.tab_widget.addTab(self.hotspot_tab_content, "Hotspots")


        # ...

        # Remove this line as the layout is already set to tab_widget which is the central widget.
        # self.setLayout(self.hotspot_layout)

        # Set tab_widget as the central widget
        self.setCentralWidget(self.tab_widget)


        # set a timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.periodic_update)
        self.timer.start(5000) # 5 seconds


    def on_hotspot_button_click(self):
      # color hotspots
      current_model = self.tab_hot_content.model()
      df = current_model._df
      # TODO: Move this to a generic color function
      if "Color" in df.columns:
        for i,row in df.iterrows():
          color = row["Color"]
          sel_string = row["Sel Str Phenix"]
          #print(sel_string)
          response = self.view_selection(sel_string)
          cmd = "color sel "+color
          #print(cmd)
          response = self.viewer.send_command(cmd)
          #print(response)


    def get_selected_rows(self,content):
        # Get a list of all selected cell indexes.
        indexes = content.selectedIndexes()
    
        # Extract the row numbers for each index.
        rows = [index.row() for index in indexes]
    
        # Get a list of unique rows.
        unique_rows = list(set(rows))
    
        return unique_rows
    @Slot()
    def on_selection_changed_atoms(self, selected, deselected):
        selected_rows = self.get_selected_rows(self.tab2_content)
        selection_int = flex.int(selected_rows)
        phenix_selection_string_modified = convert_selection_int_to_str(self.model,selection_int)
        success = self.view_selection(phenix_selection_string_modified)
        self.focus_selection()

      
    @Slot()
    def on_mouse_released_atoms(self):
        selected = self.tab2_content.selectionModel().selection()
        deselected = QtCore.QItemSelection()
        self.on_selection_changed_atoms(selected, deselected)
      
    @Slot()
    def on_selection_changed_comp(self, selected, deselected):
        selected_rows = self.get_selected_rows(self.tab_comp_content)
        model = self.tab_comp_content.model()
        df = model._df
        selected_rows = df["Atom Index"].iloc[selected_rows]
        def flatten(df):
          def is_iterable(x):
              try:
                  iter(x)
              except TypeError:
                  return False
              return True

          return [item for sublist in df.values.flatten() for item in sublist if is_iterable(sublist)]

        selected_rows = flatten(selected_rows)
        selection_int = flex.int(selected_rows)
        phenix_selection_string_modified = convert_selection_int_to_str(self.model,selection_int)
        success, chimerax_selection = self.view_selection(phenix_selection_string_modified,return_chimerax_string=True)
        if self.model.selection(phenix_selection_string_modified).count(True) == \
          self.model.selection(phenix_selection_string_modified+" and protein").count(True):
          self.focus_single_residue(chimerax_selection)
        else:
           self.focus_selection()

      
    @Slot()
    def on_mouse_released_comp(self):
        selected = self.tab2_content.selectionModel().selection()
        deselected = QtCore.QItemSelection()
        self.on_selection_changed_comp(selected, deselected)

    @Slot()
    def on_selection_changed_validation(self, selected, deselected):
        content = self.validation_tab.currentWidget()
        selected_rows = self.get_selected_rows(content)
        model = content.model()
        df = model._df
        iseq_columns = []
        has_explicit = False
        for column in df.columns:
          if column in ["Integer Selection"]:
            has_explicit = True
            iseq_columns.append(column)
        if not has_explicit:
          assert "i_seqs" in df.columns
          iseq_columns.append("i_seqs")
        #print(iseq_columns)
        #print(selected_rows)
        selected_rows = df[iseq_columns].iloc[selected_rows]
        #print(selected_rows)

        def flatten(df):
          def is_iterable(x):
            if isinstance(x,(int,np.int64)):
              return False
            try:
                iter(x)
            except TypeError:
                return False
            return True
          
          return [item for sublist in df.values.flatten() for item in sublist if is_iterable(sublist)]

        # flatten the df of indices

        selected_rows = flatten(selected_rows)
        selected_rows = [int(e) for e in selected_rows]
        selection_int = flex.int(selected_rows)
        phenix_selection_string_modified = convert_selection_int_to_str(self.model,selection_int)

        success = self.view_selection(phenix_selection_string_modified)
        self.focus_selection()

    @Slot()
    def on_mouse_released_validation(self):
        content = self.validation_tab.currentWidget()
        selected = content.selectionModel().selection()
        deselected = QtCore.QItemSelection()
        self.on_selection_changed_validation(selected, deselected)

    @Slot()
    def on_mouse_released_validation_map(self):
        content = self.validation_tab_map.currentWidget()
        selected = content.selectionModel().selection()
        deselected = QtCore.QItemSelection()
        self.on_selection_changed_validation_map(selected, deselected)

    @Slot()
    def on_selection_changed_validation_map(self, selected, deselected):
        content = self.validation_tab_map.currentWidget()
        selected_rows = self.get_selected_rows(content)
        model = content.model()
        df = model._df
    
    @Slot()
    def on_selection_changed_hotspots(self, selected, deselected):
        content = self.tab_hot_content
        selected_rows = self.get_selected_rows(content)
        model = content.model()
        df = model._df
        selected_rows = [int(e) for e in selected_rows]
        phenix_selection_string = "or".join([df.iloc[i]["Sel Str Phenix"] for i in selected_rows])
        success = self.view_selection(phenix_selection_string)
        self.focus_selection()

    @Slot()
    def on_mouse_released_hotspots(self):
        content = self.tab_hot_content
        selected = content.selectionModel().selection()
        deselected = QtCore.QItemSelection()
        self.on_selection_changed_hotspots(selected, deselected)

    @Slot()
    def on_mouse_released_restraints(self):
        content = self.restraints_tab.currentWidget()
        selected = content.selectionModel().selection()
        deselected = QtCore.QItemSelection()
        self.on_selection_changed_restraints(selected, deselected)
    
    @Slot()
    def on_selection_changed_restraints(self, selected, deselected):
        content = self.restraints_tab.currentWidget()
        selected_rows = self.get_selected_rows(content)
        model = content.model()
        df = model._df
        iseq_columns = []
        has_explicit = False
        for column in df.columns:
          if column in ["i","j","k","l"]:
            has_explicit = True
            iseq_columns.append(column)
        if not has_explicit:
          assert "i_seqs" in df.columns
          iseq_columns.append("i_seqs")
        #print(iseq_columns)
        #print(selected_rows)
        selected_rows = df[iseq_columns].iloc[selected_rows]
        #print(selected_rows)

        def flatten(df):
          def is_iterable(x):
            if isinstance(x,(int,np.int64)):
              return False
            try:
                iter(x)
            except TypeError:
                return False
            return True
          
          return [item for sublist in df.values.flatten().astype(int) for item in sublist if is_iterable(sublist)]

        # flatten the df of indices
        try:
          selected_rows = flatten(selected_rows)
          selected_rows = selected_rows.astype(int)
        except:
          selected_rows = selected_rows.values.flatten()
          selected_rows = selected_rows.astype(int)
        
        #print(selected_rows)
        selected_rows = [int(e) for e in selected_rows]
        selection_int = flex.int(selected_rows)
        phenix_selection_string_modified = convert_selection_int_to_str(self.model,selection_int)
        success = self.view_selection(phenix_selection_string_modified)
        self.focus_selection()
      
    @property
    def viewer(self):   
      return self._viewer
    
    @viewer.setter
    def viewer(self,value):
      assert isinstance(value,ChimeraXViewer), "viewer must be instance of ChimeraXViewer"
      self._viewer = value
    
    @property
    def viewer_status_alive(self):
      #print("Checking status: ",type(self.viewer.url),self.viewer.url)
      if self.viewer is not None:
        if self.viewer.url != None:
          return self.viewer.is_alive()
      return False

    @property
    def data_manager(self):
      if self._data_manager is None:
        self._data_manager = DataManager()
      return self._data_manager
    
    @property
    def model(self):
      if self._current_model_path is None:
        return None
      else:
        model = self.data_manager.get_model(filename=self.current_model_path)

        if model.crystal_symmetry() is None:
          from cctbx.maptbx.box import shift_and_box_model
          model= shift_and_box_model(model)
        model = model.select(model.selection("protein"))
        return model

    def add_map(self,filepath):
      filepath = Path(filepath)
      # load in cctbx/phenix
      self.data_manager.process_real_map_file(str(filepath))
      cmd = "volume #2 style solid opacity 0.5" # TODO: do not hardcode model number!
      response = self.viewer.send_command(cmd)

    def add_model(self,filepath):
      filepath = Path(filepath)
      # load in cctbx/phenix
      self.data_manager.process_model_file(str(filepath))
      
      self.current_model_path = str(filepath)
      self._model_name_to_path[filepath.name] = str(filepath)
      df = model_to_df(self.model)
      self.dfs[filepath.name] = df

      # add restraints dfs
      if self.process:
        self.dfs_restraints[self.current_model_name] = get_restraint_dfs_from_model(self.model,lazy_build=True)

      # TODO: enable changing model for validation


      # add validation dfs
      self.dfs_validation[self.current_model_name]["ramachandran"] = get_rama_dfs(self.model)
          
      
    @property
    def current_model_path(self):
      return self._current_model_path
      
    @current_model_path.setter
    def current_model_path(self,value):
      self._current_model_path = value
      self._current_model_name = Path(value).name

    @property
    def current_model_name(self):
      return self._current_model_name
      
    @current_model_name.setter
    def current_model_name(self,value):
      self._current_model_name = value
      self._current_model_path = self._model_name_to_path[value]
      
    @property
    def model_paths(self):
      return [str(Path(file)) for file in self.data_manager.get_model_names()]

    @property
    def model_names(self):
      return [str(Path(path).name) for path in self.model_paths]
  

     
    def calculate_q_score(self,state):
      pass
      # assert self._has_model and self._has_map, "Map and model are required"
      # model = self.model
      # mm = self.data_manager.get_map_manager()
      # mmm = map_model_manager(model=model,map_manager=mm)
      # q = qscore_np(mmm,n_probes=16)
      # df_atoms = model_to_df(model)
      # df_atoms["Q-score"] = q
      # groupby_labels= ["asym_id","seq_id","comp_id"]
      # mean_labels = ["Q-score"]
      # drop_labels = [label for label in df_atoms.columns if label not in set(groupby_labels+mean_labels)]
      # df_atoms = df_atoms.drop(columns=drop_labels)

      # agg_dict = {label:"mean" for label in mean_labels}
      # self.df_residues_q = df_atoms.groupby(groupby_labels).agg(agg_dict).reset_index()
      
        

    
    def load_file_in_chimera(self,filename):
      if not self.viewer_status_alive:
        self.viewer_start()
      assert self.viewer_status_alive, "Cannot load file if viewer is not running"
      response = self.viewer.send_command(["open "+filename,"set bgcolor white","color #1 #999fffff"]) 
      # TODO: better color handling right off the bat
      
      # update dropdown
      self.update_model_list()
    

    def upload_file_map(self,state):
      # get file name
      filename, _ = QFileDialog.getOpenFileName(self, "Select File")
      self.filename_label.setText(filename)
      self.load_file_in_chimera(filename)
      self.add_map(filename)

    def upload_file(self,state):
      
      # get file name
      filename, _ = QFileDialog.getOpenFileName(self, "Select File")
      self.filename_label.setText(filename)
      
      
      # load in chimera
      self.load_file_in_chimera(filename)
      
      self.add_model(filename)

      # update the combo box
      self.update_model_list()
      model_name = Path(filename).name
      self.set_active_model(model_name=model_name)

    def viewer_start(self):
      if not self.viewer_status_alive:
        self.viewer = ChimeraXViewer()
        self.viewer.start_viewer(json_response=True)
        self.checkbox_chimerax.setChecked(True)
      
    def viewer_stop(self):
      #print("running viewer stop")
      if self.viewer_status_alive:
        #response = self.viewer.send_command(["remotecontrol stop","close","exit"])
        response = self.viewer.send_command(["exit"])
        self.checkbox_chimerax.setChecked(False)
      
      
    
    def checkbox_changed(self, state):
        if state == 2: # 2 means the checkbox is checked
          if not self.viewer_status_alive:
            #self.viewer_start()
            pass
          self.checkbox_chimerax.setText("Remote ChimeraX connected: "+self.viewer.url)

            
        else:
          self.checkbox_chimerax.setText('Remote ChimeraX connected')
          if self.viewer_status_alive:
            #self.viewer_stop()
            for widget in self.output_widgets:
              widget.setText("")
    
    def response_error_rest(self,response):
      # returns True if the reponse failed
      if response.status_code == 200:
        return False
      else:
        return True

    def response_error_chimera(self,response):
      # returns True if Chimera logs a problem
      if (response.json()["error"] == None):
        return False
      else:
        return True



    def view_selection(self,phenix_selection_string,return_chimerax_string=False):
        # rename this to just "selection", decouple focusing 

        # Step 1. Validate
        validated,composition_phenix,remove_through_ok, has_partial_occ = self.validate_selection(phenix_selection_string)

        if not validated:
          self.output_notes.setText("Error: Malformed Phenix selection")

        elif not remove_through_ok:
          self.output_notes.setText("Warning!: the keyword 'through' and this selection are incompatible with visualization.")
        elif has_partial_occ and self.strict_selections:
          self.output_notes.setText("Warning!: the selection contains partial occupancy, which cannot be selected simultaneously in Chimera")
          self.phenix_selection_string.setText(phenix_selection_string)
        elif  composition_phenix is None:
          self.output_notes.setText("Composition is None")
        else:
          self.phenix_selection_string.setText(phenix_selection_string)
          
        if composition_phenix is not None and remove_through_ok:
          n_atoms, n_chains, n_residue_groups = composition_phenix["n_atoms"], composition_phenix["n_chains"], composition_phenix["n_residue_groups"]
          
          # self.phenix_selection_result.setText(
          #   "Atoms: {n_atoms}, Chains: {n_chains}, Residue groups: {n_residue_groups}".format(n_atoms,n_chains,n_residue_groups))
          self.phenix_selection_result.setText(
            str(n_atoms)+" atoms, "+str(n_chains)+" chains, "+str(n_residue_groups)+" residue groups selected")
            
        
        # Step 2. Translate
        if validated and remove_through_ok:
          try:
            chimerax_selection = translate_phenix_selection_string(phenix_selection_string,model_number=1)
            translated = True
            self.text_output_title.setText("ChimeraX selection string")
            self.text_output.setText(chimerax_selection)
          except:
            if self.debug:
              raise
            translated = False
            self.output_notes.setText("Error: the translation failed")
        else:
          translated = False
        
        # Send to ChimeraX
        success = False
        used_workaround = False
          
        if validated and translated and self.viewer_status_alive and remove_through_ok:
          ###################################################################################################
          # Here we send a selection to Chimera
          ###################################################################################################

          chimerax_command = "sel "+chimerax_selection
          success,response = self.send_command(chimerax_command)


          if success and not used_workaround:
            
            self.text_output.setText(chimerax_command)
            notes = "".join(response.json()["log messages"]["note"])
            self.text_response.setText(notes)
            self.output_notes.setText("")
          
          if n_atoms is not None:
            try:
              n_atoms_chimera = self.read_chimera_output_atoms_selected(response)
            except:
              n_atoms_chimera = None
            if n_atoms_chimera!= None and n_atoms_chimera !=n_atoms and not has_partial_occ:
              pass
              #QMessageBox.warning(self, "Warning", "Number of selected atoms differs between programs!")
              # TODO: Re-enable warning after debugging
          
          if return_chimerax_string:
             return success, chimerax_selection
          return success
        
    def send_command(self,command):
      
      if self.debug:
        print("Sending command: "+command)
      response = self.viewer.send_command(command)
      rest_error = self.response_error_rest(response)
      chimera_error = self.response_error_chimera(response)
      if rest_error== True:
        self.output_notes.setText("Error: failure in REST communication with ChimeraX")
      elif chimera_error== True:

        if not self.checkbox_workaround.isChecked():
          self.output_notes.setText("Error: failure reported from ChimeraX")
          notes = "".join(response.json()["log messages"]["note"])
          self.text_response.setText(notes)
        else:
          success = self.phenix_selection_workaround(phenix_selection_string)
          used_workaround = True
      else:
        # Everything is ok, no errors
        success = True

      return success, response
    
    def phenix_selection_workaround(self,phenix_selection_string):      
      # round trip integer to string selection
      selection_int = convert_selection_str_to_int(self.model,phenix_selection_string)
      phenix_selection_string_modified = convert_selection_int_to_str(self.model,selection_int)
      success = self.view_selection(phenix_selection_string_modified)
      if success:
        self.output_notes.setText("Note: the translation failed initially, but was successful by reducing it to a more basic vocabulary")
      else:
        self.output_notes.setText("Error: the translation failed once initially, and also failed after reducing to a more basic vocabulary")
      return success
        
      
    def read_chimera_output_atoms_selected(self,response):
      notes = "".join(response.json()["log messages"]["note"])
      note_split = notes.split("atom")
      if "atom" in notes:
        note_split = notes.split("atom")
        assert len(note_split)==2, "Unable to determine how many atoms ChimeraX selected"
        n_atoms_chimera = int(note_split[0])
      elif "Nothing selected" in notes:
        n_atoms_chimera = 0
      else:
        assert False,  "Unable to determine how many atoms ChimeraX selected"
      return n_atoms_chimera

    
    def validate_selection(self,phenix_selection_string):
      if self.model == None:
        model_str = """
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N  
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C  
"""
        dm = DataManager()
        dm.process_model_str("dummy_model",model_str)
        model = dm.get_model()
        
      else:
        model = self.model
        
      # make sure selection doesn't return error
      try:
        selection_int = convert_selection_str_to_int(model,phenix_selection_string)
        if self.model == None:
          composition = None
          remove_through_ok = None
          has_partial_occ = None
        else:
          model_selected = self.model.select(self.model.selection(phenix_selection_string))
          n_atoms = len(list(model_selected.get_hierarchy().atoms()))
          n_chains = len(list(model_selected.get_hierarchy().chains()))
          n_residues = len(list(model_selected.get_hierarchy().residue_groups()))
          composition = {
            "n_atoms":n_atoms,
            "n_chains":n_chains,
            "n_residue_groups":n_residues,
          }
          # check that through to colon is ok
          sel_unchanged = model.selection(phenix_selection_string)
          sel_only_colon = model.selection(phenix_selection_string.replace("through",":"))
          remove_through_ok = (sum(~(sel_unchanged==sel_only_colon))==0) and (self.model is not None)

          # check if partial occupancy
          has_partial_occ = int(sum(model.get_occ()))!=model.get_number_of_atoms()

        return True, composition, remove_through_ok, has_partial_occ
      except:
        if self.debug:
          raise
        return False,None, None, None
 
    
    def button_clicked(self):
        phenix_selection_string = self.text_input.text()
        success = self.view_selection(phenix_selection_string)
      

  
    def on_combobox_models_changed(self,index):
      if index == -1:
        index = self.combo_box_models.count()
      self.set_active_model(index=index)
      
        
    def focus_selection(self):
      if self.viewer_status_alive:
        response = self.viewer.send_command([f"view sel"])
        response = self.viewer.send_command([f"show sel atoms"])
        response = self.viewer.send_command([f"hide sel ribbon"])

    def focus_single_residue(self,sel_residue_cmd):
        # Verify a single residue is selected
        response = self.viewer.send_command(f"sel {sel_residue_cmd}") # get a response about the composition of the selection
        composition = chimerax_composition_from_selection_response(response)
        assert composition["residue"]==1, "Selection did not return a single residue:"+sel_residue_cmd
        
        # show atoms


        # get the xyz coords of the anchor atoms N,CA,C
        cmd = f"sel {sel_residue_cmd}; show sel atoms; sel intersect @N|@CA|@C"
        response = self.viewer.send_command(cmd)
        cmd = f"getcrd sel"
        response = self.viewer.send_command(cmd)
        content =  response.json()["log messages"]["note"][0]
        content = [e for e in content.split("Atom") if len(e)>0]
        coords = {}
        for line,name in zip(content,["N","CA","C"]):
            coord = np.array([float(s) for s in line.split()[-3:]])
            coords[name] = coord
            
        # Get a rotation for the residue
        R_residue = rotation_from_abc(coords["N"],coords["CA"],coords["C"])
        
        # convert to chimerax matrix string of where to move the camera to
        matrix_string = matrix_string_from_R(R_residue,center=coords["CA"],offset=25)
        


        # clip and send
        clip_cmd = ";clip near -3 far 3 position sel"
        cmd = f"view matrix camera {matrix_string}; sel {sel_residue_cmd}{clip_cmd}"
        response = self.viewer.send_command(cmd)
        return response

    def focus_model(self,model_name):
      spec = self._model_name_to_spec[model_name]
      if self.viewer_status_alive:
        response = self.viewer.send_command([f"view {spec}"])

    def set_active_model(self,model_name=None,index=None):
      #print("Running set_active_model()",model_name,index)
      if model_name is not None:
        assert model_name in self.model_names, f"Chimerax model name not found: {str([model_name])} of {str(self.model_names)}"
        found = False
        for i in range(self.combo_box_models.count()):
          
          current_text = self.combo_box_models.itemText(i)
          spec,filepath_name = None,None
          if "," in current_text:
            spec,filepath_name = current_text.split(",")
            assert filepath_name in self.model_names, f"Chimerax model name not found: {str([filepath_name])} of {str(self.model_names)}"
            if filepath_name == model_name:
              index = i
              found = True
              break
        if not found:
          assert False, f"Model name not found: {str([model_name])} of {str(self.model_names)}"
      elif index is not None:
        pass
      else:
        assert False, "Provide either model_name or index"

      # set new model
      self.combo_box_models.setCurrentIndex(index)
      current_text = self.combo_box_models.itemText(index)
      if "," in current_text:
        spec,model_name = current_text.split(",")
        assert model_name in self.model_names, f"Chimerax model name not found: {str([model_name])} of {str(self.model_names)}"
        self.current_model_name = model_name
      else:
        #print("set_active_model() function was called, but setting a model failed")
        #print("model_name",model_name,"index",index)
       # print(current_text)
        pass
        

      # focus
      self.focus_model(self.current_model_name)
        
      # load atom data
      df= self.dfs[self.current_model_name]
      model = PandasModel(df)
      self.tab2_content.setModel(model)

      # residue data
      df_comp = group_by_columns(df,["asym_id","seq_id","comp_id"])
      
      model_comp = PandasModel(df_comp,suppress_columns=["Atom Index"])
      self.tab_comp_content.setModel(model_comp)

      # load restraint data
      for restraint_type in self.restraint_types:
        restraint_dict = self.dfs_restraints[self.current_model_name]
        if restraint_type in restraint_dict:
          df = restraint_dict[restraint_type]
          #df.drop(columns=["labels"],inplace=True)
          model = PandasModel(df,suppress_columns=["i","j","k","l","labels","slack","weight","sym_op_j","rt_mx"])
          content = getattr(self,f"tab_content_{restraint_type}")
          content.setModel(model)

      # load validation data
      for validation_type in self.validation_types:
        validation_dict = self.dfs_validation[self.current_model_name]
        if validation_type in validation_dict:
          df = validation_dict[validation_type]
          #df.drop(columns=["labels"],inplace=True)
          model = PandasModel(df,suppress_columns=["Residue Type","Integer Selection"])
          content = getattr(self,f"tab_content_{validation_type}")
          content.setModel(model)

      # load hotspots data
      params = phil.parse(HotSpots.master_phil_str).extract()
      hotspots = call_program(program_class=HotSpots,data_manager=self.data_manager,params=params)
      df_hotspots = pd.read_json(hotspots)
      df_hotspots = df_hotspots.reset_index()
      model_hotspots = PandasModel(df_hotspots,suppress_columns=[
        "index","Cluster Label","x","y","z","Oversample Count","Sel Str Phenix","Sel Str Chimera"
      ])
      self.tab_hot_content.setModel(model_hotspots) 

          


    def update_model_list(self):
      if self.current_model_name is None:
        return None
      # first poll chimerax to get the model specs
      if self.viewer_status_alive:
        
        response= self.viewer.send_command("info models")
        model_list = response.json()["json values"]


        for entry in model_list:
          models = json.loads(entry)
          for model in models:

            if model["class"]=="AtomicStructure":

              spec = model["spec"]
              name_value = model["value"]
              if name_value in self.model_names:
                self._model_name_to_spec[name_value] = spec
                
      # update the dropdown menu
      combo_list = []  
      for name_value in self.model_names:
        spec = self._model_name_to_spec[name_value]
        combo_list.append((spec,name_value))

      if len(combo_list)==0:
        combo_list.append(("",""))
      
      self.combo_box_models.clear()
      for spec,name in combo_list:
        if len(spec)>0:
          self.combo_box_models.addItem(str(spec)+","+str(name))
          
      self.set_active_model(model_name=self.current_model_name)

    def periodic_update(self):
      pass
      
if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = SelectionViewerWindow()
  window.resize(800, 600) 
  window.show()

  app.exec_()
  window.viewer_stop()