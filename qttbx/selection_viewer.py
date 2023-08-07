import sys
import types
import requests
import json
from PySide2.QtCore import QAbstractTableModel, Qt, Slot
from PySide2.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QTabWidget, QTableView
import pandas as pd
from urllib.parse import urlencode
from pathlib import Path
import numpy as np
from PySide2 import QtCore
from PySide2.QtCore import Qt, QTimer
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
from iotbx.data_manager import DataManager
from qttbx.viewers.chimerax import ChimeraXViewer
from cctbx.array_family import flex
from iotbx.pdb.atom_selection import selection_string_from_selection

from qttbx.sel_convert_chimera import (
  translate_phenix_selection_string,
  convert_selection_str_to_int,
  convert_selection_int_to_str
)

from cctbx_utils import model_to_df
from cctbx_utils import get_restraint_dfs_from_model

# A pyside2 windows which connects to ChimeraX to view Phenix selections


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

    def __init__(self, parent=None):
        super(FastTableView, self).__init__(parent)

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        self.mouseReleased.emit()

      
class PandasModel(QAbstractTableModel):
    def __init__(self, df=pd.DataFrame(), parent=None):
        QAbstractTableModel.__init__(self, parent=parent)
        self._df = df

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            return self._df.columns[section]
        else:
            return self._df.index[section]

    def data(self, index, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        return str(self._df.iloc[index.row(), index.column()])

    def setData(self, index, value, role):
        if not index.isValid() or role != Qt.EditRole:
            return False
        self._df.iloc[index.row(), index.column()] = value
        self.dataChanged.emit(index, index, (Qt.DisplayRole,))
        return True

    def rowCount(self, parent=None):
        return len(self._df.values)

    def columnCount(self, _, parent=None):
        return self._df.columns.size

    def flags(self, index):
        return Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable









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
        self._data_manager = None # Stores models
        self._current_model_path = None
        self._current_model_name = None
        self._current_model_spec = None
      
        self._model_name_to_path = {}
        self._model_name_to_spec = {}
        self.dfs = {}
        self.dfs_restraints = {}
        self.restraint_types = ["bond","nonbonded","angle","dihedral","chirality","planarity","parallelity"]
        self.process = True # make restraints?
        self.debug = True
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
        
        # file
        self.filename = None
        self.upload_button = QPushButton("Load File")
        self.upload_button.setFixedWidth(100) 
        self.upload_button.clicked.connect(self.upload_file)
        layout.addWidget(self.upload_button)
        
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
            
        # workaround checkbox
        self.checkbox_workaround = QCheckBox("If selection fails, try with simpler vocabulary.")
        self.checkbox_workaround.setChecked(True)
        layout.addWidget(self.checkbox_workaround)
        
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

        # Create a QTableView for the second tab.
        self.tab2_content = FastTableView()
        #self.tab2_content.clicked.connect(self.on_atom_select)
        self.tab2_content.setSelectionBehavior(QTableView.SelectRows)
        self.tab2_content.mouseReleased.connect(self.on_mouse_released_atoms)
        self.tab_widget.addTab(self.tab2_content, "Atoms")
      
        # Third tab
        self.tab_comp_content = FastTableView()
        self.tab_comp_content.setSelectionBehavior(QTableView.SelectRows)
        self.tab_comp_content.mouseReleased.connect(self.on_mouse_released_comp)
        self.tab_widget.addTab(self.tab_comp_content, "Residues")

        #restraints
        for restraint_type in self.restraint_types:
          setattr(self,f"tab_content_{restraint_type}",FastTableView())
          content = getattr(self,f"tab_content_{restraint_type}")
          content.setSelectionBehavior(QTableView.SelectRows)
          content.mouseReleased.connect(self.on_mouse_released_restraints)
          self.tab_widget.addTab(content, restraint_type.capitalize())
          
        # Set the central widget of the Window. Widget will expand
        # to take up all the space in the window by default.
        self.setCentralWidget(self.tab_widget)

        
        
        # set a timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.periodic_update)
        self.timer.start(5000) # 5 seconds



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
        success = self.view_selection(phenix_selection_string_modified)
        self.focus_selection()

      
    @Slot()
    def on_mouse_released_comp(self):
        selected = self.tab2_content.selectionModel().selection()
        deselected = QtCore.QItemSelection()
        self.on_selection_changed_comp(selected, deselected)
      
    @Slot()
    def on_mouse_released_restraints(self):
        content = self.tab_widget.currentWidget()
        selected = content.selectionModel().selection()
        deselected = QtCore.QItemSelection()
        self.on_selection_changed_restraints(selected, deselected)
    
    @Slot()
    def on_selection_changed_restraints(self, selected, deselected):
        content = self.tab_widget.currentWidget()
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

        
        selected_rows = df[iseq_columns].iloc[selected_rows]
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
        try:
          selected_rows = flatten(selected_rows)
        except:
          selected_rows = df.values.flatten()
        print(selected_rows)
        selected_rows = [int(e) for e in selected_rows]
        selection_int = flex.int(selected_rows)
        phenix_selection_string_modified = convert_selection_int_to_str(self.model,selection_int)
        success = self.view_selection(phenix_selection_string_modified)
        self.focus_selection()
      
    @property
    def viewer(self):
      if not hasattr(self,"_viewer"):
        self._viewer = ChimeraXViewer()    
      return self._viewer
    
    @viewer.setter
    def viewer(self,value):
      assert isinstance(value,ChimeraXViewer), "viewer must be instance of ChimeraXViewer"
      self._viewer = value
    
    @property
    def viewer_status_alive(self):
      #print("Checking status: ",type(self.viewer.url),self.viewer.url)
      if self.viewer.url == None:
        return False
      else:
        return self.viewer.is_alive()

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
        
        return model

    def add_model(self,filepath):
      filepath = Path(filepath)
      # load in cctbx/phenix
      self.data_manager.process_model_file(str(filepath))
      self.current_model_path = str(filepath)
      self._model_name_to_path[filepath.name] = str(filepath)
      df = model_to_df(self.model)
      self.dfs[filepath.name] = df

      if self.process:
        self.dfs_restraints[self.current_model_name] = get_restraint_dfs_from_model(self.model,lazy_build=True)

      
          
      
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
  
    def viewer_start(self):
      if not self.viewer_status_alive:
        self.viewer.start_viewer(json_response=True)
        self.checkbox_chimerax.setChecked(True)
        
        
    def viewer_stop(self):
      if self.viewer_status_alive:
        response = self.viewer.send_command(["remotecontrol stop","close","exit"])
        self.checkbox_chimerax.setChecked(False)
    
    def load_file_in_chimera(self,filename):
      if not self.viewer_status_alive:
        self.viewer_start()
      assert self.viewer_status_alive, "Cannot load file if viewer is not running"
      response = self.viewer.send_command(["open "+filename,"set bgcolor white"])
      
      # update dropdown
      self.update_model_list()
    
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

    
      

      
      
    
    def checkbox_changed(self, state):
        if state == 2: # 2 means the checkbox is checked
          if not self.viewer_status_alive:
            self.viewer_start()
          self.checkbox_chimerax.setText("Remote ChimeraX connected: "+self.viewer.url)

            
        else:
          self.checkbox_chimerax.setText('Remote ChimeraX connected')
          if self.viewer_status_alive:
            self.viewer_stop()
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


    def view_selection(self,phenix_selection_string):
        # check if malformed phenix selection
        validated,composition_phenix,remove_through_ok, has_partial_occ = self.validate_selection(phenix_selection_string)

        if not validated:
          self.output_notes.setText("Error: Malformed Phenix selection")

        elif not remove_through_ok:
          self.output_notes.setText("Warning!: the keyword 'through' and this selection are incompatible with visualization.")
        elif has_partial_occ:
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
            
        
        
        # check if translation fails
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
          # try to send translated selection to chimerax          
          
          chimerax_command = "sel "+chimerax_selection
          response = self.viewer.send_command(chimerax_command)
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
              QMessageBox.warning(self, "Warning", "Number of selected atoms differs between programs!")
          
          return success
        
    
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

    def focus_model(self,model_name):
      spec = self._model_name_to_spec[model_name]
      if self.viewer_status_alive:
        response = self.viewer.send_command([f"view {spec}"])

    def set_active_model(self,model_name=None,index=None):
      print("Running set_active_model()",model_name,index)
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
        print("set_active_model() function was called, but setting a model failed")
        print("model_name",model_name,"index",index)
        print(current_text)
        pass
        

      # focus
      self.focus_model(self.current_model_name)
        
      # load atom data
      df= self.dfs[self.current_model_name]
      model = PandasModel(df)
      self.tab2_content.setModel(model)

      df_comp = group_by_columns(df,["asym_id","seq_id","comp_id"])
      model_comp = PandasModel(df_comp)
      self.tab_comp_content.setModel(model_comp)

      # load restraint data
      for restraint_type,df in self.dfs_restraints[self.current_model_name].items():
        model = PandasModel(df)
        content = getattr(self,f"tab_content_{restraint_type}")
        content.setModel(model)


                                 
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