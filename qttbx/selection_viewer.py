import sys
import types
import requests
import json

from urllib.parse import urlencode
from pathlib import Path
import numpy as np

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




# A pyside2 windows which connects to ChimeraX to view Phenix selections


class SelectionViewerWindow(QMainWindow):
    todo = """

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
        self.notes = []
        self.truncation_limit = 100
        self.setWindowTitle("Phenix Selection Viewer")
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
        

        
        
        
        widget = QWidget()
        widget.setLayout(layout)

        # Set the central widget of the Window. Widget will expand
        # to take up all the space in the window by default.
        self.setCentralWidget(widget)
        
        # set a timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.periodic_update)
        self.timer.start(5000) # 5 seconds
        
    
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
      if not hasattr(self,"_data_manager"):
        self._data_manager = DataManager()
      return self._data_manager
    
    @property
    def model(self):
      if not hasattr(self,"_model"):
        self._model = None
      return self._model
    
    @model.setter
    def model(self,value):
      self._model = value
      
    @property
    def model_filename(self):
      if not hasattr(self,"_model_filename"):
        self._model_filename = None
      return self._model_filename
    
    @model_filename.setter
    def model_filename(self,value):
      self._model_filename = value
    
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
      self.model_filename = filename
      self.filename_label.setText(self.model_filename)
      
      
      # load in chimera
      self.load_file_in_chimera(self.model_filename)
      
      
      # load in cctbx/phenix
      self.data_manager.process_model_file(filename)
      model = self.data_manager.get_model()
      if model.crystal_symmetry() is None:
        from cctbx.maptbx.box import shift_and_box_model
        model= shift_and_box_model(model)
      self.model = model
      # print(model)
      # print(model_p1.crystal_symmetry())

      
      
    
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

            
    def update_model_list(self):
      combo_list = []
      if self.viewer_status_alive:
        
        response= self.viewer.send_command("info models")
        model_list = response.json()["json values"]


        for entry in model_list:
          models = json.loads(entry)
          for model in models:

            if model["class"]=="AtomicStructure":

              spec = model["spec"]
              name_value = model["value"]
              combo_list.append((spec,name_value))

      if len(combo_list)==0:
        combo_list.append(("#1",""))
      
      self.combo_box_models.clear()
      for spec,name in combo_list:
        self.combo_box_models.addItem(str(spec)+" "+str(name))
    
    
    def periodic_update(self):
      #print("Periodic update")
      # update
      self.update_model_list()
      
if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = SelectionViewerWindow()
  window.show()

  app.exec_()
  window.viewer_stop()