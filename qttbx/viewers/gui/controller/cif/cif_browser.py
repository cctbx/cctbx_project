from pathlib import Path
from collections import UserDict

import pandas as pd
from PySide2.QtWidgets import (
    QFileDialog,
    QLabel
)

from ...state.table import  PandasTableModel
from ....core.cif_io import write_dataframes_to_cif_file
from ..table import TableController
from ...state.cif import CifFileData
from ...state.ref import CifFileRef


def remove_values_for_key(d, key_to_remove):
  """
  Recursively remove any key-value pairs for a given key in a nested dictionary.
  """
  if isinstance(d, (dict,UserDict)):
    # First, process the current dictionary level
    keys_to_delete = [key for key in d if key == key_to_remove]
    for key in keys_to_delete:
      del d[key]
    
    # Then, recursively process nested dictionaries and lists
    for key, value in d.items():
      remove_values_for_key(value, key_to_remove)
  elif isinstance(d, list):
    for item in d:
      remove_values_for_key(item, key_to_remove)

def rename_key_in_nested_dict(d, old_key, new_key):
  """
  Recursively rename a key in a nested dictionary.
  """
  if isinstance(d, (dict,UserDict)):
    # First, process the current dictionary level
    if old_key in d:
      d[new_key] = d.pop(old_key)
    
    # Then, recursively process nested dictionaries and lists
    for key, value in d.items():
      rename_key_in_nested_dict(value, old_key, new_key)
  elif isinstance(d, list):
    for item in d:
      rename_key_in_nested_dict(item, old_key, new_key)

class CifBrowserController(TableController):

  # Modifications to files
  rename_files = {} # More likely to use an object with getitem than dict
  suppress_files = {}

  # Modifications to data blocks and items
  rename_data_blocks = {}
  suppress_data_blocks = []

  rename_data_items = {}
  suppress_data_items = []


  # Modifications to columns
  display_columns = [] 
  column_display_names = {}
  editable_columns = [] # if not empty, enforced as exclusive
  non_editable_columns = [] # excluded regardless of presence in editable list



  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._df_dict = None
    self._table_model_dict = {}
    self._cif_ref = None
    self._table_model = PandasTableModel()



    # Signals
    self.view.save_button.clicked.connect(self.save)
    self.state.signals.ciffile_change.connect(self.update_file)
    self.view.combobox_files.currentIndexChanged.connect(lambda: self.update_file_from_label(self.view.combobox_files.currentText()))
    self.view.combobox_data_block.currentIndexChanged.connect(lambda: self.update_data_block(data_block_key=self.view.combobox_data_block.currentText()))
    self.view.combobox_data_item.currentIndexChanged.connect(lambda: self.update_data_block_item(data_item_key=self.view.combobox_data_item.currentText()))
    self.table_model.dataChanged.connect(self.on_data_changed)
    if hasattr(self.view,"load_button"):
      self.view.load_button.clicked.connect(self.showFileDialog)
    if hasattr(self.view,"edit_button"):
      #pass
      self.view.edit_button.clicked.connect(self.edit_button_clicked)
    if hasattr(self.view,"next_button"):
      self.view.next_button.clicked.connect(self.on_mouse_released)

  def showFileDialog(self):
    home_dir = Path.home()  # Cross-platform home directory
    fname = QFileDialog.getOpenFileName(self.view, 'Open file', str(home_dir))
    if fname[0]:
      filename = fname[0]
      filepath = Path(filename).absolute()

      # add to state
      data = CifFileData(filepath=filepath)
      ref = CifFileRef(data=data,show=True)
      self.state.add_ref(ref)

      # update file combobox
      self.update_file(ref)
      self.view.setComboBoxToValue(self.view.combobox_files,ref.data.filename)


  def edit_button_clicked(self):
    indices = self.view.table_view.selectedIndexes()
    if not len(indices)==1:
        return
    else:
      self.view.table_view.editRequested.emit(indices[0])

  @property
  def cif_refs(self):
    return self.state.references_ciffile 
  @property
  def filenames(self):
    # implemented to enable modification in subclasses
    filenames = [ref.data.filename for ref in self.cif_refs]
    return filenames
  @property
  def df_dict(self):
    return self._df_dict

  @df_dict.setter
  def df_dict(self,value):
    self.log("Changing df_dict")
    assert value is not None
    self._df_dict = value
    self.rename_df_dict()

    # make table models
    self._table_model_dict = self.transform_nested_dict(value,self._make_table_model)
  
  def _make_table_model(self,df):
    if isinstance(df,pd.DataFrame):
      model = PandasTableModel(df=df,display_columns=self.display_columns,column_display_names=self.column_display_names)
      return model
    else:
      return df

  def transform_nested_dict(self,d, func):
    if isinstance(d, (dict,UserDict)):
      return {k: self.transform_nested_dict(v, func) for k, v in d.items()}
    else:
      return func(d)
  @property
  def table_model_dict(self):
    # follows from df_dict but table model objects
    return self._table_model_dict

  @table_model_dict.setter
  def table_model_dict(self,value):
    self._table_model_dict = value

  @property
  def cif_ref(self):
    return self._cif_ref

  @cif_ref.setter
  def cif_ref(self,value):
    self._cif_ref = value

  def rename_df_dict(self):
    value = self.df_dict
    for old_name,new_name in self.rename_data_blocks.items():
      rename_key_in_nested_dict(value,old_name,new_name)
    for old_name,new_name in self.rename_data_items.items():
      rename_key_in_nested_dict(value,old_name,new_name)

  def rename_df_dict_undo(self):
    value = self.df_dict
    reverse_d = {v:k for k,v in self.rename_data_blocks.items()}
    for old_name,new_name in reverse_d.items():
      rename_key_in_nested_dict(value,old_name,new_name)

    reverse_d = {v:k for k,v in self.rename_data_items.items()}
    for old_name,new_name in reverse_d.items():
      rename_key_in_nested_dict(value,old_name,new_name)

  def save(self,*args):
    # Opens a save file dialog and returns the selected file path and filter

    path = Path(self.cif_ref.data.filepath)
    suggested_path = str(Path(path.parent,path.stem+"_edited"+"".join(path.suffixes)))
    filepath, _ = QFileDialog.getSaveFileName(self.view, "Save File", suggested_path, "All Files (*);")

    self.rename_df_dict_undo()
    if filepath:
        self.log(f"File selected for saving: {filepath}")
        write_dataframes_to_cif_file(self.df_dict,str(Path(filepath)))




  def on_data_changed(self,*args):
    self.log("** Data has been manually changed **")
    notification = QLabel("** Data has been manually changed **")
    self.view.layout.insertWidget(0,notification)

  def update_file_from_label(self,label):
    refs = [ref for ref in self.cif_refs if ref.data.filename==label]
    if len(refs)>0:
      self.update_file(refs[0])

  def update_file(self,ref):
    self.log("updating file: ",ref)
    self.view.combobox_files.blockSignals(True)

    if ref and ref != self.cif_ref:
      # load existing cif files
      current_filename = self.view.combobox_files.currentText()


      self.view.combobox_files.clear()
    

      self.view.combobox_files.addItems(self.filenames)
      if current_filename in self.filenames:
        self.view.setComboBoxToValue(self.view.combobox_files,current_filename)

      self.cif_ref = ref
      df_dict = ref.data.dataframes
      self.df_dict = df_dict
      current_data_block_key = self.view.combobox_data_block.currentText()
      self.view.combobox_data_block.clear()
      data_block_list = [item for item in self.df_dict.keys() if item not in self.suppress_data_blocks]
      self.view.combobox_data_block.addItems(data_block_list)
      if current_data_block_key in  data_block_list:
        self.view.setComboBoxToValue(self.view.combobox_data_block,current_data_block_key)
      data_block_key = self.view.combobox_data_block.currentText()

      # trigger update_data_block
      self.update_data_block(data_block_key=data_block_key)
      #self.view.combobox_data_block.setCurrentIndex(0)

      # switch to browser tab
      #self.parent.view.setCurrentIndex(1)

    self.view.combobox_files.blockSignals(False)

  def update_data_block(self,data_block_key=None):

    self.log("cif browser update_data_block: data_block_key="+data_block_key)
    # if data_block_key is None:
    #   data_block_key = self.view.combobox_data_block.itemText(index)
    if data_block_key not in  [None,""]:
      assert data_block_key in self.df_dict, f"{type(data_block_key)}, {data_block_key}"
      self.view.combobox_data_item.clear()
      data_block_item_list = [item for item in self.df_dict[data_block_key].keys() if item not in self.suppress_data_items]
      self.view.combobox_data_item.addItems(data_block_item_list)


      current_data_item_key = self.view.combobox_data_item.currentText()
      self.view.combobox_data_item.clear()
      self.view.combobox_data_item.addItems(data_block_item_list) # need to do twice?
      if current_data_item_key in data_block_item_list:
        self.view.setComboBoxToValue(self.view.combobox_data_item,current_data_item_key)
      data_item_key = self.view.combobox_data_item.currentText()
      self.update_data_block_item(data_block_key=data_block_key,data_item_key=data_item_key)



  def update_data_block_item(self,data_block_key=None,data_item_key=None):
    self.log(f"cif browser update block: data_block_key={data_block_key} data_item_key={data_item_key}")
    if "" not in [data_block_key,data_item_key]:
      if data_block_key is None:
        data_block_key = self.view.combobox_data_block.currentText()

      if data_item_key is None:
        data_item_key = self.view.combobox_data_item.currentText()

      assert None not in [data_block_key,data_item_key]

      data = self.df_dict[data_block_key]
      #was_modified = self.was_modified_dict[(data_block_key,data_item_key)]
      #assert data_item_key in list(data, f"{data.keys()}, {data_item_key}"
      df = data[data_item_key]
      if isinstance(df,pd.DataFrame):
        # import pdb
        # pdb.set_trace()
        #self.view.table_view = CifTableView()
        #self.view.layout.addWidget(self.view.table_view)
        #model = PandasTableModel(df=df,display_columns=self.display_columns,column_display_names=self.column_display_names)
        model = self.table_model_dict[data_block_key][data_item_key]
        #self.table_model.dataChanged.connect(self.on_data_changed)

        # Set the new model on the view (For Pyside2, not used directly), and the controller (self))
        self.table_model = model
        
        self.log("Set new table model")

