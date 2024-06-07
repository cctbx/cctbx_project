from pathlib import Path

import pandas as pd
from PySide2.QtWidgets import QFileDialog, QLabel

from ...state.table import  PandasTableModel
from ....core.pandas_utils import write_cif_file
from ..table import TableController


class CifBrowserController(TableController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._df_dict = None
    self._cif_ref = None
    self.table_model = PandasTableModel()
    #self.parent.view.toggle_tab_visible("Browser",show=True)
    self.data_has_changed = False



    # Signals
    #self.view.save_button.clicked.connect(self.save)
    self.state.signals.ciffile_change.connect(self.update_file)
    # Problem here:
    self.view.combobox_files.currentIndexChanged.connect(lambda: self.update_file_from_label(self.view.combobox_files.currentText()))
    self.view.combobox_data.currentIndexChanged.connect(lambda: self.update_data(data_key=self.view.combobox_data.currentText()))
    self.view.combobox_block.currentIndexChanged.connect(lambda: self.update_block(block_key=self.view.combobox_block.currentText()))
    self.table_model.dataChanged.connect(self.on_data_changed)


  @property
  def cif_refs(self):
    return self.state.references_ciffile 

  @property
  def df_dict(self):
    return self._df_dict

  @df_dict.setter
  def df_dict(self,value):
    print("Changing df_dict")
    assert value is not None
    key1,value1 = list(value.items())[0]
    if isinstance(value1,pd.DataFrame):
      value = {"data":value}
    self._df_dict = value

  @property
  def cif_ref(self):
    return self._cif_ref

  @cif_ref.setter
  def cif_ref(self,value):
    self._cif_ref = value

  def save(self,*args):
    # Opens a save file dialog and returns the selected file path and filter

    path = Path(self.cif_ref.data.filepath)
    suggested_path = str(Path(path.parent,path.stem+"_edited"+"".join(path.suffixes)))
    filepath, _ = QFileDialog.getSaveFileName(self.view, "Save File", suggested_path, "All Files (*);")

    if filepath:
        print(f"File selected for saving: {filepath}")

        write_cif_file(self.df_dict,str(Path(filepath)),inp_type='pandas',method='iotbx')




  def on_data_changed(self):
    notification = QLabel("** Data has been manually changed **")
    self.view.layout.insertWidget(0,notification)

  def update_file_from_label(self,label):
    refs = [ref for ref in self.cif_refs if ref.data.filename==label]
    if len(refs)>0:
      self.update_file(refs[0])

  def update_file(self,ref):
    self.view.combobox_files.blockSignals(True)

    if ref and ref != self.cif_ref:
      # load existing cif files
      current_filename = self.view.combobox_files.currentText()
      self.view.combobox_files.clear()
      
      filenames = [ref.data.filename for ref in self.cif_refs]
      self.view.combobox_files.addItems(filenames)
      if current_filename in filenames:
        self.view.setComboBoxToValue(self.view.combobox_files,current_filename)

      self.cif_ref = ref
      df_dict = ref.data.dataframes
      self.df_dict = df_dict
      current_data_key = self.view.combobox_data.currentText()
      self.view.combobox_data.clear()
      keys = list(self.df_dict.keys())
      self.view.combobox_data.addItems(keys)
      if current_data_key in keys:
        self.view.setComboBoxToValue(self.view.combobox_data,current_data_key)
      data_key = self.view.combobox_data.currentText()

      # trigger update_data
      self.update_data(data_key=data_key)
      #self.view.combobox_data.setCurrentIndex(0)

      # switch to browser tab
      #self.parent.view.setCurrentIndex(1)

    self.view.combobox_files.blockSignals(False)

  def update_data(self,data_key=None):

    print("cif browser update_data: data_key="+data_key)
    # if data_key is None:
    #   data_key = self.view.combobox_data.itemText(index)
    if data_key not in  [None,""]:
      assert data_key in self.df_dict, f"{type(data_key)}, {data_key}"
      self.view.combobox_block.clear()
      self.view.combobox_block.addItems(list(self.df_dict[data_key].keys()))


      current_block_key = self.view.combobox_block.currentText()
      keys = list(self.df_dict[data_key].keys())
      self.view.combobox_block.clear()
      self.view.combobox_block.addItems(keys)
      if current_block_key in keys:
        self.view.setComboBoxToValue(self.view.combobox_block,current_block_key)
      block_key = self.view.combobox_block.currentText()
      self.update_block(data_key=data_key,block_key=block_key)



  def update_block(self,data_key=None,block_key=None):
    print(f"cif browser update block: data_key={data_key} block_key={block_key}")
    if "" not in [data_key,block_key]:
      if data_key is None:
        data_key = self.view.combobox_data.currentText()

      if block_key is None:
        block_key = self.view.combobox_block.currentText()

      assert None not in [data_key,block_key]

      data = self.df_dict[data_key]
      #assert block_key in list(data, f"{data.keys()}, {block_key}"
      df = data[block_key]
      if isinstance(df,pd.DataFrame):
        # import pdb
        # pdb.set_trace()
        #self.view.table_view = CifTableView()
        #self.view.layout.addWidget(self.view.table_view)
        model = PandasTableModel(df)
        #self.table_model.dataChanged.connect(self.on_data_changed)
        self.view.table_view.setModel(model)
        print("Set new table model")

