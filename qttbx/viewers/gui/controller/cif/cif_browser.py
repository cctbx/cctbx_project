from pathlib import Path
import pandas as pd
from PySide2 import QtCore
from PySide2.QtWidgets import  QVBoxLayout, QWidget, QMessageBox, QLabel, QFileDialog

from iotbx.pdb.mmcif import cif_input

from ...view.widgets import  PandasTable, PandasTableView
from ...state.ref import SelectionRef
from ..controller import Controller
#from ..view.tabs.cif import add_tabs
from ....last.pandas_utils import write_cif_file



class CifBrowserController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._df_dict = None
    self._cif_ref = None
    self.table_model = PandasTable()
    self.data_has_changed = False
    # Signals
    self.view.save_button.clicked.connect(self.save)
    self.state.signals.ciffile_change.connect(self.update_file)
    self.view.combobox_data.currentIndexChanged.connect(self.update_data)
    self.view.combobox_block.currentIndexChanged.connect(self.update_block)
    self.table_model.dataChanged.connect(self.on_data_changed)


  @property
  def df_dict(self):
    return self._df_dict

  @df_dict.setter
  def df_dict(self,value):
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

  def update_file(self,ref):
    if ref:
      self.cif_ref = ref
      df_dict = ref.data.dataframes
      self.df_dict = df_dict
      #self.view.combobox_data.clear()
      self.view.combobox_data.addItems(list(self.df_dict.keys()))

      # trigger update_data
      self.update_data(0,data_key=list(self.df_dict.keys())[0])
      #self.view.combobox_data.setCurrentIndex(0)

      # switch to browser tab
      self.parent.view.setCurrentIndex(1)
    else:
      # reset
      self.cif_ref = None
      self.df_dict = None

  def update_data(self,index,data_key=None):
      if data_key is None:
        data_key = self.view.combobox_data.itemText(index)
      if data_key in self.df_dict:
        self.view.combobox_block.clear()
        self.view.combobox_block.addItems(list(self.df_dict[data_key].keys()))

        # trigger update block
        block_key = list(self.df_dict[data_key].keys())[0]
        self.update_block(0,block_key=block_key)
        #self.view.combobox_block.setCurrentIndex(0)



  def update_block(self, index,data_key=None,block_key=None):
    # first delete old table
    if self.view.layout.count() > 0:
        # Get the last item in the layout
        last_item = self.view.layout.itemAt(self.view.layout.count() - 1)
        if isinstance(last_item.widget(),PandasTableView):


          # If the item is a widget, delete it
          widget = last_item.widget()
          if widget:
              widget.deleteLater()

          # Remove the item from the layout
          self.view.layout.removeItem(last_item)
        else:
          print("type of last widget: ",type(last_item.widget()))

    # add new table
    #print("called update_block with block_key: ",block_key,type(block_key))

    if data_key is None:
      data_key = self.view.combobox_data.currentText()
    if block_key is None:
      block_key = self.view.combobox_block.itemText(index)
    #print("using data_key: ",data_key,type(data_key))
    #print("using block_key: ",block_key,type(block_key))

    if "" not in [data_key,block_key]:
      data = self.df_dict[data_key]
      if block_key in data:
        df = data[block_key]
        if isinstance(df,pd.DataFrame):
          self.view.table = PandasTableView()
          self.view.layout.addWidget(self.view.table)
          self.table_model = PandasTable(df)
          self.table_model.dataChanged.connect(self.on_data_changed)
          self.view.table.setModel(self.table_model)




      # # Clear previous tabs
      # old_layout = self.tabs_container.layout()
      # if old_layout:
      #     QWidget().setLayout(old_layout)

      # # Create new tabs based on the selected top-level key
      # layout = QVBoxLayout(self.view.tabs_container)
      # selected_key = self.view.combobox.itemText(index)
      # tabs = add_tabs(self.view.tabs_container, self.d[selected_key])
      # layout.addWidget(tabs)
