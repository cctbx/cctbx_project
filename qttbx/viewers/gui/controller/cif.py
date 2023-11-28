import pandas as pd
from PySide2 import QtCore
from PySide2.QtWidgets import  QVBoxLayout, QWidget, QMessageBox
from iotbx.pdb.mmcif import cif_input

from ..view.widgets import  PandasTableModel, FastTableView
from ..state.ref import SelectionRef
from .controller import Controller
#from ..view.tabs.cif import add_tabs
from ...last.cif_io import nest_dict, convert_to_dataframes




class CifTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.d = None
    # Signals
    self.view.load_button.clicked.connect(self.save)
    self.state.signals.model_change.connect(self.update)
    self.view.combobox.currentIndexChanged.connect(self.update_tabs)
  
  def save(self):
    QMessageBox.information(self.view, 'Stopped', 'Saving modified CIF files not yet supported')


  def update(self,*args):
    
    # make nested dict
    if self.state.active_model_ref is not None:
      input = self.state.active_model_ref.model.get_model_input()
      if isinstance(input,cif_input):
        d = nest_dict(input.cif_model)
        d = convert_to_dataframes(d)
        if len(d.keys())==1:
          head_key = list(d.keys())[0]
          self.d = d[head_key]
        else:
          self.d = d

        self.view.combobox.addItems(self.d.keys())


    
  def update_tabs(self, index):
    # first delete old table
    if self.view.layout.count() > 0:
        # Get the last item in the layout
        last_item = self.view.layout.itemAt(self.view.layout.count() - 1)
        if isinstance(last_item.widget(),FastTableView):
          

          # If the item is a widget, delete it
          widget = last_item.widget()
          if widget:
              widget.deleteLater()

          # Remove the item from the layout
          self.view.layout.removeItem(last_item)
        else:
          print("type of last widget: ",type(last_item.widget()))

    # add new table
    key = self.view.combobox.itemText(index)
    if isinstance(self.d[key],pd.DataFrame):
      df = self.d[key]
      self.view.table = FastTableView()
      self.view.layout.addWidget(self.view.table)
      model = PandasTableModel(df)
      self.view.table.setModel(model)
      


      
      # # Clear previous tabs
      # old_layout = self.tabs_container.layout()
      # if old_layout:
      #     QWidget().setLayout(old_layout)

      # # Create new tabs based on the selected top-level key
      # layout = QVBoxLayout(self.view.tabs_container)
      # selected_key = self.view.combobox.itemText(index)
      # tabs = add_tabs(self.view.tabs_container, self.d[selected_key])
      # layout.addWidget(tabs)