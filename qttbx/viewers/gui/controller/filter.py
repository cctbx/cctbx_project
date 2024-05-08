from .controller import Controller

class TableFilterController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._table = None
    self.view.combobox_comp.currentIndexChanged[str].connect(self.update_loud)
    self.view.combobox_percentile.currentIndexChanged[str].connect(self.update_loud)
    self.view.combobox_percentile_direction.currentIndexChanged[str].connect(self.update_loud)


  def update_loud(self):
    self.state.signals.filter_update.emit(self)


  def update_quiet(self):
    table = self.table

  def _init_comps(self):
    comps = self.table.df["Res"].unique()
    for comp in sorted(list(comps)):
      self.add_item_if_not_exists(self.view.combobox_comp,comp)

  @property
  def current_comp(self):
    return self.view.combobox_comp.currentText()

  @property
  def current_percentile(self):
    return float(self.view.combobox_percentile.currentText())
  
  @property
  def current_percentile_direction(self):
    return self.view.combobox_percentile_direction.currentText()

  @property
  def table(self):
    if self._table is not None:
      return self._table
    elif self.parent.table_model is not None:
      self._table = self.parent.table_model
      self._init_comps()
    
    return self._table


  def add_item_if_not_exists(self,combo_box,item):
    if combo_box.findText(item) == -1:  # Item not found
      combo_box.addItem(item)
      #print(f"Added: {item}")
    
  def setComboBoxToValue(self,comboBox, value):
    index = comboBox.findText(value)
    if index >= 0:  # The value was found
      comboBox.setCurrentIndex(index)
    else:
      print(f"Value '{value}' not found in ComboBox.")