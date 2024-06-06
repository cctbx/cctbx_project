from .controller import Controller
from ...core.parameters import params

"""
Filters are analogs of selections, but applied to non- atom-sites tables like bonds
"""

class FilterObj:
  def __init__(self,name="GenericFilter"):
    self.name = name

  def filter_df(self,df):
    raise NotImplementedError

  def _get_filter_col(self,df):
    if "comp_id" in df:
      return "comp_id"
    elif "Res" in df:
      return "Res"


class ComponentFilterObj(FilterObj):
  def __init__(self,name):
    super().__init__(name)


  @property
  def comp_id(self):
    return self.name

  def filter_df(self,df):
    col = self._get_filter_col(df)
    return df[df[col]==self.name]

class FilterProtein(ComponentFilterObj):
  def __init__(self):
    super().__init__(name="Protein")

  def filter_df(self,df):
    col = self._get_filter_col(df)

    return df[df[col].isin(params.protein_comp_ids)]

class FilterSolvent(ComponentFilterObj):
  def __init__(self):
    super().__init__(name="Solvent")

  def filter_df(self,df):
    col = self._get_filter_col(df)
    return df[df[col].isin(params.solvent_comp_ids)]

class FilterLigands(ComponentFilterObj):
  def __init__(self):
    super().__init__(name="Ligands")

  def filter_df(self,df):
    col = self._get_filter_col(df)
    return df[(~df[col].isin(params.solvent_comp_ids) & ~df[col].isin(params.protein_comp_ids))]

class CompositeFilter:
  def __init__(self,filters):
    self.filters = filters

  def filter_df(self,df):
    for filter_obj in self.filters:
      df = filter_obj.filter_df(df)
    return df

class FilterAll(FilterObj):
  def __init__(self,name="All"):
    super().__init__(name)

  def filter_df(self,df):
    return df

class FilterPercentile(FilterObj):
  def __init__(self,name):
    super().__init__(name)

  def filter_df(self,df):
    return df
    
class TableFilterController(Controller):
  other_filters = [FilterAll(),FilterPercentile("5% Worst"),FilterPercentile("5% Best"),FilterLigands(),FilterSolvent(),FilterProtein()]
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._table = None
    self.view.combobox_comp.currentIndexChanged[str].connect(self.update_loud)
    self.view.combobox_filter.currentIndexChanged[str].connect(self.update_loud)
    self.comp_filters = []
  
  @property
  def comp_filter_dict(self):
    return {filter.name:filter for filter in self.comp_filters}


  @property
  def other_filter_dict(self):
    return {filter.name:filter for filter in self.other_filters}

  def update_loud(self):
    comp = self.current_comp
    filter_obj_comp = self.comp_filter_dict[comp]
    filter_obj_other = self.other_filter_dict[self.current_filter_other]
    filter_obj = CompositeFilter([filter_obj_comp,filter_obj_other])
    self.state.signals.filter_update.emit(filter_obj,False)


  def update_quiet(self):
    self.table

  def _init_comps(self):
    # Do all
    filter_obj = FilterAll("All")
    self.comp_filters.append(filter_obj)
    self.add_item_if_not_exists(self.view.combobox_comp,filter_obj.name)

    comps = self.state.mol.sites["comp_id"].unique()
    for comp in sorted(list(comps)):
      filter_obj = ComponentFilterObj(name=comp)
      self.comp_filters.append(filter_obj)
      self.add_item_if_not_exists(self.view.combobox_comp,filter_obj.name)

  def _init_filters(self):
    for filter_obj in self.other_filters:
      self.add_item_if_not_exists(self.view.combobox_filter,filter_obj.name)

  @property
  def current_comp(self):
    return self.view.combobox_comp.currentText()

  @property
  def current_filter_other(self):
    return self.view.combobox_filter.currentText()

  # @property
  # def current_percentile(self):
  #   return float(self.view.combobox_percentile.currentText())
  
  # @property
  # def current_percentile_direction(self):
  #   return self.view.combobox_percentile_direction.currentText()

  @property
  def table(self):
    if self._table is not None:
      return self._table
    elif self.parent.table_model is not None:
      self._table = self.parent.table_model
      self._init_comps()
      self._init_filters()
    
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