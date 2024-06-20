
from .controller import Controller
from ..state.table import PandasTableModel
from ...core.parameters import params

"""
Filters are analogs of selections, but applied to non- atom-sites tables like bonds
"""

class FilterObj:
  def __init__(self,name="GenericFilter",geometry_type="any",mmcif_prefix="comp_id"):
    self.name = name
    self.geometry_type = geometry_type
    self.mmcif_prefix = mmcif_prefix

  def filter_df(self,df):
    raise NotImplementedError

  def _get_filter_cols(self,df):
    cols = [col for col in df.columns if self.mmcif_prefix in col]
    return cols

class AtomFilterObj(FilterObj):
  def __init__(self,name,geometry_type="any",mmcif_prefix="atom_id"):
    super().__init__(name,geometry_type=geometry_type,mmcif_prefix=mmcif_prefix)

  def filter_df(self,df):
    cols = self._get_filter_cols(df)
    return df[(df[cols]==self.name).all(axis=1)]

class ComponentFilterObj(FilterObj):
  def __init__(self,name,geometry_type="any",mmcif_prefix="comp_id"):
    super().__init__(name,geometry_type=geometry_type,mmcif_prefix=mmcif_prefix)

  @property
  def comp_id(self):
    return self.name

  def filter_df(self,df):
    cols = self._get_filter_cols(df)
    return df[(df[cols]==self.name).all(axis=1)]

class FilterProtein(ComponentFilterObj):
  def __init__(self):
    super().__init__(name="Protein")

  def filter_df(self,df):
    cols = self._get_filter_cols(df)
    return df[df[cols].isin(params.protein_comp_ids).all(axis=1)]

class FilterSolvent(ComponentFilterObj):
  def __init__(self):
    super().__init__(name="Solvent")

  def filter_df(self,df):
    cols = self._get_filter_cols(df)
    return df[df[cols].isin(params.solvent_comp_ids).all(axis=1)]

class FilterLigands(ComponentFilterObj):
  def __init__(self):
    super().__init__(name="Ligands")

  def filter_df(self,df):
    cols = self._get_filter_cols(df)
    return df[(~df[cols].isin(params.solvent_comp_ids) & ~df[cols].isin(params.protein_comp_ids)).all(axis=1)]

class CompositeFilter:
  def __init__(self,filters):
    self.filters = filters

  @property
  def geometry_type(self):
    types = [filter_obj.geometry_type for filter_obj in self.filters]
    assert len(set(types))==1, "Cannot mix geometry types in a composite filter"
    return types[0]

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
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._table = None
    self.view.combobox_comp.currentIndexChanged[str].connect(self.update_filter)
    self.view.combobox_filter.currentIndexChanged[str].connect(self.update_filter)
    self.state.signals.geometry_filter_from_restraint.connect(self.add_restraint_filter)

    self.comp_filters = []
    self.other_filters = [FilterAll(),FilterPercentile("5% Worst"),FilterPercentile("5% Best"),FilterLigands(),FilterSolvent(),FilterProtein()]
    self.restraint_filters = [FilterAll()]
  
  @property
  def comp_filter_dict(self):
    return {filter.name:filter for filter in self.comp_filters}


  @property
  def other_filter_dict(self):
    return {filter.name:filter for filter in self.other_filters}

  def add_restraint_filter(self,filter_obj):
    if filter_obj.geometry_type == "All" or filter_obj.geometry_type == self.parent.geometry_type:
      self.restraint_filters = [filter_obj]
      self.update_filter()

  def update_filter(self):
    df = self.parent.dataframe
    if df is None:
      return

    comp = self.current_comp
    filter_obj_comp = self.comp_filter_dict[comp]
    filter_obj_other = self.other_filter_dict[self.current_filter_other]
    filter_obj_restraint = self.restraint_filters[0]
    filter_obj = CompositeFilter([filter_obj_comp,filter_obj_other,filter_obj_restraint])
    #self.state.signals.filter_update.emit(filter_obj,False)

    df_filtered = filter_obj.filter_df(df)
    df_filtered = df_filtered.rename(columns=self.parent.rename_columns)
    suppress_cols = [c.lower() for c in self.parent.suppress_columns]
    suppress_cols+= [c.capitalize() for c in self.parent.suppress_columns]
    table_model = PandasTableModel(df_filtered,suppress_columns=suppress_cols)
    self.parent.table_model = table_model

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
    text =  self.view.combobox_comp.currentText()
    if text == "":
      text = "All"
    return text

  @property
  def current_filter_other(self):
    text =  self.view.combobox_filter.currentText()
    if text == "":
      text = "All"
    return text

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