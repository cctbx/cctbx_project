from itertools import chain

import pandas as pd

from PySide2.QtCore import QAbstractTableModel, Qt
from PySide2.QtWidgets import QApplication, QTableView, QStyledItemDelegate, QAbstractItemView, QPushButton, QVBoxLayout, QWidget
from PySide2.QtGui import QBrush, QColor

"""
Model for Pandas Table generally
"""
class PandasTableModel(QAbstractTableModel):
  def __init__(self, df=pd.DataFrame(), display_columns=[], column_display_names={}, capitalize=False, remove_underscores=False, transpose=False):
    super().__init__()
    self._df = df
    self._transpose = transpose
    
    if len(display_columns) == 0:
      display_columns = list(df.columns)
    self.display_columns = display_columns
    self.visible_columns = self.generate_visible_columns(df, self.display_columns)
    self.col_map = {i: list(self.df.columns).index(col) for i, col in enumerate(self.visible_columns)}
    self.sort_order = Qt.AscendingOrder
    self.editable_rows = set(range(len(df)))  # All rows are editable by default
    self.was_modified = False
    self.color_map = {}
    self.font_map = {}
    
    # Initialize column display names mapping
    self.column_display_names = {}
    for col in self.visible_columns:
      if col in column_display_names:
        new_col = column_display_names[col]
      else:
        new_col = col
      self.column_display_names[col] = new_col

    if capitalize:
      self.column_display_names = {col: new_col.capitalize() for col, new_col in self.column_display_names.items()}
    if remove_underscores:
      self.column_display_names = {col: new_col.replace("_", " ") for col, new_col in self.column_display_names.items()}

  def __len__(self):
    return self.df.shape[0]

  @property
  def df(self):
    return self._df

  @staticmethod
  def split_cif_key_prefix(cif_key):
    # return everything prior to the last _ separated suffix
    return "_".join(cif_key.split("_")[:-1])

  def generate_visible_columns(self, df, display_columns):
    # Takes display columns as a list of column names OR prefixes
    # Returns 'real' column names in the df that match
    visible_columns = []
    for col in df.columns:
      if col in display_columns:
        visible_columns.append(col)
      else:
        prefix = self.split_cif_key_prefix(col)
        if prefix in display_columns:
          visible_columns.append(col)
    return visible_columns

  def sort(self, column, order=Qt.AscendingOrder):
    if len(self) > 0:
      col_name = self.visible_columns[column]
      self.layoutAboutToBeChanged.emit()
      self.df.sort_values(by=col_name, ascending=(order == Qt.AscendingOrder), inplace=True)
      self.sort_order = Qt.DescendingOrder if self.sort_order == Qt.AscendingOrder else Qt.AscendingOrder
      self.layoutChanged.emit()

  def headerData(self, section, orientation, role=Qt.DisplayRole):
    if len(self) > 0:
      if role != Qt.DisplayRole:
        return None
      if self._transpose:
        if orientation == Qt.Horizontal:
          return str(self.df.index[section])
        else:
          if section < len(self.visible_columns):
            return self.column_display_names.get(self.visible_columns[section], self.visible_columns[section])
          else:
            return None
      else:
        if orientation == Qt.Horizontal:
          if section < len(self.visible_columns):
            return self.column_display_names.get(self.visible_columns[section], self.visible_columns[section])
          else:
            return None
        else:
          return str(self.df.index[section])

  def data(self, index, role=Qt.DisplayRole):
    if not index.isValid():
      return None

    if role == Qt.DisplayRole:
      if self._transpose:
        value = self.df.iloc[index.column(), self.col_map[index.row()]]
      else:
        value = self.df.iloc[index.row(), self.col_map[index.column()]]
        
      if isinstance(value, list):
        return "" if pd.isna(value).any() else str(value)

      return "" if pd.isna(value) else str(value)

  def flags(self, index):
    # Make items selectable and enabled, but not editable
    return Qt.ItemIsSelectable | Qt.ItemIsEnabled

  def setData(self, index, value, role):
    if not index.isValid() or role != Qt.EditRole:
      return False

    if value == "":
      if self._transpose:
        value = self.df.iloc[index.column(), self.col_map[index.row()]]  # Keep original value
      else:
        value = self.df.iloc[index.row(), self.col_map[index.column()]]  # Keep original value

    if self._transpose:
      self.df.iloc[index.column(), self.col_map[index.row()]] = value
    else:
      self.df.loc[index.row(), self.visible_columns[index.column()]] = value
    self.dataChanged.emit(index, index, (Qt.DisplayRole,))

    if "action" not in self.df.columns:
      self.df["action"] = pd.NA
    if self._transpose:
      self.df.iloc[index.column(), self.df.columns.get_loc("action")] = "edit"
    else:
      self.df.loc[index.row(), "action"] = "edit"

    return True

  def rowCount(self, parent=None):
    return len(self.visible_columns) if self._transpose else len(self.df)

  def columnCount(self, parent=None):
    return len(self.df) if self._transpose else len(self.visible_columns)

  def toggleRowEditability(self, row):
    if row in self.editable_rows:
      self.editable_rows.remove(row)
    else:
      self.editable_rows.add(row)
    self.dataChanged.emit(self.index(row, 0), self.index(row, self.columnCount() - 1), [Qt.EditRole])



  def set_column_display_names(self, column_display_names):
    """Update the column display names mapping."""
    self.column_display_names.update(column_display_names)
    self.headerDataChanged.emit(Qt.Horizontal, 0, self.columnCount())