
import pandas as pd

from PySide2.QtCore import QAbstractTableModel, Qt


"""
Model for Pandas Table generally
"""
class PandasTableModel(QAbstractTableModel):
  def __init__(self, df=pd.DataFrame(), parent=None, suppress_columns=[]):
    super().__init__(parent=parent)
    self._df = df
    self.suppress_columns = set(suppress_columns)
    self.visible_columns = [col for col in self._df.columns if col not in self.suppress_columns]
    self.col_map = {i: list(self._df.columns).index(col) for i, col in enumerate(self.visible_columns)}
    self.sort_order = Qt.AscendingOrder
    self.editable_rows = set(range(len(df)))  # All rows are editable by default

  def __len__(self):
    return self.df.shape[0]

  @property
  def df(self):
    return self._df

  def sort(self, column, order=Qt.AscendingOrder):
    if len(self) > 0:
      col_name = self.visible_columns[column]
      self.layoutAboutToBeChanged.emit()
      self._df.sort_values(by=col_name, ascending=(order == Qt.AscendingOrder), inplace=True)
      self.sort_order = Qt.DescendingOrder if self.sort_order == Qt.AscendingOrder else Qt.AscendingOrder
      self.layoutChanged.emit()

  def headerData(self, section, orientation, role=Qt.DisplayRole):
    if len(self) > 0:
      if role != Qt.DisplayRole:
        return None
      if orientation == Qt.Horizontal:
        return self.visible_columns[section]
      else:
        return str(self._df.index[section])

  def data(self, index, role=Qt.DisplayRole):
    if not index.isValid():
      return None

    if role == Qt.DisplayRole:
      value = self.df.iloc[index.row(), self.col_map[index.column()]]
      if isinstance(value, list):
        return "" if pd.isna(value).any() else str(value)

      return "" if pd.isna(value) else str(value)
    elif role == Qt.BackgroundRole:
      pass
    return None

  def flags(self, index):
    flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable
    if index.row() in self.editable_rows:
      flags |= Qt.ItemIsEditable
    return flags

  def setData(self, index, value, role):
    if not index.isValid() or role != Qt.EditRole:
        return False

    if value == "":
        value = self._df.iloc[index.row(), self.col_map[index.column()]]  # Keep original value

    self._df.loc[index.row(), self.visible_columns[index.column()]] = value
    self.dataChanged.emit(index, index, (Qt.DisplayRole,))

    if "action" not in self.df.columns:
        self.df["action"] = pd.NA
    self._df.loc[index.row(), "action"] = "edit"

    return True
  def rowCount(self, parent=None):
    return len(self._df)

  def columnCount(self, parent=None):
    return len(self.visible_columns)

  def toggleRowEditability(self, row):
    if row in self.editable_rows:
      self.editable_rows.remove(row)
    else:
      self.editable_rows.add(row)
    self.dataChanged.emit(self.index(row, 0), self.index(row, self.columnCount() - 1), [Qt.EditRole])
