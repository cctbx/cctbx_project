from __future__ import absolute_import, division, print_function

from PySide2.QtCore import Qt, QAbstractItemModel, QModelIndex

import libtbx.phil

# =============================================================================
def check_phil(phil, scope=True, definition=True, raise_error=True):
  """
  Convenience function for checking if the input is a libtbx.phil.scope
  only or a libtbx.phil.definition only or either.

  Parameters
  ----------
    phil: object
      The object to be tested
    scope: bool
      Flag to check if phil is a libtbx.phil.scope
    definition: bool
      Flag to check if phil is a libtbx.phil.definition
    raise_error: bool
      If true, a RuntimeError is raised if the check(s) fail

  Returns
  -------
    value: bool
  """
  value = False
  if scope:  # check for only libtbx.phil.scope
    value = isinstance(phil, libtbx.phil.scope)
  if definition:  # check for only libtbx.phil.definition
    value = isinstance(phil, libtbx.phil.definition)
  if scope and definition:  # check for either
    value = isinstance(phil, libtbx.phil.scope) or isinstance(phil, libtbx.phil.definition)
  if (scope and definition) and not value and raise_error:
    raise RuntimeError('A libtbx.phil.scope or libtbx.phil.definition is expected.')
  elif scope and not value and raise_error:
    raise RuntimeError('A libtbx.phil.scope is expected.')
  elif definition and not value and raise_error:
    raise RuntimeError('A libtbx.phil.definition is expected.')
  return value

# =============================================================================
class PhilItem(object):
  """
  """

  # ---------------------------------------------------------------------------
  def __init__(self, parent=None):

    self._parent = parent
    self._children = list()
    self._type = Qt.UserRole + 1

    # PHIL information
    self.definition = None
    self.full_path = None
    self.value = None

    # widget information
    self.label_text = None
    self.tooltip = None

  # ---------------------------------------------------------------------------
  def set_phil(self, phil, scope=True, definition=True):
    check_phil(phil, scope=scope, definition=definition)

    self.definition = phil
    self.full_path = phil.full_path()
    if phil.is_definition:
      self.value = phil.extract()
    else:
      self.value = ''

    # set label text for widget
    self.label_text = ' '.join(self.definition.name.split('_')).capitalize()
    if self.definition.short_caption is not None:
      self.label_text = self.definition.short_caption

    # set tooltip text to caption or help (preferred)
    self.tooltip = self.full_path
    if self.definition.caption is not None:
      self.tooltip = self.full_path + '\n\n' + self.definition.caption
    if self.definition.help is not None:
      self.tooltip = self.full_path + '\n\n' + self.definition.help

  # ---------------------------------------------------------------------------
  def parent(self):
    return self._parent

  # ---------------------------------------------------------------------------
  def appendChild(self, item):
    self._children.append(item)

  # ---------------------------------------------------------------------------
  def childCount(self):
    return len(self._children)

  # ---------------------------------------------------------------------------
  def child(self, row):
    if row < self.childCount():
      return self._children[row]
    else:
      raise RuntimeError('There is no child at row {row}.'.format(row=row))

  # ---------------------------------------------------------------------------
  def row(self):
    if self._parent is None:
      return 0
    else:
      return self._parent._children.index(self)

  # ---------------------------------------------------------------------------
  def type(self):
    return self._type

  # ---------------------------------------------------------------------------
  def data(self, role):
    if role == Qt.DisplayRole or role == Qt.EditRole:
      return str(self.value)

  # ---------------------------------------------------------------------------
  def setData(self, value, role):
    if role == Qt.EditRole:
      self.value = value

# =============================================================================
class PhilModel(QAbstractItemModel):
  """
  The model component of a PHIL scope. This is used to map widgets to
  the underlying PHIL scope and to update the data.
  """

  # ---------------------------------------------------------------------------
  def __init__(self, parent=None):
    """
    """
    QAbstractItemModel.__init__(self, parent)

    self._header_labels = ['parameter', 'value']
    self._root = PhilItem()

    # PHIL
    self._scope = None
    self._extract = None

  # ---------------------------------------------------------------------------
  def initialize_model(self, phil_scope):

    check_phil(phil_scope, scope=True, definition=False)

    self._scope = phil_scope
    self._extract = phil_scope.extract()

    # fill out model data with values from PHIL scope
    self.beginResetModel()
    self._populate_tree(self._scope, self._root)
    self.endResetModel()

  # ---------------------------------------------------------------------------
  def get_phil_extract_value(self, full_path):
    """
    Return the value given the full path

    Parameters
    ----------
      full_path: str
        The full PHIL path

    Returns
    -------
      value: object
        The value stored in full_path. Can be a list if .multiple is True.
    """
    value = self._extract
    for subpath in full_path.split('.'):
      if isinstance(value, libtbx.phil.scope_extract_list):
        break
      value = getattr(value, subpath)
    return value

  # ---------------------------------------------------------------------------
  def set_phil_extract_value(self, full_path, value):
    """
    Set the value given the full path

    Parameters
    ----------
      full_path: str
        The full PHIL path
      value: object
        object to be stored at full_path

    Returns
    -------
      Nothing
    """
    paths = full_path.split('.')
    extract = self._extract
    for subpath in paths[:-1]:
      extract = getattr(extract, subpath)
    setattr(extract, paths[-1], value)

  # ---------------------------------------------------------------------------
  def get_master_phil(self):
    """
    Function for getting the full master PHIL scope

    Parameters
    ---------
      None

    Returns
    -------
      master_phil: libtbx.phil.scope
        The original master PHIL scope
    """
    return self._scope

  # ---------------------------------------------------------------------------
  def get_working_phil(self, diff_only=True):
    """
    Function for getting the working PHIL scope

    Parameters
    ----------
      diff_only: bool
        If True, only the differences are returned

    Returns
    -------
      working_phil: libtbx.phil.scope
        The working PHIL scope based on non-default settings
    """
    working_phil = self._scope.format(python_object=self._extract)
    if diff_only:
      working_phil = self._scope.fetch_diff(working_phil)
    return working_phil

  # ---------------------------------------------------------------------------
  def get_phil_extract(self):
    """
    Function for getting the PHIL extract

    Parameters
    ----------
      None

    Returns
    -------
      extract: libtbx.phil.extract
    """
    return self._extract

  # ---------------------------------------------------------------------------
  def _populate_tree(self, leaf, branch):
    """
    Recursively creates a tree of PhilItem objects for each PHIL
    definition

    Parameters
    ----------
      leaf: libtbx.phil.scope or libtbx.phil.definition
        If leaf is a PHIL definition, a PhilItem is created to represent
        the definition, otherwise, new root is created with the name.
      branch: PhilItem
        A model item created from a leaf.

    Returns
    -------
      Nothing
    """
    if leaf.is_definition:
      item = PhilItemFactory.get_phil_item(leaf, branch)
      branch.appendChild(item)
    else:
      new_root = PhilItemFactory.get_phil_item(leaf, branch)
      branch.appendChild(new_root)
      for sub_branch in leaf.objects:
        self._populate_tree(sub_branch, new_root)

  # ---------------------------------------------------------------------------
  def headerData(self, section, orientation, role):
    if role == Qt.DisplayRole:
      if orientation == Qt.Horizontal:
        return self._header_labels[section]
    else:
      return None

  # ---------------------------------------------------------------------------
  def parent(self, index):
    if not index.isValid():
      return QModelIndex()

    child = index.internalPointer()
    parent = child.parent()

    if parent == self._root:
      return QModelIndex()

    return self.createIndex(parent.row(), 0, parent)

  # ---------------------------------------------------------------------------
  def index(self, row, column, parent=QModelIndex()):
    if not self.hasIndex(row, column, parent):
      return QModelIndex()

    if not parent.isValid():
      parent_item = self._root
    else:
      parent_item = parent.internalPointer()

    child = parent_item.child(row)
    if child is not None:
      return self.createIndex(row, column, child)
    else:
      return QModelIndex()

  # ---------------------------------------------------------------------------
  def rowCount(self, parent=QModelIndex()):
    if parent.column() > 0:
      return 0

    if not parent.isValid():
      parent_item = self._root
    else:
      parent_item = parent.internalPointer()

    return parent_item.childCount()

  # ---------------------------------------------------------------------------
  def columnCount(self, parent=QModelIndex()):
    return len(self._header_labels)

  # ---------------------------------------------------------------------------
  def flags(self, index):
    flags = QAbstractItemModel.flags(self, index)

    if index.column() == 1:
      flags = Qt.ItemIsEditable | flags

    return flags

  # ---------------------------------------------------------------------------
  def data(self, index, role):
    if not index.isValid():
      return None

    item = index.internalPointer()
    value = item.data(role)
    if role == Qt.DisplayRole:
      if index.column() == 0:
        return item.label_text
      elif index.column() == 1:
        return value
    elif role == Qt.EditRole:
      if index.column() == 1:
        return value
    else:
      return None

  # ---------------------------------------------------------------------------
  def setData(self, index, value, role):
    if role == Qt.EditRole:
      if index.column() == 1:
        item = index.internalPointer()
        item.setData(value, role)
        self.dataChanged.emit(index, index, [Qt.EditRole])
        return True
    return False

# =============================================================================
class PhilItemFactory(object):

  mapping = {
    'choice': PhilItem,
  }

  @classmethod
  def get_phil_item(self, phil, parent):
    if check_phil(phil, scope=False, definition=True, raise_error=False):
      item_type = self.mapping.get(phil.type.phil_type, PhilItem)
      item = item_type(parent=parent)
      item.set_phil(phil)
    else:
      item = PhilItem(parent=parent)
      item.set_phil(phil)
    return item

# =============================================================================
