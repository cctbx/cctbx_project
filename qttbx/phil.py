from qttbx.qt.QtCore import Qt, QAbstractItemModel, QModelIndex, QThread, QCoreApplication
from qttbx.qt.QtGui import QBrush

from libtbx import Auto
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
def label_for_definition(definition):
  """Return a human-readable label for a PHIL definition or scope.

  Uses ``definition.short_caption`` when set; otherwise falls back to the
  parameter name with underscores replaced by spaces and the first letter
  capitalized (e.g., ``"atom_selection"`` becomes ``"Atom selection"``).

  Used uniformly by widget consumers (``PhilItem.label_text`` for tree-view
  column 0, ``PhilField`` for the form-view editor label, and
  ``RepeatableScopeWidget`` for multi-scope tab text) so all three label
  sources agree.

  Parameters
  ----------
    definition: libtbx.phil.definition or libtbx.phil.scope

  Returns
  -------
    str
  """
  if definition.short_caption is not None:
    return definition.short_caption
  return ' '.join(definition.name.split('_')).capitalize()

# =============================================================================
class PhilItem(object):
  """One node in the PhilModel tree.

  Carries the PHIL definition (or scope) plus per-instance value and
  display-only metadata (label text, tooltip). For scopes with
  ``.multiple = True`` there is one PhilItem per scope-extract instance,
  distinguished by ``instance_index``.
  """

  # ---------------------------------------------------------------------------
  def __init__(self, parent=None):
    """Initialize an empty PhilItem.

    Parameters
    ----------
      parent: PhilItem, optional
        The item's parent in the tree (None only for the artificial root).
    """
    self._parent = parent
    self._children = list()
    self._type = Qt.UserRole + 1

    # PHIL information
    self.definition = None
    self.full_path = None
    self.value = None

    # Multi-scope instance index. None for non-multiple items; an integer
    # for each scope_extract instance when the parent scope has multiple=True.
    self.instance_index = None

    # Definition-intrinsic validation result.
    # None when valid (or unchecked); a non-empty str describes why invalid.
    # Computed at PhilModel.setData / initialize_model time.
    self._validation_error = None

    # widget information
    self.label_text = None
    self.tooltip = None

  # ---------------------------------------------------------------------------
  def set_phil(self, phil, scope=True, definition=True):
    """Bind this item to a PHIL definition or scope.

    Stores the definition/scope object, its full dotted path, an initial
    value (extract for definitions, ``''`` for scopes), a display label
    and a tooltip (using ``short_caption`` / ``caption`` / ``help`` when
    available, falling back to the full path).

    Parameters
    ----------
      phil: libtbx.phil.scope or libtbx.phil.definition
      scope: bool
        Forwarded to :func:`check_phil`.
      definition: bool
        Forwarded to :func:`check_phil`.
    """
    check_phil(phil, scope=scope, definition=definition)

    self.definition = phil
    self.full_path = phil.full_path()
    if phil.is_definition:
      self.value = phil.extract()
    else:
      self.value = ''

    # set label text for widget
    self.label_text = label_for_definition(self.definition)

    # set tooltip text to caption or help (preferred)
    self.tooltip = self.full_path
    if self.definition.caption is not None:
      self.tooltip = self.full_path + '\n\n' + self.definition.caption
    if self.definition.help is not None:
      self.tooltip = self.full_path + '\n\n' + self.definition.help

  # ---------------------------------------------------------------------------
  def parent(self):
    """Return the parent PhilItem, or None for the root.

    Returns
    -------
      PhilItem or None
    """
    return self._parent

  # ---------------------------------------------------------------------------
  def appendChild(self, item):
    """Append ``item`` to the end of this item's children.

    Parameters
    ----------
      item: PhilItem
    """
    self._children.append(item)

  # ---------------------------------------------------------------------------
  def childCount(self):
    """Return the number of children.

    Returns
    -------
      int
    """
    return len(self._children)

  # ---------------------------------------------------------------------------
  def child(self, row):
    """Return the child at the given row.

    Parameters
    ----------
      row: int

    Returns
    -------
      PhilItem

    Raises
    ------
      RuntimeError
        If ``row`` is out of range.
    """
    if row < self.childCount():
      return self._children[row]
    else:
      raise RuntimeError('There is no child at row {row}.'.format(row=row))

  # ---------------------------------------------------------------------------
  def row(self):
    """Return this item's index among its parent's children, or 0 for the root.

    Returns
    -------
      int
    """
    if self._parent is None:
      return 0
    else:
      return self._parent._children.index(self)

  # ---------------------------------------------------------------------------
  def type(self):
    """Return the Qt user-role type tag for this item.

    Returns
    -------
      int
        ``Qt.UserRole + 1``.
    """
    return self._type

  # ---------------------------------------------------------------------------
  def data(self, role):
    """Return data for the given Qt role.

    Parameters
    ----------
      role: int
        ``Qt.DisplayRole`` returns ``str(self.value)`` (for cell rendering);
        ``Qt.EditRole`` returns the native Python value (for editor seeding).

    Returns
    -------
      object or None
    """
    if role == Qt.DisplayRole:
      return str(self.value)
    if role == Qt.EditRole:
      return self.value

  # ---------------------------------------------------------------------------
  def setData(self, value, role):
    """Store ``value`` when ``role`` is ``Qt.EditRole``.

    Parameters
    ----------
      value: object
      role: int
    """
    if role == Qt.EditRole:
      self.value = value

  # ---------------------------------------------------------------------------
  def extract_path(self):
    """Return the indexed path to this item's value in ``PhilModel._extract``.

    Returns
    -------
      path: list
        A list of segments where each is either a ``str`` (regular attribute
        access on a scope_extract) or a ``(name: str, index: int)`` tuple
        (access into a ``scope_extract_list`` element). Walk this with
        ``PhilModel.get_phil_extract_value`` / ``set_phil_extract_value``
        (Task 2 updates them to accept this form).

        Returns ``None`` for the root and the master-scope item (they are
        structural, not value-bearing).
    """
    # Build by walking up to (but not including) the master-scope item.
    # _root is the artificial root holding the master scope at child(0);
    # the master scope itself has the empty name '' and is structural.
    segments = []
    walker = self
    while walker is not None and walker.definition is not None:
      name = walker.definition.name
      if walker.instance_index is not None:
        segments.append((name, walker.instance_index))
      else:
        segments.append(name)
      walker = walker.parent()
    if not segments:
      return None
    segments.reverse()
    # Drop the leading master-scope segment (it has no name in the extract).
    # The master scope's PhilItem has definition.name == '' (top-level
    # libtbx.phil.parse() yields a scope with empty name); other items
    # always have a real name.
    if segments and segments[0] == '':
      segments = segments[1:]
    if not segments:
      return None
    return segments

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
    """Bind the model to a master PHIL scope.

    Stores the master scope and a fresh extract, then rebuilds the
    internal PhilItem tree (paired with ``beginResetModel`` /
    ``endResetModel`` so attached views refresh).

    Must be called once before any other operation on this model.

    Parameters
    ----------
      phil_scope: libtbx.phil.scope
        The master scope (template). Required to be a scope, not a
        definition; :func:`check_phil` enforces this.
    """
    check_phil(phil_scope, scope=True, definition=False)

    self._scope = phil_scope
    self._extract = phil_scope.extract()

    # fill out model data with values from PHIL scope
    self.beginResetModel()
    self._populate_tree(self._scope, self._root)
    self.endResetModel()

  # ---------------------------------------------------------------------------
  def restore_defaults(self):
    """Reset the working extract to the master scope's defaults.

    Equivalent to ``self.initialize_model(self.get_master_phil())``: pairs
    with ``beginResetModel`` / ``endResetModel`` so attached views refresh,
    and rebuilds the multi-instance tree to the master's default instance
    counts.

    Useful for a "Reset to defaults" button in a consumer dialog.
    """
    self.initialize_model(self._scope)

  # ---------------------------------------------------------------------------
  def get_phil_extract_value(self, full_path):
    """Return the value at ``full_path``.

    Parameters
    ----------
      full_path: str or list
        A dotted path (e.g. ``"refinement.macro_cycles"``) for the legacy
        single-instance case, OR an indexed path (list of str and
        ``(name, index)`` tuples) for paths that traverse a ``.multiple = True``
        scope.

    Returns
    -------
      value: object
        The value stored at ``full_path``. Can be a list if ``.multiple`` is True
        and the path doesn't index into the list.
    """
    if hasattr(full_path, "split"):
      # Legacy dotted-string form.
      value = self._extract
      for subpath in full_path.split("."):
        if isinstance(value, libtbx.phil.scope_extract_list):
          break
        value = getattr(value, subpath)
      return value
    # Indexed-path form.
    parent, last = self._walk_extract(full_path)
    if isinstance(last, tuple):
      name, index = last
      return getattr(parent, name)[index]
    return getattr(parent, last)

  # ---------------------------------------------------------------------------
  def set_phil_extract_value(self, full_path, value):
    """Set the value at ``full_path``.

    Parameters
    ----------
      full_path: str or list
        Same forms as :py:meth:`get_phil_extract_value`.
      value: object
        Object to be stored at ``full_path``.
    """
    if hasattr(full_path, "split"):
      paths = full_path.split(".")
      extract = self._extract
      for subpath in paths[:-1]:
        extract = getattr(extract, subpath)
      setattr(extract, paths[-1], value)
      return
    parent, last = self._walk_extract(full_path)
    if isinstance(last, tuple):
      name, index = last
      getattr(parent, name)[index] = value
    else:
      setattr(parent, last, value)

  # ---------------------------------------------------------------------------
  def _walk_extract(self, path):
    """Walk ``self._extract`` by an indexed path (list of str / (str, int)).

    Used internally by ``get_phil_extract_value`` and ``set_phil_extract_value``
    when invoked with the indexed-path form returned by ``PhilItem.extract_path()``.

    Parameters
    ----------
      path: list
        Each element is either a ``str`` (attribute access) or a ``(name, index)``
        tuple (element access into a ``scope_extract_list``).

    Returns
    -------
      parent, last_segment: (object, str or tuple)
        The container in which the addressed value lives, and the last segment
        to use for the get/set operation.
    """
    parent = self._extract
    for segment in path[:-1]:
      if isinstance(segment, tuple):
        name, index = segment
        parent = getattr(parent, name)[index]
      else:
        parent = getattr(parent, segment)
    return parent, path[-1]

  # ---------------------------------------------------------------------------
  def _refresh_intrinsic_validity(self, item):
    """Compute and cache definition-intrinsic validation on ``item``.

    Checks only constraints encoded on the PHIL definition itself
    (``value_min`` / ``value_max`` for int/float; choice membership for
    ``choice``; ``size_min`` / ``size_max`` for list-native types).
    Widget-author-level constraints (``must_exist`` on PathWidget,
    ``min_length`` on StrWidget, etc.) are NOT enforced here; consumers
    that need them must use the editor-time validation path.

    Stores the result on ``item._validation_error`` (None when valid,
    otherwise a human-readable error string).

    Parameters
    ----------
      item: PhilItem
        The tree node to validate. Non-definition items are skipped.
    """
    if item.definition is None or not item.definition.is_definition:
      item._validation_error = None
      return
    d = item.definition
    v = item.value
    if v is None or v is Auto:
      item._validation_error = None
      return
    pt = d.type.phil_type if d.type is not None else None
    err = None
    if pt in ("int", "float"):
      vmin = getattr(d.type, "value_min", None)
      vmax = getattr(d.type, "value_max", None)
      if vmin is not None and v < vmin:
        err = "value {v} below minimum {m}".format(v=v, m=vmin)
      elif vmax is not None and v > vmax:
        err = "value {v} exceeds maximum {m}".format(v=v, m=vmax)
    elif pt == "choice":
      multi = getattr(d.type, "multi", False)
      valid_choices = [w.value.lstrip("*") for w in d.words]
      if multi:
        bad = [c for c in (v or []) if c not in valid_choices]
        if bad:
          err = "unknown choice(s): {b}".format(b=bad)
      else:
        if v not in valid_choices:
          err = "unknown choice: {v!r}".format(v=v)
    elif pt in ("ints", "floats", "strings", "words"):
      n = len(v) if hasattr(v, "__len__") else 0
      smin = getattr(d.type, "size_min", None)
      smax = getattr(d.type, "size_max", None)
      if smin is not None and n < smin:
        err = "{n} items below minimum {m}".format(n=n, m=smin)
      elif smax is not None and n > smax:
        err = "{n} items exceed maximum {m}".format(n=n, m=smax)
      elif pt in ("ints", "floats"):
        vmin = getattr(d.type, "value_min", None)
        vmax = getattr(d.type, "value_max", None)
        if vmin is not None or vmax is not None:
          for i, x in enumerate(v):
            if x is None:
              continue
            if vmin is not None and x < vmin:
              err = "row {i}: value {x} below minimum {m}".format(i=i, x=x, m=vmin)
              break
            if vmax is not None and x > vmax:
              err = "row {i}: value {x} exceeds maximum {m}".format(i=i, x=x, m=vmax)
              break
    item._validation_error = err

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
    """Recursively build the PhilItem tree.

    For a definition: build one PhilItem.
    For a scope with multiple=True: build one PhilItem per instance present in
    the current ``_extract`` (instance_index = 0..N-1), each with its own
    subtree of value-bearing items addressing that instance.
    For a non-multi scope: build one PhilItem with no instance_index.

    Parameters
    ----------
      leaf: libtbx.phil.scope or libtbx.phil.definition
      branch: PhilItem
    """
    if leaf.is_definition:
      item = PhilItem(parent=branch)
      item.set_phil(leaf)
      branch.appendChild(item)
      self._refresh_intrinsic_validity(item)
      return
    # Scope. Check multi-ness.
    if getattr(leaf, "multiple", False):
      instances = self._instance_count_at(branch, leaf.name)
      for idx in range(instances):
        sub_root = PhilItem(parent=branch)
        sub_root.set_phil(leaf, scope=True, definition=False)
        sub_root.instance_index = idx
        branch.appendChild(sub_root)
        for sub_branch in leaf.objects:
          self._populate_tree(sub_branch, sub_root)
      return
    # Non-multi scope (existing behavior).
    new_root = PhilItem(parent=branch)
    new_root.set_phil(leaf, scope=True, definition=False)
    branch.appendChild(new_root)
    for sub_branch in leaf.objects:
      self._populate_tree(sub_branch, new_root)

  # ---------------------------------------------------------------------------
  def _instance_count_at(self, branch, name):
    """Return the number of scope_extract instances of ``name`` at ``branch``.

    Walks ``self._extract`` from the root down to ``branch``'s extract path,
    then accesses the named attribute (which must be a ``scope_extract_list``)
    and returns its length.

    Returns 1 if the path doesn't yet exist in the extract (i.e. a fresh
    master with one default instance per multi-scope).
    """
    path = []
    walker = branch
    while walker is not None and walker.definition is not None:
      if walker.instance_index is not None:
        path.append((walker.definition.name, walker.instance_index))
      else:
        path.append(walker.definition.name)
      walker = walker.parent()
    path.reverse()
    # Drop the master-scope segment (its name is empty in the extract).
    if path and path[0] == "":
      path = path[1:]
    ex = self._extract
    for seg in path:
      if isinstance(seg, tuple):
        n, i = seg
        ex = getattr(ex, n)[i]
      else:
        ex = getattr(ex, seg)
    attr = getattr(ex, name, None)
    if attr is None:
      return 1
    return len(attr)

  # ---------------------------------------------------------------------------
  def _item_for_path(self, full_path):
    """
    Walk the tree by dotted path, returning the matching PhilItem.

    Parameters
    ----------
      full_path: str
        The full dotted PHIL path (e.g. "refinement.macro_cycles").

    Returns
    -------
      item: PhilItem
        The model item at full_path.

    Raises
    ------
      ValueError
        If any segment of the path does not match an existing PhilItem.
    """
    # _root.child(0) is the master scope (no name); descend from there.
    item = self._root.child(0)
    for segment in full_path.split('.'):
      found = None
      for i in range(item.childCount()):
        child = item.child(i)
        if child.definition.name == segment:
          found = child
          break
      if found is None:
        raise ValueError(
          "PHIL path not found: {p}".format(p=full_path))
      item = found
    return item

  # ---------------------------------------------------------------------------
  def iter_definitions(self):
    """Iterate over every PHIL definition in the master scope.

    Walks the master ``libtbx.phil.scope`` (set by
    :py:meth:`initialize_model`) in document order, descending into nested
    scopes. Yields a ``(phil_path, definition)`` pair for each
    ``libtbx.phil.definition`` encountered; scopes are traversed but not
    themselves yielded.

    Each definition (including ones inside ``.multiple = True`` scopes)
    is yielded exactly once, addressed by its template ``full_path()``;
    per-instance addresses are not enumerated. Returns nothing if
    :py:meth:`initialize_model` has not been called.

    Yields
    ------
    (str, libtbx.phil.definition)
      The definition's dotted path (via ``definition.full_path()``) and
      the definition object itself.
    """
    if self._scope is None:
      return
    stack = [self._scope]
    while stack:
      node = stack.pop(0)
      if node.is_definition:
        yield node.full_path(), node
        continue
      # Scope: descend into its objects in document order.
      stack[:0] = list(node.objects)

  # ---------------------------------------------------------------------------
  def value_at_path(self, phil_path):
    """Return the current value at ``phil_path``, or None.

    Thin wrapper around :py:meth:`get_phil_extract_value` that adds a
    "not a definition" guard: scopes (including ``.multiple = True``
    scopes) yield ``None`` rather than a ``scope_extract`` /
    ``scope_extract_list``. Useful when callers want a scalar PHIL
    value and want to silently ignore non-leaves.

    Parameters
    ----------
    phil_path: str
      The full dotted PHIL path (e.g. ``"input_model"``,
      ``"refinement.macro_cycles"``).

    Returns
    -------
    object or None
      The scalar value of the definition (``None`` for unset PHIL
      ``None`` defaults), or ``None`` when ``phil_path`` does not exist
      or addresses a scope rather than a definition.
    """
    try:
      item = self._item_for_path(phil_path)
    except ValueError:
      return None
    if item.definition is None or not item.definition.is_definition:
      return None
    return self.get_phil_extract_value(phil_path)

  # ---------------------------------------------------------------------------
  def path_for_index(self, index):
    """Return the dotted PHIL path for ``index``, or ``None``.

    Inverse of :py:meth:`index_for_path` / :py:meth:`persistent_index_for_path`:
    looks up the :class:`PhilItem` stored on the index and returns its
    ``full_path``. Returns ``None`` when ``index`` is invalid or when the
    item has no associated definition (the artificial root or the master
    scope item).

    Used by :class:`qttbx.widgets.data_manager.DataManagerTableModel` to
    map ``dataChanged`` notifications from PhilModel back to PHIL paths
    while updating its binding caches.

    Parameters
    ----------
      index: QModelIndex
        Any index from this model, regardless of column.

    Returns
    -------
      str or None
        The PHIL ``full_path`` (e.g. ``"input_model"``,
        ``"refinement.macro_cycles"``) when ``index`` points at a real
        item, otherwise ``None``.
    """
    if not index.isValid():
      return None
    item = index.internalPointer()
    if item is None or item.definition is None:
      return None
    return item.full_path

  # ---------------------------------------------------------------------------
  def set_value_at_path(self, phil_path, value):
    """Set the value at ``phil_path`` and notify attached views.

    Thin wrapper around :py:meth:`setData` that resolves ``phil_path``
    to its value-column QModelIndex via :py:meth:`index_for_path` and
    issues a ``Qt.EditRole`` write. The underlying ``setData`` updates
    the PhilItem, propagates into ``self._extract``, refreshes
    intrinsic validity, and emits ``dataChanged`` -- so external
    callers using this helper benefit from the same plumbing as a
    delegate-driven edit.

    Parameters
    ----------
      phil_path: str
        The full dotted PHIL path (e.g. ``"input_model"``). Must not
        traverse any ``.multiple = True`` scope; use
        :py:meth:`persistent_index_for_path` + :py:meth:`setData` for
        repeated-scope paths.
      value: object
        The new Python value (matching the definition's PHIL type).

    Returns
    -------
      bool
        ``True`` when the write was accepted (mirrors
        :py:meth:`setData`'s return).
    """
    index = self.index_for_path(phil_path)
    return self.setData(index, value, Qt.EditRole)

  # ---------------------------------------------------------------------------
  def get_definition(self, full_path):
    """
    Return the libtbx.phil definition or scope at the given dotted path.

    Parameters
    ----------
      full_path: str
        The full dotted PHIL path.

    Returns
    -------
      definition: libtbx.phil.definition or libtbx.phil.scope
        The PHIL object at full_path. Includes type information
        (definition.type.phil_type, value_min, value_max, choices, ...)
        used by widgets to enforce limits.

    Raises
    ------
      ValueError
        If the path does not exist.
    """
    return self._item_for_path(full_path).definition

  # ---------------------------------------------------------------------------
  def definition_for_path(self, full_path):
    """Return the ``libtbx.phil.definition`` at ``full_path``.

    Thin alias of :py:meth:`get_definition` named so consumers reading
    binding code (e.g.
    :class:`qttbx.widgets.data_manager.DataManagerWidget`) can phrase
    their lookup in domain terms (``defn = pm.definition_for_path(p)``)
    without inheriting the broader ``get_definition`` name shape.

    Parameters
    ----------
      full_path: str
        The dotted PHIL path.

    Returns
    -------
      libtbx.phil.definition or libtbx.phil.scope
        The PHIL object at ``full_path``.

    Raises
    ------
      ValueError
        If ``full_path`` does not exist.
    """
    return self.get_definition(full_path)

  # ---------------------------------------------------------------------------
  def instances_for_path(self, full_path):
    """Return the list of current values at a ``.multiple = True`` definition.

    For a ``.multiple = True`` path/scalar definition (not a multi-*scope*),
    the PHIL extract stores a list of values. This helper returns that
    list as a plain Python list, dropping the singleton ``[None]``
    default that a freshly-extracted multi-definition holds.

    Parameters
    ----------
      full_path: str
        The dotted PHIL path of a ``.multiple = True`` definition.

    Returns
    -------
      list
        A fresh list of current instance values. Empty list when the
        extract is unset, equals ``[None]`` (the default state), or the
        path does not exist.
    """
    try:
      value = self.value_at_path(full_path)
    except Exception:
      return []
    if value is None:
      return []
    if hasattr(value, "encode"):
      # A scalar string -- caller likely passed a non-multiple path; treat
      # as a single-element list.
      return [value]
    out = [v for v in value if v is not None]
    return out

  # ---------------------------------------------------------------------------
  def append_instance(self, full_path, value):
    """Append ``value`` to a ``.multiple = True`` definition's list.

    For a ``.multiple = True`` path/scalar definition, the underlying
    PHIL extract stores a list of values. This helper appends an entry
    and propagates the update through ``set_value_at_path`` so attached
    views receive a ``dataChanged`` notification.

    A freshly-initialized multi-definition extract starts as
    ``[None]`` (PHIL's "no instances yet" placeholder); this helper
    replaces that placeholder with ``[value]`` rather than leaving a
    trailing ``None`` in the list.

    Parameters
    ----------
      full_path: str
        The dotted PHIL path of a ``.multiple = True`` definition.
      value: object
        The new instance value to append.

    Returns
    -------
      bool
        ``True`` when the write was accepted (mirrors
        :py:meth:`set_value_at_path`).
    """
    current = self.value_at_path(full_path)
    if current is None or hasattr(current, "encode"):
      new_list = [value]
    else:
      filtered = [v for v in current if v is not None]
      filtered.append(value)
      new_list = filtered
    return self.set_value_at_path(full_path, new_list)

  # ---------------------------------------------------------------------------
  def remove_instance_with_value(self, full_path, value):
    """Remove the instance equal to ``value`` from a ``.multiple = True``
    definition's list.

    Locates the first list entry equal to ``value`` and removes it,
    leaving the surrounding entries undisturbed. When the resulting
    list is empty, the extract is set to ``[None]`` to match PHIL's
    "no instances yet" placeholder convention used by
    :py:meth:`initialize_model`.

    Parameters
    ----------
      full_path: str
        The dotted PHIL path of a ``.multiple = True`` definition.
      value: object
        The value of the instance to remove (compared with ``==``).

    Returns
    -------
      bool
        ``True`` when an instance was found and removed; ``False`` when
        no matching instance exists (the extract is left unchanged).
    """
    current = self.value_at_path(full_path)
    if current is None or hasattr(current, "encode"):
      return False
    filtered = [v for v in current if v is not None]
    if value not in filtered:
      return False
    filtered.remove(value)
    if not filtered:
      new_list = [None]
    else:
      new_list = filtered
    self.set_value_at_path(full_path, new_list)
    return True

  # ---------------------------------------------------------------------------
  def index_for_path(self, full_path):
    """
    Return a QModelIndex for the value column at the given dotted path.

    Parameters
    ----------
      full_path: str
        The full dotted PHIL path.

    Returns
    -------
      index: QModelIndex
        Pointing at column 1 (the value column) of the PhilItem at
        full_path. The index's internalPointer() is the PhilItem.

    Raises
    ------
      ValueError
        If the path does not exist, or if it traverses a scope with
        .multiple = True (caller must use persistent_index_for_path() in
        that case -- added in v3).
    """
    item = self._item_for_path(full_path)
    # Walk up the ancestor scopes (skip the leaf itself: a definition with
    # ``.multiple = True`` is a list-valued leaf handled by MultipleWidget,
    # not a path-traversal hazard).
    walker = item.parent()
    while walker is not None and walker is not self._root:
      d = walker.definition
      if d is not None and getattr(d, 'multiple', False):
        raise ValueError(
          "Path traverses a repeated scope: {p}".format(p=full_path))
      walker = walker.parent()
    return self.createIndex(item.row(), 1, item)

  # ---------------------------------------------------------------------------
  def persistent_index_for_path(self, full_path, scope_indices=None):
    """Return a ``QPersistentModelIndex`` for the value column at ``full_path``.

    Unlike :py:meth:`index_for_path`, this method handles paths that traverse
    scopes with ``.multiple = True``.

    Parameters
    ----------
      full_path: str
        The full dotted PHIL path.
      scope_indices: list of int, optional
        One index per multi-scope segment encountered in the path, in path
        order. Required when ``full_path`` traverses any ``.multiple = True``
        scope; ignored otherwise.

    Returns
    -------
      qpi: QPersistentModelIndex
        Stable index that survives insertions/removals of sibling instances.

    Raises
    ------
      ValueError
        If the path doesn't exist or ``scope_indices`` doesn't match the number
        of multi-scope segments.
    """
    from qttbx.qt.QtCore import QPersistentModelIndex
    scope_indices = list(scope_indices) if scope_indices is not None else []
    # Walk down the tree by full_path segments, consuming a scope_index at
    # each multi-scope encountered.
    item = self._root.child(0)
    for segment in full_path.split("."):
      found = None
      for i in range(item.childCount()):
        child = item.child(i)
        if child.definition.name == segment:
          found = child
          break
      if found is None:
        raise ValueError("PHIL path not found: {p}".format(p=full_path))
      if found.instance_index is not None:
        # This is a multi-scope segment; consume an index.
        if not scope_indices:
          raise ValueError(
            "scope_indices required for multi-scope segment {s} in {p}".format(
              s=segment, p=full_path))
        idx = scope_indices.pop(0)
        # Re-find the right instance.
        target = None
        for i in range(item.childCount()):
          child = item.child(i)
          if (child.definition.name == segment
              and child.instance_index == idx):
            target = child
            break
        if target is None:
          raise ValueError(
            "no instance {i} of {s} in path {p}".format(
              i=idx, s=segment, p=full_path))
        item = target
      else:
        item = found
    if scope_indices:
      raise ValueError(
        "extra scope_indices for non-multi path {p}".format(p=full_path))
    return QPersistentModelIndex(self.createIndex(item.row(), 1, item))

  # ---------------------------------------------------------------------------
  def _master_scope_at_path(self, full_path):
    """Return the master ``libtbx.phil.scope`` object at ``full_path``.

    Walks ``self._scope.objects`` by name to reach the scope template.
    """
    walker = self._scope
    for segment in full_path.split("."):
      found = None
      for obj in walker.objects:
        if obj.name == segment:
          found = obj
          break
      if found is None:
        raise ValueError("PHIL path not found: {p}".format(p=full_path))
      walker = found
    return walker

  # ---------------------------------------------------------------------------
  def _parent_item_for_path(self, full_path, scope_indices=None):
    """Return the PhilItem that is the parent of the multi-scope at ``full_path``.

    Walks the tree by name segments, consuming ``scope_indices`` for any
    intermediate multi-scope. Returns the item ABOVE the final segment.
    """
    scope_indices = list(scope_indices) if scope_indices is not None else []
    item = self._root.child(0)              # master-scope item
    segments = full_path.split(".")
    for segment in segments[:-1]:
      target = None
      for i in range(item.childCount()):
        child = item.child(i)
        if child.definition.name != segment:
          continue
        if child.instance_index is not None:
          if not scope_indices:
            raise ValueError(
              "scope_indices required for multi-scope segment {s}".format(s=segment))
          if child.instance_index == scope_indices[0]:
            target = child
            scope_indices.pop(0)
            break
        else:
          target = child
          break
      if target is None:
        raise ValueError("PHIL path not found: {p}".format(p=full_path))
      item = target
    return item, segments[-1]

  # ---------------------------------------------------------------------------
  def add_scope_instance(self, full_path, scope_indices=None):
    """Add a new instance of the multi-scope at ``full_path``.

    Pairs ``beginInsertRows``/``endInsertRows`` so attached views update
    consistently.

    Parameters
    ----------
      full_path: str
        The dotted PHIL path of a scope with ``.multiple = True``.
      scope_indices: list of int, optional
        Required when intermediate scopes on the path are also ``.multiple``.

    Returns
    -------
      qpi: QPersistentModelIndex
        Persistent index of the new instance's PhilItem.

    Raises
    ------
      ValueError
        If the target path is not a ``.multiple = True`` scope.
    """
    app = QCoreApplication.instance()
    assert app is not None and QThread.currentThread() == app.thread(), \
        "PhilModel mutated from non-GUI thread"
    from qttbx.qt.QtCore import QPersistentModelIndex
    master_scope = self._master_scope_at_path(full_path)
    if not getattr(master_scope, "multiple", False):
      raise ValueError(
        "add_scope_instance requires a multi-scope path; "
        "{p} is not multi".format(p=full_path))
    # 1. Find the parent PhilItem (above the multi-scope) and the multi-scope name.
    parent_item, scope_name = self._parent_item_for_path(
      full_path, scope_indices=scope_indices)
    # 2. Compute the parent QModelIndex. Always createIndex of the parent
    # PhilItem -- this works for the master-scope item (its row() is 0) as
    # well as for any deeper non-multi parent.
    parent_qmi = self.createIndex(parent_item.row(), 0, parent_item)
    # 3. Find the parent's extract container. The master-scope item's
    # extract_path() is None (the master is structural; extract starts
    # below it), so dispatch to self._extract directly in that case.
    if parent_item is self._root.child(0):
      parent_extract = self._extract
    else:
      parent_extract = self.get_phil_extract_value(parent_item.extract_path())
    fresh = master_scope.extract()
    getattr(parent_extract, scope_name).append(fresh)
    new_index = len(getattr(parent_extract, scope_name)) - 1
    # 4. Find the row position to insert the new PhilItem.
    # Multi-scope instances are consecutive in the parent's children. Count
    # existing siblings with the same name and determine the insertion row.
    insert_row = 0
    for i in range(parent_item.childCount()):
      if parent_item.child(i).definition.name == scope_name:
        insert_row = i + 1                  # append after existing siblings
    # 5. Begin row insertion, build the subtree, end row insertion.
    self.beginInsertRows(parent_qmi, insert_row, insert_row)
    new_item = PhilItem(parent=parent_item)
    new_item.set_phil(master_scope, scope=True, definition=False)
    new_item.instance_index = new_index
    parent_item._children.insert(insert_row, new_item)
    for sub_branch in master_scope.objects:
      self._populate_tree(sub_branch, new_item)
    self.endInsertRows()
    return QPersistentModelIndex(
      self.createIndex(insert_row, 0, new_item))

  # ---------------------------------------------------------------------------
  def remove_scope_instance(self, full_path, instance_index, scope_indices=None):
    """Remove the ``instance_index``-th instance of the multi-scope at ``full_path``.

    Pairs ``beginRemoveRows``/``endRemoveRows`` so attached views and any
    outstanding ``QPersistentModelIndex`` references invalidate cleanly.

    Parameters
    ----------
      full_path: str
        The dotted PHIL path of a scope with ``.multiple = True``.
      instance_index: int
        The zero-based instance to remove.
      scope_indices: list of int, optional
        Required when intermediate scopes on the path are also ``.multiple``.
    """
    app = QCoreApplication.instance()
    assert app is not None and QThread.currentThread() == app.thread(), \
        "PhilModel mutated from non-GUI thread"
    master_scope = self._master_scope_at_path(full_path)
    if not getattr(master_scope, "multiple", False):
      raise ValueError(
        "remove_scope_instance requires a multi-scope path; "
        "{p} is not multi".format(p=full_path))
    parent_item, scope_name = self._parent_item_for_path(
      full_path, scope_indices=scope_indices)
    parent_qmi = self.createIndex(parent_item.row(), 0, parent_item)
    if parent_item is self._root.child(0):
      parent_extract = self._extract
    else:
      parent_extract = self.get_phil_extract_value(parent_item.extract_path())
    # Find the PhilItem with the matching instance_index AND its row position.
    target_row = None
    for i in range(parent_item.childCount()):
      child = parent_item.child(i)
      if (child.definition.name == scope_name
          and child.instance_index == instance_index):
        target_row = i
        break
    if target_row is None:
      raise ValueError(
        "no instance {i} of {n} at {p}".format(
          i=instance_index, n=scope_name, p=full_path))
    # Remove from extract list and from tree.
    self.beginRemoveRows(parent_qmi, target_row, target_row)
    del getattr(parent_extract, scope_name)[instance_index]
    del parent_item._children[target_row]
    # Re-number the instance_index of remaining siblings with the same name.
    next_index = 0
    for i in range(parent_item.childCount()):
      child = parent_item.child(i)
      if child.definition.name == scope_name:
        child.instance_index = next_index
        next_index += 1
    self.endRemoveRows()

  # ---------------------------------------------------------------------------
  def headerData(self, section, orientation, role):
    """QAbstractItemModel override: header label for ``parameter`` and ``value`` columns.

    Parameters
    ----------
      section: int
      orientation: Qt.Orientation
      role: int

    Returns
    -------
      str or None
        The header label for ``Qt.DisplayRole`` + horizontal orientation;
        ``None`` otherwise.
    """
    if role == Qt.DisplayRole:
      if orientation == Qt.Horizontal:
        return self._header_labels[section]
    else:
      return None

  # ---------------------------------------------------------------------------
  def parent(self, index):
    """QAbstractItemModel override: return the parent QModelIndex.

    Parameters
    ----------
      index: QModelIndex

    Returns
    -------
      QModelIndex
        Invalid (root) for top-level items; otherwise the parent's index
        at column 0.
    """
    if not index.isValid():
      return QModelIndex()

    child = index.internalPointer()
    parent = child.parent()

    if parent == self._root:
      return QModelIndex()

    return self.createIndex(parent.row(), 0, parent)

  # ---------------------------------------------------------------------------
  def index(self, row, column, parent=QModelIndex()):
    """QAbstractItemModel override: build a QModelIndex for the given row/column.

    Parameters
    ----------
      row: int
      column: int
      parent: QModelIndex
        Invalid means "child of the artificial root".

    Returns
    -------
      QModelIndex
        Invalid when ``hasIndex`` rejects, or when no child exists at
        ``row``.
    """
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
    """QAbstractItemModel override: number of children under ``parent``.

    Parameters
    ----------
      parent: QModelIndex

    Returns
    -------
      int
        0 when ``parent.column() > 0`` (children only attach at column 0);
        otherwise the parent item's ``childCount()``.
    """
    if parent.column() > 0:
      return 0

    if not parent.isValid():
      parent_item = self._root
    else:
      parent_item = parent.internalPointer()

    return parent_item.childCount()

  # ---------------------------------------------------------------------------
  def columnCount(self, parent=QModelIndex()):
    """QAbstractItemModel override: always 2 (parameter / value).

    Returns
    -------
      int
    """
    return len(self._header_labels)

  # ---------------------------------------------------------------------------
  def flags(self, index):
    """QAbstractItemModel override: mark the value column editable.

    Parameters
    ----------
      index: QModelIndex

    Returns
    -------
      Qt.ItemFlags
    """
    flags = QAbstractItemModel.flags(self, index)

    if index.column() == 1:
      flags = Qt.ItemIsEditable | flags

    return flags

  # ---------------------------------------------------------------------------
  def data(self, index, role):
    """QAbstractItemModel override: per-cell data per Qt role.

    Returns the item's ``label_text`` for column 0 + ``Qt.DisplayRole``;
    the native value for column 1 + ``Qt.DisplayRole`` or ``Qt.EditRole``.

    Surfaces definition-intrinsic validation results on column 1 via
    ``Qt.BackgroundRole`` (faint pink when invalid) and ``Qt.ToolTipRole``
    (the error message), in addition to the standard label / value roles.

    Parameters
    ----------
      index: QModelIndex
      role: int

    Returns
    -------
      object or None
    """
    if not index.isValid():
      return None

    item = index.internalPointer()
    if role == Qt.BackgroundRole and index.column() == 1:
      if item._validation_error:
        from qttbx.widgets.phil._colors import invalid_background
        return QBrush(invalid_background())
      return None
    if role == Qt.ToolTipRole and index.column() == 1:
      if item._validation_error:
        return item._validation_error
      return None
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
    """QAbstractItemModel override: write a value to the value column.

    Asserts the call is on the GUI thread (silent cross-thread mutations
    would corrupt model+view state). On a successful write to column 1
    with ``Qt.EditRole``, updates the PhilItem, propagates the value
    into the underlying PHIL extract via ``set_phil_extract_value`` (using
    the item's :py:meth:`PhilItem.extract_path`), and emits ``dataChanged``.

    Parameters
    ----------
      index: QModelIndex
      value: object
      role: int

    Returns
    -------
      bool
        True on a write that took effect; False otherwise.
    """
    app = QCoreApplication.instance()
    assert app is not None and QThread.currentThread() == app.thread(), \
        "PhilModel mutated from non-GUI thread"
    if role == Qt.EditRole:
      if index.column() == 1:
        item = index.internalPointer()
        item.setData(value, role)
        if item.definition is not None and item.definition.is_definition:
          ep = item.extract_path()
          if ep is not None:
            self.set_phil_extract_value(ep, value)
        self._refresh_intrinsic_validity(item)
        self.dataChanged.emit(
          index, index, [Qt.EditRole, Qt.BackgroundRole, Qt.ToolTipRole])
        return True
    return False
