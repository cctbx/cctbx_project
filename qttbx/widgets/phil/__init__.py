"""qttbx.widgets.phil -- PySide2 widgets for editing PHIL parameters.

This package provides three layers of an MVC-style widget system over the
existing ``qttbx.phil.PhilModel``:

1. ``PhilWidget`` (base) and its per-type subclasses -- pure view widgets
   that wrap a Qt primitive and expose a uniform value / validity contract.
   Scalar widgets (``IntWidget``, ``BoolWidget``, ``FloatWidget``,
   ``KeyWidget``, ``ChoiceWidget``) and text/list widgets (``StrWidget``,
   ``QstrWidget``, ``PathWidget``, ``StringsWidget``, ``WordsWidget``,
   ``IntsWidget``, ``FloatsWidget``, ``ChoiceMultiWidget``) are registered
   as the default for their PHIL type by ``register_builtin_widgets``;
   long-text variants (``StrTextWidget``, ``QstrTextWidget``,
   ``PathsTextWidget``, ``StringsTextWidget``, ``WordsTextWidget``,
   ``IntsTextWidget``, ``FloatsTextWidget``) are opt-in via
   ``PhilField(model, path, widget=...)``.
2. ``PhilField`` -- a controller binding a PhilWidget to a PhilModel for
   form-style consumers; performs bidirectional sync and validation
   feedback.
3. ``PhilItemDelegate`` (in ``qttbx.widgets.phil.delegate``) -- the tree-view
   counterpart, used as a QStyledItemDelegate.

The per-type widget classes are dispatched through a registry: callers use
``widget_for_definition(definition)`` and the right widget class is
instantiated based on ``definition.type.phil_type``. ``choice`` definitions
with ``definition.type.multi == True`` dispatch to ``ChoiceMultiWidget``
instead of ``ChoiceWidget``.
"""

from PySide2.QtCore import Qt, QModelIndex, QSignalBlocker, Signal
from PySide2.QtWidgets import QHBoxLayout, QLabel, QWidget


class PhilWidget(QWidget):
  """Base class for PHIL-type widgets.

  A PhilWidget is a pure view: it knows how to display and edit a value of
  a particular PHIL type, but it does not own the data. The value is
  derived from the underlying Qt primitive on demand (no separate stored
  field), so the widget state and the reported value cannot diverge.

  Subclasses override the abstract methods (``value``, ``setValue``, and
  optionally ``_parse`` / ``_format``) and place a Qt primitive
  (QLineEdit, QSpinBox, QCheckBox, QComboBox, ...) inside themselves.

  Attributes
  ----------
    valueChanged: Signal(object)
      Emitted on commit (focus loss, Enter, choice change). The payload is
      the widget's current value.
    validityChanged: Signal(bool)
      Emitted when validity transitions True <-> False.
    definition: libtbx.phil.definition
      The PHIL definition the widget edits.
  """

  # Emitted when the widget's value commits (focus loss, Enter, choice change).
  valueChanged = Signal(object)
  # Emitted when validity transitions True <-> False.
  validityChanged = Signal(bool)

  def __init__(self, definition, parent=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
        The PHIL definition this widget edits. Stored on ``self.definition``.
      parent: QWidget, optional
    """
    QWidget.__init__(self, parent)
    self.definition = definition
    self._allow_none = True
    self._use_auto = False
    self._error = ""

  # -- value contract -- subclasses must override -----------------------------

  def value(self):
    """Return the current Python value of the widget. Abstract.

    Subclasses must override.

    Raises
    ------
      NotImplementedError
        Always (this is the abstract base implementation).
    """
    raise NotImplementedError("PhilWidget subclasses must override value()")

  def setValue(self, value):
    """Set the widget to display the given value. Abstract.

    Subclasses must override and must make this idempotent: when the new
    value matches the current value (read from the underlying Qt primitive),
    return early without firing signals.

    Parameters
    ----------
      value: object

    Raises
    ------
      NotImplementedError
        Always (this is the abstract base implementation).
    """
    raise NotImplementedError("PhilWidget subclasses must override setValue()")

  # -- validation -------------------------------------------------------------

  def isValid(self):
    """Return True iff the current widget state parses to a legal value.

    Returns
    -------
      bool
    """
    return self._error == ""

  def errorString(self):
    """Return the current error message, or '' if valid.

    Returns
    -------
      str
    """
    return self._error

  def validate(self):
    """Re-run validation; emit ``validityChanged`` if state changed.

    Calls ``self.value()`` in a try/except, stores ``str(exception)`` in
    ``_error`` on failure, clears it on success, and emits
    ``validityChanged(bool)`` only on transition.
    """
    was_valid = (self._error == "")
    try:
      self.value()
      self._error = ""
    except Exception as e:
      self._error = str(e)
    is_valid = (self._error == "")
    if is_valid != was_valid:
      self.validityChanged.emit(is_valid)

  # -- Auto / None ------------------------------------------------------------

  def setAllowNone(self, enable=True):
    """Set whether None is a permitted value.

    Parameters
    ----------
      enable: bool
    """
    self._allow_none = bool(enable)

  def setUseAuto(self, enable=True):
    """Set whether libtbx.Auto is a permitted value.

    Parameters
    ----------
      enable: bool
    """
    self._use_auto = bool(enable)

  # -- enable / read-only -----------------------------------------------------

  def setReadOnly(self, ro=True):
    """Set the widget read-only.

    The base implementation is a no-op. Subclasses with text-shaped
    primitives override to forward to ``QLineEdit.setReadOnly``.

    Parameters
    ----------
      ro: bool
    """
    pass

  # -- subclass hooks ---------------------------------------------------------

  def _parse(self, text):
    """Convert displayed text/state to a Python value or raise.

    Used by text-shaped subclasses (Int/Float/Key) as the parse callable
    given to ``ValidatedLineEdit``. Non-text-shaped subclasses (Bool,
    Choice) need not override.

    Parameters
    ----------
      text: str

    Raises
    ------
      ValueError
        If the text cannot be interpreted as a value of this widget's type.
    """
    raise NotImplementedError

  def _format(self, value):
    """Convert a Python value to displayed text.

    Used by text-shaped subclasses; non-text-shaped subclasses need not
    override.

    Parameters
    ----------
      value: object

    Returns
    -------
      str
    """
    raise NotImplementedError


# ---------------------------------------------------------------------------
# Widget registry
# ---------------------------------------------------------------------------

_widget_registry = {}    # phil_type (str) -> PhilWidget subclass

# PHIL types whose values are already lists; wrapping such a definition with
# ``.multiple = True`` in MultipleWidget would mean "list of list", which is
# not a real PHIL pattern.
_LIST_NATIVE_TYPES = {"ints", "floats", "strings", "words"}


def register_widget(phil_type, widget_cls):
  """Register a PhilWidget subclass for a PHIL type name.

  Multiple calls with the same ``phil_type`` overwrite the previous entry;
  last-write-wins lets sites customize the dispatch table by importing
  their own module after the builtin registration runs.

  Parameters
  ----------
    phil_type: str
      A PHIL type name as it appears on ``definition.type.phil_type``
      (e.g. ``'int'``, ``'float'``, ``'bool'``, ``'key'``, ``'choice'``).
    widget_cls: type
      A PhilWidget subclass with a ``(definition, parent=None)``
      constructor.
  """
  _widget_registry[phil_type] = widget_cls


def widget_for_definition(definition, parent=None, fallback=None):
  """Build the right PhilWidget instance for the given PHIL definition.

  Special cases:
    - ``choice`` definitions with ``definition.type.multi == True`` dispatch
      to ``ChoiceMultiWidget`` (a non-registry-based dispatch).
    - Definitions with ``.multiple = True`` whose type is NOT in
      ``_LIST_NATIVE_TYPES`` are wrapped in ``MultipleWidget``.

  Parameters
  ----------
    definition: libtbx.phil.definition
      A PHIL definition. ``definition.type.phil_type`` selects the widget
      class via the registry.
    parent: QWidget, optional
    fallback: type, optional
      A ``PhilWidget`` subclass used when ``definition.type.phil_type`` is
      not registered. When None (the default), an unregistered type raises
      ``ValueError`` with the list of registered types and a suggestion
      to import the relevant ``qttbx.widgets.phil.<module>``.

  Returns
  -------
    widget: PhilWidget
      A fresh instance of the registered widget class.

  Raises
  ------
    ValueError
      If ``definition.type.phil_type`` is unregistered AND no ``fallback``
      was given.
  """
  phil_type = definition.type.phil_type
  if phil_type == "choice" and getattr(definition.type, "multi", False):
    from qttbx.widgets.phil.choice_widget import ChoiceMultiWidget
    return ChoiceMultiWidget(definition, parent=parent)
  if (getattr(definition, "multiple", False)
      and phil_type not in _LIST_NATIVE_TYPES):
    from qttbx.widgets.phil.multiple import MultipleWidget
    inner_cls = _widget_registry.get(phil_type, fallback)
    if inner_cls is None:
      raise ValueError(_unregistered_widget_message(phil_type))
    return MultipleWidget(definition, inner_cls=inner_cls, parent=parent)
  cls = _widget_registry.get(phil_type)
  if cls is None:
    if fallback is not None:
      return fallback(definition, parent=parent)
    raise ValueError(_unregistered_widget_message(phil_type))
  return cls(definition, parent=parent)


def _unregistered_widget_message(phil_type):
  """Return a friendly message naming the missing type and registered ones.

  Parameters
  ----------
    phil_type: str
      The PHIL type name that was not found in the registry.

  Returns
  -------
    str
  """
  registered = sorted(_widget_registry)
  return (
    "no widget registered for PHIL type {pt!r}; "
    "registered types: {r}; "
    "did you forget to import the relevant qttbx.widgets.phil.<module>?"
    .format(pt=phil_type, r=registered))


def register_builtin_widgets():
  """Register all standard PHIL-type widgets shipped in this package.

  Idempotent. Per-type modules are imported inside the function body to
  break the otherwise-circular import (each per-type module imports
  ``register_widget`` from this module). Called once at module load.
  """
  from qttbx.widgets.phil.int_widget import IntWidget
  from qttbx.widgets.phil.bool_widget import BoolWidget
  from qttbx.widgets.phil.float_widget import FloatWidget
  from qttbx.widgets.phil.key_widget import KeyWidget
  from qttbx.widgets.phil.choice_widget import ChoiceWidget
  register_widget("int", IntWidget)
  register_widget("bool", BoolWidget)
  register_widget("float", FloatWidget)
  register_widget("key", KeyWidget)
  register_widget("choice", ChoiceWidget)
  from qttbx.widgets.phil.str_widget import StrWidget
  register_widget("str", StrWidget)
  from qttbx.widgets.phil.qstr_widget import QstrWidget
  register_widget("qstr", QstrWidget)
  from qttbx.widgets.phil.path_widget import PathWidget
  register_widget("path", PathWidget)
  from qttbx.widgets.phil.strings_widget import StringsWidget
  register_widget("strings", StringsWidget)
  from qttbx.widgets.phil.words_widget import WordsWidget
  register_widget("words", WordsWidget)
  from qttbx.widgets.phil.ints_widget import IntsWidget
  register_widget("ints", IntsWidget)
  from qttbx.widgets.phil.floats_widget import FloatsWidget
  register_widget("floats", FloatsWidget)


register_builtin_widgets()


# ---------------------------------------------------------------------------
# PhilField -- form-view controller
# ---------------------------------------------------------------------------


class PhilField(QWidget):
  """Form-view controller binding a PhilWidget to a PhilModel.

  Layout
  ------
    [ label ] [ inner widget ] [ error icon (hidden when valid) ]

  Bidirectional sync
  ------------------
    inner widget.valueChanged -> model.setData(index_for_path(full_path), v)
    model.dataChanged          -> inner.setValue(...) when our index falls
                                  within the changed range AND Qt.EditRole
                                  is in the roles list

  Both directions use ``QSignalBlocker`` around the programmatic
  ``setValue`` call to avoid recursion (model<->widget signal loop).

  Attributes
  ----------
    validityChanged: Signal(bool)
      Forwarded from the inner widget; consumer dialogs can wire this to a
      QDialogButtonBox OK button's ``setEnabled`` to gate submission.
  """

  validityChanged = Signal(bool)

  def __init__(self, model, full_path=None, widget=None, parent=None,
               persistent_index=None):
    """
    Parameters
    ----------
      model: qttbx.phil.PhilModel
        The model owning the value at ``full_path`` or ``persistent_index``.
        Must outlive this field.
      full_path: str, optional
        Dotted path to the PHIL definition. Mutually exclusive with
        ``persistent_index``; one must be given.
      widget: type, optional
        Override of the widget class to instantiate. Defaults to whatever
        the registry returns for the definition's PHIL type.
      parent: QWidget, optional
      persistent_index: QPersistentModelIndex, optional
        An alternative to ``full_path`` for binding to a specific node that
        may move under sibling insertions/removals (used by
        ``RepeatableScopeWidget``). Mutually exclusive with ``full_path``.
    """
    QWidget.__init__(self, parent)
    if (full_path is None) == (persistent_index is None):
      raise ValueError(
        "PhilField requires exactly one of full_path or persistent_index")
    self._model = model
    if persistent_index is not None:
      self._persistent_index = persistent_index
      self._index = QModelIndex(persistent_index)
      self._full_path = None
      item = self._index.internalPointer()
      definition = item.definition
    else:
      self._persistent_index = None
      self._full_path = full_path
      self._index = model.index_for_path(full_path)
      definition = model.get_definition(full_path)

    # Build inner widget via registry, unless an override class is given.
    if widget is None:
      self._widget = widget_for_definition(definition, parent=self)
    else:
      self._widget = widget(definition, parent=self)

    from qttbx.phil import label_for_definition
    self._label = QLabel(
      label_for_definition(self._widget.definition), parent=self)
    base_tooltip = self._widget.definition.help or (full_path or
      self._widget.definition.full_path())
    self._label.setToolTip(base_tooltip)
    self._widget.setToolTip(base_tooltip)
    self._base_tooltip = base_tooltip

    # Error icon: a small label that is shown only when the widget is invalid.
    from qttbx.widgets.phil._colors import error_emphasis
    self._error_icon = QLabel("⚠", parent=self)   # WARNING SIGN
    self._error_icon.setStyleSheet(
      "QLabel { color: %s; font-weight: bold; }" % error_emphasis().name())
    self._error_icon.hide()

    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._label)
    layout.addWidget(self._widget, stretch=1)
    layout.addWidget(self._error_icon)

    # Initial value: model -> widget (signals blocked).
    initial = model.data(self._index, Qt.EditRole)
    blocker = QSignalBlocker(self._widget)
    self._widget.setValue(initial)
    del blocker

    # Wire bidirectional sync.
    self._widget.valueChanged.connect(self._on_widget_value_changed)
    self._widget.validityChanged.connect(self._on_validity_changed)
    self._widget.validityChanged.connect(self.validityChanged)
    self._model.dataChanged.connect(self._on_data_changed)

    # Initial visual state.
    self._on_validity_changed(self._widget.isValid())

  # -- public API -------------------------------------------------------------

  def widget(self):
    """Return the inner PhilWidget instance.

    Returns
    -------
      widget: PhilWidget
    """
    return self._widget

  def isValid(self):
    """Forward to the inner widget.

    Returns
    -------
      bool
    """
    return self._widget.isValid()

  def errorString(self):
    """Forward to the inner widget.

    Returns
    -------
      str
    """
    return self._widget.errorString()

  def commit(self):
    """Force commit of any pending edit.

    Consumers (dialogs, "Run" buttons) call this before reading the model
    so that a value the user has typed but not yet committed (focus still
    in the widget) lands in the model.

    Returns
    -------
      bool
        True if the inner widget was valid and the value was written to the
        model. False if the inner widget is currently invalid (no write
        occurred); consumers should typically refuse to proceed.
    """
    if self._widget.isValid():
      self._on_widget_value_changed(self._widget.value())
      return True
    return False

  def setEnabled(self, enable=True):
    """Enable or disable the field.

    Parameters
    ----------
      enable: bool
    """
    QWidget.setEnabled(self, enable)
    self._widget.setEnabled(enable)

  def setReadOnly(self, ro=True):
    """Forward read-only state to the inner widget.

    Parameters
    ----------
      ro: bool
    """
    self._widget.setReadOnly(bool(ro))

  # -- internal slots ---------------------------------------------------------

  def _on_widget_value_changed(self, value):
    """Slot for inner widget's valueChanged; writes to the model.

    When this field was constructed with a persistent_index, refresh the
    QModelIndex from the persistent one so we write to the right row even
    after sibling insertions/removals.

    Parameters
    ----------
      value: object
        The widget's committed value.
    """
    if self._persistent_index is not None:
      self._index = QModelIndex(self._persistent_index)
    self._model.setData(self._index, value, Qt.EditRole)

  def _on_validity_changed(self, is_valid):
    """Slot for inner widget's validityChanged; updates visual feedback.

    Toggles the error-icon visibility, updates its tooltip with the error
    message, and appends the error message to the widget's tooltip when
    invalid.

    Parameters
    ----------
      is_valid: bool
    """
    if is_valid:
      self._error_icon.hide()
      self._error_icon.setToolTip("")
      self._widget.setToolTip(self._base_tooltip)
    else:
      err = self._widget.errorString()
      self._error_icon.show()
      self._error_icon.setToolTip(err)
      self._widget.setToolTip(self._base_tooltip + "\n\n" + err)

  def _on_data_changed(self, top_left, bottom_right, roles):
    """Slot for ``model.dataChanged``; pushes external updates into the widget.

    Performs a bounding-box check so we only re-sync when our own index
    falls in the changed range and Qt.EditRole is among the changed roles
    (or no roles list was given). Programmatic ``setValue`` is wrapped in
    ``QSignalBlocker`` to break the model<->widget signal loop.

    Parameters
    ----------
      top_left: QModelIndex
      bottom_right: QModelIndex
      roles: list of int
        The roles changed; an empty list means "all roles".
    """
    if roles and Qt.EditRole not in roles:
      return
    if self._persistent_index is not None:
      self._index = QModelIndex(self._persistent_index)
    if top_left.parent() != self._index.parent():
      return
    if not (top_left.row() <= self._index.row() <= bottom_right.row()):
      return
    if not (top_left.column() <= self._index.column() <= bottom_right.column()):
      return
    new_value = self._model.data(self._index, Qt.EditRole)
    blocker = QSignalBlocker(self._widget)
    self._widget.setValue(new_value)
    del blocker
