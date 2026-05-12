"""Base text-input widget shared by PHIL line-edit-shaped widgets."""

from qttbx.qt.QtCore import Signal
from qttbx.qt.QtWidgets import QLineEdit, QPlainTextEdit


def _invalid_le_style():
  """Return the QSS for an invalid QLineEdit, derived from the active palette."""
  from qttbx.widgets.phil._colors import invalid_background, invalid_border
  return "QLineEdit { border: 1px solid %s; background: %s; }" % (
    invalid_border().name(), invalid_background().name())


def _invalid_te_style():
  """Return the QSS for an invalid QPlainTextEdit, derived from the active palette."""
  from qttbx.widgets.phil._colors import invalid_background, invalid_border
  return "QPlainTextEdit { border: 1px solid %s; background: %s; }" % (
    invalid_border().name(), invalid_background().name())


class ValidatedLineEdit(QLineEdit):
  """A QLineEdit that runs a user-supplied parse function on every keystroke.

  On every textEdited signal, the parse callable is invoked. If it raises,
  the widget paints a red border / pink background and stores the exception
  message; otherwise the widget restores its baseline stylesheet. On focus
  loss or Enter (editingFinished), valueChanged(value) is emitted iff the
  current text parses successfully. validityChanged(bool) is emitted on
  every transition between valid and invalid.

  Attributes
  ----------
    valueChanged: Signal(object)
      Emitted on commit (focus loss or Enter) when the current text parses.
    validityChanged: Signal(bool)
      Emitted on every transition between valid and invalid state.
  """

  valueChanged = Signal(object)
  validityChanged = Signal(bool)

  def __init__(self, parse, parent=None):
    """
    Parameters
    ----------
      parse: callable
        A function (str -> object) that converts the line edit's text to a
        Python value or raises Exception with a human-readable message.
      parent: QWidget, optional
        Qt parent widget.
    """
    QLineEdit.__init__(self, parent)
    self._parse = parse
    self._error = "no value"     # empty text starts invalid; consumer overrides
    self._was_valid = False
    self._normal_style = self.styleSheet()
    self.textEdited.connect(self._on_text_edited)
    self.editingFinished.connect(self._commit)

  # -- public API -------------------------------------------------------------

  def value(self):
    """Return the parsed Python value of the current text.

    Returns
    -------
      value: object
        Whatever the parse callable returns for the current text. Re-raises
        any exception the parse callable raises.
    """
    return self._parse(self.text())

  def isValid(self):
    """Return True iff the current text parses without error.

    Returns
    -------
      bool
    """
    return self._error == ""

  def errorString(self):
    """Return the current parse error message, or '' if valid.

    Returns
    -------
      str
    """
    return self._error

  def validate(self):
    """Re-run the parse callable on the current text.

    Updates the widget's stylesheet (red on invalid, baseline on valid),
    stores the error message if any, and emits validityChanged(bool) only
    when the validity state transitions.
    """
    try:
      self._parse(self.text())
      self._error = ""
    except Exception as e:
      self._error = str(e)
    is_valid = (self._error == "")
    if is_valid:
      self.setStyleSheet(self._normal_style)
    else:
      self.setStyleSheet(_invalid_le_style())
    if is_valid != self._was_valid:
      self._was_valid = is_valid
      self.validityChanged.emit(is_valid)

  # -- internal slots ---------------------------------------------------------

  def _on_text_edited(self, _text):
    """Slot for QLineEdit.textEdited; re-runs validation."""
    self.validate()

  def _commit(self):
    """Slot for QLineEdit.editingFinished; emits valueChanged when valid."""
    self.validate()
    if self.isValid():
      self.valueChanged.emit(self.value())


class ValidatedTextEdit(QPlainTextEdit):
  """A QPlainTextEdit that runs a user-supplied parse function on every keystroke.

  Sister class of ``ValidatedLineEdit`` (above) for long-text widgets.

  On every textChanged signal, the parse callable is invoked. If it raises,
  the widget paints a red border / pink background and stores the exception
  message; otherwise the widget restores its baseline stylesheet. On focus
  loss, ``valueChanged(value)`` is emitted iff the current text parses
  successfully. ``validityChanged(bool)`` is emitted on every transition
  between valid and invalid.

  ``QPlainTextEdit`` has no ``editingFinished`` signal, so commit fires on
  focus-out only (Enter inserts a newline, which is the expected behavior
  for a multi-line editor).

  Attributes
  ----------
    valueChanged: Signal(object)
      Emitted on commit (focus loss) when the current text parses.
    validityChanged: Signal(bool)
      Emitted on every transition between valid and invalid state.
  """

  valueChanged = Signal(object)
  validityChanged = Signal(bool)

  def __init__(self, parse, parent=None):
    """
    Parameters
    ----------
      parse: callable
        A function (str -> object) converting the text-edit's plainText to a
        Python value or raising Exception with a human-readable message.
      parent: QWidget, optional
    """
    QPlainTextEdit.__init__(self, parent)
    self._parse = parse
    self._error = "no value"
    self._was_valid = False
    self._normal_style = self.styleSheet()
    self.textChanged.connect(self._on_text_changed)

  # -- public API -------------------------------------------------------------

  def value(self):
    """Return the parsed Python value of the current text.

    Returns
    -------
      value: object
        Whatever ``parse`` returns for ``self.toPlainText()``. Re-raises any
        exception the parse callable raises.
    """
    return self._parse(self.toPlainText())

  def isValid(self):
    """Return True iff the current text parses without error.

    Returns
    -------
      bool
    """
    return self._error == ""

  def errorString(self):
    """Return the current parse error message, or '' if valid.

    Returns
    -------
      str
    """
    return self._error

  def validate(self):
    """Re-run the parse callable on the current text.

    Updates the widget's stylesheet (red on invalid, baseline on valid),
    stores the error message if any, and emits ``validityChanged(bool)`` only
    when the validity state transitions.
    """
    try:
      self._parse(self.toPlainText())
      self._error = ""
    except Exception as e:
      self._error = str(e)
    is_valid = (self._error == "")
    if is_valid:
      self.setStyleSheet(self._normal_style)
    else:
      self.setStyleSheet(_invalid_te_style())
    if is_valid != self._was_valid:
      self._was_valid = is_valid
      self.validityChanged.emit(is_valid)

  # -- Qt overrides -----------------------------------------------------------

  def focusOutEvent(self, event):
    """Commit on focus loss; emit ``valueChanged`` if valid."""
    self.validate()
    if self.isValid():
      self.valueChanged.emit(self.value())
    QPlainTextEdit.focusOutEvent(self, event)

  # -- internal slots ---------------------------------------------------------

  def _on_text_changed(self):
    """Slot for QPlainTextEdit.textChanged; re-runs validation."""
    self.validate()
