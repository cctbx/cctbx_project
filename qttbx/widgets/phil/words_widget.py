"""WordsWidget, WordsTextWidget -- PHIL words list editors."""

from PySide2.QtCore import QSignalBlocker
from PySide2.QtWidgets import QHBoxLayout, QVBoxLayout

from libtbx import Auto
from libtbx.phil import tokenizer
from libtbx.phil.tokenizer import quote_python_str
from qttbx.widgets.phil import PhilWidget
from qttbx.widgets.phil.text_base import ValidatedLineEdit, ValidatedTextEdit


def _tokenize_words(text):
  """Parse ``text`` into a list of ``tokenizer.word`` instances.

  Mirrors how PHIL files are tokenized: whitespace-separated, double-quoted
  segments preserve embedded spaces.

  Parameters
  ----------
    text: str

  Returns
  -------
    words: list of libtbx.phil.tokenizer.word
  """
  return list(tokenizer.word_iterator(
    input_string=text,
    file_name=None,
    list_of_settings=[tokenizer.settings(contiguous_word_characters="")]))


def _format_words(value):
  """Format a list of ``tokenizer.word`` (or strings) for display.

  Re-emits ``"..."`` quoting for tokens that originally had a quote_token or
  contain whitespace; embedded ``"`` characters are escaped.

  Parameters
  ----------
    value: list, None, or libtbx.Auto

  Returns
  -------
    text: str
  """
  if value is None or value is Auto:
    return ""
  out = []
  for w in value:
    s = str(w.value)
    quote = getattr(w, "quote_token", None)
    if quote or " " in s or "\t" in s:
      out.append(quote_python_str(quote_token='"', string=s))
    else:
      out.append(s)
  return " ".join(out)


def _check_size_words(items, size_min, size_max):
  """Raise ValueError if ``items`` violates the size bounds.

  Parameters
  ----------
    items: sequence
    size_min: int or None
    size_max: int or None
  """
  n = len(items)
  if size_min is not None and n < size_min:
    raise ValueError("{n} items below minimum {m}".format(n=n, m=size_min))
  if size_max is not None and n > size_max:
    raise ValueError("{n} items exceed maximum {m}".format(n=n, m=size_max))


class WordsWidget(PhilWidget):
  """PHIL words list editor (whitespace-separated single line).

  The PHIL ``words`` type is a ``list[tokenizer.word]``. Tokens may be
  double-quoted to preserve embedded whitespace, mirroring how PHIL files
  are tokenized. For multi-line entry use ``WordsTextWidget``.
  """

  def __init__(self, definition, parent=None,
               size_min=None, size_max=None, width=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      size_min, size_max: int, optional
        Inclusive bounds on the number of items.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._size_min = size_min
    self._size_max = size_max
    self._line_edit = ValidatedLineEdit(parse=self._parse, parent=self)
    if width is not None:
      self._line_edit.setMinimumWidth(int(width))
    layout = QHBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._line_edit)
    self._line_edit.valueChanged.connect(self._emit_value_changed)
    self._line_edit.validityChanged.connect(self.validityChanged)
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)
    self._line_edit.validate()

  def value(self):
    """Return the current value (list of tokenizer.word / None / Auto)."""
    return self._parse(self._line_edit.text())

  def setValue(self, value):
    """Display ``value`` in the line edit. Idempotent and signal-blocked."""
    new_text = _format_words(value)
    if new_text == self._line_edit.text():
      return
    blocker = QSignalBlocker(self._line_edit)
    self._line_edit.setText(new_text)
    del blocker
    self._line_edit.validate()

  def isValid(self):
    """Forward to the inner ValidatedLineEdit."""
    return self._line_edit.isValid()

  def errorString(self):
    """Forward to the inner ValidatedLineEdit."""
    return self._line_edit.errorString()

  def setReadOnly(self, ro=True):
    """Forward read-only state to the inner ValidatedLineEdit."""
    self._line_edit.setReadOnly(bool(ro))

  def _parse(self, text):
    """Parse callable used by the inner ValidatedLineEdit."""
    if text.strip() == "":
      if self._use_auto:
        return Auto
      if self._allow_none:
        return None
      raise ValueError("at least one item required")
    words = _tokenize_words(text)
    _check_size_words(words, self._size_min, self._size_max)
    return words

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)


class WordsTextWidget(PhilWidget):
  """PHIL words list editor (multi-line; whitespace-separated within lines).

  Long-form variant of ``WordsWidget``. Tokens are parsed via the PHIL
  tokenizer across the whole text buffer, so newlines simply act as
  whitespace separators and double-quoted segments preserve embedded
  spaces.
  """

  def __init__(self, definition, parent=None,
               size_min=None, size_max=None, width=None, lines=None):
    """
    Parameters
    ----------
      definition: libtbx.phil.definition
      parent: QWidget, optional
      size_min, size_max: int, optional
        Inclusive bounds on the number of items.
      width: int, optional
        Pixel hint passed to ``setMinimumWidth``.
      lines: int, optional
        Visible-row hint; converted to a fixed height via ``fontMetrics``.
    """
    PhilWidget.__init__(self, definition, parent)
    self._size_min = size_min
    self._size_max = size_max
    self._text_edit = ValidatedTextEdit(parse=self._parse, parent=self)
    if width is not None:
      self._text_edit.setMinimumWidth(int(width))
    if lines is not None:
      h = self._text_edit.fontMetrics().lineSpacing() * int(lines) + 8
      self._text_edit.setFixedHeight(h)
    layout = QVBoxLayout(self)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(self._text_edit)
    self._text_edit.valueChanged.connect(self._emit_value_changed)
    self._text_edit.validityChanged.connect(self.validityChanged)
    initial = None
    try:
      initial = self.definition.extract()
    except Exception:
      initial = None
    self.setValue(initial)
    self._text_edit.validate()

  def value(self):
    """Return the current value (list of tokenizer.word / None / Auto)."""
    return self._parse(self._text_edit.toPlainText())

  def setValue(self, value):
    """Display ``value`` in the text edit. Idempotent and signal-blocked."""
    new_text = _format_words(value)
    if new_text == self._text_edit.toPlainText():
      return
    blocker = QSignalBlocker(self._text_edit)
    self._text_edit.setPlainText(new_text)
    del blocker
    self._text_edit.validate()

  def isValid(self):
    """Forward to the inner ValidatedTextEdit."""
    return self._text_edit.isValid()

  def errorString(self):
    """Forward to the inner ValidatedTextEdit."""
    return self._text_edit.errorString()

  def setReadOnly(self, ro=True):
    """Forward read-only state to the inner ValidatedTextEdit."""
    self._text_edit.setReadOnly(bool(ro))

  def _parse(self, text):
    """Parse callable used by the inner ValidatedTextEdit."""
    if text.strip() == "":
      if self._use_auto:
        return Auto
      if self._allow_none:
        return None
      raise ValueError("at least one item required")
    words = _tokenize_words(text)
    _check_size_words(words, self._size_min, self._size_max)
    return words

  def _emit_value_changed(self, _v):
    """Re-derive current value (resolving Auto/None) and emit."""
    try:
      v = self.value()
    except Exception:
      return
    self.valueChanged.emit(v)
