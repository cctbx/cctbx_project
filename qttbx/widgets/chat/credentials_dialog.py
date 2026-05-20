"""Credentials dialog widgets for chat agents.

The base ``CredentialsDialog`` provides the common UX (instructions,
masked input, show/hide toggle, save-to-keychain checkbox, save and
cancel buttons). Each agent subclasses to set its own title and
instructions.
"""

from qttbx.qt import QtWidgets


class CredentialsDialog(QtWidgets.QDialog):
  """Parent-agnostic credentials prompt. Subclass and configure via
  constructor args."""

  def __init__(self, title, instructions, field_label, parent=None):
    super().__init__(parent)
    self._title = title
    self._instructions_text = instructions
    self._field_label = field_label
    self._result_value = None

    self.setWindowTitle(title)
    self._build_ui()

  def _build_ui(self):
    layout = QtWidgets.QVBoxLayout(self)

    self._instructions_label = QtWidgets.QLabel(self._instructions_text)
    self._instructions_label.setWordWrap(True)
    self._instructions_label.setOpenExternalLinks(True)
    layout.addWidget(self._instructions_label)

    row = QtWidgets.QHBoxLayout()
    row.addWidget(QtWidgets.QLabel(self._field_label + ":"))
    self._value_edit = QtWidgets.QLineEdit()
    self._value_edit.setEchoMode(QtWidgets.QLineEdit.EchoMode.Password)
    row.addWidget(self._value_edit)
    self._toggle_btn = QtWidgets.QPushButton("Show")
    self._toggle_btn.setCheckable(True)
    self._toggle_btn.clicked.connect(self._on_toggle_show)
    row.addWidget(self._toggle_btn)
    layout.addLayout(row)

    self._save_to_keychain = QtWidgets.QCheckBox(
      "Save to system keychain (recommended)")
    self._save_to_keychain.setChecked(True)
    layout.addWidget(self._save_to_keychain)

    btns = QtWidgets.QHBoxLayout()
    self._cancel_btn = QtWidgets.QPushButton("Cancel")
    self._cancel_btn.clicked.connect(self._on_cancel)
    self._save_btn = QtWidgets.QPushButton("Save && Start")
    self._save_btn.clicked.connect(self._on_save)
    btns.addStretch(1)
    btns.addWidget(self._cancel_btn)
    btns.addWidget(self._save_btn)
    layout.addLayout(btns)

  # ---- public API used by callers ------------------------------------------

  def set_value(self, value):
    """Pre-populate (e.g., for an 'Update key' flow)."""
    self._value_edit.setText(value)

  def result_value(self):
    """Returns the saved value, or None if cancelled or save not invoked."""
    return self._result_value

  def save_to_keychain_checked(self):
    return self._save_to_keychain.isChecked()

  # ---- slots ---------------------------------------------------------------

  def _on_toggle_show(self):
    # Toggle based on the current echo mode so the slot works whether the
    # button's checked state has been flipped by Qt (real click) or not
    # (direct programmatic call from tests).
    if self._value_edit.echoMode() == QtWidgets.QLineEdit.EchoMode.Password:
      self._value_edit.setEchoMode(QtWidgets.QLineEdit.EchoMode.Normal)
      self._toggle_btn.setText("Hide")
      self._toggle_btn.setChecked(True)
    else:
      self._value_edit.setEchoMode(QtWidgets.QLineEdit.EchoMode.Password)
      self._toggle_btn.setText("Show")
      self._toggle_btn.setChecked(False)

  def _on_cancel(self):
    self._result_value = None
    self.reject()

  def _on_save(self):
    self._result_value = self._value_edit.text()
    self.accept()


class AnthropicCredentialsDialog(CredentialsDialog):
  """Anthropic-specific credential prompt."""

  TITLE = "PhenixChat needs an Anthropic API key"
  INSTRUCTIONS = (
    "PhenixChat uses Claude (Anthropic). You need an API key to start "
    "chatting.<br><br>"
    "Get a key at "
    "<a href='https://console.anthropic.com'>console.anthropic.com</a>.<br><br>"
    "Tip: you can also set the env var <code>PHENIX_ANTHROPIC_API_KEY</code> "
    "to skip this dialog.")

  def __init__(self, parent=None):
    super().__init__(title=self.TITLE,
                     instructions=self.INSTRUCTIONS,
                     field_label="API key",
                     parent=parent)
