import os
import sys

from libtbx.utils import format_cpu_times

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
try:
  from qttbx.qt.QtWidgets import QApplication
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)


def _app():
  from qttbx.widgets.font_init import init_default_app_font
  app = QApplication.instance() or QApplication(sys.argv)
  init_default_app_font(app)
  return app


def exercise_credentials_dialog_base_returns_none_on_cancel():
  from qttbx.widgets.chat.credentials_dialog import CredentialsDialog
  app = _app()
  dlg = CredentialsDialog(title="Test", instructions="Enter a key.",
                          field_label="API key")
  # No interaction; programmatic cancel returns None.
  dlg._on_cancel()
  assert dlg.result_value() is None


def exercise_credentials_dialog_save_returns_entered_value():
  from qttbx.widgets.chat.credentials_dialog import CredentialsDialog
  app = _app()
  dlg = CredentialsDialog(title="Test", instructions="Enter a key.",
                          field_label="API key")
  dlg.set_value("sk-test-1234")
  dlg._on_save()
  assert dlg.result_value() == "sk-test-1234"


def exercise_credentials_dialog_save_strips_whitespace():
  """A pasted key often carries surrounding whitespace / a trailing newline
  ('  sk-ant-abc\n'). Unstripped it reaches the SDK and httpx rejects the
  Authorization header -> a spurious auth failure. _on_save must .strip()."""
  from qttbx.widgets.chat.credentials_dialog import CredentialsDialog
  app = _app()
  dlg = CredentialsDialog(title="Test", instructions="Enter a key.",
                          field_label="API key")
  dlg.set_value("  sk-ant-abc\n")
  dlg._on_save()
  assert dlg.result_value() == "sk-ant-abc", repr(dlg.result_value())


def exercise_credentials_dialog_show_hide_toggle():
  from qttbx.qt.QtWidgets import QLineEdit
  from qttbx.widgets.chat.credentials_dialog import CredentialsDialog
  app = _app()
  dlg = CredentialsDialog(title="Test", instructions="Enter a key.",
                          field_label="API key")
  # Default: masked
  assert dlg._value_edit.echoMode() == QLineEdit.EchoMode.Password
  dlg._on_toggle_show()
  assert dlg._value_edit.echoMode() == QLineEdit.EchoMode.Normal
  dlg._on_toggle_show()
  assert dlg._value_edit.echoMode() == QLineEdit.EchoMode.Password


def exercise_anthropic_credentials_dialog_uses_anthropic_metadata():
  from qttbx.widgets.chat.credentials_dialog import (
    AnthropicCredentialsDialog)
  app = _app()
  dlg = AnthropicCredentialsDialog()
  # Anthropic-specific title and console link in instructions.
  assert "Anthropic" in dlg._title or "Anthropic" in dlg._instructions_text
  assert "console.anthropic.com" in dlg._instructions_text


def exercise():
  exercise_credentials_dialog_base_returns_none_on_cancel()
  exercise_credentials_dialog_save_returns_entered_value()
  exercise_credentials_dialog_save_strips_whitespace()
  exercise_credentials_dialog_show_hide_toggle()
  exercise_anthropic_credentials_dialog_uses_anthropic_metadata()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
