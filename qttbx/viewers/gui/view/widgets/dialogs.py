from PySide2.QtWidgets import (
    QDialog,
    QTextEdit,
    QVBoxLayout,
)


class InfoDialog(QDialog):
  def __init__(self, text, title="Info",parent=None):
    super().__init__(parent)
    self.setWindowTitle(title)

    self.text_edit = QTextEdit(self)
    self.text_edit.setReadOnly(True)
    self.text_edit.setText(text)

    layout = QVBoxLayout(self)
    layout.addWidget(self.text_edit)