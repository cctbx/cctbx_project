from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QWidget
)

class DefocusLineEdit(QLineEdit):
  def __init__(self,parent=None,text=''):
    super().__init__(parent,text=text)

  def focusOutEvent(self, event):
    self.parent().switch_to_label()

class EditableLabel(QWidget):
  def __init__(self, parent=None,text=''):
    self.parent_explicit = parent
    super(EditableLabel, self).__init__(parent)
    self._text = text

    self.layout = QVBoxLayout()
    # self.layout.setAlignment(Qt.AlignLeft)
    # self.layout.setContentsMargins(0, 0, 0, 0)

    #self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    self.setMinimumHeight(self.parent_explicit.height())
    self.layout.setContentsMargins(0,0,0,10) # left, top, right, bottom
    self.layout.setSpacing(0)
    self.layout.setAlignment(Qt.AlignCenter)
    self.setLayout(self.layout)

    #spacer_top= QSpacerItem(0, 0, QSizePolicy.Expanding, QSizePolicy.Expanding)
    #self.layout.addItem(spacer_top)
    self.label = QLabel(self.text)
    #self.label.setMinimumHeight(self.height())
    self.layout.addWidget(self.label)
    self.label.show()

    self.lineEdit = DefocusLineEdit(parent=self,text=self.text)
    #self.lineEdit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    #self.setMinimumHeight(self.parent_explicit.height())
    #self.lineEdit.setParent(self)
    self.layout.addWidget(self.lineEdit)
    #self.lineEdit.setMinimumHeight(self.height())
    #spacer_bottom= QSpacerItem(0, self.parent_explicit.height()*0.4, QSizePolicy.Expanding, QSizePolicy.Expanding)
    #self.layout.addItem(spacer_bottom)
    self.lineEdit.hide()
    self.lineEdit.editingFinished.connect(self.switch_to_label)

    self.label.mousePressEvent = self.switch_to_edit

  @property
  def text(self):
    return self._text

  @text.setter
  def text(self,value):
    self._text = value
    self.label.setText(value)
    self.lineEdit.setText(value)

  def setText(self,value):
    self.text = value

  def switch_to_edit(self, event):
    self.label.hide()
    self.lineEdit.show()
    self.lineEdit.setFocus()

  def switch_to_label(self):
    self.text = self.lineEdit.text()
    self.label.setText(self.text)
    self.parent_explicit.name = self.text
    self.lineEdit.hide()
    self.label.show()

# class EditableLabel(QWidget):
#   def __init__(self, parent,text):
#     super().__init__()
#     self.parent = parent
#     self.text = text

#     self.layout = QVBoxLayout()


#     self.setLayout(self.layout)

#     self.label = QLabel(self.text)
#     self.layout.addWidget(self.label)
#     self.label.show()

#     self.lineEdit = DefocusLineEdit(self.text)
#     self.layout.addWidget(self.lineEdit)
#     self.lineEdit.hide()
#     self.lineEdit.editingFinished.connect(self.switch_to_label)

#     self.label.mousePressEvent = self.switch_to_edit

#   def switch_to_edit(self, event):
#     self.log("switching to edit")
#     self.label.hide()
#     self.lineEdit.show()
#     self.lineEdit.setFocus()

#   def switch_to_label(self):
#     self.log("switching to label")
#     self.text = self.lineEdit.text()
#     self.label.setText(self.text)
#     self.lineEdit.hide()
#     self.label.show()

#   def focusOutEvent(self, event):
#     self.switch_to_label()


# class EditableLabel(QWidget):
#   def __init__(self,parent, text):
#     super(EditableLabel, self).__init__()
#     self.parent = parent
#     self.text = text
#     self.label = QLabel(text)
#     self.lineEdit = QLineEdit(text)
#     self.lineEdit.hide()
#     self.lineEdit.returnPressed.connect(self.switchToLabel)

#     layout = QVBoxLayout()

#     layout.addWidget(self.label)
#     layout.addWidget(self.lineEdit)
#     self.setLayout(layout)

#     self.label.mousePressEvent = self.switchToLineEdit

#   def switchToLineEdit(self, event):
#     self.label.hide()
#     self.lineEdit.show()
#     self.lineEdit.setFocus()

#   def switchToLabel(self):
#     self.label.setText(self.lineEdit.text())
#     self.label.show()
#     self.lineEdit.hide()
