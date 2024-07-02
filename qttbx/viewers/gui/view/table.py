from collections import defaultdict

from PySide2.QtCore import Qt, Signal
from PySide2.QtGui import QFontMetrics, QCursor, QBrush, QColor
from PySide2.QtCore import QModelIndex
from PySide2.QtWidgets import (
    QApplication,
    QHeaderView,
    QVBoxLayout,
    QHBoxLayout,
    QLineEdit,
    QLabel,
    QPushButton,
    QTableView,
    QMenu,
    QDialog, 
    QFormLayout, 
    QStyledItemDelegate
)
class CustomDelegate(QStyledItemDelegate):
  def __init__(self, color_map={}, font_map={}, parent=None):
    super().__init__(parent)
    self.color_map = color_map # (row,col): state.color.Color object
    self.font_map = font_map   # (row,col): {'bold': True and/or 'italic': True} 

  def paint(self, painter, option, index):
    row = index.row()
    col = index.column()

    # Set background color
    if (row, col) in self.color_map:
      color = self.color_map[(row, col)]
      painter.save()
      painter.fillRect(option.rect, QBrush(QColor(*color.RGBA)))
      painter.restore()

    # Set font style
    painter.save()
    if (row, col) in self.font_map:
      font = option.font
      font.setBold(self.font_map[(row, col)].get('bold', False))
      font.setItalic(self.font_map[(row, col)].get('italic', False))
      painter.setFont(font)
    else:
      painter.setFont(option.font)
    
    # Draw the text
    painter.drawText(option.rect, Qt.AlignVCenter | Qt.AlignHCenter, index.data())
    painter.restore()

  def update_color_map(self, row, col, color):
    self.color_map[(row, col)] = color

  def update_font_map(self, row, col, bold=False, italic=False):
    self.font_map[(row, col)] = {'bold': bold, 'italic': italic}

class GenericEditDialog(QDialog):

  def __init__(self, row_dict=None,edit_allowed_items=[],edit_denied_items=[]):
    super().__init__()
    input_names = list(row_dict.keys())
    defaults_dict = defaultdict(None)
    defaults_dict.update(row_dict)

    self.setWindowTitle(f'Edit Row')
    mainLayout = QVBoxLayout()
    self.inputsLayout = QFormLayout()

    if input_names is None:
      input_names = self.input_names
    if defaults_dict is None:
      defaults_dict = {name: "" for name in input_names}

     # Dictionary to store QLineEdit references and initial values
    self.inputFields = {}
    self.initialValues = {}

    for name in input_names:
      #inputLayout = QHBoxLayout()
      inputLabel = QLabel(name)
      default_value = str(defaults_dict.get(name, ""))
      added = False
      if (len(edit_allowed_items)>0 and name in edit_allowed_items) or (len(edit_allowed_items)==0):
        if name not in edit_denied_items:
            inputField = QLineEdit(self)
            inputField.setText(default_value)
            added = True
      if not added:
        inputField = QLabel(self)
        inputField.setText(default_value)
    
      self.inputsLayout.addRow(inputLabel, inputField)
    #   inputLayout.addWidget(inputLabel)
    #   inputLayout.addWidget(inputField)
    #   self.inputsLayout.addLayout(inputLayout)

      # Store QLineEdit reference with its name as key
      self.inputFields[name] = inputField
      # Store the initial value
      self.initialValues[name] = default_value

    mainLayout.addLayout(self.inputsLayout)
    self.setLayout(mainLayout)

    self.buttonsLayout = QHBoxLayout()
    self.cancelButton = QPushButton('Cancel', self)
    self.cancelButton.clicked.connect(self.reject)
    self.buttonsLayout.addWidget(self.cancelButton)

    self.acceptButton = QPushButton('Accept', self)
    self.acceptButton.clicked.connect(self.accept)
    self.buttonsLayout.addWidget(self.acceptButton)

    mainLayout.addLayout(self.buttonsLayout)

  def collectInputValues(self):
    # Collect values from QLineEdit widgets
    values = {name: field.text() for name, field in self.inputFields.items()}
    return values

  def getModifiedFields(self):
    # Compare initial values with current values and return modified fields
    modified_fields = {name: field.text() for name, field in self.inputFields.items() if field.text() != self.initialValues[name]}
    return modified_fields

  def exec_(self):
    # Get the cursor's current position
    cursorPos = QCursor.pos()

    # Move the dialog to the cursor's position initially
    self.move(cursorPos.x(), cursorPos.y())

    # Adjust position to ensure the dialog is fully visible on the screen
    screen = QApplication.desktop().screenNumber(cursorPos)
    screen_geom = QApplication.desktop().screenGeometry(screen)

    # Calculate the dialog's geometry after the initial move
    dialog_geom = self.geometry()
    dialog_x = cursorPos.x()
    dialog_y = cursorPos.y()
    dialog_width = dialog_geom.width()
    dialog_height = dialog_geom.height()

    # Check right boundary
    if dialog_x + dialog_width > screen_geom.right():
        dialog_x = screen_geom.right() - dialog_width

    # Check bottom boundary
    if dialog_y + dialog_height > screen_geom.bottom():
        dialog_y = screen_geom.bottom() - dialog_height

    # Check left boundary
    if dialog_x < screen_geom.left():
        dialog_x = screen_geom.left()

    # Check top boundary
    if dialog_y < screen_geom.top():
        dialog_y = screen_geom.top()

    # Move dialog to the adjusted position
    self.move(dialog_x, dialog_y)
    
    result = super().exec_()
    # Collect values if dialog was accepted
    if result == QDialog.Accepted:
      return self.collectInputValues()
    else:
      return None


"""
View for a Pandas Table generally
"""

class PandasTableView(QTableView):
    """
    TODO: Refactor to more MVC-like architecture.
    """
    mouseReleased = Signal()
    editRequested = Signal(object) # QModelIndex
    edit_dialog = GenericEditDialog
    MAX_COLUMN_WIDTH = 200

    def __init__(self, parent=None, default_col_width=75):
        super().__init__(parent)
        self.horizontalHeader().setDefaultSectionSize(default_col_width)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showContextMenu)

        # make sortable
        self.setSortingEnabled(True)

        # Set table selection granularity
        self.setSelectionBehavior(QTableView.SelectRows)

    def showContextMenu(self, position):
        index = self.indexAt(position)
        if not index.isValid():
            return

        menu = QMenu(self)
        toggleEditAction = menu.addAction("Edit")
        action = menu.exec_(self.viewport().mapToGlobal(position))

        if action == toggleEditAction:
            self.editRequested.emit(index)

    def toggleRowEditability(self, row):
        # Assuming your model has a 'toggleRowEditability' method to handle this
        self.model().toggleRowEditability(row)

    def deleteRow(self, row):
        if row < 0 or row >= self.model().rowCount():
            return

        self.model().beginRemoveRows(QModelIndex(), row, row)
        self.model()._df.drop(self.model()._df.index[row], inplace=True)
        self.model().endRemoveRows()

    def selected_indices(self):
        indexes = self.selectedIndexes()
        rows = [index.row() for index in indexes]
        unique_rows = sorted(list(set(rows)))
        return unique_rows

    def selected_rows(self):
        return self.model().df.iloc[self.selected_indices()]
        
    def mouseReleaseEvent(self, event):
      super(PandasTableView, self).mouseReleaseEvent(event)
      # Emit the signal whenever the mouse is released
      self.mouseReleased.emit()


    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Space:
            selection = self.selectionModel()
            current_index = selection.currentIndex()
            next_row = current_index.row() + 1
            if next_row < self.model().rowCount():
                next_index = self.model().index(next_row, 0)  # Move to the first column of the next row
                selection.setCurrentIndex(next_index, selection.ClearAndSelect | selection.Rows)
                event.accept()
                self.mouseReleased.emit()
                return
        super(PandasTableView, self).keyPressEvent(event)
        

    def adjustColumnWidths(self, margin=20):
        tableView = self
        model = tableView.model()
        if not model:
            return

        for column in range(model.columnCount(None)):
            maxWidth = 0
            # Ensure using the correct font for metrics
            fontMetrics = QFontMetrics(tableView.font())

            # Calculate header width with margin
            headerText = model.headerData(column, Qt.Horizontal)
            headerWidth = fontMetrics.boundingRect(headerText).width() + margin

            # Iterate over each row to find the maximum required width
            for row in range(model.rowCount(None)):
                index = model.index(row, column)
                text = model.data(index, Qt.DisplayRole)
                if text:
                    # Calculate text width and consider margin
                    textWidth = fontMetrics.boundingRect(str(text)).width() + margin
                    maxWidth = max(maxWidth, textWidth)

            # Set column width to the maximum of header and cell width, constrained by MAX_COLUMN_WIDTH
            finalWidth = min(max(headerWidth, maxWidth), self.MAX_COLUMN_WIDTH)
            tableView.setColumnWidth(column, finalWidth)

