
from PySide2.QtCore import Qt, Signal
from PySide2.QtWidgets import QHeaderView, QTableView
from PySide2.QtGui import QFontMetrics
from PySide2.QtCore import QModelIndex

from PySide2.QtWidgets import QTableView
from PySide2.QtWidgets import QHeaderView, QTableView

"""
View for a Pandas Table generally
"""

class PandasTableView(QTableView):
    """
    TODO: Refactor to more MVC-like architecture.
    """
    mouseReleased = Signal()
    MAX_COLUMN_WIDTH = 200

    def __init__(self, parent=None, default_col_width=75):
        super().__init__(parent)
        self.horizontalHeader().setDefaultSectionSize(default_col_width)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)

        # make sortable
        self.setSortingEnabled(True)

        # Set table selection granularity
        self.setSelectionBehavior(QTableView.SelectRows)


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
