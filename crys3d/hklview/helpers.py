from __future__ import absolute_import, division, print_function
# helper module for our own classes and widgets

import numpy as np
import matplotlib.pyplot as plt
from PySide2 import QtWidgets
from PySide2.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PySide2.QtCore import Qt, QEvent, QAbstractTableModel, QModelIndex
from PySide2.QtGui import QCursor
from PySide2.QtWidgets import ( QCheckBox, QTableWidget, QAction, QMenu,
      QTableView, QDialog,  QSpinBox, QLabel, QComboBox, QGridLayout, QGroupBox,
      QScrollArea, QVBoxLayout
     )
import math, csv
from io import StringIO


class HeaderDataTableWidget(QTableWidget):
  def __init__(self, *args, **kwargs):
    QTableWidget.__init__(self, *args, **kwargs)
    self.mousebutton = None
    self.selectedrows = []
  def mousePressEvent(self, event):
    if event.type() == QEvent.MouseButtonPress:
      self.mousebutton = None
      if event.button() == Qt.RightButton:
        self.mousebutton = Qt.RightButton
    QTableWidget.mousePressEvent(self, event)
  def mouseDoubleClickEvent(self, event):
    if event.type() == QEvent.MouseButtonDblClick:
      if event.button() == Qt.LeftButton:
        self.mousebutton = QEvent.MouseButtonDblClick
    QTableWidget.mouseDoubleClickEvent(self, event)


class MillerArrayTableForm(QDialog):
  def __init__(self, parent=None):
    super(MillerArrayTableForm, self).__init__(parent.window)
    self.setWindowTitle("Tabulated Reflection Data")
    self.precision_spinBox = QSpinBox()
    self.precision_spinBox.setSingleStep(1)
    self.precision_spinBox.setRange(1, 20)
    self.precision_spinBox.setValue(3)
    self.precision_spinBox.valueChanged.connect(parent.onPrecisionChanged)
    precision_labeltxt = QLabel()
    precision_labeltxt.setText("Precision:")
    self.SortComboBox = QComboBox()
    self.SortComboBox.activated.connect(parent.onSortComboBoxSelchange)
    sort_labeltxt = QLabel()
    sort_labeltxt.setText("Sort according to:")
    self.sortChkbox = QCheckBox()
    self.sortChkbox.setCheckState(Qt.Unchecked)
    self.sortChkbox.setText("Ascending order")
    self.sortChkbox.clicked.connect(parent.onSortChkbox)
    self.myGroupBox = QGroupBox()
    self.layout = QGridLayout()
    self.layout.addWidget(precision_labeltxt,       0, 0, 1, 1)
    self.layout.addWidget(self.precision_spinBox,   0, 1, 1, 1)
    self.layout.addWidget(sort_labeltxt,            0, 2, 1, 1)
    self.layout.addWidget(self.SortComboBox,        0, 3, 1, 1)
    self.layout.addWidget(self.sortChkbox,          0, 4, 1, 1)
    self.layout.addWidget(parent.millerarraytable,  1, 0, 1, 5)
    self.layout.setColumnStretch (0 ,0)
    self.layout.setColumnStretch (1 ,0)
    self.layout.setColumnStretch (2 ,0)
    self.layout.setColumnStretch (3 ,0)
    self.layout.setColumnStretch (4 ,1)
    self.myGroupBox.setLayout(self.layout)
    self.mainLayout = QGridLayout()
    self.mainLayout.addWidget(self.myGroupBox,     0, 0)
    self.setLayout(self.mainLayout)
  def eventFilter(self, source, event):
    if (event.type() == QEvent.KeyPress and
      event.matches(QKeySequence.Copy)):
      self.parent().parent.millerarraytable.copySelection()
      return True
    return super(MillerArrayTableForm, self).eventFilter(source, event)


class MillerArrayTableView(QTableView):
  def __init__(self, *args, **kwargs):
    QTableView.__init__(self, *args, **kwargs)
    myqa = QAction("Copy selected cells...", self)
    myqa.setData( ("Copying selection" ))
    self.tablemenu = QMenu(self)
    self.tablemenu.addAction(myqa)
    self.tablemenu.triggered.connect(self.onTableMenuAction)
    self.setContextMenuPolicy(Qt.CustomContextMenu)
    self.customContextMenuRequested.connect(self.onRightClick)
  def onRightClick(self, QPos=None):
    parent=self.sender()
    self.tablemenu.move(QCursor.pos())
    self.tablemenu.show()
  def onTableMenuAction(self, action):
    data = action.data()
    if data == "Copying selection":
      self.copySelection()
  def copySelection(self):
    # from https://stackoverflow.com/questions/40225270/copy-paste-multiple-items-from-qtableview-in-pyqt4
    selection = self.selectedIndexes()
    if selection:
      rows = sorted(index.row() for index in selection)
      columns = sorted(index.column() for index in selection)
      rowcount = rows[-1] - rows[0] + 1
      colcount = columns[-1] - columns[0] + 1
      table = [[''] * colcount for _ in range(rowcount)]
      for index in selection:
        row = index.row() - rows[0]
        column = index.column() - columns[0]
        table[row][column] = index.data()
      stream = StringIO()
      csv.writer(stream, delimiter='\t').writerows(table)
      self.parent().parent().parent().parent.app.clipboard().setText(stream.getvalue())


class MillerArrayTableModel(QAbstractTableModel):
  def __init__(self, data, headerdata, parent=None):
    super(MillerArrayTableModel, self).__init__(parent.window)
    # input data are a list of column arrays organised as:
    # [[list of H], [list of K], [list of L], [list of millerdata1], [list of millersigmas1], [list of millerdata2]... ]
    # Use zip to transpose it from a list of columns data to matching rows of data values for the table
    self._data = list(zip(*data))
    # values strictly smaller than minima in case of any NaN present in the columns
    # used for lambda function when sorting table and encountering NaN values in sorting column
    if len(data[0]):
      self.minvals = [ (min(col)-0.1) if type(col[0]) is not str else 1 for col in data ]
    self.columnheaderdata = headerdata
    self.precision = 4
  def rowCount(self, parent=None):
    return len(self._data)
  def columnCount(self, parent=None):
    return len(self._data[0]) if self.rowCount() else 0
  def data(self, index, role=Qt.DisplayRole):
    if role == Qt.DisplayRole:
      row = index.row()
      if 0 <= row < self.rowCount():
        column = index.column()
        if 0 <= column < self.columnCount():
          val = self._data[row][column]
          if not (type(val) is float or type(val) is int):
            return val
          if math.isnan(val):
            return None
          if abs(val) < float("1e-%d" %self.precision):
            fstr = "%" + "%d" %self.precision
            fstr += ".%de" %self.precision
            return float(fstr %val)
          p = 10 ** self.precision
          if val > 0:
            return float(math.floor((val * p) + 0.5))/p
          return float(math.ceil((val * p) - 0.5))/p
  def clear(self):
    rows = self.rowCount()
    self.beginRemoveRows(QModelIndex(), 0, rows-1 )
    self.endRemoveRows()
    cols = self.columnCount()
    self.beginRemoveColumns(QModelIndex(), 0, cols-1 )
    self.endRemoveColumns()
  def headerData(self, index, orientation, role):
    if orientation == Qt.Horizontal and role == Qt.DisplayRole:
      return self.columnheaderdata[index]
    if orientation == Qt.Vertical and role == Qt.DisplayRole:
      return index + 1 # just return the row number
  def sort(self, col, order):
    """Sort table by given column number.
    """
    self.layoutAboutToBeChanged.emit()
    if order == Qt.AscendingOrder:
      #print(self.columnheaderdata[col] + " sort AscendingOrder")
      self._data = sorted(self._data, key= lambda data: self.minvals[col] if math.isnan(data[col]) else data[col])
    if order == Qt.DescendingOrder:
      #print(self.columnheaderdata[col] + " sort DescendingOrder")
      self._data = sorted(self._data, key= lambda data: self.minvals[col] if math.isnan(data[col]) else data[col], reverse=True)
    self.layoutChanged.emit()


# Dialog box with MatPlotLib colour gradient charts from 
# http://matplotlib.org/examples/color/colormaps_reference.html

import numpy as np
import sys, time
import matplotlib.pyplot as plt
from PySide2 import QtWidgets
from PySide2.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# list all colour maps except their reverse which end with "_r"
cmaps = [ c for c in plt.colormaps() if "_r" not in c]
gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

dpi=50


class MplCanvas(FigureCanvas):
  def __init__(self, parent=None, width=5, height=4, dpi=dpi):
    self.parent=parent
    self.fig = Figure(figsize=(10, 20), dpi=dpi, facecolor=(1, 1, 1), edgecolor=(0.5, 0, 0))
    self.axes = self.fig.subplots(nrows=len(cmaps), ncols=1)
    # alignment of each subplot.
    self.fig.subplots_adjust(top=0.995, 
                             bottom=0.01, 
                             left=0.22, # leave room for the label naming this colour map
                             right=0.98 # leave a small gap before edge of canvas
                             )
    # total size of canvas
    self.fig.set_size_inches(7,40) # total size of canvas
    super(MplCanvas, self).__init__(self.fig)
    cid = self.fig.canvas.mpl_connect('button_press_event', self.on_press)

  def on_press(self, event):
    if event.inaxes is not None:
      self.parent.selcolmap = cmaps[event.inaxes.get_subplotspec().rowspan.start]
      self.parent.labeltxt.setText('Colour gradient map is: %s, Click a map to select a different one' %self.parent.selcolmap )

# TODO work out scaling of canvas to match QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
# and
class MPLColourSchemes(QtWidgets.QDialog):
  def __init__(self, parent=None):
    super(MPLColourSchemes, self).__init__(parent.window)
    self.parent = parent
    self.isOK = False
    self.selcolmap = ""
    #self.setWindowFlags(Qt.Tool)
    # Create the maptlotlib FigureCanvas object, 
    # which defines a single set of axes as self.axes.
    self.labeltxt = QtWidgets.QLabel()
    self.labeltxt.setText("Click on a gradient map for colouring data values")
    sc = MplCanvas(self, dpi=dpi)
    for ax, name in zip(sc.axes, cmaps):
      ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
      ax.set_axis_off()
      pos = list(ax.get_position().bounds)
      x_text = pos[0] - 0.21
      y_text = pos[1] + pos[3]/2.
      sc.fig.text(x_text, y_text, name, va='center', ha='left', fontsize=15)
    #sc.fig.tight_layout()
    scroll = QtWidgets.QScrollArea()
    scroll.setWidget(sc)
    scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    self.reversecheckbox = QtWidgets.QCheckBox()
    self.reversecheckbox.setText("Reverse colour mapping")
    self.reversecheckbox.clicked.connect(self.onReverseMap)
    self.OKbtn = QtWidgets.QPushButton("OK")
    self.OKbtn.clicked.connect(self.onOK)
    self.Cancelbtn =  QtWidgets.QPushButton("Cancel")
    self.Cancelbtn.clicked.connect(self.onCancel)
    gridlayout = QtWidgets.QGridLayout() 
    gridlayout.addWidget(self.labeltxt,          0, 0, 1, 2)
    gridlayout.addWidget(self.reversecheckbox,   1, 0, 1, 2)
    gridlayout.addWidget(scroll,                 2, 0, 1, 2)
    gridlayout.addWidget(self.OKbtn,             3, 0, 1, 1)
    gridlayout.addWidget(self.Cancelbtn,         3, 1, 1, 1)
    self.setLayout(gridlayout)
    scw = self.parent.app.style().pixelMetric(QtWidgets.QStyle.PM_ScrollBarExtent)
    self.setFixedWidth(sc.width() + scw*2 ) # fudge
    #self.setFixedWidth(sc.width() +scroll.verticalScrollBar().width() )
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )

  def onOK(self):
    self.isOK = True
    if hasattr(self.parent,"onColourChartSelect"):
      self.parent.onColourChartSelect(self.selcolmap)
    self.hide()

  def onCancel(self):
    self.isOK = False
    self.hide()

  def onReverseMap(self):
    if self.reversecheckbox.isChecked():
      if not self.selcolmap.endswith( "_r"):
        self.selcolmap = self.selcolmap + "_r"
    else:
      if self.selcolmap.endswith( "_r"):
        self.selcolmap = self.selcolmap[:-2]
    self.labeltxt.setText('Selected colour gradient map: %s' %self.selcolmap )




