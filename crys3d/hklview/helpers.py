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

cmaps = [  'viridis', 'plasma', 'inferno', 'magma',
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c',
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'
          ]

import numpy as np
import sys, time
import matplotlib.pyplot as plt
from PySide2 import QtWidgets
from PySide2.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

dpi=100


class MplCanvas(FigureCanvas):
  def __init__(self, parent=None, width=5, height=4, dpi=dpi):
    self.selcolmap = ""
    self.parent=parent
    self.fig = Figure(figsize=(10, 20), dpi=dpi, facecolor=(1, 1, 1), edgecolor=(0.5, 0, 0))
    self.axes = self.fig.subplots(nrows=len(cmaps), ncols=1)
    self.fig.subplots_adjust(top=0.995, bottom=0.01, left=0.20, right=0.98)
    self.fig.set_size_inches(7,40)
    super(MplCanvas, self).__init__(self.fig)
    cid = self.fig.canvas.mpl_connect('button_press_event', self.on_press)

  def on_press(self, event):
    print("waffle")
    if event.inaxes is not None:
      self.selcolmap = cmaps[event.inaxes.get_subplotspec().rowspan.start]
      self.parent.labeltxt.setText('Selected colour gradient map: %s' %self.selcolmap )

app = QtWidgets.QApplication(sys.argv)


class MPLColourSchemes(QtWidgets.QDialog):
  def __init__(self, *args, **kwargs):
    super(MPLColourSchemes, self).__init__(*args, **kwargs)
    #self.setWindowFlags(Qt.Tool)
    # Create the maptlotlib FigureCanvas object, 
    # which defines a single set of axes as self.axes.
    self.labeltxt = QtWidgets.QLabel()
    self.labeltxt.setText("Select a colour gradient mapping for data values")
    sc = MplCanvas(self, dpi=dpi)
    for ax, name in zip(sc.axes, cmaps):
      ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
      ax.set_axis_off()
      pos = list(ax.get_position().bounds)
      x_text = pos[0] - 0.19
      y_text = pos[1] + pos[3]/2.
      sc.fig.text(x_text, y_text, name, va='center', ha='left', fontsize=15)
    #sc.fig.tight_layout()
    scroll = QtWidgets.QScrollArea()
    scroll.setWidget(sc)
    scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    OKbtn = QtWidgets.QPushButton("OK")
    Cancelbtn =  QtWidgets.QPushButton("Cancel")
    vbox = QtWidgets.QVBoxLayout() 
    vbox.addWidget(self.labeltxt)
    vbox.addWidget(scroll)
    vbox.addWidget(OKbtn)
    vbox.addWidget(Cancelbtn)

    self.setLayout(vbox)
    scw = app.style().pixelMetric(QtWidgets.QStyle.PM_ScrollBarExtent)
    self.setFixedWidth(sc.width() + scw*2 )
    #self.setFixedWidth(sc.width() +scroll.verticalScrollBar().width() )
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.show()







