from __future__ import absolute_import, division, print_function
# helper module for our own classes and widgets

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from .qt import QtWidgets
from .qt import Qt, QEvent, QAbstractTableModel, QModelIndex
from .qt import QCursor, QKeySequence, QLineEdit, QTextDocument
from .qt import ( QAbstractItemView, QCheckBox, QTableWidget, QAction, QDoubleSpinBox,
      QMenu, QTableView, QDialog, QSpinBox, QLabel, QComboBox, QGridLayout, QGroupBox,
      QScrollArea, QVBoxLayout, QHeaderView, QTableWidgetItem, QSizePolicy, QPlainTextEdit,
      QPushButton
     )
import math, csv, sys
from io import StringIO


class FindDialog(QDialog):
  def __init__(self, parent=None):
    super(FindDialog, self).__init__(parent)
    self.setWindowTitle("Find a string")
    self.setWindowFlags(Qt.Tool)
    mainLayout = QGridLayout()
    self.label = QLabel()
    self.label.setText("Enter search string")
    self.label.setWordWrap(True)
    self.findtextbox = QLineEdit()
    self.texteditclient = None
    sp = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
    self.findtextbox.setSizePolicy(sp)
    self.label.setSizePolicy(sp)
    self.cancelbtn = QPushButton()
    self.cancelbtn.setText("Cancel")
    self.cancelbtn.clicked.connect(self.onCancel)
    self.findbtn = QPushButton()
    self.findbtn.setText("Find next")
    self.findbtn.clicked.connect(self.onFind)
    mainLayout.addWidget(self.label,  0, 0, 1, 3)
    mainLayout.addWidget(self.findtextbox,  1, 0, 1, 3)
    mainLayout.addWidget(self.findbtn,  2, 0, 1, 1)
    mainLayout.addWidget(self.cancelbtn,  2, 2, 1, 1)
    self.setLayout(mainLayout)
    self.setSizePolicy(sp)
    #self.setFixedSize( self.sizeHint() )
  def onFind(self):
    if self.texteditclient:
      cursor = self.texteditclient.textCursor()
      if not self.texteditclient.find(self.findtextbox.text(), QTextDocument.FindCaseSensitively):
        cursor.setPosition(0)
        self.texteditclient.setTextCursor(cursor) # wrap around if reached end of text
      self.texteditclient.activateWindow()
  def onCancel(self):
    self.hide()



class MyQPlainTextEdit(QPlainTextEdit):
  def __init__(self, parent=None):
    super(MyQPlainTextEdit,self).__init__(parent)
  def keyPressEvent(self, event):
    cursor = self.textCursor()
    if event.key() == Qt.Key_F and event.modifiers() == Qt.ControlModifier:
      cursor.setPosition(0)
      self.setTextCursor(cursor)
      # In order to use the same FindDialog for all MyQPlainTextEdit instances
      # self.finddlg is assigned to a FindDialog instance by NGL_HKLViewer.
      self.finddlg.texteditclient = self
      self.finddlg.show()
      self.finddlg.activateWindow()
    elif event.key() == Qt.Key_F3 and self.finddlg.findtextbox.text() != "":
      if not self.find(self.finddlg.findtextbox.text(), QTextDocument.FindCaseSensitively):
        cursor.setPosition(0)
        self.setTextCursor(cursor) # wrap around if reached end of text
    else:
      QPlainTextEdit.keyPressEvent(self, event)


class MyQDoubleSpinBox(QDoubleSpinBox):
  def __init__(self, parent=None):
    super(MyQDoubleSpinBox,self).__init__(parent)
  def stepBy(self, steps):
    QDoubleSpinBox.stepBy(self, steps)
    if hasattr(self, "onStepBy"):
      self.onStepBy()
  def mouseReleaseEvent(self, event):
    # used by NGL_HKLViewer.fontspinBox since changing fonts and resizing windows is slow and
    # valueChanged() will then mistakenly be fired twice when clicking once on the spin control
    QDoubleSpinBox.mouseReleaseEvent(self, event)
    if hasattr(self, "onMouseRelease"):
      self.onMouseRelease()


class MyhorizontalHeader(QHeaderView):
# Assigned to HeaderDataTableWidget (NGL_HKLViewer.millertable) but in
# NGL_HKLViewer.createFileInfoBox() as to avoid a very long chain of parents when
# accessing select_millertable_column_dlg()
# Display the labels of the columns in the millertable
# Right-clicking will invoke NGL_HKLViewer.select_millertable_column_dlg
  def __init__(self, parent=None):
    super(MyhorizontalHeader,self).__init__(Qt.Horizontal, parent)
    self.mousebutton = None
    self.parent = parent
    self.setSectionsClickable(True)
    self.setSectionsMovable(True)
    self.setToolTip("Right-click to specify which columns to show")
  def mousePressEvent(self, event):
    if event.type() == QEvent.MouseButtonPress:
      if event.button() == Qt.RightButton:
        self.parent.parent.onSelect_millertable_column_dlg()
    QHeaderView.mousePressEvent(self, event)


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
  def keyPressEvent(self, event):
    super().keyPressEvent(event)
    if event.key() == Qt.Key.Key_C and (event.modifiers() & Qt.KeyboardModifier.ControlModifier):
      copied_cells = sorted(self.selectedIndexes())
      copy_text = ''
      max_column = copied_cells[-1].column()
      for c in copied_cells:
        copy_text += self.item(c.row(), c.column()).text()
        if c.column() == max_column:
          copy_text += '\n'
        else:
          copy_text += '\t'
      QtWidgets.QApplication.instance().clipboard().setText(copy_text)


class MillerTableColumnHeaderDialog(QDialog):
  # Assigned to NGL_HKLViewer.select_millertable_column_dlg allowing the user to
  # specify with checkboxes which column of the millertable should be visible.
  def __init__(self, parent=None):
    super(MillerTableColumnHeaderDialog, self).__init__(parent.window)
    self.setWindowFlags(Qt.Tool)
    self.setSizeGripEnabled(True)
    self.setWindowTitle("Dataset properties to display")
    self.parent = parent
    self.selectcolumnstable = QTableWidget()
    self.layout = QGridLayout()
    self.layout.addWidget(self.selectcolumnstable,  0, 0)
    self.setLayout(self.layout)
    self.selectcolumnstable.itemChanged.connect(self.onSelectColumnsTableItemChanged  )
    self.selectcolumnstable.itemPressed.connect(self.onSelectColumnsTableItemChanged  )
    self.unfeedback = False
  def make_new_selection_table(self):
    self.unfeedback = True
    nrows = len(self.parent.colnames_select_lst)
    self.selectcolumnstable.clearContents()
    self.selectcolumnstable.setColumnCount(1)
    self.selectcolumnstable.setRowCount(nrows)
    for row,(philname, (short_caption, helpstr), is_selected) in enumerate(self.parent.colnames_select_lst):
      item = QTableWidgetItem(short_caption)
      item.setData(Qt.UserRole, philname) # associated phil parameter name is stored with the setData function
      item.setFlags((Qt.ItemIsUserCheckable | Qt.ItemIsEnabled) )
      item.setToolTip(helpstr)
      if is_selected:
        item.setCheckState(Qt.Checked)
      else:
        item.setCheckState(Qt.Unchecked)
      self.selectcolumnstable.setItem(row, 0, item)
    self.unfeedback = False
    self.resize()
  def onSelectColumnsTableItemChanged(self, item):
    if self.unfeedback:
      return
    philname = item.data(Qt.UserRole) # get the phil parameter name stored as data
    if item.checkState()==Qt.Unchecked:
      is_selected = False
    else:
      is_selected = True
    current_philstr = "selected_info.%s = %s" %(philname, is_selected)
    self.parent.send_message(current_philstr)
  def resize(self):
    self.selectcolumnstable.resizeColumnsToContents()
    self.setFixedWidth(
      self.selectcolumnstable.verticalHeader().width() +
      self.selectcolumnstable.horizontalHeader().length() +
      self.selectcolumnstable.frameWidth()*2 +
      self.selectcolumnstable.verticalScrollBar().width()*2
      )



class MillerArrayTableForm(QDialog):
  # dialog for showing miller array data, sigmas and hkl values
  def __init__(self, parent=None):
    super(MillerArrayTableForm, self).__init__(parent.window)
    self.setWindowFlag(Qt.WindowContextHelpButtonHint,False);
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
    self.doubleClicked.connect(self.onDoubleClick)
  def onDoubleClick(self, index):
    hkl = (int(index.siblingAtColumn(0).data()),
           int(index.siblingAtColumn(1).data()),
           int(index.siblingAtColumn(2).data()))
    self.parent().parent().parent().parent.HighlightReflection(hkl)
  def onRightClick(self, QPos=None):
    self.parent().parent().parent().parent.HighlightReflection("deselect") # deselects any highlighted hkl
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
            return str(val)
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
    # when sorting the list list define any possible NaN values as -sys.maxsize so these will come last in the sorted list
    if order == Qt.AscendingOrder:
      #print(self.columnheaderdata[col] + " sort AscendingOrder")
      self._data = sorted(self._data, key= lambda data: -sys.maxsize if math.isnan(data[col]) else data[col])
    if order == Qt.DescendingOrder:
      #print(self.columnheaderdata[col] + " sort DescendingOrder")
      self._data = sorted(self._data, key= lambda data: -sys.maxsize if math.isnan(data[col]) else data[col], reverse=True)
    self.layoutChanged.emit()


# Dialog box with MatPlotLib colour gradient charts from
# http://matplotlib.org/examples/color/colormaps_reference.html


# list all colour maps except their reverse which end with "_r"
cmaps = [ c for c in plt.colormaps() if not c.endswith("_r")]
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
    self.cid1 = self.fig.canvas.mpl_connect('button_press_event', self.on_press)

  def on_press(self, event):
    if event.inaxes is not None:
      self.parent.selcolmap = cmaps[event.inaxes.get_subplotspec().rowspan.start]
      self.parent.updatelabel()
      self.parent.EnactColourMapSelection()


# TODO work out scaling of canvas to match QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
# and
class MPLColourSchemes(QtWidgets.QDialog):
  def __init__(self, parent=None):
    super(MPLColourSchemes, self).__init__(parent.window)
    self.setWindowFlags(Qt.Tool)
    self.parent = parent
    self.selcolmap = ""
    self.datatypestr = ""
    self.powscale = 1
    #self.setWindowFlags(Qt.Tool)
    # Create the maptlotlib FigureCanvas object,
    # which defines a single set of axes as self.axes.
    self.labeltxt = QtWidgets.QLabel()
    self.labeltxt.setText("Click on a gradient map for colouring data values")
    self.mycanvas = MplCanvas(self, dpi=dpi)
    self.draw_axes_and_text()
    scroll = QtWidgets.QScrollArea()
    scroll.setWidget(self.mycanvas)
    scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    self.reversecheckbox = QtWidgets.QCheckBox()
    self.reversecheckbox.setText("Reverse colour mapping")
    self.reversecheckbox.clicked.connect(self.onReverseMap)
    self.powscale_label = QtWidgets.QLabel()
    self.powscale_label.setText("Power factor for map scaling:")
    self.powscaleslider = QtWidgets.QSlider(Qt.Horizontal)
    self.powscaleslider.setMinimum(-20)
    self.powscaleslider.setMaximum(20)
    self.powscaleslider.setTickPosition(QtWidgets.QSlider.TicksAbove)
    self.powscaleslider.setTickInterval(1)
    self.powscaleslider.sliderReleased.connect(self.onReleasePowscaleslider)
    self.powscaleslider.valueChanged.connect(self.onValueChangedPowscaleslider)
    self.powscaletxtbox = QtWidgets.QLineEdit('')
    self.powscaletxtbox.setReadOnly(True)
    self.OKbtn = QtWidgets.QPushButton("OK")
    self.OKbtn.clicked.connect(self.onOK)
    gridlayout = QtWidgets.QGridLayout()
    gridlayout.addWidget(self.labeltxt,          0, 0, 1, 2)
    gridlayout.addWidget(self.reversecheckbox,   1, 0, 1, 1)
    gridlayout.addWidget(self.powscale_label,    2, 0, 1, 1)
    gridlayout.addWidget(self.powscaleslider,    2, 1, 1, 1)
    gridlayout.addWidget(self.powscaletxtbox,    2, 2, 1, 1)
    gridlayout.addWidget(scroll,                 3, 0, 1, 3)
    gridlayout.addWidget(self.OKbtn,             4, 1, 1, 1)
    self.setLayout(gridlayout)

  def draw_axes_and_text(self):
    for ax, name in zip(self.mycanvas.axes, cmaps):
      ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
      ax.set_axis_off()
      pos = list(ax.get_position().bounds)
      x_text = pos[0] - 0.21
      y_text = pos[1] + pos[3]/2.
      self.mycanvas.fig.text(x_text, y_text, name, va='center', ha='left',
                             fontsize= self.parent.app.font().pointSize()*1.5 )

  def resizeEvent(self, event):
    # MplCanvas doesn't resize with the rest of QtWidgets whenever
    # triggered by a font size changes from the main GUI. So resize here instead
    if event.type() == QEvent.Type.Resize:
      ltxt = len(self.mycanvas.fig.texts)
      for i in range(ltxt): # delete from the end as the lists is changed
        self.mycanvas.fig.texts[ltxt-(i+1)].remove()
      self.draw_axes_and_text()
      qs = self.sizeHint()
      self.mycanvas.resize(qs.width()*0.93, qs.height()*4)
    QtWidgets.QDialog.resizeEvent(self, event)

  def EnactColourMapSelection(self):
    if hasattr(self.parent,"onColourChartSelect"):
      self.parent.onColourChartSelect(self.selcolmap, self.powscale)

  def showEvent(self, event):
    self.updatelabel()
    self.powscaletxtbox.setText("%2.2f" %self.powscale )

  def onOK(self):
    self.hide()

  def updatelabel(self):
    self.labeltxt.setText('Selected colour gradient map: %s for %s data' %(self.selcolmap, self.datatypestr) )

  def onReverseMap(self):
    if self.reversecheckbox.isChecked():
      if not self.selcolmap.endswith( "_r"):
        self.selcolmap = self.selcolmap + "_r"
    else:
      if self.selcolmap.endswith( "_r"):
        self.selcolmap = self.selcolmap[:-2]
    self.updatelabel()
    self.EnactColourMapSelection()

  def onReleasePowscaleslider(self):
    self.EnactColourMapSelection()

  def onValueChangedPowscaleslider(self):
    val= self.powscaleslider.value()
    # want to raise the colour scaling to a power bigger than 0
    # so compute powscale from an exponential of the slider value
    self.powscale = math.pow(1.1, val) # 1.1 varies sufficiently slowly for the slider range [-10,10]
    self.powscaletxtbox.setText("%2.2f" %self.powscale )

  def setPowerScaleSliderVal(self, power):
    self.powscale = power
    val = math.log(power)/math.log(1.1)
    self.powscaleslider.setValue(int(val))
    self.updatelabel()

  def setDataType(self, datatypestr):
    self.datatypestr = datatypestr
    self.updatelabel()
    if datatypestr == "Map coeffs" or datatypestr == "Phases":
      self.powscaleslider.setDisabled(True)
      self.setPowerScaleSliderVal(1.0)
    else:
      self.powscaleslider.setEnabled(True)
