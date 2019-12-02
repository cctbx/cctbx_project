from __future__ import absolute_import, division, print_function

#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      oeffner
#
# Created:     09/05/2019
# Copyright:   (c) oeffner 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from PySide2.QtCore import Qt, QAbstractTableModel, QEvent, QModelIndex, QSize, QTimer
from PySide2.QtWidgets import (  QAction, QApplication, QCheckBox,
        QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMenu, QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QVBoxLayout, QWidget )

from PySide2.QtGui import QColor, QFont, QCursor, QKeySequence
from PySide2.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )
import sys, zmq, subprocess, time, traceback, shutil, os.path, zlib, math, csv, io



class MakeNewDataForm(QDialog):
  def __init__(self, parent=None):
    super(MakeNewDataForm, self).__init__(parent)
    self.setWindowTitle("Create New Reflection Data")
    myGroupBox = QGroupBox("Python expression for newdata")
    layout = QGridLayout()
    layout.addWidget(parent.operationlabeltxt,     0, 0, 1, 2)
    layout.addWidget(parent.MillerLabel1,           1, 0, 1, 2)
    layout.addWidget(parent.MillerComboBox,        2, 0, 1, 1)
    layout.addWidget(parent.MillerLabel2,          2, 1, 1, 1)
    layout.addWidget(parent.MillerLabel3,          3, 0, 1, 2)
    layout.addWidget(parent.operationtxtbox,       4, 0, 1, 2)
    layout.addWidget(parent.newlabelLabel,          5, 0, 1, 1)
    layout.addWidget(parent.newlabeltxtbox,         5, 1, 1, 1)
    layout.addWidget(parent.operationbutton,       6, 0, 1, 2)
    layout.setRowStretch (0, 1)
    layout.setRowStretch (1 ,0)
    myGroupBox.setLayout(layout)
    mainLayout = QGridLayout()
    mainLayout.addWidget(myGroupBox,     0, 0)
    self.setLayout(mainLayout)
    m = self.fontMetrics().width( "asdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdf")
    self.setMinimumWidth(m)
    #self.setFixedSize( self.sizeHint() )


class SettingsForm(QDialog):
  def __init__(self, parent=None):
    super(SettingsForm, self).__init__(parent)
    self.setWindowTitle("Settings")
    myGroupBox = QGroupBox("Stuff")
    layout = QGridLayout()
    layout.addWidget(parent.mousespeed_labeltxt,     0, 0, 1, 1)
    layout.addWidget(parent.mousemoveslider,         0, 1, 1, 3)
    layout.addWidget(parent.mousesensitxtbox,        0, 4, 1, 1)
    layout.addWidget(parent.Fontsize_labeltxt,       1, 0, 1, 1)
    layout.addWidget(parent.fontspinBox,             1, 4, 1, 1)
    layout.addWidget(parent.cameraPerspectCheckBox,  2, 0, 1, 1)
    layout.addWidget(parent.bufsize_labeltxt,        3, 0, 1, 1)
    layout.addWidget(parent.bufsizespinBox,          3, 4, 1, 1)
    layout.addWidget(parent.ttipalpha_labeltxt,      4, 0, 1, 1)
    layout.addWidget(parent.ttipalpha_spinBox,        4, 4, 1, 1)
    layout.setRowStretch (0, 1)
    layout.setRowStretch (1 ,0)
    myGroupBox.setLayout(layout)
    mainLayout = QGridLayout()
    mainLayout.addWidget(myGroupBox,     0, 0)
    self.setLayout(mainLayout)
    self.setFixedSize( self.sizeHint() )


class MillerArrayTableForm(QDialog):
  def __init__(self, parent=None):
    super(MillerArrayTableForm, self).__init__(parent)
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
      self.parent().millerarraytable.copySelection()
      return True
    return super(MillerArrayTableForm, self).eventFilter(source, event)


class MyTableView(QTableView):
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
      stream = io.StringIO()
      csv.writer(stream, delimiter='\t').writerows(table)
      app.clipboard().setText(stream.getvalue())


class MyTableWidget(QTableWidget):
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


class MillerArrayTableModel(QAbstractTableModel):
  def __init__(self, data, headerdata, parent=None):
    super(MillerArrayTableModel, self).__init__(parent)
    # input data are a list of column arrays organised as:
    # [[list of H], [list of K], [list of L], [list of millerdata1], [list of millersigmas1], [list of millerdata2]... ]
    # We use zip to transpose it from a list of columns data to matching rows of data values for the table
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
      print(self.columnheaderdata[col] + " sort AscendingOrder")
      self._data = sorted(self._data, key= lambda data: data[col])
    if order == Qt.DescendingOrder:
      print(self.columnheaderdata[col] + " sort DescendingOrder")
      self._data = sorted(self._data, key= lambda data: data[col], reverse=True)
    self.layoutChanged.emit()


class NGL_HKLViewer(QWidget):
  def __init__(self, parent=None):
    super(NGL_HKLViewer, self).__init__(parent)

    self.verbose = 0
    self.UseOSbrowser = False
    self.jscriptfname = ""
    self.devmode = False
    self.hklfname = ""
    self.handshakewait = 5
    for e in sys.argv:
      if "verbose" in e:
        self.verbose = e.split("verbose=")[1]
      if "UseOSbrowser" in e:
        self.UseOSbrowser = True
      if "jscriptfname" in e:
        self.jscriptfname = e.split("jscriptfname=")[1]
      if "devmode" in e or "debug" in e:
        self.devmode = True
      if 'htmlfname' in e:
        self.hklfname = e.split("htmlfname=")[1]
      if 'handshakewait' in e:
        self.handshakewait = e.split('handshakewait=')[1]

    self.zmq_context = None
    self.unfeedback = False

    self.originalPalette = QApplication.palette()

    self.openFileNameButton = QPushButton("Open reflection file")
    self.openFileNameButton.setDefault(True)
    self.openFileNameButton.clicked.connect(self.onOpenReflectionFile)

    self.debugbutton = QPushButton("Debug")
    self.debugbutton.clicked.connect(self.DebugInteractively)

    self.settingsbtn = QPushButton("Settings")
    self.settingsbtn.clicked.connect(self.SettingsDialog)

    self.mousespeed_labeltxt = QLabel()
    self.mousespeed_labeltxt.setText("Mouse speed:")
    self.mousemoveslider = QSlider(Qt.Horizontal)
    self.mousemoveslider.setMinimum(0)
    self.mousemoveslider.setMaximum(200)
    self.mousemoveslider.setValue(0)
    self.mousemoveslider.sliderReleased.connect(self.onFinalMouseSensitivity)
    self.mousemoveslider.valueChanged.connect(self.onMouseSensitivity)
    self.mousesensitxtbox = QLineEdit('')
    self.mousesensitxtbox.setReadOnly(True)

    self.fontspinBox = QDoubleSpinBox()
    self.fontspinBox.setSingleStep(1)
    self.fontspinBox.setRange(4, 50)
    self.font = QFont()
    self.font.setFamily(self.font.defaultFamily())
    self.fontspinBox.setValue(self.font.pointSize())
    self.fontspinBox.valueChanged.connect(self.onFontsizeChanged)
    self.Fontsize_labeltxt = QLabel()
    self.Fontsize_labeltxt.setText("Font size:")

    self.cameraPerspectCheckBox = QCheckBox()
    self.cameraPerspectCheckBox.setText("Perspective camera")
    self.cameraPerspectCheckBox.clicked.connect(self.onCameraPerspect)
    self.cameraPerspectCheckBox.setCheckState(Qt.Unchecked)

    self.bufsizespinBox = QSpinBox()
    self.bufsizespinBox.setSingleStep(1)
    self.bufsizespinBox.setRange(1, 100)
    self.bufsizespinBox.setValue(10)
    self.bufsize_labeltxt = QLabel()
    self.bufsize_labeltxt.setText("Text buffer size (Kbytes):")

    self.ttipalpha = 0.85
    self.ttipalpha_spinBox = QDoubleSpinBox()
    self.ttipalpha_spinBox.setSingleStep(0.05)
    self.ttipalpha_spinBox.setRange(0.0, 1.0)
    self.ttipalpha_spinBox.setValue(self.ttipalpha)
    self.ttipalpha_spinBox.valueChanged.connect(self.onTooltipAlphaChanged)
    self.ttipalpha_labeltxt = QLabel()
    self.ttipalpha_labeltxt.setText("Tooltip Opacity:")

    self.settingsform = SettingsForm(self)

    self.MillerComboBox = QComboBox()
    self.MillerComboBox.activated.connect(self.onMillerComboSelchange)
    #self.MillerComboBox.setSizeAdjustPolicy(QComboBox.AdjustToContents)

    self.MillerLabel1 = QLabel()
    self.MillerLabel1.setText("and 'data2' and or 'sigma2' variable from the")
    self.MillerLabel2 = QLabel()
    self.MillerLabel2.setText("column")
    self.MillerLabel3 = QLabel()
    self.MillerLabel3.setText("Example: 'newdata=data/sigmas; newsigmas= -42*sigmas' ")
    self.newlabelLabel = QLabel()
    self.newlabelLabel.setText("Column label for new data:")
    self.newlabeltxtbox = QLineEdit('')
    self.operationlabeltxt = QLabel()
    self.operationtxtbox = QLineEdit('')
    self.operationbutton = QPushButton("OK")
    self.operationbutton.clicked.connect(self.onMakeNewData)
    self.makenewdataform = MakeNewDataForm(self)
    self.makenewdataform.setModal(True)

    labels = ["Column label", "Type", "Space group", "# HKLs", "Span of HKLs",
       "Min Max data", "Min Max sigmas", "d_min, d_max", "Symmetry unique", "Anomalous"]
    self.millertable = MyTableWidget(0, len(labels))
    self.millertable.setHorizontalHeaderLabels(labels)
    self.millertable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    # don't allow editing this table
    self.millertable.setEditTriggers(QTableWidget.NoEditTriggers)
    self.millertable.cellPressed.connect(self.onMillerTableCellPressed)
    self.millertable.cellDoubleClicked.connect(self.onMillerTableCellPressed)
    self.millertable.itemSelectionChanged.connect(self.onMillerTableitemSelectionChanged)

    self.millerarraytable = MyTableView(self)
    self.millerarraytable.setSortingEnabled(False)
    self.millerarraytable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    self.millerarraytableform = MillerArrayTableForm(self)
    self.millerarraytablemodel = None
    self.millerarraytable.installEventFilter(self.millerarraytableform) # for keyboard copying to clipboard

    self.createExpansionBox()
    self.createFileInfoBox()
    self.CreateSliceTabs()
    self.createRadiiScaleGroupBox()
    self.createBinsBox()
    self.CreateFunctionTabs()

    self.mainLayout = QGridLayout()
    self.mainLayout.addWidget(self.FileInfoBox,         0, 0)
    self.mainLayout.addWidget(self.functionTabWidget,   1, 0)
    self.mainLayout.addWidget(self.settingsbtn,         2, 0, 1, 1)
    self.cpath = ""
    if self.UseOSbrowser==False:
      self.BrowserBox = QWebEngineView()
      self.BrowserBox.setAttribute(Qt.WA_DeleteOnClose)
      # omitting name for QWebEngineProfile() means it is private/off-the-record with no cache files
      self.webprofile = QWebEngineProfile(parent=self.BrowserBox)
      self.webpage = QWebEnginePage( self.webprofile, self.BrowserBox)
      if self.devmode:
        #self.webpage.setUrl("https://webglreport.com/")
        self.webpage.setUrl("chrome://gpu")
      else:
        self.webpage.setUrl("https://cctbx.github.io/")
      self.cpath = self.webprofile.cachePath()
      self.BrowserBox.setPage(self.webpage)
      self.mainLayout.addWidget(self.BrowserBox,          0, 1, 5, 3)
      self.mainLayout.setColumnStretch(2, 1)
    self.mainLayout.setRowStretch(0, 1)
    self.mainLayout.setRowStretch(1, 0)
    self.mainLayout.setRowStretch(2, 1)
    #mainLayout.setRowStretch(3, 0)
    self.mainLayout.setColumnStretch(4, 0)
    self.setLayout(self.mainLayout)

    self.setWindowTitle("HKL-Viewer")
    self.cctbxproc = None
    self.LaunchCCTBXPython()
    self.out = None
    self.err = None
    self.comboviewwidth = 0
    self.currentphilstringdict = {}
    self.hklscenes_arrays = []
    self.millerarraylabels = []
    self.scenearraylabels = []
    self.array_infotpls = []
    self.matching_arrays = []
    self.bin_infotpls = None
    self.bin_opacities= None
    self.html_url = ""
    self.spacegroups = {}
    self.info = []
    self.infostr = ""
    self.tncsvec = []
    self.fileisvalid = False
    self.NewFileLoaded = False
    self.NewMillerArray = False
    self.NewHKLscenes = False
    self.updatingNbins = False
    self.binstableitemchanges = False
    self.canexit = False
    self.closing = False
    self.indices = None
    self.datalst = []
    self.tabulate_miller_array = None
    self.binTableCheckState = None
    self.millertablemenu = QMenu(self)
    self.millertablemenu.triggered.connect(self.onMillerTableMenuAction)
    self.resize(QSize(self.FileInfoBox.size().width()*2, self.size().height() ))
    self.show()


  def closeEvent(self, event):
    self.NGL_HKL_command('NGL_HKLviewer.action = is_terminating')
    self.closing = True
    self.setVisible(False)

    #while not self.canexit:
    #  time.sleep(1)
    #self.webprofile.clearHttpCache()
    del self.webprofile
    del self.webpage
    #del self.BrowserBox

    #shutil.rmtree(cpath)
    print("HKLviewer exiting now.")
    nc = 0
    sleeptime = 0.2
    while not self.canexit and nc < 5: # until cctbx.python has finished or after 5 sec
      time.sleep(sleeptime)
      self.ProcessMessages()
      nc += sleeptime

    print("accepting close.event")
    self.cctbxproc.terminate()
    self.out, self.err = self.cctbxproc.communicate()
    #print( str(self.out) + "\n" + str(self.err) )
    self.cctbxproc.wait()
    self.BrowserBox.close()
    self.mainLayout.removeWidget(self.BrowserBox)
    self.BrowserBox.deleteLater()
    #self.BrowserBox.destroy()

    """
    self.webprofile.clearHttpCache()
    self.webpage = None
    self.BrowserBox = None
    self.webprofile = None
    """
    event.accept()


  def onOpenReflectionFile(self):
    options = QFileDialog.Options()
    fileName, filtr = QFileDialog.getOpenFileName(self,
            "Open a reflection file", "",
            "All Files (*);;MTZ Files (*.mtz);;CIF (*.cif)", "", options)
    if fileName:
      self.HKLnameedit.setText(fileName)
      self.setWindowTitle("HKL-viewer: " + fileName)
      #self.infostr = ""
      self.textInfo.setPlainText("")
      self.fileisvalid = False
      self.NGL_HKL_command('NGL_HKLviewer.filename = "%s"' %fileName )
      self.MillerComboBox.clear()
      self.BinDataComboBox.clear()
      self.tncsvec = []
      self.ClipPlaneChkBox.setChecked(False)
      self.expandP1checkbox.setChecked(False)
      self.expandAnomalouscheckbox.setChecked(False)
      self.sysabsentcheckbox.setChecked(False)
      self.missingcheckbox.setChecked(False)
      self.onlymissingcheckbox.setChecked(False)
      self.showslicecheckbox.setChecked(False)


  def SettingsDialog(self):
    self.settingsform.show()


  def ProcessMessages(self):
    """
    Deal with the messages posted to this GUI by cmdlineframes.py
    """
    if self.cctbxproc:
      if self.cctbxproc.stdout:
        self.out = self.cctbxproc.stdout.read().decode("utf-8")
      if self.cctbxproc.stderr:
        self.err = self.cctbxproc.stderr.read().decode("utf-8")
    if self.out:
      print(self.out.decode("utf-8"))
    if self.err:
      print(self.err.decode("utf-8"))

    if self.zmq_context:
      try:
        binmsg = self.socket.recv(flags=zmq.NOBLOCK) #To empty the socket from previous messages
        msg = zlib.decompress(binmsg)
        nan = float("nan") # workaround for "evaluating" any NaN values in the messages received
        msgstr = msg.decode()
        self.infodict = eval(msgstr)
        if self.infodict:
          if self.infodict.get("current_ phil_ strings"):
            philstringdict = self.infodict.get("current_ phil_ strings", {})
            for k, v in philstringdict.items():
              try:
                self.currentphilstringdict[k] = eval(v)
              except Exception as e:
                self.currentphilstringdict[k] = v
            self.UpdateGUI()

          if self.infodict.get("hklscenes_arrays"):
            self.hklscenes_arrays = self.infodict.get("hklscenes_arrays", [])
            self.scenearraylabels = [ e[3] for e in self.hklscenes_arrays ]

          if self.infodict.get("array_infotpls"):
            self.array_infotpls = self.infodict.get("array_infotpls",[])
            self.millerarraylabels = [ e[0] for e in self.array_infotpls ]

          if self.infodict.get("bin_data_label"):
            self.BinDataComboBox.setCurrentText(self.infodict["bin_data_label"])

          if self.infodict.get("bin_infotpls"):
            self.bin_infotpls = self.infodict["bin_infotpls"]

            self.nbins = len(self.bin_infotpls)
            self.updatingNbins = True
            self.Nbins_spinBox.setValue(self.nbins)
            self.updatingNbins = False
            self.binstable.clearContents()
            self.binstable.setRowCount(self.nbins)
            for row,bin_infotpl in enumerate(self.bin_infotpls):
              for col,elm in enumerate(bin_infotpl):
                # only allow changing the last column with opacity values
                if col != 3:
                  item = QTableWidgetItem(str(elm))
                else:
                  item = QTableWidgetItem()
                  item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
                  item.setCheckState(Qt.Checked)
                item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                item.setFlags(item.flags() ^ Qt.ItemIsSelectable )
                self.binstable.setItem(row, col, item)
            if self.bin_opacities:
              self.update_table_opacities()

          if self.infodict.get("bin_opacities"):
            self.bin_opacities = self.infodict["bin_opacities"]
            if self.binstable.rowCount() > 0:
              self.update_table_opacities()

          if self.infodict.get("html_url"):
            self.html_url = self.infodict["html_url"]
            if self.UseOSbrowser==False:
              self.BrowserBox.setUrl(self.html_url)
              # workaround for background colour bug in chromium
              # https://bugreports.qt.io/browse/QTBUG-41960
              self.BrowserBox.page().setBackgroundColor(QColor(255, 255, 255, 0.0) )

          if self.infodict.get("spacegroups"):
            spgs = self.infodict.get("spacegroups",[])
            self.spacegroups = { i : e for i,e in enumerate(spgs) }
            self.SpaceGroupComboBox.clear()
            self.SpaceGroupComboBox.addItems( list(self.spacegroups.values()) )

          if self.infodict.get("tabulate_miller_array"):
            print("received table")
            self.tabulate_miller_array = self.infodict["tabulate_miller_array"]
            self.indices = self.tabulate_miller_array[0]
            #labels = ["H", "K", "L"] + [ ld[0] for ld in self.tabulate_miller_array[1:] ]
            labels = [ ld[0] for ld in self.tabulate_miller_array ]
            self.datalst =  [ ld[1] for ld in self.tabulate_miller_array ]

            if self.millerarraytable.model():
              self.millerarraytable.model().clear()
            self.millerarraytablemodel = MillerArrayTableModel(self.datalst, labels, self)
            self.millerarraytable.setModel(self.millerarraytablemodel)
            self.millerarraytable.horizontalHeader().setHighlightSections(False)
            self.millerarraytable.setSortingEnabled(False)

            self.millerarraytable_sortorder = ["unsorted"] * (len(self.datalst) + 3)
            self.millerarraytable.resizeColumnsToContents()
            self.millerarraytableform.layout.setRowStretch (0, 0)
            self.millerarraytableform.mainLayout.setRowStretch (0, 0)
            tablewidth = 0
            for e in range(self.millerarraytablemodel.columnCount()):
              tablewidth +=  self.millerarraytable.columnWidth(e)

            self.millerarraytableform.SortComboBox.clear()
            self.millerarraytableform.SortComboBox.addItems(["unsorted"] + labels )
            self.millerarraytableform.SortComboBox.view().setMinimumWidth(self.comboviewwidth)
            self.millerarraytableform.resize(tablewidth, self.millerarraytableform.size().height())
            self.millerarraytableform.show()

          currentinfostr = ""
          if self.infodict.get("info"):
            currentinfostr = self.infodict.get("info",[])
            if self.closing:
              print(currentinfostr)

            if "Purging HKLViewFrame" in currentinfostr:
              self.canexit = True

          if self.infodict.get("tncsvec"):
            self.tncsvec = self.infodict.get("tncsvec",[])
          if len(self.tncsvec) == 0:
            self.clipTNCSBtn.setDisabled(True)
          else:
            self.clipTNCSBtn.setEnabled(True)

          if self.infodict.get("file_name"):
            self.HKLnameedit.setText( self.infodict.get("file_name", "") )
            self.setWindowTitle("HKL-viewer: " + self.infodict.get("file_name", "") )

          if self.infodict.get("NewFileLoaded"):
            self.NewFileLoaded = self.infodict.get("NewFileLoaded",False)

          if self.infodict.get("NewHKLscenes"):
            self.NewHKLscenes = self.infodict.get("NewHKLscenes",False)

          if self.infodict.get("NewMillerArray"):
            self.NewMillerArray = self.infodict.get("NewMillerArray",False)

          self.fileisvalid = True
          #print("ngl_hkl_infodict: " + str(ngl_hkl_infodict))

          if currentinfostr:
            #print(currentinfostr)
            self.infostr += currentinfostr + "\n"
            # display no more than self.bufsize bytes of text
            self.infostr = self.infostr[-1000*self.bufsizespinBox.value():]
            self.textInfo.setPlainText(self.infostr)
            self.textInfo.verticalScrollBar().setValue( self.textInfo.verticalScrollBar().maximum()  )

          if (self.NewFileLoaded or self.NewMillerArray) and self.NewHKLscenes:
            #print("got hklscenes: " + str(self.hklscenes_arrays))
            self.NewMillerArray = False
            if self.millerarraytablemodel:
              self.millerarraytablemodel.clear()
              self.millerarraytablemodel = MillerArrayTableModel([[]], [], self)

            self.MillerComboBox.clear()
            self.MillerComboBox.addItems( self.millerarraylabels )
            self.MillerComboBox.setCurrentIndex(-1) # unselect the first item in the list
            self.comboviewwidth = 0
            for e in self.millerarraylabels:
              self.comboviewwidth = max(self.comboviewwidth, self.MillerComboBox.fontMetrics().width( e) )
            self.MillerComboBox.view().setMinimumWidth(self.comboviewwidth)

            self.millertable.clearContents()
            self.millertable.setRowCount(len(self.array_infotpls))
            for n,millarr in enumerate(self.array_infotpls):
              for m,elm in enumerate(millarr):
                self.millertable.setItem(n, m, QTableWidgetItem(str(elm)))
            self.functionTabWidget.setDisabled(True)
            self.NewFileLoaded = False

          if self.NewHKLscenes:
            self.BinDataComboBox.clear()
            self.BinDataComboBox.addItems(["Resolution"] + self.scenearraylabels )
            self.BinDataComboBox.view().setMinimumWidth(self.comboviewwidth)
            #self.BinDataComboBox.setCurrentIndex(-1) # unselect the first item in the list
            self.NewHKLscenes = False

      except Exception as e:
        errmsg = str(e)
        if "Resource temporarily unavailable" not in errmsg:
          print( errmsg  +  traceback.format_exc(limit=10) )
        pass


  def UpdateGUI(self):
    """
    NGL_HKLviewer.viewer.scale_colors_multiplicity, False
    #NGL_HKLviewer.viewer.nth_power_scale_radii, -1
    #NGL_HKLviewer.viewer.slice_mode, False
    NGL_HKLviewer.viewer.uniform_size, False
    NGL_HKLviewer.tabulate_miller_array_ids, []
    #NGL_HKLviewer.clip_plane.l, 0
    #NGL_HKLviewer.viewer.expand_to_p1, False
    #NGL_HKLviewer.clip_plane.h, 2
    NGL_HKLviewer.viewer.sigma_radius, False
    NGL_HKLviewer.viewer.map_to_asu, False
    #NGL_HKLviewer.clip_plane.k, 0
    #NGL_HKLviewer.viewer.show_systematic_absences, False
    #NGL_HKLviewer.viewer.NGL.tooltip_alpha, 0.85
    NGL_HKLviewer.viewer.d_min, None
    NGL_HKLviewer.tooltips_in_script, False
    NGL_HKLviewer.using_space_subgroup, False
    #NGL_HKLviewer.viewer.expand_anomalous, False
    #NGL_HKLviewer.viewer.scale, 1
    #NGL_HKLviewer.filename, C:/Users/oeffner/Buser/NGL_HKLviewer/mymtz.mtz
    NGL_HKLviewer.viewer.scale_radii_multiplicity, False
    NGL_HKLviewer.viewer.data, None
    #NGL_HKLviewer.viewer.NGL.mouse_sensitivity, 0.2
    NGL_HKLviewer.action, *'is_running'
    #NGL_HKLviewer.clip_plane.angle_around_vector, 0
    #NGL_HKLviewer.viewer.slice_index, 0
    NGL_HKLviewer.shape_primitive, *'spheres'
    #NGL_HKLviewer.nbins, 6
    #NGL_HKLviewer.clip_plane.hkldist, 0
    NGL_HKLviewer.viewer.phase_color, False
    NGL_HKLviewer.viewer.show_labels, True
    NGL_HKLviewer.viewer.labels, None
    NGL_HKLviewer.viewer.show_anomalous_pairs, False
    NGL_HKLviewer.miller_array_operations,
    NGL_HKLviewer.viewer.symmetry_file, None
    #NGL_HKLviewer.spacegroup_choice, None
    #NGL_HKLviewer.clip_plane.is_parallel, False
    NGL_HKLviewer.scene_bin_thresholds, None
    NGL_HKLviewer.viewer.scene_id, 3
    NGL_HKLviewer.viewer.color_scheme, *rainbow
    NGL_HKLviewer.viewer.black_background, True
    NGL_HKLviewer.viewer.show_axes, True
    NGL_HKLviewer.bin_scene_label, Resolution
    NGL_HKLviewer.merge_data, False
    NGL_HKLviewer.mouse_moved, True
    #NGL_HKLviewer.viewer.show_missing, False
    NGL_HKLviewer.viewer.slice_axis, *h
    NGL_HKLviewer.viewer.keep_constant_scale, True
    NGL_HKLviewer.viewer.show_data_over_sigma, False
    NGL_HKLviewer.viewer.NGL.camera_type, *orthographic
    NGL_HKLviewer.clip_plane.clipwidth, None
    NGL_HKLviewer.viewer.sqrt_scale_colors, False
    NGL_HKLviewer.viewer.sigma_color, False
    NGL_HKLviewer.clip_plane.is_real_space_frac_vec, False
    NGL_HKLviewer.viewer.NGL.fixorientation, False
    NGL_HKLviewer.viewer.inbrowser, True
    NGL_HKLviewer.viewer.show_only_missing, False
    NGL_HKLviewer.viewer.NGL.bin_opacities, [(1.0, 0), (1.0, 1), (1.0, 2)]
    NGL_HKLviewer.clip_plane.bequiet, False
    """
    self.unfeedback = True
    self.power_scale_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.viewer.nth_power_scale_radii'])

    self.ManualPowerScalecheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.nth_power_scale_radii'] == -1 )
    self.radii_scale_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.viewer.scale'])
    self.showslicecheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.slice_mode'])
    self.hvec_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.clip_plane.h'])
    self.kvec_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.clip_plane.k'])
    self.lvec_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.clip_plane.l'])
    self.expandP1checkbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.expand_to_p1'])
    self.expandAnomalouscheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.expand_anomalous'])
    self.sysabsentcheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.show_systematic_absences'])
    self.ttipalpha_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.viewer.NGL.tooltip_alpha'])
    #self.HKLnameedit.setText( self.currentphilstringdict['NGL_HKLviewer.filename'])
    #self.setWindowTitle("HKL-viewer: " + self.currentphilstringdict['NGL_HKLviewer.filename'])
    self.mousemoveslider.setValue( 2000*self.currentphilstringdict['NGL_HKLviewer.viewer.NGL.mouse_sensitivity'])
    #self.rotavecangle_slider.setValue( self.currentphilstringdict['NGL_HKLviewer.clip_plane.angle_around_vector'])
    self.rotavecangle_labeltxt.setText("Angle rotated: %2.f" %self.currentphilstringdict['NGL_HKLviewer.clip_plane.angle_around_vector'])

    self.sliceindexspinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.viewer.slice_index'])
    self.Nbins_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.nbins'])
    if self.currentphilstringdict['NGL_HKLviewer.spacegroup_choice']:
      self.SpaceGroupComboBox.setCurrentIndex(  self.currentphilstringdict['NGL_HKLviewer.spacegroup_choice'] )
    self.clipParallelBtn.setChecked( self.currentphilstringdict['NGL_HKLviewer.clip_plane.is_parallel'])
    self.missingcheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.show_missing'])
    axidx = -1
    for axidx,c in enumerate(self.sliceaxis.values()):
      if c in self.currentphilstringdict['NGL_HKLviewer.viewer.slice_axis']:
        break

    self.SliceLabelComboBox.setCurrentIndex( axidx )
    self.cameraPerspectCheckBox.setChecked( "perspective" in self.currentphilstringdict['NGL_HKLviewer.viewer.NGL.camera_type'])
    self.ClipPlaneChkBox.setChecked( self.currentphilstringdict['NGL_HKLviewer.clip_plane.clipwidth'] != None )
    if self.currentphilstringdict['NGL_HKLviewer.clip_plane.clipwidth']:
      self.clipwidth_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.clip_plane.clipwidth'])
    self.hkldist_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.clip_plane.hkldist'])
    self.realspacevecBtn.setChecked( "realspace" in self.currentphilstringdict['NGL_HKLviewer.clip_plane.fractional_vector'])
    self.recipvecBtn.setChecked( "reciprocal" in self.currentphilstringdict['NGL_HKLviewer.clip_plane.fractional_vector'])
    self.clipTNCSBtn.setChecked( "tncs" in self.currentphilstringdict['NGL_HKLviewer.clip_plane.fractional_vector'])
    self.fixedorientcheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.NGL.fixorientation'])
    self.onlymissingcheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.show_only_missing'])
    self.unfeedback = False

    pass



  def onSortComboBoxSelchange(self, i):
    if i==0: # i.e. unsorted
      labels = [ ld[0] for ld in self.tabulate_miller_array ]
      self.datalst =  [ ld[1] for ld in self.tabulate_miller_array ]
      if self.millerarraytable.model():
        self.millerarraytable.model().clear()
      self.millerarraytablemodel = MillerArrayTableModel(self.datalst, labels, self)
      self.millerarraytable.setModel(self.millerarraytablemodel)
      return
    idx = i-1
    if self.millerarraytableform.sortChkbox.checkState() == Qt.Unchecked:
      self.millerarraytable.sortByColumn(idx, Qt.SortOrder.DescendingOrder)
    else:
      self.millerarraytable.sortByColumn(idx, Qt.SortOrder.AscendingOrder)


  def onSortChkbox(self):
    self.onSortComboBoxSelchange(self.millerarraytableform.SortComboBox.currentIndex() )


  def onPrecisionChanged(self, val):
    self.millerarraytablemodel.precision = val
    self.millerarraytable.resizeColumnsToContents()


  def onFinalMouseSensitivity(self):
    val = self.mousemoveslider.value()/2000.0
    self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.mouse_sensitivity = %f' %val)


  def onMouseSensitivity(self):
    val = self.mousemoveslider.value()/2000.0
    self.mousesensitxtbox.setText("%2.2f" %val )


  def onTooltipAlphaChanged(self, val):
    self.ttipalpha = val
    self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.tooltip_alpha = %f' %val)


  def onFontsizeChanged(self, val):
    font = app.font()
    font.setPointSize(val);
    app.setFont(font);
    self.settingsform.setFixedSize( self.settingsform.sizeHint() )


  def onCameraPerspect(self,val):
    if self.cameraPerspectCheckBox.isChecked():
      self.NGL_HKL_command("NGL_HKLviewer.viewer.NGL.camera_type = perspective")
    else:
      self.NGL_HKL_command("NGL_HKLviewer.viewer.NGL.camera_type = orthographic")


  def ExpandToP1(self):
    if self.expandP1checkbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.expand_to_p1 = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.expand_to_p1 = False')


  def ExpandAnomalous(self):
    if self.expandAnomalouscheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.expand_anomalous = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.expand_anomalous = False')


  def showSysAbsent(self):
    if self.sysabsentcheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_systematic_absences = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_systematic_absences = False')


  def showMissing(self):
    if self.missingcheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_missing = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_missing = False')


  def showOnlyMissing(self):
    if self.onlymissingcheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_only_missing = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_only_missing = False')


  def showSlice(self):
    if self.showslicecheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.slice_mode = True')
      if self.expandP1checkbox.isChecked():
        self.NGL_HKL_command("""NGL_HKLviewer.viewer {
                                                       expand_to_p1 = True
                                                       inbrowser = False
                                                    }
                             """)
      if self.expandAnomalouscheckbox.isChecked():
        self.NGL_HKL_command("""NGL_HKLviewer.viewer {
                                                       expand_anomalous = True
                                                       inbrowser = False
                                                     }
                             """)
      self.SliceLabelComboBox.setEnabled(True)
      self.sliceindexspinBox.setEnabled(True)
    else:
      self.NGL_HKL_command("""NGL_HKLviewer.viewer {
                                                      slice_mode = False
                                                      inbrowser = True
                                                    }
                            """)
      self.SliceLabelComboBox.setDisabled(True)
      self.sliceindexspinBox.setDisabled(True)


  def onSliceComboSelchange(self,i):
    # 4th element in each table row is the min-max span of hkls.
    rmin = self.array_infotpls[self.MillerComboBox.currentIndex()][4][0][i]
    rmax = self.array_infotpls[self.MillerComboBox.currentIndex()][4][1][i]
    self.sliceindexspinBox.setRange(rmin, rmax)
    self.NGL_HKL_command("NGL_HKLviewer.viewer.slice_axis = %s" % self.sliceaxis.values()[i] )


  def onSliceIndexChanged(self, val):
    self.sliceindex = val
    self.NGL_HKL_command("NGL_HKLviewer.viewer.slice_index = %d" %self.sliceindex)


  def onBindataComboSelchange(self,i):
    if self.BinDataComboBox.currentText():
      if self.BinDataComboBox.currentIndex() > 0:
        bin_scene_label = str(self.BinDataComboBox.currentIndex()-1)
      else:
        bin_scene_label = "Resolution"
      self.NGL_HKL_command("NGL_HKLviewer.bin_scene_label = %s" % bin_scene_label )
      bin_opacitieslst = []
      for i in range(self.nbins):
        bin_opacitieslst.append((1.0, i)) #   ("1.0, %d" %i)
      self.bin_opacities = str(bin_opacitieslst)
      self.OpaqueAllCheckbox.setCheckState(Qt.Checked)
      #self.OpaqueAllCheckbox.setTristate(false)
      #self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.bin_opacities = "%s"' %self.bin_opacities)


  def update_table_opacities(self, allalpha=None):
    bin_opacitieslst = eval(self.bin_opacities)
    self.binstable_isready = False
    for binopacity in bin_opacitieslst:
      if not allalpha:
        alpha = binopacity[0]  #float(binopacity.split(",")[0])
      else:
        alpha = allalpha
      bin = binopacity[1]  #int(binopacity.split(",")[1])
      item = QTableWidgetItem()
      item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
      if alpha == 0.0:
        item.setCheckState(Qt.Unchecked)
      if alpha > 0.0 and alpha < 1.0:
        item.setCheckState(Qt.PartiallyChecked)
      if alpha == 1.0:
        item.setCheckState(Qt.Checked)
      item.setText(str(alpha))
      item.setFlags(item.flags() ^ Qt.ItemIsEditable)
      item.setFlags(item.flags() ^ Qt.ItemIsSelectable )
      self.binstable.setItem(bin, 3, item)
    self.binstable_isready = True


  def SetAllOpaqueCheckboxes(self):
    if self.binstableitemchanges:
      return
    bin_opacitieslst = eval(self.bin_opacities)
    nbins = len(bin_opacitieslst)
    sum = 0
    for binopacity in bin_opacitieslst:
      sum += binopacity[0]  # float(binopacity.split(",")[0])
    if sum >= nbins:
      self.OpaqueAllCheckbox.setCheckState(Qt.Checked)
    if sum == 0:
      self.OpaqueAllCheckbox.setCheckState(Qt.Unchecked)
    if sum >0.0 and sum < nbins:
      self.OpaqueAllCheckbox.setCheckState(Qt.PartiallyChecked)


  def onBinsTableitemPressed(self, item):
    #print( "in itemPressed %s,  %s" %(item.text(), str( item.checkState())) )
    self.binTableCheckState = item.checkState()
    self.bintableAlpha = float(item.text())


  def onBinsTableitemClicked(self, item):
    #print( "in itemClicked  %s,  %s" %(item.text(), str( item.checkState())) )
    pass


  def onBinsTableItemChanged(self, item):
    #print( "in itemChanged %s,  %s" %(item.text(), str( item.checkState())) )
    bin = item.row()
    column = item.column()
    try:
      if not self.bin_opacities:
        return
      bin_opacitieslst = eval(self.bin_opacities)
      alpha = max(0.0, min(1.0, float(item.text()) ) ) # between 0 and 1 only
      try:
        (oldalpha, bin) = bin_opacitieslst[bin]
        if oldalpha == float(item.text()):
          if item.checkState()==Qt.Unchecked:
            alpha = 0.0
          else:
            alpha = 1.0
      except Exception as e:
        pass

      if column==3 and self.binstable_isready: # changing opacity
        bin_opacitieslst[bin] = (alpha, bin)
        self.bin_opacities = str(bin_opacitieslst)
        self.SetAllOpaqueCheckboxes()
        self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.bin_opacities = "%s"' %self.bin_opacities )
    except Exception as e:
      print( str(e)  +  traceback.format_exc(limit=10) )

  """
  def onBinsTableItemSelectionChanged(self):
    item = self.binstable.currentItem()
    #print( "in SelectionChanged %s,  %s" %(item.text(), str( item.checkState())) )
    try:
      self.currentSelectedBinsTableVal = float(item.text())
    except Exception as e:
      pass
  """

  def onOpaqueAll(self):
    self.binstableitemchanges = True
    bin_opacitieslst = eval(self.bin_opacities)
    nbins = len(bin_opacitieslst)
    bin_opacitieslst = []
    self.binstable_isready = False
    if self.OpaqueAllCheckbox.isChecked():
      for i in range(nbins):
        bin_opacitieslst.append((1.0, i))  #  ("1.0, %d" %i)
    else:
      for i in range(nbins):
        bin_opacitieslst.append((0.0, i))  #   ("0.0, %d" %i)
    self.bin_opacities = str(bin_opacitieslst)
    self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.bin_opacities = "%s"' %self.bin_opacities)
    self.binstableitemchanges = False
    self.binstable_isready = True


  """
  def onLoadFinished(self, val):
    pass
    #print("web page finished loading now")


  def onBinsTableCellentered(self, row, col):
    pass
    #print( "in Cellentered " + self.binstable.currentItem().text() )

  """

  def onNbinsChanged(self, val):
    self.nbins = val
    if not self.updatingNbins: # avoid possible endless loop to cctbx
      self.NGL_HKL_command("NGL_HKLviewer.nbins = %d" %self.nbins)


  def onRadiiScaleChanged(self, val):
    if self.unfeedback:
      return
    self.NGL_HKL_command("""
      NGL_HKLviewer.viewer {
        nth_power_scale_radii = %f
        scale = %f
      }
      """ %(self.power_scale_spinBox.value(), self.radii_scale_spinBox.value() )
    )


  def onPowerScaleChanged(self, val):
    if self.unfeedback:
      return
    self.NGL_HKL_command("""
      NGL_HKLviewer.viewer {
        nth_power_scale_radii = %f
        scale = %f
      }
      """ %(self.power_scale_spinBox.value(), self.radii_scale_spinBox.value() )
    )


  def onManualPowerScale(self):
    if self.unfeedback:
      return
    if self.ManualPowerScalecheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.nth_power_scale_radii = %f' %self.power_scale_spinBox.value())
      self.power_scale_spinBox.setEnabled(True)
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.nth_power_scale_radii = -1.0')
      self.power_scale_spinBox.setEnabled(False)
      #self.power_scale_spinBox.setValue(-1.0)


  def createExpansionBox(self):
    self.SpaceGroupComboBox = QComboBox()
    self.SpaceGroupComboBox.activated.connect(self.SpacegroupSelchange)

    self.SpacegroupLabel = QLabel()
    self.SpacegroupLabel.setText("Space Subgroups")

    self.expandP1checkbox = QCheckBox()
    self.expandP1checkbox.setText("Expand to P1")
    self.expandP1checkbox.clicked.connect(self.ExpandToP1)

    self.expandAnomalouscheckbox = QCheckBox()
    self.expandAnomalouscheckbox.setText("Show Friedel pairs")
    self.expandAnomalouscheckbox.clicked.connect(self.ExpandAnomalous)

    self.sysabsentcheckbox = QCheckBox()
    self.sysabsentcheckbox.setText("Show Systematic Absences")
    self.sysabsentcheckbox.clicked.connect(self.showSysAbsent)

    self.missingcheckbox = QCheckBox()
    self.missingcheckbox.setText("Show Missing")
    self.missingcheckbox.clicked.connect(self.showMissing)

    self.onlymissingcheckbox = QCheckBox()
    self.onlymissingcheckbox.setText("Only Show Missing")
    self.onlymissingcheckbox.clicked.connect(self.showOnlyMissing)

    self.ExpansionBox = QGroupBox("Expansions")
    layout = QGridLayout()
    layout.addWidget(self.SpacegroupLabel,           0, 0)
    layout.addWidget(self.SpaceGroupComboBox,        0, 1)
    layout.addWidget(self.expandP1checkbox,          1, 0)
    layout.addWidget(self.expandAnomalouscheckbox,   1, 1)
    layout.addWidget(self.sysabsentcheckbox,         2, 0)
    layout.addWidget(self.missingcheckbox,           3, 0)
    layout.addWidget(self.onlymissingcheckbox,       3, 1)
    layout.setRowStretch(0,0)
    layout.setRowStretch(1,0)
    layout.setRowStretch(2,0)
    layout.setRowStretch(3,1)
    self.ExpansionBox.setLayout(layout)


  def CreateSliceTabs(self):
    self.showslicecheckbox = QCheckBox()
    self.showslicecheckbox.setText("Show Slice")
    self.showslicecheckbox.clicked.connect(self.showSlice)

    self.sliceindexspinBox = QDoubleSpinBox()
    self.sliceindex = 0
    self.sliceindexspinBox.setValue(self.sliceindex)
    self.sliceindexspinBox.setDecimals(0)
    self.sliceindexspinBox.setSingleStep(1)
    self.sliceindexspinBox.setRange(0, 20)
    self.sliceindexspinBox.valueChanged.connect(self.onSliceIndexChanged)

    self.SliceLabelComboBox = QComboBox()
    self.SliceLabelComboBox.activated.connect(self.onSliceComboSelchange)
    self.sliceaxis = { 0:"h", 1:"k", 2:"l" }
    self.SliceLabelComboBox.addItems( list( self.sliceaxis.values()) )
    self.SliceLabelComboBox.setDisabled(True)
    self.sliceindexspinBox.setDisabled(True)

    self.sliceTabWidget = QTabWidget()
    tab1 = QWidget()
    layout1 = QGridLayout()
    layout1.addWidget(self.showslicecheckbox,         0, 0, 1, 1)
    layout1.addWidget(self.SliceLabelComboBox,        0, 1, 1, 1)
    layout1.addWidget(self.sliceindexspinBox,         0, 2, 1, 1)
    tab1.setLayout(layout1)

    tab2 = QWidget()
    layout2 = QGridLayout()

    self.recipvecBtn = QRadioButton(tab2)
    self.recipvecBtn.setText("as fractional values in reciprocal space")
    self.recipvecBtn.setChecked(False)
    self.recipvecBtn.clicked.connect(self.onClipPlaneChkBox)
    layout2.addWidget(self.recipvecBtn,       0, 0, 1, 2)

    self.realspacevecBtn = QRadioButton(tab2)
    self.realspacevecBtn.setText("as fractional values in real space")
    self.realspacevecBtn.setChecked(True)
    self.realspacevecBtn.clicked.connect(self.onClipPlaneChkBox)
    layout2.addWidget(self.realspacevecBtn,    1, 0, 1, 2)

    self.clipTNCSBtn = QRadioButton(tab2)
    self.clipTNCSBtn.setText("as tNCS vector")
    self.clipTNCSBtn.setChecked(False)
    self.clipTNCSBtn.clicked.connect(self.onClipPlaneChkBox)
    layout2.addWidget(self.clipTNCSBtn,    1, 2, 1, 1)

    self.hvec_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.hvec_spinBox.setValue(2.0)
    self.hvec_spinBox.setDecimals(3)
    #self.hvec_spinBox.setSingleStep(0.5)
    self.hvec_spinBox.setRange(-100.0, 100.0)
    self.hvec_spinBox.valueChanged.connect(self.onHvecChanged)
    self.hvec_Label = QLabel()
    self.hvec_Label.setText("R1")
    layout2.addWidget(self.hvec_Label,      2, 0, 1, 1)
    layout2.addWidget(self.hvec_spinBox,    3, 0, 1, 1)

    self.kvec_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.kvec_spinBox.setValue(0.0)
    self.kvec_spinBox.setDecimals(3)
    self.kvec_spinBox.setSingleStep(0.5)
    self.kvec_spinBox.setRange(-100.0, 100.0)
    self.kvec_spinBox.valueChanged.connect(self.onKvecChanged)
    self.kvec_Label = QLabel()
    self.kvec_Label.setText("R2")
    layout2.addWidget(self.kvec_Label,      2, 1, 1, 1)
    layout2.addWidget(self.kvec_spinBox,    3, 1, 1, 1)

    self.lvec_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.lvec_spinBox.setValue(0.0)
    self.lvec_spinBox.setDecimals(3)
    self.lvec_spinBox.setSingleStep(0.5)
    self.lvec_spinBox.setRange(-100.0, 100.0)
    self.lvec_spinBox.valueChanged.connect(self.onLvecChanged)
    self.lvec_Label = QLabel()
    self.lvec_Label.setText("R3")
    layout2.addWidget(self.lvec_Label,      2, 2, 1, 1)
    layout2.addWidget(self.lvec_spinBox,    3, 2, 1, 1)

    self.hkldist_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.hkldistval = 0.0
    self.hkldist_spinBox.setValue(self.hkldistval)
    self.hkldist_spinBox.setDecimals(3)
    self.hkldist_spinBox.setSingleStep(0.5)
    self.hkldist_spinBox.setRange(-100.0, 100.0)
    self.hkldist_spinBox.valueChanged.connect(self.onHKLdistChanged)
    self.hkldist_Label = QLabel()
    self.hkldist_Label.setText("Distance from Origin")
    layout2.addWidget(self.hkldist_Label,      4, 0, 1, 1)
    layout2.addWidget(self.hkldist_spinBox,    4, 1, 1, 1)

    self.clipwidth_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.clipwidth_spinBox.setValue(0.5 )
    self.clipwidth_spinBox.setDecimals(3)
    self.clipwidth_spinBox.setSingleStep(0.05)
    self.clipwidth_spinBox.setRange(0.0, 100.0)
    self.clipwidth_spinBox.valueChanged.connect(self.onClipwidthChanged)
    self.clipwidth_Label = QLabel()
    self.clipwidth_Label.setText("Clip Plane Width")
    layout2.addWidget(self.clipwidth_Label,      5, 0, 1, 1)
    layout2.addWidget(self.clipwidth_spinBox,    5, 1, 1, 1)

    self.rotavecangle_labeltxt = QLabel()
    self.rotavecangle_labeltxt.setText("Angle rotated: 0.0")
    self.rotavecangle_slider = QSlider(Qt.Horizontal)
    self.rotavecangle_slider.setMinimum(0)
    self.rotavecangle_slider.setMaximum(360)
    self.rotavecangle_slider.setValue(0)
    self.rotavecangle_slider.setSingleStep(5)
    self.rotavecangle_slider.sliderReleased.connect(self.onFinalRotaVecAngle)
    self.rotavecangle_slider.valueChanged.connect(self.onRotaVecAngle)
    layout2.addWidget(self.rotavecangle_labeltxt,  6, 0, 1, 1)
    layout2.addWidget(self.rotavecangle_slider,    6, 1, 1, 2)


    self.ClipBox = QGroupBox("Specify vector components (R1,R2,R3)")
    self.ClipBox.setLayout(layout2)

    layout3 = QGridLayout()
    self.ClipPlaneChkBox = QCheckBox(self.sliceTabWidget)
    self.ClipPlaneChkBox.setText("Slice reflections with a clip plane oriented")
    #self.ClipPlaneChkBox.clicked.connect(self.onClipPlaneChkBox)
    self.ClipPlaneChkBox.stateChanged.connect(self.onClipPlaneChkBox)
    self.clipLabel = QLabel()

    self.clipParallelBtn = QRadioButton(tab2)
    self.clipParallelBtn.setText("parallel to vector below")
    self.clipParallelBtn.setChecked(False)
    self.clipParallelBtn.clicked.connect(self.onClipPlaneChkBox)

    self.clipNormalBtn = QRadioButton(tab2)
    self.clipNormalBtn.setText("perpendicular to vector below")
    self.clipNormalBtn.setChecked(True)
    self.clipNormalBtn.clicked.connect(self.onClipPlaneChkBox)

    layout3.addWidget(self.ClipPlaneChkBox,  0, 0)
    layout3.addWidget(self.clipParallelBtn,  1, 0)
    layout3.addWidget(self.clipNormalBtn,    1, 1)
    layout3.addWidget(self.ClipBox,          2, 0, 1, 2)
    tab2.setLayout(layout3)

    self.sliceTabWidget.addTab(tab1, "Explicit Slicing")
    self.sliceTabWidget.addTab(tab2, "Clip Plane Slicing")
    self.ClipBox.setDisabled(True)
    self.clipParallelBtn.setChecked(True)
    self.clipNormalBtn.setDisabled(True)
    self.clipParallelBtn.setDisabled(True)


  def onClipPlaneChkBox(self):
    if self.unfeedback:
      return
    if self.ClipPlaneChkBox.isChecked():
      self.clipNormalBtn.setDisabled(False)
      self.clipParallelBtn.setDisabled(False)
      self.ClipBox.setDisabled(False)
      if len(self.tncsvec):
        self.clipTNCSBtn.setDisabled(False)
        self.clipwidth_spinBox.setValue(4)
      if self.realspacevecBtn.isChecked(): fracvectype = "realspace"
      if self.recipvecBtn.isChecked(): fracvectype = "reciprocal"
      if self.clipTNCSBtn.isChecked(): fracvectype = "tncs"
      if self.clipTNCSBtn.isChecked() and len(self.tncsvec):
        fracvectype = "tncs"
        self.hvec_spinBox.setValue(self.tncsvec[0])
        self.kvec_spinBox.setValue(self.tncsvec[1])
        self.lvec_spinBox.setValue(self.tncsvec[2])
        self.hvec_spinBox.setDisabled(True)
        self.kvec_spinBox.setDisabled(True)
        self.lvec_spinBox.setDisabled(True)
      else:
        self.hvec_spinBox.setDisabled(False)
        self.kvec_spinBox.setDisabled(False)
        self.lvec_spinBox.setDisabled(False)
      philstr = """NGL_HKLviewer.clip_plane {
  h = %s
  k = %s
  l = %s
  hkldist = %s
  clipwidth = %s
  is_parallel = %s
  fractional_vector = %s
}
  NGL_HKLviewer.viewer.NGL.fixorientation = %s
        """ %(self.hvec_spinBox.value(), self.kvec_spinBox.value(), self.lvec_spinBox.value(),\
              self.hkldistval, self.clipwidth_spinBox.value(), \
              str(self.clipParallelBtn.isChecked()), fracvectype, \
              str(self.fixedorientcheckbox.isChecked()) )

      self.NGL_HKL_command(philstr)
    else:
      self.ClipBox.setDisabled(True)
      self.clipNormalBtn.setDisabled(True)
      self.clipParallelBtn.setDisabled(True)
      self.clipTNCSBtn.setDisabled(True)
      self.NGL_HKL_command("NGL_HKLviewer.clip_plane.clipwidth = None")



  def onFinalRotaVecAngle(self):
    if self.unfeedback:
      return
    val = self.rotavecangle_slider.value()
    self.NGL_HKL_command("""NGL_HKLviewer.clip_plane {
    angle_around_vector = %f
    bequiet = False
}""" %val)


  def onRotaVecAngle(self):
    if self.unfeedback:
      return
    val = self.rotavecangle_slider.value()
    self.NGL_HKL_command("""NGL_HKLviewer.clip_plane {
    angle_around_vector = %f
    bequiet = True
}""" %val)


  def onClipwidthChanged(self, val):
    if not self.unfeedback:
      self.NGL_HKL_command("NGL_HKLviewer.clip_plane.clipwidth = %f" %self.clipwidth_spinBox.value())


  def onHKLdistChanged(self, val):
    self.hkldistval = val
    if not self.unfeedback:
      self.NGL_HKL_command("NGL_HKLviewer.clip_plane.hkldist = %f" %self.hkldistval)


  def onHvecChanged(self, val):
    if not self.unfeedback:
      self.NGL_HKL_command("NGL_HKLviewer.clip_plane.h = %f" %self.hvec_spinBox.value())


  def onKvecChanged(self, val):
    if not self.unfeedback:
      self.NGL_HKL_command("NGL_HKLviewer.clip_plane.k = %f" %self.kvec_spinBox.value())


  def onLvecChanged(self, val):
    if not self.unfeedback:
      self.NGL_HKL_command("NGL_HKLviewer.clip_plane.l = %f" %self.lvec_spinBox.value())


  def onFixedorient(self):
    if not self.unfeedback:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.fixorientation = %s' \
                                    %str(self.fixedorientcheckbox.isChecked()))


  def onMillerTableCellPressed(self, row, col):
    #print( "in millertable CellPressed " + self.millertable.currentItem().text() )
    if self.millertable.mousebutton == Qt.RightButton:
      self.MillerTableContextMenuHandler(QCursor.pos(), row)
    if self.millertable.mousebutton == QEvent.MouseButtonDblClick:
      # quickly display data with a double click
      for i,scenelabel in enumerate(self.scenearraylabels):
        if self.millerarraylabels[row] == scenelabel:
          self.DisplayData(i)


  def onMillerTableitemSelectionChanged(self):
    self.millertable.selectedrows = list(set([ e.row() for e in self.millertable.selectedItems() ]))


  def MillerTableContextMenuHandler(self, pos, row):
    self.millertablemenu.clear()
    # Tag menu items with data being int or a (string, int) tuple.
    # These are being checked for in onMillerTableMenuAction() and appropriate
    # action taken
    for i,scenelabel in enumerate(self.scenearraylabels):
      if self.millerarraylabels[row] == scenelabel or self.millerarraylabels[row] + " + " in scenelabel:
        #print(i, scenelabel)
        myqa = QAction("Display %s data" %scenelabel, self, triggered=self.testaction)
        myqa.setData(i)
        self.millertablemenu.addAction(myqa)
    myqa = QAction("Make new data as a function of this data...", self, triggered=self.testaction)
    myqa.setData( ("newdata_1", row ))
    self.millertablemenu.addAction(myqa)
    myqa = QAction("Make new data as a function of this data and another data set...", self, triggered=self.testaction)
    myqa.setData( ("newdata_2", row ))
    self.millertablemenu.addAction(myqa)
    #myqa = QAction("Show a table of this data set...", self, triggered=self.testaction)
    #myqa.setData( ("tabulate_data", row ))
    if len(self.millertable.selectedrows) > 0:
      arraystr = ""
      for i,r in enumerate(self.millertable.selectedrows):
        arraystr += self.millerarraylabels[r]
        if i < len(self.millertable.selectedrows)-1:
          arraystr += " and "
      myqa = QAction("Show a table of %s data ..." %arraystr, self, triggered=self.testaction)
      myqa.setData( ("tabulate_data", self.millertable.selectedrows ))
      self.millertablemenu.addAction(myqa)
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.millertablemenu.exec_(QCursor.pos())


  def onMillerTableMenuAction(self, action):
    data = action.data()
    # depending on what menu item the user clicked data is either an int or a (string, int) tuple
    if data is not None:
      if type(data) is int:
        idx = data
        self.DisplayData(idx)
      else:
        (strval, idx) = data
        self.operate_arrayidx1 = idx
        if strval=="newdata_1":
          self.operationlabeltxt.setText("Enter a python expression of " + self.millerarraylabels[idx] + " 'data' and or 'sigma' variable")
          self.MillerLabel1.setDisabled(True)
          self.MillerLabel2.setDisabled(True)
          self.MillerComboBox.setDisabled(True)
          self.MillerLabel3.setText("Example: 'newdata = data / flex.sqrt(sigmas); newsigmas= - 42*sigmas' ")
          self.operate_arrayidx2 = None
          self.makenewdataform.show()
        if strval=="newdata_2":
          self.operationlabeltxt.setText("Enter a python expression of " + self.millerarraylabels[idx] + " 'data1' and or 'sigma1' variable")
          self.MillerLabel1.setEnabled(True)
          self.MillerLabel2.setEnabled(True)
          self.MillerComboBox.setEnabled(True)
          self.MillerLabel3.setText("Example: 'newdata = data1 - flex.pow(data2); newsigmas= sigmas1 - data2 / sigmas1' ")
          self.makenewdataform.show()
        if strval=="tabulate_data":
          self.NGL_HKL_command('NGL_HKLviewer.tabulate_miller_array_ids = "%s"' %str(idx))


  def DisplayData(self, idx):
    self.NGL_HKL_command("NGL_HKLviewer.viewer.scene_id = %d" %idx)
    if self.fileisvalid:
      self.functionTabWidget.setEnabled(True)
      self.expandAnomalouscheckbox.setEnabled(True)
      self.expandP1checkbox.setEnabled(True)
      # don' allow anomalous expansion for data that's already anomalous
      arrayinfo = self.array_infotpls[idx]
      isanomalous = arrayinfo[-1]
      spacegroup = arrayinfo[2]
      label = arrayinfo[0]
      if isanomalous:
        self.expandAnomalouscheckbox.setDisabled(True)
      if spacegroup=='P 1 (No. 1)':
        self.expandP1checkbox.setDisabled(True)
    else:
      self.functionTabWidget.setDisabled(True)
    self.SpaceGroupComboBox.clear()
    self.SpaceGroupComboBox.addItems( list(self.spacegroups.values() ))


  def onMakeNewData(self):
    mtpl = (self.operationtxtbox.text(), self.newlabeltxtbox.text() ,
              self.operate_arrayidx1, self.operate_arrayidx2 )
    self.NGL_HKL_command('NGL_HKLviewer.miller_array_operations = "[ %s ]"' %str(mtpl) )
    self.makenewdataform.accept()


  def onMillerComboSelchange(self, i):
    self.operate_arrayidx2 = i


  def testaction(self):
    pass


  def createFileInfoBox(self):
    self.datasetLabel = QLabel()
    self.datasetLabel.setText("Display a data set with a double-click or right-click it for more options.")
    self.FileInfoBox = QGroupBox("Reflection File Information")
    self.HKLnameedit = QLineEdit('')
    self.HKLnameedit.setReadOnly(True)
    self.textInfo = QTextEdit()
    self.textInfo.setLineWrapMode(QTextEdit.NoWrap)
    self.textInfo.setReadOnly(True)
    layout = QGridLayout()
    layout.addWidget(self.openFileNameButton,     0, 0, 1, 2)
    if self.devmode:
      layout.addWidget(self.debugbutton,            0, 2, 1, 1)
    layout.addWidget(self.HKLnameedit,            1, 0, 1, 3)
    layout.addWidget(self.datasetLabel,           2, 0, 1, 3)
    layout.addWidget(self.millertable,            3, 0, 1, 3)
    layout.addWidget(self.textInfo,               4, 0, 1, 3)
    #layout.setColumnStretch(1, 2)
    self.FileInfoBox.setLayout(layout)


  def createRadiiScaleGroupBox(self):
    self.RadiiScaleGroupBox = QGroupBox("Radii Size of HKL Spheres")

    self.ManualPowerScalecheckbox = QCheckBox()
    self.ManualPowerScalecheckbox.setText("Manual Power Scaling of Sphere Radii")
    self.ManualPowerScalecheckbox.clicked.connect(self.onManualPowerScale)

    self.power_scale_spinBox = QDoubleSpinBox(self.RadiiScaleGroupBox)
    self.power_scale_spinBox.setValue(0.33)
    self.power_scale_spinBox.setDecimals(2)
    self.power_scale_spinBox.setSingleStep(0.05)
    self.power_scale_spinBox.setRange(0.0, 1.0)
    self.power_scale_spinBox.valueChanged.connect(self.onPowerScaleChanged)
    self.power_scale_spinBox.setEnabled(False)
    self.powerscaleLabel = QLabel()
    self.powerscaleLabel.setText("Power scale Factor")

    self.radii_scale_spinBox = QDoubleSpinBox(self.RadiiScaleGroupBox)
    self.radii_scale_spinBox.setValue(1.0)
    self.radii_scale_spinBox.setDecimals(1)
    self.radii_scale_spinBox.setSingleStep(0.1)
    self.radii_scale_spinBox.setRange(0.2, 2.0)
    self.radii_scale_spinBox.valueChanged.connect(self.onRadiiScaleChanged)
    self.radiiscaleLabel = QLabel()
    self.radiiscaleLabel.setText("Linear Scale Factor")

    layout = QGridLayout()
    layout.addWidget(self.ManualPowerScalecheckbox, 1, 0, 1, 2)
    layout.addWidget(self.powerscaleLabel,          2, 0, 1, 2)
    layout.addWidget(self.power_scale_spinBox,      2, 1, 1, 2)
    layout.addWidget(self.radiiscaleLabel,          3, 0, 1, 2)
    layout.addWidget(self.radii_scale_spinBox,      3, 1, 1, 2)
    layout.setColumnStretch (0, 1)
    layout.setColumnStretch (1 ,0)
    self.RadiiScaleGroupBox.setLayout(layout)


  def createBinsBox(self):
    self.binstable = QTableWidget(0, 4)
    self.binstable_isready = False
    labels = ["no. of HKLs", "lower bin value", "upper bin value", "opacity"]
    self.binstable.setHorizontalHeaderLabels(labels)
    self.binstable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    self.bindata_labeltxt = QLabel()
    self.bindata_labeltxt.setText("Binning according to:")
    self.Nbins_spinBox = QSpinBox()
    self.Nbins_spinBox.setSingleStep(1)
    self.Nbins_spinBox.setRange(1, 40)
    self.Nbins_spinBox.valueChanged.connect(self.onNbinsChanged)
    self.Nbins_labeltxt = QLabel()
    self.Nbins_labeltxt.setText("Number of bins:")

    self.OpaqueAllCheckbox = QCheckBox()
    #self.OpaqueAllCheckbox.setTristate()
    self.OpaqueAllCheckbox.setCheckState(Qt.Checked)
    self.OpaqueAllCheckbox.setText("Show all data in bins")
    self.OpaqueAllCheckbox.clicked.connect(self.onOpaqueAll)

    self.binstable.itemChanged.connect(self.onBinsTableItemChanged  )
    self.binstable.itemClicked.connect(self.onBinsTableitemClicked  )
    self.binstable.itemPressed.connect(self.onBinsTableitemPressed  )
    #self.binstable.itemSelectionChanged.connect(self.onBinsTableItemSelectionChanged  )

    self.BinDataComboBox = QComboBox()
    self.BinDataComboBox.activated.connect(self.onBindataComboSelchange)
    self.BinsGroupBox = QGroupBox("Bins")
    layout = QGridLayout()
    layout.addWidget(self.bindata_labeltxt, 0, 0)
    layout.addWidget(self.BinDataComboBox, 0, 1)
    layout.addWidget(self.Nbins_labeltxt, 0, 2)
    layout.addWidget(self.Nbins_spinBox, 0, 3)
    layout.addWidget(self.OpaqueAllCheckbox, 1, 2)
    layout.addWidget(self.binstable, 2, 0, 1, 4)
    layout.setColumnStretch(0, 0)
    layout.setColumnStretch(1, 2)
    layout.setColumnStretch(3, 1)
    self.BinsGroupBox.setLayout(layout)


  def DebugInteractively(self):
    import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )


  def CreateFunctionTabs(self):
    self.functionTabWidget = QTabWidget()
    tab1 = QWidget()
    layout1 = QGridLayout()
    layout1.addWidget(self.ExpansionBox,     0, 0)
    layout1.setRowStretch (0 ,0)
    tab1.setLayout(layout1)

    tab2 = QWidget()
    layout2 = QGridLayout()

    self.fixedorientcheckbox = QCheckBox(self.sliceTabWidget)
    self.fixedorientcheckbox.setText("Fix orientation but allow zoom and translation")
    self.fixedorientcheckbox.clicked.connect(self.onFixedorient)
    layout2.addWidget(self.fixedorientcheckbox,   0, 0)

    layout2.addWidget(self.sliceTabWidget,     1, 0)
    tab2.setLayout(layout2)

    tab3 = QWidget()
    layout3 = QGridLayout()
    layout3.addWidget(self.RadiiScaleGroupBox,     0, 0)
    tab3.setLayout(layout3)

    tab4 = QWidget()
    layout4 = QGridLayout()
    layout4.addWidget(self.BinsGroupBox,     0, 0)
    tab4.setLayout(layout4)

    self.functionTabWidget.addTab(tab1, "Expand")
    self.functionTabWidget.addTab(tab2, "Slice")
    self.functionTabWidget.addTab(tab3, "Size")
    self.functionTabWidget.addTab(tab4, "Bins")
    self.functionTabWidget.setDisabled(True)


  def SpacegroupSelchange(self,i):
    self.NGL_HKL_command("NGL_HKLviewer.spacegroup_choice = %d" %i)


  def find_free_port(self):
    import socket
    s = socket.socket()
    s.bind(('', 0))      # Bind to a free port provided by the host.
    port = s.getsockname()[1]
    s.close()
    return port


  def LaunchCCTBXPython(self):
    self.sockport = self.find_free_port()
    self.zmq_context = zmq.Context()
    self.socket = self.zmq_context.socket(zmq.PAIR)
    self.socket.bind("tcp://127.0.0.1:%s" %self.sockport)
    try: msg = self.socket.recv(flags=zmq.NOBLOCK) #To empty the socket from previous messages
    except Exception as e: pass
    cmdargs = 'cctbx.python -i -c "from crys3d.hklview import cmdlineframes;' \
     + ' myHKLview = cmdlineframes.HKLViewFrame(useGuiSocket=%s, high_quality=True,' %self.sockport \
     + ' jscriptfname = \'%s\', ' %self.jscriptfname \
     + ' verbose=%s, UseOSBrowser=%s, htmlfname=\'%s\', handshakewait=%s )"\n'\
       %(self.verbose, str(self.UseOSbrowser), self.hklfname, self.handshakewait)
    self.cctbxproc = subprocess.Popen( cmdargs, shell=True, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr)
    """
    input = 'from crys3d.hklview import cmdlineframes;' \
     + ' myHKLview = cmdlineframes.HKLViewFrame(useGuiSocket=%s, high_quality=True,' %self.sockport \
     + ' jscriptfname = \'%s\', ' %self.jscriptfname \
     + ' verbose=%s, UseOSBrowser=%s, htmlfname=\'%s\', handshakewait=%s )\n'\
       %(self.verbose, str(self.UseOSbrowser), self.hklfname, self.handshakewait)
    self.cctbxproc = subprocess.Popen( 'cctbx.python.bat', shell=False, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr)
    self.out, self.err = self.cctbxproc.communicate(bytes(input,"utf-8") )
    """
    #self.out, self.err = self.cctbxproc.communicate()
    #time.sleep(1)


  def NGL_HKL_command(self, cmdstr):
    if sys.version_info.major==3:
      self.socket.send(bytes(cmdstr,"utf-8"))
    else:
      self.socket.send(bytes(cmdstr))



if __name__ == '__main__':
  try:
    """
    If chromium webgl error on MacOS try using commandline arguments:
    --enable-webgl-software-rendering and --ignore-gpu-blacklist
    """
    debugtrue = False
    for e in sys.argv:
      if "devmode" in e or "debug" in e:
        debugtrue = True
    if debugtrue:
      # some useful flags as per https://doc.qt.io/qt-5/qtwebengine-debugging.html
      sys.argv.append("--remote-debugging-port=9742")
      sys.argv.append("--single-process")
      sys.argv.append("--js-flags='--expose_gc'")

    app = QApplication(sys.argv)
    guiobj = NGL_HKLViewer()

    timer = QTimer()
    timer.setInterval(10)
    timer.timeout.connect(guiobj.ProcessMessages)
    timer.start()

    ret = app.exec_()
    #app.quit()

    """    # no need to delete cache for private/off-the-record profile
    guiobj.webpage.deleteLater()
    guiobj.webpage = None

    guiobj.BrowserBox.deleteLater()
    guiobj.BrowserBox.destroy()
    guiobj.deleteLater()
    guiobj.destroy()

    del QWebEngineView
    del QWebEngineProfile
    del QWebEnginePage
    #del sys.modules["QWebEngineView"]
    #del sys.modules["QWebEngineProfile"]
    #del sys.modules["QWebEnginePage"]

    del guiobj
    del app
    gc.collect()

    present = True
    while present:
      present = False
      try:
        if os.path.exists(guiobj.cpath):
          shutil.rmtree(guiobj.cpath)
        cpath2 = guiobj.cpath.replace("cache/","")
        if os.path.exists(cpath2):
          shutil.rmtree(cpath2)
      except Exception as delerr:
        time.sleep(0.2)
        present = True

    """
    sys.exit(ret)
  except Exception as e:
    print( str(e)  +  traceback.format_exc(limit=10) )
