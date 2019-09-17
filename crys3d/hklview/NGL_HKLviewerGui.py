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

from PySide2.QtCore import Qt, QTimer, QEvent
from PySide2.QtWidgets import ( QAction, QApplication, QCheckBox, QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QMenu, QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QVBoxLayout, QWidget )

from PySide2.QtGui import QColor, QFont, QCursor
from PySide2.QtWebEngineWidgets import QWebEngineView
import sys, zmq, subprocess, time, traceback



class MakeNewDataForm(QDialog):
  def __init__(self, parent=None):
    super(MakeNewDataForm, self).__init__(parent)
    self.setWindowTitle("Create New Reflection Data")
    myGroupBox = QGroupBox("Stuff")
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


class MyTableWidget(QTableWidget):
  def __init__(self, *args, **kwargs):
    QTableWidget.__init__(self, *args, **kwargs)
    self.mousebutton = None

  def mousePressEvent(self, event):
    if event.type() == QEvent.MouseButtonPress:
      self.mousebutton = None
      if event.button() == Qt.RightButton:
        print( "Right clicked")
        self.mousebutton = Qt.RightButton
      else:
        print("Left clicked")
    QTableWidget.mousePressEvent(self, event)


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
        self.UseOSbrowser = e.split("UseOSbrowser=")[1]
      if "jscriptfname" in e:
        self.jscriptfname = e.split("jscriptfname=")[1]
      if "devmode" in e:
        self.devmode = True
      if 'htmlfname' in e:
        self.hklfname = e.split("htmlfname=")[1]
      if 'handshakewait' in e:
        self.handshakewait = e.split('handshakewait=')[1]

    self.zmq_context = None

    self.originalPalette = QApplication.palette()

    self.openFileNameButton = QPushButton("Load reflection file")
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
    #self.fontspinBox.setValue(self.font.pixelSize())
    self.fontspinBox.valueChanged.connect(self.onFontsizeChanged)
    self.Fontsize_labeltxt = QLabel()
    self.Fontsize_labeltxt.setText("Font size:")

    self.cameraPerspectCheckBox = QCheckBox()
    self.cameraPerspectCheckBox.setText("Perspective camera")
    self.cameraPerspectCheckBox.clicked.connect(self.onCameraPerspect)
    self.cameraPerspectCheckBox.setCheckState(Qt.Unchecked)

    self.bufsize = 5000
    self.bufsizespinBox = QSpinBox()
    self.bufsizespinBox.setSingleStep(500)
    self.bufsizespinBox.setRange(500, 100000)
    self.bufsizespinBox.setValue(self.bufsize)
    self.bufsizespinBox.valueChanged.connect(self.onTextbuffersizeChanged)
    self.bufsize_labeltxt = QLabel()
    self.bufsize_labeltxt.setText("Text buffer size (bytes):")

    self.ttipalpha = 0.75
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


    self.HKLnameedit = QLineEdit('')
    self.HKLnameedit.setReadOnly(True)
    self.textInfo = QTextEdit()
    self.textInfo.setLineWrapMode(QTextEdit.NoWrap)
    self.textInfo.setReadOnly(True)

    labels = ["Column label", "Type", "Space group", "# HKLs", "Span of HKLs",
       "Min Max data", "Min Max sigmas", "d_min, d_max", "Symmetry unique", "Anomalous"]
    self.millertable = MyTableWidget(0, len(labels))
    self.millertable.setHorizontalHeaderLabels(labels)
    self.millertable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    # don't allow editing this table
    self.millertable.setEditTriggers(QTableWidget.NoEditTriggers)
    self.millertable.cellPressed.connect(self.onMillerTableCellPressed)
    #self.millertable.installEventFilter(self)

    self.createExpansionBox()
    self.createFileInfoBox()
    self.CreateSliceTabs()
    self.createRadiiScaleGroupBox()
    self.createBinsBox()
    self.CreateFunctionTabs()

    mainLayout = QGridLayout()
    mainLayout.addWidget(self.FileInfoBox,         0, 0)
    mainLayout.addWidget(self.functionTabWidget,   1, 0)
    mainLayout.addWidget(self.settingsbtn,         2, 0, 1, 1)

    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    if self.UseOSbrowser==False:
      self.BrowserBox = QWebEngineView()
      mainLayout.addWidget(self.BrowserBox,          0, 1, 5, 3)
      self.BrowserBox.setUrl("https://cctbx.github.io/")
      #self.BrowserBox.setUrl("https://webglreport.com/")
      #self.BrowserBox.loadFinished.connect(self.onLoadFinished)
      mainLayout.setColumnStretch(2, 1)

    mainLayout.setRowStretch(0, 1)
    mainLayout.setRowStretch(1, 0)
    mainLayout.setRowStretch(2, 1)
    #mainLayout.setRowStretch(3, 0)
    mainLayout.setColumnStretch(4, 0)
    self.setLayout(mainLayout)

    self.setWindowTitle("HKL-Viewer")
    self.cctbxproc = None
    self.LaunchCCTBXPython()
    self.out = None
    self.err = None
    self.comboviewwidth = 0
    self.hklscenes_arrays = []
    self.millerarraylabels = []
    self.scenearraylabels = []
    self.array_infotpls = []
    self.matching_arrays = []
    self.bin_infotpls = None
    self.bin_opacities= None
    self.html_url = ""
    self.spacegroups = []
    self.info = []
    self.infostr = ""
    self.fileisvalid = False
    self.NewFileLoaded = False
    self.NewMillerArray = False
    self.NewHKLscenes = False
    self.updatingNbins = False
    self.binstableitemchanges = False
    self.millertablemenu = QMenu(self)
    self.millertablemenu.triggered.connect(self.onMillerTableMenuAction)

    self.show()


  def SettingsDialog(self):
    self.settingsform.show()


  """
  def eventFilter(self, source, event):
    if (event.type() == QEvent.MouseButtonPress and
     source is self.millertable):
      if event.button() == Qt.RightButton:
        print("Right button clicked")
    #return False

    if (event.type() == QEvent.ContextMenu and
     source is self.millertable):
      self.CustomContextMenuHandler(event.pos())

      #menu = QMenu()
      #menu.addAction('Open Window')
      #if menu.exec_(event.globalPos()):
      #  item = source.itemAt(event.pos())
      #  print(item.text())

      return True
    return super(NGL_HKLViewer, self).eventFilter(source, event)
  """


  def update(self):
    if self.cctbxproc:
      if self.cctbxproc.stdout:
        print(self.cctbxproc.stdout.read().decode("utf-8"))
      if self.cctbxproc.stderr:
        print(self.cctbxproc.stderr.read().decode("utf-8"))
    if self.out:
      print(self.out.decode("utf-8"))
    if self.err:
      print(self.err.decode("utf-8"))
    if self.zmq_context:
      try:
        msg = self.socket.recv(flags=zmq.NOBLOCK) #To empty the socket from previous messages
        nan = float("nan") # workaround for "evaluating" any NaN values in the messages received
        msgstr = msg.decode()
        self.infodict = eval(msgstr)
        #print("received from cctbx: " + str(self.infodict))
        if self.infodict:

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
            self.spacegroups = self.infodict.get("spacegroups",[])
            self.SpaceGroupComboBox.clear()
            self.SpaceGroupComboBox.addItems( self.spacegroups )

          if self.infodict.get("merge_data"):
            self.mergedata = self.infodict["merge_data"]

          currentinfostr = ""
          if self.infodict.get("info"):
            currentinfostr = self.infodict.get("info",[])

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
            self.infostr = self.infostr[-self.bufsize:]
            self.textInfo.setPlainText(self.infostr)
            self.textInfo.verticalScrollBar().setValue( self.textInfo.verticalScrollBar().maximum()  )

          if (self.NewFileLoaded or self.NewMillerArray) and self.NewHKLscenes:
            #if self.mergedata == True : val = Qt.CheckState.Checked
            #if self.mergedata == None : val = Qt.CheckState.PartiallyChecked
            #if self.mergedata == False : val = Qt.CheckState.Unchecked
            #self.mergecheckbox.setCheckState(val )
            #print("got hklscenes: " + str(self.hklscenes_arrays))
            self.NewMillerArray = False

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


  def onFinalMouseSensitivity(self):
    val = self.mousemoveslider.value()/100.0
    self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.mouse_sensitivity = %f' %val)


  def onMouseSensitivity(self):
    val = self.mousemoveslider.value()/100.0
    self.mousesensitxtbox.setText("%2.2f" %val )


  def onTextbuffersizeChanged(self, val):
    self.bufsize = val


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


  def MergeData(self):
    if self.mergecheckbox.checkState()== Qt.CheckState.Checked:
      self.NGL_HKL_command('NGL_HKLviewer.mergedata = True')
    if self.mergecheckbox.checkState()== Qt.CheckState.PartiallyChecked:
      self.NGL_HKL_command('NGL_HKLviewer.mergedata = None')
    if self.mergecheckbox.checkState()== Qt.CheckState.Unchecked:
      self.NGL_HKL_command('NGL_HKLviewer.mergedata = False')


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
    self.NGL_HKL_command("NGL_HKLviewer.viewer.slice_axis = %s" % self.sliceaxis[i] )


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
      if alpha < 0.5:
        item.setCheckState(Qt.Unchecked)
      else:
        item.setCheckState(Qt.Checked)
      item.setFlags(item.flags() ^ Qt.ItemIsEditable)
      self.binstable.setItem(bin, 3, item)
    self.binstable_isready = True


  def SetOpaqueAll(self):
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


  def onBinsTableItemChanged(self, item):
    bin = item.row()
    column = item.column()
    try:
      if item.checkState()==Qt.Unchecked:
        alpha = 0.0
      else:
        alpha = 1.0
      if column==3 and self.binstable_isready: # changing opacity
        #assert (alpha <= 1.0 and alpha >= 0.0)
        bin_opacitieslst = eval(self.bin_opacities)
        bin_opacitieslst[bin] = (alpha, bin)
        self.bin_opacities = str(bin_opacitieslst)
        self.SetOpaqueAll()
        self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.bin_opacities = "%s"' %self.bin_opacities )
    except Exception as e:
      print( str(e)  +  traceback.format_exc(limit=10) )
      #self.binstable.currentItem().setText( self.currentSelectedBinsTableVal)


  def onBinsTableItemSelectionChanged(self):
    row = self.binstable.currentItem().row()
    column = self.binstable.currentItem().column()
    self.currentSelectedBinsTableVal = self.binstable.currentItem().text()
    #print( "in itemSelectionChanged " + self.currentSelectedBinsTableVal)


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


  def onBinsTableitemActivated(self, item):
    row = item.row()
    column = item.column()
    currentval = item.text()
    #print( "in itemActivated " + currentval)


  def onBinsTableCellentered(self, row, col):
    pass
    #print( "in Cellentered " + self.binstable.currentItem().text() )

  """

  def onNbinsChanged(self, val):
    self.nbins = val
    if not self.updatingNbins: # avoid possible endless loop to cctbx
      self.NGL_HKL_command("NGL_HKLviewer.nbins = %d" %self.nbins)


  def onRadiiScaleChanged(self, val):
    self.radii_scale = val
    self.NGL_HKL_command("""
      NGL_HKLviewer.viewer {
        nth_power_scale_radii = %f
        scale = %f
      }
      """ %(self.nth_power_scale, self.radii_scale)
    )


  def onPowerScaleChanged(self, val):
    self.nth_power_scale = val
    self.NGL_HKL_command("""
      NGL_HKLviewer.viewer {
        nth_power_scale_radii = %f
        scale = %f
      }
      """ %(self.nth_power_scale, self.radii_scale)
    )


  def onManualPowerScale(self):
    if self.ManualPowerScalecheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.nth_power_scale_radii = %f' %self.nth_power_scale)
      self.power_scale_spinBox.setEnabled(True)
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.nth_power_scale_radii = -1.0')
      self.power_scale_spinBox.setEnabled(False)
      self.nth_power_scale = -1.0


  def onOpenReflectionFile(self):
    options = QFileDialog.Options()
    fileName, filtr = QFileDialog.getOpenFileName(self,
            "Load reflections file",
            "",
            "All Files (*);;MTZ Files (*.mtz);;CIF (*.cif)", "", options)
    if fileName:
      self.HKLnameedit.setText(fileName)
      #self.infostr = ""
      self.textInfo.setPlainText("")
      self.fileisvalid = False
      self.NGL_HKL_command('NGL_HKLviewer.filename = "%s"' %fileName )
      self.MillerComboBox.clear()
      self.BinDataComboBox.clear()


  def createExpansionBox(self):
    self.SpaceGroupComboBox = QComboBox()
    self.SpaceGroupComboBox.activated.connect(self.SpacegroupSelchange)

    self.SpacegroupLabel = QLabel()
    self.SpacegroupLabel.setText("Space Subgroups")

    self.mergecheckbox = QCheckBox()
    self.mergecheckbox.setText("Merge data")
    #self.mergecheckbox.setTristate (True)
    self.mergecheckbox.clicked.connect(self.MergeData)

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
    #layout.addWidget(self.mergecheckbox,             1, 0)
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
    self.sliceaxis = [ "h", "k", "l" ]
    self.SliceLabelComboBox.addItems( self.sliceaxis )
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

    self.hvec_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.hvecval = 2.0
    self.hvec_spinBox.setValue(self.hvecval)
    self.hvec_spinBox.setDecimals(6)
    self.hvec_spinBox.setSingleStep(0.5)
    self.hvec_spinBox.setRange(-100.0, 10.0)
    self.hvec_spinBox.valueChanged.connect(self.onHvecChanged)
    self.hvec_Label = QLabel()
    self.hvec_Label.setText("H")
    layout2.addWidget(self.hvec_Label,      0, 0, 1, 1)
    layout2.addWidget(self.hvec_spinBox,    0, 1, 1, 1)

    self.kvec_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.kvecval = 0.0
    self.kvec_spinBox.setValue(self.kvecval)
    self.kvec_spinBox.setDecimals(6)
    self.kvec_spinBox.setSingleStep(0.5)
    self.kvec_spinBox.setRange(-100.0, 100.0)
    self.kvec_spinBox.valueChanged.connect(self.onKvecChanged)
    self.kvec_Label = QLabel()
    self.kvec_Label.setText("K")
    layout2.addWidget(self.kvec_Label,      1, 0, 1, 1)
    layout2.addWidget(self.kvec_spinBox,    1, 1, 1, 1)

    self.lvec_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.lvecval = 0.0
    self.lvec_spinBox.setValue(self.lvecval)
    self.lvec_spinBox.setDecimals(6)
    self.lvec_spinBox.setSingleStep(0.5)
    self.lvec_spinBox.setRange(-100.0, 100.0)
    self.lvec_spinBox.valueChanged.connect(self.onLvecChanged)
    self.lvec_Label = QLabel()
    self.lvec_Label.setText("L")
    layout2.addWidget(self.lvec_Label,      2, 0, 1, 1)
    layout2.addWidget(self.lvec_spinBox,    2, 1, 1, 1)

    self.hkldist_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.hkldistval = 0.0
    self.hkldist_spinBox.setValue(self.hkldistval)
    self.hkldist_spinBox.setDecimals(6)
    self.hkldist_spinBox.setSingleStep(0.5)
    self.hkldist_spinBox.setRange(-100.0, 100.0)
    self.hkldist_spinBox.valueChanged.connect(self.onHKLdistChanged)
    self.hkldist_Label = QLabel()
    self.hkldist_Label.setText("Distance from Origin")
    layout2.addWidget(self.hkldist_Label,      3, 0, 1, 1)
    layout2.addWidget(self.hkldist_spinBox,    3, 1, 1, 1)

    self.clipwidth_spinBox = QDoubleSpinBox(self.sliceTabWidget)
    self.clipwidthval = 0.5
    self.clipwidth_spinBox.setValue(self.clipwidthval )
    self.clipwidth_spinBox.setDecimals(6)
    self.clipwidth_spinBox.setSingleStep(0.05)
    self.clipwidth_spinBox.setRange(0.0, 100.0)
    self.clipwidth_spinBox.valueChanged.connect(self.onClipwidthChanged)
    self.clipwidth_Label = QLabel()
    self.clipwidth_Label.setText("Clip Plane Width")
    layout2.addWidget(self.clipwidth_Label,      4, 0, 1, 1)
    layout2.addWidget(self.clipwidth_spinBox,    4, 1, 1, 1)

    self.ClipBox = QGroupBox("Normal Vector to Clip Plane")
    self.ClipBox.setLayout(layout2)

    layout3 = QGridLayout()
    self.ClipPlaneChkBox = QCheckBox(self.sliceTabWidget)
    self.ClipPlaneChkBox.setText("Slice reflections with a clip plane oriented")

    self.ClipPlaneChkBox.clicked.connect(self.onClipPlaneChkBox)
    self.clipLabel = QLabel()
    self.clipLabel.setText(" perpendicularly to reciprocal vector below pointing out.")

    layout3.addWidget(self.ClipPlaneChkBox, 0, 0)
    layout3.addWidget(self.clipLabel,       1, 0)
    layout3.addWidget(self.ClipBox,         2, 0)
    tab2.setLayout(layout3)
    self.sliceTabWidget.addTab(tab1, "Explicit Slicing")
    self.sliceTabWidget.addTab(tab2, "Clip Plane Slicing")
    self.ClipBox.setDisabled(True)


  def onClipPlaneChkBox(self):
    if self.ClipPlaneChkBox.isChecked():
      self.ClipBox.setDisabled(False)
      philstr = """NGL_HKLviewer.normal_clip_plane {
  h = %s
  k = %s
  l = %s
  hkldist = %s
  clipwidth = %s
}
  NGL_HKLviewer.viewer.NGL.fixorientation = %s
      """ %(self.hvecval, self.kvecval, self.lvecval, self.hkldistval, self.clipwidthval, \
                              str(self.fixedorientcheckbox.isChecked()) )
      self.NGL_HKL_command(philstr)
    else:
      self.ClipBox.setDisabled(True)
      self.NGL_HKL_command("NGL_HKLviewer.normal_clip_plane.clipwidth = None")


  def onClipwidthChanged(self, val):
    self.clipwidthval = val
    self.NGL_HKL_command("NGL_HKLviewer.normal_clip_plane.clipwidth = %f" %self.clipwidthval)


  def onHKLdistChanged(self, val):
    self.hkldistval = val
    self.NGL_HKL_command("NGL_HKLviewer.normal_clip_plane.hkldist = %f" %self.hkldistval)


  def onHvecChanged(self, val):
    self.hvecval = val
    self.NGL_HKL_command("NGL_HKLviewer.normal_clip_plane.h = %f" %self.hvecval)


  def onKvecChanged(self, val):
    self.kvecval = val
    self.NGL_HKL_command("NGL_HKLviewer.normal_clip_plane.k = %f" %self.kvecval)


  def onLvecChanged(self, val):
    self.lvecval = val
    self.NGL_HKL_command("NGL_HKLviewer.normal_clip_plane.l = %f" %self.lvecval)


  def onFixedorient(self):
    self.NGL_HKL_command('NGL_HKLviewer.viewer.NGL.fixorientation = %s' \
                                    %str(self.fixedorientcheckbox.isChecked()))


  def onMillerTableCellPressed(self, row, col):
    #print( "in millertable CellPressed " + self.millertable.currentItem().text() )
    if self.millertable.mousebutton == Qt.RightButton:
      self.MillerTableContextMenuHandler(QCursor.pos(), row)


  def MillerTableContextMenuHandler(self, pos, row):
    self.millertablemenu.clear()

    for i,scenelabel in enumerate(self.scenearraylabels):
      #[e for e in self.millerarraylabels if e in scenelabel]
      if self.millerarraylabels[row] == scenelabel or self.millerarraylabels[row] + " + " in scenelabel:
        print(i, scenelabel)
        myqa = QAction("Display %s data" %scenelabel, self, triggered=self.testaction)
        myqa.setData(i)
        self.millertablemenu.addAction(myqa)
    myqa = QAction("Make new data as a function of this data...", self, triggered=self.testaction)
    myqa.setData( ("newdata_1", row ))
    self.millertablemenu.addAction(myqa)
    myqa = QAction("Make new data as a function of this data and another data set...", self, triggered=self.testaction)
    myqa.setData( ("newdata_2", row ))
    self.millertablemenu.addAction(myqa)
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.millertablemenu.exec_(QCursor.pos())


  def onMillerTableMenuAction(self, action):
    data = action.data()
    if data is not None:
      if type(data) is int: # as set above when matching millerarraylabels with scenearraylabels
        idx = data
        self.NGL_HKL_command("NGL_HKLviewer.viewer.scene_id = %d" %idx)
        if self.fileisvalid:
          self.functionTabWidget.setEnabled(True)
          self.expandAnomalouscheckbox.setEnabled(True)
          # don' allow anomalous expansion for data that's already anomalous
          for arrayinfo in self.array_infotpls:
            isanomalous = arrayinfo[-1]
            label = arrayinfo[0]
            if isanomalous and label == self.MillerComboBox.currentText()[: len(label) ]:
              self.expandAnomalouscheckbox.setDisabled(True)
        else:
          self.functionTabWidget.setDisabled(True)
        self.SpaceGroupComboBox.clear()
        self.SpaceGroupComboBox.addItems( self.spacegroups )

      else:
        (strval, idx) = data
        self.operate_arrayidx1 = idx
        if strval=="newdata_1":
          self.operationlabeltxt.setText("Enter a python expression of " + self.millerarraylabels[idx] + " 'data' and or 'sigma' variable")
          self.MillerLabel1.setDisabled(True)
          self.MillerLabel2.setDisabled(True)
          self.MillerComboBox.setDisabled(True)
          self.MillerLabel3.setText("Example: 'newdata=data/sigmas; newsigmas= -42*sigmas' ")
          self.operate_arrayidx2 = None
        else:
          self.operationlabeltxt.setText("Enter a python expression of " + self.millerarraylabels[idx] + " 'data1' and or 'sigma1' variable")
          self.MillerLabel1.setEnabled(True)
          self.MillerLabel2.setEnabled(True)
          self.MillerComboBox.setEnabled(True)
          self.MillerLabel3.setText("Example: 'newdata=data1+data2; newsigmas= sigmas1 - data2/sigmas1' ")

        self.makenewdataform.show()


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
    self.datasetLabel.setText("Data sets:")

    self.FileInfoBox = QGroupBox("Reflection File Information")
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
    self.nth_power_scale = 0.33
    self.power_scale_spinBox.setValue(self.nth_power_scale)
    self.power_scale_spinBox.setDecimals(2)
    self.power_scale_spinBox.setSingleStep(0.05)
    self.power_scale_spinBox.setRange(0.0, 1.0)
    self.power_scale_spinBox.valueChanged.connect(self.onPowerScaleChanged)
    self.power_scale_spinBox.setEnabled(False)
    self.powerscaleLabel = QLabel()
    self.powerscaleLabel.setText("Power scale Factor")

    self.radii_scale_spinBox = QDoubleSpinBox(self.RadiiScaleGroupBox)
    self.radii_scale = 1.5
    self.radii_scale_spinBox.setValue(self.radii_scale)
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
    self.binstable.itemSelectionChanged.connect(self.onBinsTableItemSelectionChanged  )

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
    #time.sleep(1)


  def NGL_HKL_command(self, cmdstr):
    #print("sending:\n" + cmdstr)
    self.socket.send(bytes(cmdstr,"utf-8"))



if __name__ == '__main__':
  try:
    app = QApplication(sys.argv)

    guiobj = NGL_HKLViewer()

    timer = QTimer()
    timer.setInterval(20)
    timer.timeout.connect(guiobj.update)
    timer.start()

    if guiobj.cctbxproc:
      guiobj.cctbxproc.terminate()
    sys.exit(app.exec_())
  except Exception as e:
    print( str(e)  +  traceback.format_exc(limit=10) )
