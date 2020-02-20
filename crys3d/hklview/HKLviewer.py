# -*- coding: utf-8 -*-
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
from __future__ import absolute_import, division, print_function

from PySide2.QtCore import Qt, QEvent, QSize, QSettings, QTimer
from PySide2.QtWidgets import (  QAction, QApplication, QCheckBox,
        QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMainWindow, QMenu, QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QVBoxLayout, QWidget )

from PySide2.QtGui import QColor, QFont, QCursor
from PySide2.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )
import sys, zmq, subprocess, time, traceback, zlib, io, os


try: # if invoked by cctbx.python or some such
  from crys3d.hklview import HKLviewerGui, QtChromiumCheck
  from crys3d.hklview.helpers import MillerArrayTableView, MillerArrayTableForm, MillerArrayTableModel
except Exception as e: # if invoked by a generic python that doesn't know cctbx modules
  import HKLviewerGui, QtChromiumCheck
  from helpers import MillerArrayTableView, MillerArrayTableForm, MillerArrayTableModel

class MakeNewDataForm(QDialog):
  def __init__(self, parent=None):
    super(MakeNewDataForm, self).__init__(parent.window)
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
    super(SettingsForm, self).__init__(parent.window)
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

    layout.addWidget(parent.ttiplabeltxt,            4, 0, 1, 1)
    layout.addWidget(parent.ttipClickradio,          4, 1, 1, 1)
    layout.addWidget(parent.ttipHoverradio,          4, 2, 1, 1)
    layout.addWidget(parent.ttipalphalabeltxt,       4, 3, 1, 1)
    layout.addWidget(parent.ttipalpha_spinBox,       4, 4, 1, 1)

    layout.setRowStretch (0, 1)
    layout.setRowStretch (1 ,0)
    myGroupBox.setLayout(layout)
    mainLayout = QGridLayout()
    mainLayout.addWidget(myGroupBox,     0, 0)
    self.setLayout(mainLayout)
    self.setFixedSize( self.sizeHint() )


class MyQMainWindow(QMainWindow):
  def __init__(self, parent):
    super(MyQMainWindow, self).__init__()
    self.parent = parent

  def closeEvent(self, event):
    self.parent.closeEvent(event)
    event.accept()



class NGL_HKLViewer(HKLviewerGui.Ui_MainWindow):
  def __init__(self, thisapp):
    self.window = MyQMainWindow(self)
    self.setupUi(self.window)
    self.app = thisapp
    self.actionOpen_reflection_file.triggered.connect(self.onOpenReflectionFile)
    self.actiondebug.triggered.connect(self.DebugInteractively)
    self.actionSettings.triggered.connect(self.SettingsDialog)
    self.actionExit.triggered.connect(self.window.close)
    self.actionSave_reflection_file.triggered.connect(self.onSaveReflectionFile)
    self.functionTabWidget.setCurrentIndex(0) # if accidentally set to a different tab in the Qtdesigner

    self.UseOSBrowser = False
    self.devmode = False
    for e in sys.argv:
      if "UseOSBrowser" in e:
        self.UseOSBrowser = True
      if "devmode" in e or "debug" in e:
        self.devmode = True

    self.zmq_context = None
    self.unfeedback = False

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
    self.ttiplabeltxt = QLabel()
    self.ttiplabeltxt.setText("Tooltips")
    self.ttipalphalabeltxt = QLabel()
    self.ttipalphalabeltxt.setText("Opacity:")
    self.ttipHoverradio = QRadioButton()
    self.ttipHoverradio.setText( "Hovering")
    self.ttipHoverradio.clicked.connect(self.onShowTooltips)
    self.ttipClickradio = QRadioButton()
    self.ttipClickradio.setText( "Clicked")
    self.ttipClickradio.clicked.connect(self.onShowTooltips)


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

    self.millerarraytable = MillerArrayTableView(self.window)
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
    self.CreateOtherBox()
    self.functionTabWidget.setDisabled(True)

    self.cpath = ""
    if self.UseOSBrowser==False:
      self.InitBrowser()

    self.window.setWindowTitle("HKL-Viewer")
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
    self.millertablemenu = QMenu(self.window)
    self.millertablemenu.triggered.connect(self.onMillerTableMenuAction)
    #self.resize(QSize(self.FileInfoBox.size().width()*2, self.size().height() ))
    self.functionTabWidget.setDisabled(True)
    self.window.show()


  def closeEvent(self, event):
    self.PhilToJsRender('NGL_HKLviewer.action = is_terminating')
    self.closing = True
    self.window.setVisible(False)
    self.webpage.deleteLater() # avoid "Release of profile requested but WebEnginePage still not deleted. Expect troubles !"
    print("HKLviewer closing down...")
    nc = 0
    sleeptime = 0.2
    while not self.canexit: #and nc < 5: # until cctbx.python has finished or after 5 sec
      time.sleep(sleeptime)
      self.ProcessMessages()
      nc += sleeptime
    self.cctbxproc.terminate()
    self.out, self.err = self.cctbxproc.communicate()
    self.cctbxproc.wait()
    self.BrowserBox.close()
    self.BrowserBox.deleteLater()
    event.accept()


  def InitBrowser(self):
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
    self.BrowserBox.setAttribute(Qt.WA_DeleteOnClose)


  def onOpenReflectionFile(self):
    options = QFileDialog.Options()
    fileName, filtr = QFileDialog.getOpenFileName(self.window,
            "Open a reflection file", "",
            "All Files (*);;MTZ Files (*.mtz);;CIF (*.cif)", "", options)
    if fileName:
      #self.HKLnameedit.setText(fileName)
      self.window.setWindowTitle("HKL-viewer: " + fileName)
      #self.infostr = ""
      self.textInfo.setPlainText("")
      self.fileisvalid = False
      self.PhilToJsRender('NGL_HKLviewer.openfilename = "%s"' %fileName )
      self.MillerComboBox.clear()
      self.BinDataComboBox.clear()
      self.tncsvec = []
      #self.ClipPlaneChkGroupBox.setChecked(False)
      self.expandP1checkbox.setChecked(False)
      self.expandAnomalouscheckbox.setChecked(False)
      self.sysabsentcheckbox.setChecked(False)
      self.missingcheckbox.setChecked(False)
      self.onlymissingcheckbox.setChecked(False)
      self.showsliceGroupCheckbox.setChecked(False)


  def onSaveReflectionFile(self):
    options = QFileDialog.Options()
    fileName, filtr = QFileDialog.getSaveFileName(self.window,
            "Save to a new reflection file", "",
            "All Files (*);;MTZ Files (*.mtz)", "", options)
    if fileName:
      self.PhilToJsRender('NGL_HKLviewer.savefilename = "%s"' %fileName )



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
    currentinfostr = None
    if self.out:
      currentinfostr = self.out.decode("utf-8")
      print(self.out.decode("utf-8"))
    if self.err:
      currentinfostr += self.err.decode("utf-8")
      print(self.err.decode("utf-8"))


    if self.zmq_context:
      try:
        binmsg = self.socket.recv(flags=zmq.NOBLOCK) #To empty the socket from previous messages
        msg = zlib.decompress(binmsg)
        nan = float("nan") # workaround for "evaluating" any NaN values in the messages received
        msgstr = msg.decode()
        self.infodict = eval(msgstr)
        if self.infodict:
          if self.infodict.get("WebGL_error"):
            self.BrowserBox.close()
            sys.argv.append("--enable-webgl-software-rendering")
            self.InitBrowser()

          if self.infodict.get("current_phil_strings"):
            philstringdict = self.infodict.get("current_phil_strings", {})
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
            if self.UseOSBrowser==False:
              self.BrowserBox.setUrl(self.html_url)
              # workaround for background colour bug in chromium
              # https://bugreports.qt.io/browse/QTBUG-41960
              self.BrowserBox.page().setBackgroundColor(QColor(127, 127, 127, 0.0) )

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
            self.millerarraytableform.resize(tablewidth, self.millerarraytable.rowHeight(0)*15)
            self.millerarraytableform.show()

          currentinfostr = ""
          if self.infodict.get("info"):
            currentinfostr = self.infodict.get("info",[])
            if self.closing:
              print(currentinfostr)

            if "Destroying HKLViewFrame" in currentinfostr:
              self.canexit = True

          if self.infodict.get("tncsvec"):
            self.tncsvec = self.infodict.get("tncsvec",[])
          if len(self.tncsvec) == 0:
            self.clipTNCSBtn.setDisabled(True)
          else:
            self.clipTNCSBtn.setEnabled(True)

          if self.infodict.get("file_name"):
            self.window.setWindowTitle("HKL-viewer: " + self.infodict.get("file_name", "") )

          if self.infodict.get("NewFileLoaded"):
            self.NewFileLoaded = self.infodict.get("NewFileLoaded",False)

          if self.infodict.get("NewHKLscenes"):
            self.NewHKLscenes = self.infodict.get("NewHKLscenes",False)

          if self.infodict.get("NewMillerArray"):
            self.NewMillerArray = self.infodict.get("NewMillerArray",False)

          if self.infodict.get("used_nth_power_scale_radii", None) is not None:
            self.unfeedback = True
            self.power_scale_spinBox.setValue( self.infodict.get("used_nth_power_scale_radii", 0.0))
            self.unfeedback = False

          self.fileisvalid = True
          #print("ngl_hkl_infodict: " + str(ngl_hkl_infodict))

          if currentinfostr:
            #print(currentinfostr)
            self.infostr += currentinfostr + "\n"
            # display no more than self.bufsize bytes of text
            self.infostr = self.infostr[-1000*self.bufsizespinBox.value():]
            self.textInfo.setPlainText(self.infostr)
            self.textInfo.verticalScrollBar().setValue( self.textInfo.verticalScrollBar().maximum()  )
            currentinfostr = ""

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
            #self.functionTabWidget.setDisabled(True)
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
    self.unfeedback = True
    self.power_scale_spinBox.setEnabled( self.currentphilstringdict['NGL_HKLviewer.viewer.nth_power_scale_radii'] >= 0.0 )
    self.ManualPowerScalecheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.nth_power_scale_radii'] >= 0.0 )
    self.radii_scale_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.viewer.scale'])
    self.showsliceGroupCheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.slice_mode'])
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
    self.rotavecangle_labeltxt.setText("Angle rotated: %2.fº" %self.currentphilstringdict['NGL_HKLviewer.clip_plane.angle_around_vector'])


    self.sliceindexspinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.viewer.slice_index'])
    self.Nbins_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.nbins'])
    if self.currentphilstringdict['NGL_HKLviewer.spacegroup_choice'] is not None:
      self.SpaceGroupComboBox.setCurrentIndex(  self.currentphilstringdict['NGL_HKLviewer.spacegroup_choice'] )
    self.clipParallelBtn.setChecked( self.currentphilstringdict['NGL_HKLviewer.clip_plane.is_parallel'])
    self.missingcheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.show_missing'])
    self.onlymissingcheckbox.setEnabled( self.currentphilstringdict['NGL_HKLviewer.viewer.show_missing'] )
    axidx = -1
    for axidx,c in enumerate(self.sliceaxis.values()):
      if c in self.currentphilstringdict['NGL_HKLviewer.viewer.slice_axis']:
        break

    #self.SliceLabelComboBox.setCurrentIndex( axidx )
    self.cameraPerspectCheckBox.setChecked( "perspective" in self.currentphilstringdict['NGL_HKLviewer.viewer.NGL.camera_type'])
    self.ClipPlaneChkGroupBox.setChecked( self.currentphilstringdict['NGL_HKLviewer.clip_plane.clipwidth'] != None )
    if self.currentphilstringdict['NGL_HKLviewer.clip_plane.clipwidth']:
      self.clipwidth_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.clip_plane.clipwidth'])
    self.hkldist_spinBox.setValue( self.currentphilstringdict['NGL_HKLviewer.clip_plane.hkldist'])
    self.realspacevecBtn.setChecked( "realspace" in self.currentphilstringdict['NGL_HKLviewer.clip_plane.fractional_vector'])
    self.recipvecBtn.setChecked( "reciprocal" in self.currentphilstringdict['NGL_HKLviewer.clip_plane.fractional_vector'])
    self.clipTNCSBtn.setChecked( "tncs" in self.currentphilstringdict['NGL_HKLviewer.clip_plane.fractional_vector'])
    self.fixedorientcheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.NGL.fixorientation'])
    self.onlymissingcheckbox.setChecked( self.currentphilstringdict['NGL_HKLviewer.viewer.show_only_missing'])
    if self.currentphilstringdict['NGL_HKLviewer.real_space_unit_cell_scale_fraction'] is not None:
      self.DrawRealUnitCellBox.setChecked(True)
      self.unitcellslider.setValue( self.currentphilstringdict['NGL_HKLviewer.real_space_unit_cell_scale_fraction'] * self.unitcellslider.maximum())
    else:
      self.DrawRealUnitCellBox.setChecked(False)
    if self.currentphilstringdict['NGL_HKLviewer.reciprocal_unit_cell_scale_fraction'] is not None:
      self.DrawReciprocUnitCellBox.setChecked(True)
      self.reciprocunitcellslider.setValue( self.currentphilstringdict['NGL_HKLviewer.reciprocal_unit_cell_scale_fraction'] * self.reciprocunitcellslider.maximum())
    else:
      self.DrawReciprocUnitCellBox.setChecked(False)
    self.unfeedback = False



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
    if type(self.millerarraytablemodel._data[0][idx]) is str:
      print("Cannot sort this column.")
      return
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
    self.PhilToJsRender('NGL_HKLviewer.viewer.NGL.mouse_sensitivity = %f' %val)


  def onMouseSensitivity(self):
    val = self.mousemoveslider.value()/2000.0
    self.mousesensitxtbox.setText("%2.2f" %val )


  def onTooltipAlphaChanged(self, val):
    if self.unfeedback:
      return
    self.ttipalpha = val
    self.PhilToJsRender('NGL_HKLviewer.viewer.NGL.tooltip_alpha = %f' %val)


  def onShowTooltips(self,val):
    if self.ttipClickradio.isChecked():
      self.PhilToJsRender("NGL_HKLviewer.viewer.NGL.show_tooltips = click")
    if self.ttipHoverradio.isChecked():
      self.PhilToJsRender("NGL_HKLviewer.viewer.NGL.show_tooltips = hover")


  def onFontsizeChanged(self, val):
    font = self.app.font()
    font.setPointSize(val);
    self.app.setFont(font);
    self.settingsform.setFixedSize( self.settingsform.sizeHint() )


  def onCameraPerspect(self,val):
    if self.cameraPerspectCheckBox.isChecked():
      self.PhilToJsRender("NGL_HKLviewer.viewer.NGL.camera_type = perspective")
    else:
      self.PhilToJsRender("NGL_HKLviewer.viewer.NGL.camera_type = orthographic")


  def ExpandRefls(self):
    if self.ExpandReflsGroupBox.isChecked():
      self.PhilToJsRender('NGL_HKLviewer.viewer.expand_to_p1 = True')
      self.PhilToJsRender('NGL_HKLviewer.viewer.expand_anomalous = True')
    else:
      self.PhilToJsRender('NGL_HKLviewer.viewer.expand_to_p1 = False')
      self.PhilToJsRender('NGL_HKLviewer.viewer.expand_anomalous = False')


  def ExpandToP1(self):
    if self.expandP1checkbox.isChecked():
      self.PhilToJsRender('NGL_HKLviewer.viewer.expand_to_p1 = True')
    else:
      self.PhilToJsRender('NGL_HKLviewer.viewer.expand_to_p1 = False')


  def ExpandAnomalous(self):
    if self.expandAnomalouscheckbox.isChecked():
      self.PhilToJsRender('NGL_HKLviewer.viewer.expand_anomalous = True')
    else:
      self.PhilToJsRender('NGL_HKLviewer.viewer.expand_anomalous = False')


  def showSysAbsent(self):
    if self.sysabsentcheckbox.isChecked():
      self.PhilToJsRender('NGL_HKLviewer.viewer.show_systematic_absences = True')
    else:
      self.PhilToJsRender('NGL_HKLviewer.viewer.show_systematic_absences = False')


  def showMissing(self):
    if self.missingcheckbox.isChecked():
      self.PhilToJsRender('NGL_HKLviewer.viewer.show_missing = True')
      self.onlymissingcheckbox.setEnabled(True)
    else:
      self.PhilToJsRender("""NGL_HKLviewer.viewer {
                                                     show_missing = False
                                                     show_only_missing = False
                                                   }
                          """)
      self.onlymissingcheckbox.setEnabled(False)


  def showOnlyMissing(self):
    if self.onlymissingcheckbox.isChecked():
      self.PhilToJsRender('NGL_HKLviewer.viewer.show_only_missing = True')
    else:
      self.PhilToJsRender('NGL_HKLviewer.viewer.show_only_missing = False')


  def showSlice(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      self.ClipPlaneChkGroupBox.setChecked(False)
      self.PhilToJsRender('NGL_HKLviewer.viewer.slice_mode = True')
      if self.expandP1checkbox.isChecked():
        self.PhilToJsRender("""NGL_HKLviewer.viewer {
                                                       expand_to_p1 = True
                                                       inbrowser = False
                                                    }
                             """)
      if self.expandAnomalouscheckbox.isChecked():
        self.PhilToJsRender("""NGL_HKLviewer.viewer {
                                                       expand_anomalous = True
                                                       inbrowser = False
                                                     }
                             """)
      self.ClipPlaneChkGroupBox.setChecked(False)
      self.PhilToJsRender("NGL_HKLviewer.clip_plane.clipwidth = None")
    else:
      self.ClipPlaneChkGroupBox.setChecked(True)
      self.PhilToJsRender("""NGL_HKLviewer.viewer {
                                                      slice_mode = False
                                                      inbrowser = True
                                                    }
                            """)


  def onSliceComboSelchange(self,i):
    if self.unfeedback:
      return
    # 4th element in each table row is the min-max span of hkls.
    rmin = self.array_infotpls[self.MillerComboBox.currentIndex()][4][0][i]
    rmax = self.array_infotpls[self.MillerComboBox.currentIndex()][4][1][i]
    self.sliceindexspinBox.setRange(rmin, rmax)
    self.PhilToJsRender("NGL_HKLviewer.viewer.slice_axis = %s" % self.sliceaxis[i] )


  def onSliceIndexChanged(self, val):
    self.sliceindex = val
    self.PhilToJsRender("NGL_HKLviewer.viewer.slice_index = %d" %self.sliceindex)


  def onBindataComboSelchange(self,i):
    if self.BinDataComboBox.currentText():
      if self.BinDataComboBox.currentIndex() > 0:
        bin_scene_label = str(self.BinDataComboBox.currentIndex()-1)
      else:
        bin_scene_label = "Resolution"
      self.PhilToJsRender("NGL_HKLviewer.bin_scene_label = %s" % bin_scene_label )
      bin_opacitieslst = []
      for i in range(self.nbins):
        bin_opacitieslst.append((1.0, i)) #   ("1.0, %d" %i)
      self.bin_opacities = str(bin_opacitieslst)
      self.OpaqueAllCheckbox.setCheckState(Qt.Checked)
      #self.OpaqueAllCheckbox.setTristate(false)
      #self.PhilToJsRender('NGL_HKLviewer.viewer.NGL.bin_opacities = "%s"' %self.bin_opacities)


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
        self.PhilToJsRender('NGL_HKLviewer.viewer.NGL.bin_opacities = "%s"' %self.bin_opacities )
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
    self.PhilToJsRender('NGL_HKLviewer.viewer.NGL.bin_opacities = "%s"' %self.bin_opacities)
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
      self.PhilToJsRender("NGL_HKLviewer.nbins = %d" %self.nbins)


  def onRadiiScaleChanged(self, val):
    if self.unfeedback:
      return
    self.PhilToJsRender("""
      NGL_HKLviewer.viewer {
        nth_power_scale_radii = %f
        scale = %f
      }
      """ %(self.power_scale_spinBox.value(), self.radii_scale_spinBox.value() )
    )


  def onPowerScaleChanged(self, val):
    if self.unfeedback:
      return
    self.PhilToJsRender("""
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
      self.PhilToJsRender('NGL_HKLviewer.viewer.nth_power_scale_radii = %f' %self.power_scale_spinBox.value())
      #self.power_scale_spinBox.setEnabled(True)
    else:
      self.PhilToJsRender('NGL_HKLviewer.viewer.nth_power_scale_radii = -1.0')
      #self.power_scale_spinBox.setEnabled(False)


  def createExpansionBox(self):
    self.SpaceGroupComboBox.activated.connect(self.SpacegroupSelchange)
    self.expandP1checkbox.clicked.connect(self.ExpandToP1)
    self.expandAnomalouscheckbox.clicked.connect(self.ExpandAnomalous)
    self.ExpandReflsGroupBox.clicked.connect(self.ExpandRefls)
    self.sysabsentcheckbox.clicked.connect(self.showSysAbsent)
    self.missingcheckbox.clicked.connect(self.showMissing)
    self.onlymissingcheckbox.clicked.connect(self.showOnlyMissing)


  def CreateSliceTabs(self):
    self.SliceReflectionsBox.clicked.connect(self.onSliceReflectionsBoxclicked)
    self.showsliceGroupCheckbox.clicked.connect(self.showSlice)

    self.sliceindex = 0
    self.sliceindexspinBox.setValue(self.sliceindex)
    self.sliceindexspinBox.setSingleStep(1)
    self.sliceindexspinBox.setRange(-100, 100)
    self.sliceindexspinBox.valueChanged.connect(self.onSliceIndexChanged)

    self.SliceLabelComboBox.activated.connect(self.onSliceComboSelchange)
    self.sliceaxis = { 0:"h", 1:"k", 2:"l" }
    self.SliceLabelComboBox.addItems( list( self.sliceaxis.values()) )
    self.fixedorientcheckbox.clicked.connect(self.onFixedorient)
    self.recipvecBtn.setText("as fractional values in reciprocal space")
    self.recipvecBtn.setChecked(False)
    self.recipvecBtn.clicked.connect(self.onClipPlaneChkBox)
    self.realspacevecBtn.setText("as fractional values in real space")
    self.realspacevecBtn.setChecked(True)
    self.realspacevecBtn.clicked.connect(self.onClipPlaneChkBox)
    self.clipTNCSBtn.setText("as tNCS vector")
    self.clipTNCSBtn.setChecked(False)
    self.clipTNCSBtn.clicked.connect(self.onClipPlaneChkBox)

    vprec = 2
    self.hvec_spinBox.setValue(2.0)
    self.hvec_spinBox.setDecimals(vprec)
    self.hvec_spinBox.setRange(-100.0, 100.0)
    self.hvec_spinBox.valueChanged.connect(self.onHvecChanged)
    self.kvec_spinBox.setValue(0.0)
    self.kvec_spinBox.setDecimals(vprec)
    self.kvec_spinBox.setSingleStep(0.5)
    self.kvec_spinBox.setRange(-100.0, 100.0)
    self.kvec_spinBox.valueChanged.connect(self.onKvecChanged)
    self.lvec_spinBox.setValue(0.0)
    self.lvec_spinBox.setDecimals(vprec)
    self.lvec_spinBox.setSingleStep(0.5)
    self.lvec_spinBox.setRange(-100.0, 100.0)
    self.lvec_spinBox.valueChanged.connect(self.onLvecChanged)
    self.hkldistval = 0.0
    self.hkldist_spinBox.setValue(self.hkldistval)
    self.hkldist_spinBox.setDecimals(vprec)
    self.hkldist_spinBox.setSingleStep(0.5)
    self.hkldist_spinBox.setRange(-100.0, 100.0)
    self.hkldist_spinBox.valueChanged.connect(self.onHKLdistChanged)
    self.clipwidth_spinBox.setValue(0.5 )
    self.clipwidth_spinBox.setDecimals(vprec)
    self.clipwidth_spinBox.setSingleStep(0.05)
    self.clipwidth_spinBox.setRange(0.0, 100.0)
    self.clipwidth_spinBox.valueChanged.connect(self.onClipwidthChanged)
    self.rotavecangle_labeltxt.setText("Angle rotated: 0 º")
    self.rotavecangle_slider.setValue(0)
    self.rotavecangle_slider.setSingleStep(3)
    self.rotavecangle_slider.sliderReleased.connect(self.onFinalRotaVecAngle)
    self.rotavecangle_slider.valueChanged.connect(self.onRotaVecAngleChanged)
    self.ClipPlaneChkGroupBox.clicked.connect(self.onClipPlaneChkBox)
    self.clipParallelBtn.setText("parallel to vector below")
    self.clipParallelBtn.setChecked(False)
    self.clipParallelBtn.clicked.connect(self.onClipPlaneChkBox)

    self.clipNormalBtn.setText("perpendicular to vector below")
    self.clipNormalBtn.setChecked(True)
    self.clipNormalBtn.clicked.connect(self.onClipPlaneChkBox)
    self.clipParallelBtn.setChecked(True)


  def onSliceReflectionsBoxclicked(self):
    if self.SliceReflectionsBox.isChecked():
      self.showsliceGroupCheckbox.setEnabled(True)
      self.ClipPlaneChkGroupBox.setEnabled(True)
      self.fixedorientcheckbox.setEnabled(True)
      self.onClipPlaneChkBox()
    else:
      self.showsliceGroupCheckbox.setEnabled(False)
      self.ClipPlaneChkGroupBox.setEnabled(False)
      self.fixedorientcheckbox.setEnabled(False)
      self.PhilToJsRender("""NGL_HKLviewer {
                                              clip_plane.clipwidth = None
                                              viewer.slice_mode = False
                                              viewer.NGL.fixorientation = False
                                           }
                          """)


  def onClipPlaneChkBox(self):
    if self.unfeedback:
      return
    if self.ClipPlaneChkGroupBox.isChecked():
      self.showsliceGroupCheckbox.setChecked(False)
      self.PhilToJsRender("""NGL_HKLviewer.viewer {
                                                      slice_mode = False
                                                      inbrowser = True
                                                    }
                            """)
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

      self.PhilToJsRender(philstr)
    else:
      self.showsliceGroupCheckbox.setChecked(True)
      self.PhilToJsRender("""NGL_HKLviewer.viewer.slice_mode = True
                             NGL_HKLviewer.clip_plane.clipwidth = None
                           """)


  def onRotaVecAngleChanged(self, val):
    if self.unfeedback:
      return
    self.PhilToJsRender("""NGL_HKLviewer.clip_plane {
    angle_around_vector = %f
    bequiet = False
}""" %val)


  def onFinalRotaVecAngle(self):
    if self.unfeedback:
      return
    val = self.rotavecangle_slider.value()
    self.PhilToJsRender("""NGL_HKLviewer.clip_plane {
    angle_around_vector = %f
    bequiet = False
}""" %val)


  def onClipwidthChanged(self, val):
    if not self.unfeedback:
      self.PhilToJsRender("NGL_HKLviewer.clip_plane.clipwidth = %f" %self.clipwidth_spinBox.value())


  def onHKLdistChanged(self, val):
    self.hkldistval = val
    if not self.unfeedback:
      self.PhilToJsRender("NGL_HKLviewer.clip_plane.hkldist = %f" %self.hkldistval)


  def onHvecChanged(self, val):
    if not self.unfeedback:
      self.PhilToJsRender("NGL_HKLviewer.clip_plane.h = %f" %self.hvec_spinBox.value())


  def onKvecChanged(self, val):
    if not self.unfeedback:
      self.PhilToJsRender("NGL_HKLviewer.clip_plane.k = %f" %self.kvec_spinBox.value())


  def onLvecChanged(self, val):
    if not self.unfeedback:
      self.PhilToJsRender("NGL_HKLviewer.clip_plane.l = %f" %self.lvec_spinBox.value())


  def onFixedorient(self):
    if not self.unfeedback:
      self.PhilToJsRender('NGL_HKLviewer.viewer.NGL.fixorientation = %s' \
                                    %str(self.fixedorientcheckbox.isChecked()))


  def onMillerTableCellPressed(self, row, col):
    #print( "in millertable CellPressed " + self.millertable.currentItem().text() )
    if self.millertable.mousebutton == Qt.RightButton:
      self.MillerTableContextMenuHandler(QCursor.pos(), row)
    if self.millertable.mousebutton == QEvent.MouseButtonDblClick:
      # quickly display data with a double click
      for i,scenelabel in enumerate(self.scenearraylabels):
        if self.millerarraylabels[row] == scenelabel:
          self.DisplayData(i, row)


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
        myqa = QAction("Display %s data" %scenelabel, self.window, triggered=self.testaction)
        myqa.setData((i, row))
        self.millertablemenu.addAction(myqa)
    myqa = QAction("Make new data as a function of this data...", self.window, triggered=self.testaction)
    myqa.setData( ("newdata_1", row ))
    self.millertablemenu.addAction(myqa)
    myqa = QAction("Make new data as a function of this data and another data set...", self.window, triggered=self.testaction)
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
      myqa = QAction("Show a table of %s data ..." %arraystr, self.window, triggered=self.testaction)
      myqa.setData( ("tabulate_data", self.millertable.selectedrows ))
      self.millertablemenu.addAction(myqa)
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.millertablemenu.exec_(QCursor.pos())


  def onMillerTableMenuAction(self, action):
    data = action.data()
    # depending on what menu item the user clicked data is either an int or a (string, int) tuple
    if data is not None:
      if type(data[0]) is int:
        idx,row = data
        self.DisplayData(idx,row)
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
          self.PhilToJsRender('NGL_HKLviewer.tabulate_miller_array_ids = "%s"' %str(idx))


  def DisplayData(self, idx, row):
    self.PhilToJsRender("NGL_HKLviewer.viewer.scene_id = %d" %idx)
    if self.fileisvalid:
      self.functionTabWidget.setEnabled(True)
      self.expandAnomalouscheckbox.setEnabled(True)
      self.expandP1checkbox.setEnabled(True)
      # don' allow anomalous expansion for data that's already anomalous
      arrayinfo = self.array_infotpls[row]
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
    self.PhilToJsRender('NGL_HKLviewer.miller_array_operations = "[ %s ]"' %str(mtpl) )
    self.makenewdataform.accept()


  def onMillerComboSelchange(self, i):
    self.operate_arrayidx2 = i


  def testaction(self):
    pass


  def createFileInfoBox(self):
    labels = ["Column label", "Type", "Space group", "# HKLs", "Span of HKLs",
       "Min Max data", "Min Max sigmas", "d_min, d_max", "Symmetry unique", "Anomalous"]
    self.millertable.setHorizontalHeaderLabels(labels)
    self.millertable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    # don't allow editing this table
    self.millertable.setEditTriggers(QTableWidget.NoEditTriggers)
    self.millertable.cellPressed.connect(self.onMillerTableCellPressed)
    self.millertable.cellDoubleClicked.connect(self.onMillerTableCellPressed)
    self.millertable.itemSelectionChanged.connect(self.onMillerTableitemSelectionChanged)


  def createRadiiScaleGroupBox(self):
    self.ManualPowerScalecheckbox.clicked.connect(self.onManualPowerScale)
    self.power_scale_spinBox.valueChanged.connect(self.onPowerScaleChanged)
    self.radii_scale_spinBox.valueChanged.connect(self.onRadiiScaleChanged)


  def createBinsBox(self):
    self.binstable_isready = False
    labels = ["no. of HKLs", "lower bin value", "upper bin value", "opacity"]
    self.binstable.setHorizontalHeaderLabels(labels)
    self.binstable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    self.Nbins_spinBox.valueChanged.connect(self.onNbinsChanged)
    self.OpaqueAllCheckbox.clicked.connect(self.onOpaqueAll)

    self.binstable.itemChanged.connect(self.onBinsTableItemChanged  )
    self.binstable.itemClicked.connect(self.onBinsTableitemClicked  )
    self.binstable.itemPressed.connect(self.onBinsTableitemPressed  )
    self.BinDataComboBox.activated.connect(self.onBindataComboSelchange)


  def CreateOtherBox(self):
    self.DrawRealUnitCellBox.clicked.connect(self.onDrawUnitCellBoxClick)
    self.DrawReciprocUnitCellBox.clicked.connect(self.onDrawReciprocUnitCellBoxClick)
    self.unitcellslider.sliderReleased.connect(self.onUnitcellScale)
    self.reciprocunitcellslider.sliderReleased.connect(self.onReciprocUnitcellScale)
    self.ResetViewBtn.clicked.connect(self.onResetViewBtn)


  def onResetViewBtn(self):
    self.PhilToJsRender('NGL_HKLviewer.action = reset_view')


  def onDrawReciprocUnitCellBoxClick(self):
    if not self.unfeedback:
      if self.DrawReciprocUnitCellBox.isChecked():
        val = self.reciprocunitcellslider.value()/self.reciprocunitcellslider.maximum()
        self.PhilToJsRender("NGL_HKLviewer.reciprocal_unit_cell_scale_fraction = %f" %val)
      else:
        self.PhilToJsRender("NGL_HKLviewer.reciprocal_unit_cell_scale_fraction = None")


  def onDrawUnitCellBoxClick(self):
    if not self.unfeedback:
      if self.DrawRealUnitCellBox.isChecked():
        val = self.unitcellslider.value()/self.unitcellslider.maximum()
        self.PhilToJsRender("NGL_HKLviewer.real_space_unit_cell_scale_fraction = %f" %val)
      else:
        self.PhilToJsRender("NGL_HKLviewer.real_space_unit_cell_scale_fraction = None")


  def onUnitcellScale(self):
    if self.unfeedback:
      return
    val = self.unitcellslider.value()/self.unitcellslider.maximum()
    self.PhilToJsRender("NGL_HKLviewer.real_space_unit_cell_scale_fraction = %f" %val)


  def onReciprocUnitcellScale(self):
    if self.unfeedback:
      return
    val = self.reciprocunitcellslider.value()/self.reciprocunitcellslider.maximum()
    self.PhilToJsRender("NGL_HKLviewer.reciprocal_unit_cell_scale_fraction = %f" %val)


  def DebugInteractively(self):
    import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )


  def SpacegroupSelchange(self,i):
    self.PhilToJsRender("NGL_HKLviewer.spacegroup_choice = %d" %i)


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
    """
    cmdargs = 'cctbx.python -i -c "from crys3d.hklview import cmdlineframes;' \
     + ' cmdlineframes.HKLViewFrame(useGuiSocket=%s, high_quality=True,' %self.sockport \
     + ' jscriptfname = \'%s\', ' %self.jscriptfname \
     + ' verbose=%s, UseOSBrowser=%s, htmlfname=\'%s\', handshakewait=%s )"\n'\
       %(self.verbose, str(self.UseOSbrowser), self.hklfname, self.handshakewait)
    """

    guiargs = [ 'useGuiSocket=' + str(self.sockport),
               'high_quality=True',
               'UseOSBrowser=' + str(self.UseOSBrowser)
              ]
    cmdargs =  'cctbx.python -i -c "from crys3d.hklview import cmdlineframes;' \
     + ' cmdlineframes.run()" ' + ' '.join( guiargs + sys.argv[1:])
    self.cctbxproc = subprocess.Popen( cmdargs, shell=True, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr)


  def PhilToJsRender(self, cmdstr):
    if sys.version_info.major==3:
      self.socket.send(bytes(cmdstr,"utf-8"))
    else:
      self.socket.send(bytes(cmdstr))


def run():
  try:
    debugtrue = False
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = " "
    for e in sys.argv:
      if "devmode" in e or "debug" in e:
        debugtrue = True
        # some useful flags as per https://doc.qt.io/qt-5/qtwebengine-debugging.html
        os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--remote-debugging-port=9742 --single-process --js-flags='--expose_gc'"

    settings = QSettings("CCTBX", "HKLviewer" )
    settings.beginGroup("SomeSettings")
    QWebEngineViewFlags = settings.value("QWebEngineViewFlags", None)
    settings.endGroup()

    if QWebEngineViewFlags is None: # avoid doing this test over and over again on the same PC
      QWebEngineViewFlags = ""
      print("testing if WebGL works in QWebEngineView....")
      cmdargs = [ sys.executable, QtChromiumCheck.__file__ ]
      webglproc = subprocess.Popen( cmdargs, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      procout, procerr = webglproc.communicate()
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      if not "WebGL works" in procout.decode():
        QWebEngineViewFlags = " --enable-webgl-software-rendering --ignore-gpu-blacklist"
    if "verbose" in sys.argv[1:]:
      print("using flags for QWebEngineView: " + QWebEngineViewFlags)
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] += QWebEngineViewFlags

    app = QApplication(sys.argv)
    guiobj = NGL_HKLViewer(app)
    timer = QTimer()
    timer.setInterval(20)
    timer.timeout.connect(guiobj.ProcessMessages)
    timer.start()
    ret = app.exec_()

    settings = QSettings("CCTBX", "HKLviewer" )
    settings.beginGroup("SomeSettings")
    QWebEngineViewFlags = settings.setValue("QWebEngineViewFlags", QWebEngineViewFlags)
    settings.endGroup()

    sys.exit(ret)
  except Exception as e:
    print( str(e)  +  traceback.format_exc(limit=10) )


if (__name__ == "__main__") :
  run()
