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

from PySide2.QtCore import Qt, QEvent, QItemSelectionModel, QSize, QSettings, QTimer
from PySide2.QtWidgets import (  QAction, QCheckBox, QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMainWindow, QMenu, QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QPlainTextEdit, QVBoxLayout, QWidget )

from PySide2.QtGui import QColor, QFont, QCursor
from PySide2.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )
import sys, zmq, subprocess, time, traceback, zlib, io, os


try: # if invoked by cctbx.python or some such
  from crys3d.hklview import HKLviewerGui, QtChromiumCheck
  from crys3d.hklview.helpers import ( MillerArrayTableView, MillerArrayTableForm,
                                     MillerArrayTableModel, MPLColourSchemes )
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
    #layout.addWidget(parent.MillerLabel1,           1, 0, 1, 2)
    layout.addWidget(parent.MillerComboBox,        1, 0, 1, 1)
    layout.addWidget(parent.MillerLabel2,          1, 1, 1, 1)
    layout.addWidget(parent.MillerLabel3,          2, 0, 1, 2)
    layout.addWidget(parent.operationtxtbox,       3, 0, 1, 2)
    layout.addWidget(parent.newlabelLabel,         4, 0, 1, 1)
    layout.addWidget(parent.newlabeltxtbox,        4, 1, 1, 1)
    layout.addWidget(parent.operationbutton,       5, 0, 1, 2)
    layout.setRowStretch (0, 1)
    layout.setRowStretch (1 ,0)
    myGroupBox.setLayout(layout)
    mainLayout = QGridLayout()
    mainLayout.addWidget(myGroupBox,     0, 0)
    self.setLayout(mainLayout)
    m = self.fontMetrics().width( "asdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdf")


class AboutForm(QDialog):
  def __init__(self, parent=None):
    super(AboutForm, self).__init__(parent.window)
    self.setWindowTitle("About HKLviewer")
    self.setWindowFlags(Qt.Tool)
    mainLayout = QGridLayout()
    self.aboutlabel = QLabel()
    self.aboutlabel.setWordWrap(True)
    aboutstr = """<html><head/><body><p><span style=" font-weight:600;">
    HKLviewer, </span>A reflection data viewer for crystallography
    <span style=" font-weight:600;"><br/>Developers: </span>Dr. Robert D. Oeffner<br/>
    Cambridge Institute for Medical Research, University of Cambridge.<br/>
    HKLviewer is part of the <a href="http://cci.lbl.gov/docs/cctbx/"> CCTBX library</a>
    as well as derived software thereof.<br/>
    Please direct any bug reports to cctbx@cci.lbl.gov or rdo20@cam.ac.uk

    </p></body></html>"""
    self.aboutlabel.setText(aboutstr)
    self.copyrightstxt = QTextEdit()
    self.copyrightstxt.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    self.copyrightstxt.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    self.copyrightstxt.setReadOnly(True)
    self.copyrightstxt.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
    self.copyrightstxt.setMinimumSize(QSize(400, 100))
    self.copyrightstxt.setMaximumSize(QSize(16777215, 16777215))
    self.OKbtn = QPushButton()
    self.OKbtn.setText("OK")
    self.OKbtn.clicked.connect(self.onOK)
    mainLayout.addWidget(self.aboutlabel,  0, 0, 1, 3)
    mainLayout.addWidget(self.copyrightstxt,  1, 0, 1, 3)
    mainLayout.addWidget(self.OKbtn,  2, 1, 1, 1)
    self.setLayout(mainLayout)
    self.setMinimumSize(QSize(350, 200))
    self.setFixedSize( self.sizeHint() )
  def onOK(self):
    self.hide()


class SettingsForm(QDialog):
  def __init__(self, parent=None):
    super(SettingsForm, self).__init__(parent.window)
    self.setWindowTitle("HKLviewer Settings")
    self.setWindowFlags(Qt.Tool)
    myGroupBox = QGroupBox("Stuff")
    layout = QGridLayout()
    layout.addWidget(parent.mousespeed_labeltxt,     0, 0, 1, 1)
    layout.addWidget(parent.mousemoveslider,         0, 1, 1, 3)
    layout.addWidget(parent.mousesensitxtbox,        0, 4, 1, 1)
    layout.addWidget(parent.Fontsize_labeltxt,       1, 0, 1, 1)
    layout.addWidget(parent.fontspinBox,             1, 4, 1, 1)
    layout.addWidget(parent.BrowserFontsize_labeltxt, 2, 0, 1, 1)
    layout.addWidget(parent.browserfontspinBox,      2, 4, 1, 1)
    layout.addWidget(parent.cameraPerspectCheckBox,  3, 0, 1, 1)
    layout.addWidget(parent.bufsize_labeltxt,        4, 0, 1, 1)
    layout.addWidget(parent.clearbufbtn,             4, 2, 1, 2)
    layout.addWidget(parent.bufsizespinBox,          4, 4, 1, 1)

    layout.addWidget(parent.ttiplabeltxt,            5, 0, 1, 1)
    layout.addWidget(parent.ttipClickradio,          5, 1, 1, 1)
    layout.addWidget(parent.ttipHoverradio,          5, 2, 1, 1)
    layout.addWidget(parent.ttipalphalabeltxt,       5, 3, 1, 1)
    layout.addWidget(parent.ttipalpha_spinBox,       5, 4, 1, 1)

    layout.setRowStretch (0, 1)
    layout.setRowStretch (1 ,0)
    myGroupBox.setLayout(layout)
    mainLayout = QGridLayout()
    mainLayout.addWidget(myGroupBox,     0, 0)
    self.setLayout(mainLayout)
    self.setFixedSize( self.sizeHint() )


class WebEngineDebugForm(QDialog):
  def __init__(self, parent=None):
    super(WebEngineDebugForm, self).__init__(None,
                        Qt.WindowMinimizeButtonHint | # Want minimise and maximise buttons on window.
                        Qt.WindowMaximizeButtonHint | # As they are not the default for QDialog we must
                        Qt.WindowCloseButtonHint |    # add them with flags at creation
                        Qt.CustomizeWindowHint |
                        Qt.WindowTitleHint |
                        Qt.WindowSystemMenuHint
                        )
    self.setWindowTitle("Chrome QWebEngineDebug")
    browser = QWebEngineView()
    mainLayout = QGridLayout()
    mainLayout.addWidget(browser, 0, 0)
    self.setLayout(mainLayout)
    webpage = QWebEnginePage( parent.webprofile, browser)
    browser.setPage(webpage)
    browser.page().setInspectedPage(parent.webpage )
    self.show()


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
    self.actiondebug.setVisible(False)
    self.UseOSBrowser = False
    self.devmode = False
    for e in sys.argv:
      if "UseOSBrowser" in e:
        self.UseOSBrowser = True
      if "devmode" in e:
        self.devmode = True
        self.actiondebug.setVisible(True)
      if  "debug" in e:
        self.actiondebug.setVisible(True)

    self.zmq_context = None
    self.unfeedback = False
    self.cctbxpythonversion = None

    self.mousespeed_labeltxt = QLabel()
    self.mousespeed_labeltxt.setText("Mouse speed:")
    self.mousemoveslider = QSlider(Qt.Horizontal)
    self.mousemoveslider.setMinimum(2)
    self.mousemoveslider.setMaximum(200)
    self.mousemoveslider.setValue(2)
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
    self.fontsize = self.font.pointSize()

    self.browserfontspinBox = QDoubleSpinBox()
    self.browserfontspinBox.setSingleStep(1)
    self.browserfontspinBox.setRange(4, 50)
    self.browserfontspinBox.setValue(self.font.pointSize())
    self.browserfontspinBox.valueChanged.connect(self.onBrowserFontsizeChanged)
    self.BrowserFontsize_labeltxt = QLabel()
    self.BrowserFontsize_labeltxt.setText("Browser font size:")
    self.browserfontsize = None

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
    self.clearbufbtn = QPushButton()
    self.clearbufbtn.setText("Clear all text")
    self.clearbufbtn.clicked.connect(self.onClearTextBuffer)
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
    self.ttip_click_invoke = "hover"

    self.ColourMapSelectDlg = MPLColourSchemes(self)
    self.ColourMapSelectDlg.setWindowTitle("HKLviewer Colour Gradient Maps")

    self.settingsform = SettingsForm(self)
    self.aboutform = AboutForm(self)
    self.webpagedebugform = None

    self.MillerComboBox = QComboBox()
    self.MillerComboBox.activated.connect(self.onMillerComboSelchange)
    self.operationlabeltxt = QLabel()
    self.operationlabeltxt.setWordWrap(True)
    self.MillerLabel2 = QLabel()
    self.MillerLabel2.setText("column")
    self.MillerLabel3 = QLabel()
    self.MillerLabel3.setText("Example: 'newdata=data/sigmas; newsigmas= -42*sigmas' ")
    self.newlabelLabel = QLabel()
    self.newlabelLabel.setText("Column label for new data:")
    self.newlabeltxtbox = QLineEdit('')
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
    self.CreateVectorsBox()
    self.functionTabWidget.setDisabled(True)

    self.cpath = ""
    if self.UseOSBrowser==False:
      self.InitBrowser()
    else:
      self.BrowserBox.setMaximumWidth(0)

    self.window.setWindowTitle("HKLviewer")
    self.cctbxproc = None
    self.LaunchCCTBXPython()
    self.out = None
    self.err = None
    self.comboviewwidth = 0
    self.currentphilstringdict = {}
    self.hklscenes_arrays = []
    self.millerarraylabels = []
    self.scenearraylabeltypes = []
    self.array_infotpls = []
    self.currentmillarray_idx = None
    self.matching_arrays = []
    self.bin_infotpls = None
    self.bin_opacities= None
    self.lowerbinvals = []
    self.upperbinvals = []
    self.html_url = None
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
    self.isfirsttime = False
    self.closing = False
    self.indices = None
    self.datalst = []
    self.tabulate_miller_array = None
    self.binTableCheckState = None
    self.millertablemenu = QMenu(self.window)
    self.millertablemenu.triggered.connect(self.onMillerTableMenuAction)
    #self.resize(QSize(self.FileInfoBox.size().width()*2, self.size().height() ))
    self.functionTabWidget.setDisabled(True)
    self.window.statusBar().showMessage("")
    self.actionOpen_reflection_file.triggered.connect(self.onOpenReflectionFile)
    self.actiondebug.triggered.connect(self.DebugInteractively)
    self.actionSave_Current_Image.triggered.connect(self.onSaveImage)
    self.actionSettings.triggered.connect(self.SettingsDialog)
    self.actionAbout.triggered.connect(self.aboutform.show)
    self.actionExit.triggered.connect(self.window.close)
    self.actionSave_reflection_file.triggered.connect(self.onSaveReflectionFile)
    self.functionTabWidget.setCurrentIndex(0) # if accidentally set to a different tab in the Qtdesigner
    self.window.show()


  def AppAboutToQuit(self):
    print("in AppAboutToQuit")


  def closeEvent(self, event):
    self.PhilToJsRender('action = is_terminating')
    self.closing = True
    #self.window.setVisible(False)
    if self.UseOSBrowser == False:
      self.webpage.deleteLater() # avoid "Release of profile requested but WebEnginePage still not deleted. Expect troubles !"
    print("HKLviewer closing down...")
    nc = 0
    sleeptime = 0.2
    timeout = 10
    while not self.canexit and nc < timeout: # until cctbx.python has finished or after 5 sec
      time.sleep(sleeptime)
      self.ProcessMessages()
      nc += sleeptime
    if nc>= timeout:
      print("Terminating hanging cctbx.python process...")
    self.cctbxproc.terminate()
    self.out, self.err = self.cctbxproc.communicate()
    self.cctbxproc.wait()
    if self.UseOSBrowser == False:
      if self.webpagedebugform and self.devmode:
        self.webpagedebugform.close()
        self.webpagedebugform.deleteLater()
      self.BrowserBox.close()
      self.BrowserBox.deleteLater()
    event.accept()


  def InitBrowser(self):
    # omitting name for QWebEngineProfile() means it is private/off-the-record with no cache files
    self.webprofile = QWebEngineProfile(parent=self.BrowserBox)
    self.webpage = QWebEnginePage( self.webprofile, self.BrowserBox)
    if self.devmode:
      if hasattr(self.webpage, "setInspectedPage"): # older versions of Qt5 hasn't got chromium debug kit
        self.webpage.setUrl("chrome://gpu")
        self.webpagedebugform = WebEngineDebugForm(self)
        self.webpagedebugform.resize( self.window.size())
      else:
        self.webpage.setUrl("https://webglreport.com/")
    else:
      self.webpage.setUrl("https://cctbx.github.io/")
    self.cpath = self.webprofile.cachePath()
    self.BrowserBox.setPage(self.webpage)
    self.BrowserBox.setAttribute(Qt.WA_DeleteOnClose)


  def Browser_download_requested(self, download_item):
    options = QFileDialog.Options()
    fileName, filtr = QFileDialog.getSaveFileName(self.window,
            "Save screenshot to file", download_item.path(),
            "PNG Files (*.png);;All Files (*)", "", options)
    if fileName:
      download_item.setPath(fileName)
      download_item.accept()
      self.download_item = download_item
      download_item.finished.connect(self._download_finished)


  def _download_finished(self):
    file_path = self.download_item.path()
    print("File Downloaded to: %s" %file_path)


  def onOpenReflectionFile(self):
    options = QFileDialog.Options()
    fileName, filtr = QFileDialog.getOpenFileName(self.window,
            "Open a reflection file", "",
            "MTZ Files (*.mtz);;HKL Files (*.hkl);;CIF Files (*.cif);;SCA Files (*.sca);;All Files (*)", "", options)
    if fileName:
      #self.HKLnameedit.setText(fileName)
      self.window.setWindowTitle("HKLviewer: " + fileName)
      #self.infostr = ""
      self.textInfo.setPlainText("")
      self.fileisvalid = False
      self.PhilToJsRender('openfilename = "%s"' %fileName )
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
            "All Files (*);;CIF Files (*.cif);;MTZ Files (*.mtz)", "", options)
    if fileName:
      self.PhilToJsRender('savefilename = "%s"' %fileName )


  def SettingsDialog(self):
    self.settingsform.show()
    # don't know why valueChanged.connect() method only takes effect from here on
    self.fontspinBox.valueChanged.connect(self.onFontsizeChanged)


  def onColourChartSelect(self, selcolmap, powscale): # called when user clicks OK in ColourMapSelectDlg
    if selcolmap != "":
      self.PhilToJsRender("""viewer.color_scheme = %s
viewer.color_powscale = %s""" %(selcolmap, powscale) )


  def ProcessMessages(self):
    """
    Deal with the messages posted to this GUI by cmdlineframes.py
    """
    if self.webpagedebugform is not None:
      self.webpagedebugform.update()
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

        if "cctbx.python.version:" in msgstr:
          self.cctbxpythonversion = msgstr
          self.PhilToJsRender("""NGL {
  fontsize = %d
  show_tooltips = %s
}
""" %(self.browserfontsize, self.ttip_click_invoke) )

          if self.cctbxpythonversion == 'cctbx.python.version: 2':
            # use NGL's download feature for images since websocket_server fails to handle large streams
            self.webpage.profile().downloadRequested.connect(self.Browser_download_requested)
          return

        self.infodict = eval(msgstr)
        if self.infodict:
          if self.infodict.get("WebGL_error"):
            self.BrowserBox.close()
            sys.argv.append("--enable-webgl-software-rendering")
            self.InitBrowser()

          if self.infodict.get("current_phil_strings"):
            philstringdict = self.infodict.get("current_phil_strings", {})
            for unlikely_dict_keyname, val in philstringdict.items():
              try:
                self.currentphilstringdict[unlikely_dict_keyname] = eval(val)
              except Exception as e:
                self.currentphilstringdict[unlikely_dict_keyname] = val
            self.UpdateGUI()

          if self.infodict.get("copyrights"):
            self.copyrightpaths = self.infodict["copyrights"]
            txts = ""
            for copyrighttitle, fname in self.copyrightpaths:
              with open(fname,"r") as f:
                txts = txts + copyrighttitle + ":\n\n"
                txts = txts + f.read()
                txts = txts + "\n" + "#" * 80  + "\n"
            self.aboutform.copyrightstxt.setText(txts)
            self.aboutform.setFixedSize( self.aboutform.sizeHint() )

          if self.infodict.get("scene_array_label_types"):
            self.scenearraylabeltypes = self.infodict.get("scene_array_label_types", [])

          if self.infodict.get("array_infotpls"):
            self.array_infotpls = self.infodict.get("array_infotpls",[])
            #self.millerarraylabels = [ ",".join(e[0]) for e in self.array_infotpls ]
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
            self.lowerbinvals = []
            self.upperbinvals = []
            self.binstable_isready = False
            for row,bin_infotpl in enumerate(self.bin_infotpls):
              for col,elm in enumerate(bin_infotpl):
                # only allow changing the last column with opacity values
                if col == 0:
                  item = QTableWidgetItem(str(elm))
                  item.setFlags(item.flags() ^ Qt.ItemIsEnabled)
                if col==1:
                  self.lowerbinvals.append(elm)
                if col==2:
                  item = QTableWidgetItem(str(elm))
                  item.setFlags(item.flags() ^ Qt.ItemIsEnabled)
                  self.upperbinvals.append(elm)
                if col == 1: # allow changing bin thresholds
                  item = QTableWidgetItem(str(elm))
                  item.setFlags(item.flags() | Qt.ItemIsEditable)
                  item.setToolTip("Change value for adjusting the number of reflections in this bin")
                if col == 3:
                  item = QTableWidgetItem()
                  item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
                  item.setCheckState(Qt.Checked)
                  item.setToolTip("Change visibility of this bin either by ticking or unticking the " +
                                   "check box or by entering an opacity value between 0 and 1")
                self.binstable.setItem(row, col, item)
            self.binstable_isready = True
            if self.bin_opacities:
              self.update_table_opacities()

          if self.infodict.get("bin_opacities"):
            self.bin_opacities = self.infodict["bin_opacities"]
            if self.binstable.rowCount() > 0:
              self.update_table_opacities()

          if self.infodict.get("html_url") and self.html_url is None:
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

          currentinfostr = ""
          if self.infodict.get("info"):
            currentinfostr = self.infodict.get("info",[])
            if self.closing:
              print(currentinfostr)

            if "Destroying HKLViewFrame" in currentinfostr:
              self.canexit = True

          if self.infodict.get("tabulate_miller_array"):
            currentinfostr = "Received table data"
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
            #for e in range(self.millerarraytablemodel.columnCount()):
            #  tablewidth +=  self.millerarraytable.columnWidth(e)
            self.millerarraytable.resizeColumnsToContents()

            self.millerarraytableform.SortComboBox.clear()
            self.millerarraytableform.SortComboBox.addItems(["unsorted"] + labels )
            self.millerarraytableform.SortComboBox.view().setMinimumWidth(self.comboviewwidth)
            self.millerarraytableform.resize(tablewidth, self.millerarraytable.rowHeight(0)*15)
            self.millerarraytableform.show()

          if self.infodict.get("tncsvec"):
            self.tncsvec = self.infodict.get("tncsvec",[])
          self.unfeedback = True
          if self.infodict.get("all_vectors"):
            self.rotvec = None
            self.all_vectors = self.infodict.get("all_vectors",[])

            self.vectortable2.clearContents()
            self.vectortable2.setRowCount(len(self.all_vectors)+1)
            for row, (opnr, label, order, v, hklop, hkls, abcs) in enumerate(self.all_vectors):
              for col,elm in enumerate((label, hklop, hkls, abcs)):
                item = QTableWidgetItem(str(elm))
                if col == 0:
                  item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled ^ Qt.ItemIsEditable)
                  item.setCheckState(Qt.Unchecked)
                item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                self.vectortable2.setItem(row, col, item)

            rc = self.vectortable2.rowCount()-1 # last row is for user defined vector
            item = QTableWidgetItem("new vector")
            item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled ^ Qt.ItemIsEditable)
            item.setCheckState(Qt.Unchecked)
            self.vectortable2.setItem(rc, 0, item)

            item = QTableWidgetItem()
            item.setFlags(item.flags() | Qt.ItemIsEditable)
            self.vectortable2.setItem(rc, 1, item)

            item = QTableWidgetItem()
            item.setFlags(item.flags() | Qt.ItemIsEditable)
            self.vectortable2.setItem(rc, 2, item)

            item = QTableWidgetItem()
            item.setFlags(item.flags() | Qt.ItemIsEditable)
            self.vectortable2.setItem(rc, 3, item)
            self.vectortable2.resizeColumnsToContents()
          self.unfeedback = False
          if self.infodict.get("file_name"):
            self.window.setWindowTitle("HKLviewer: " + self.infodict.get("file_name", "") )

          if self.infodict.get("NewFileLoaded"):
            self.NewFileLoaded = self.infodict.get("NewFileLoaded",False)
            self.vectortable2.clearContents()
            self.functionTabWidget.setDisabled(True)

          if self.infodict.get("NewHKLscenes"):
            self.NewHKLscenes = self.infodict.get("NewHKLscenes",False)

          if self.infodict.get("NewMillerArray"):
            self.NewMillerArray = self.infodict.get("NewMillerArray",False)

          if self.infodict.get("StatusBar"):
            self.window.statusBar().showMessage( self.infodict.get("StatusBar", "") )

          if self.infodict.get("clicked_HKL"):
            (h,k,l) = self.infodict.get("clicked_HKL", ( ) )
          if self.infodict.get("orig_hkl_ids"):
            if self.millerarraytablemodel is not None \
             and self.millerarraytableform.SortComboBox.currentIndex() == 0:
               # can only match hkls in the unsorted table
              orig_hkl_ids = self.infodict.get("orig_hkl_ids", [])
              mode = QItemSelectionModel.Select | QItemSelectionModel.Rows
              self.millerarraytableform.setFocus()
              self.millerarraytable.setFocus()
              for ids in orig_hkl_ids:
                self.millerarraytable.selectRow(ids)
                #self.millerarraytable.selectionModel().select(ids, mode)


          if self.infodict.get("ColourChart") and self.infodict.get("ColourPowerScale"):
            self.ColourMapSelectDlg.selcolmap = self.infodict.get("ColourChart",False)
            self.ColourMapSelectDlg.powscale = self.infodict.get("ColourPowerScale",False)
            self.ColourMapSelectDlg.show()

          if self.infodict.get("bin_labels_type_idxs"):
            bin_labels_type_idxs = self.infodict.get("bin_labels_type_idxs",False)
            self.BinDataComboBox.clear()
            # fill combobox with labels of data that can be used for binning
            for label,labeltype,idx in bin_labels_type_idxs:
              self.BinDataComboBox.addItem(label, (labeltype, idx) )
            self.BinDataComboBox.view().setMinimumWidth(self.comboviewwidth)

          if self.infodict.get("used_nth_power_scale_radii", None) is not None:
            self.unfeedback = True
            self.power_scale_spinBox.setValue( self.infodict.get("used_nth_power_scale_radii", 0.0))
            self.unfeedback = False

          if self.infodict.get("spacegroup_info"):
            spacegroup_info = self.infodict.get("spacegroup_info",False)
            unitcell_info = self.infodict.get("unitcell_info",False)
            htmlstr = '''<html><head/><body><p><span style=" font-weight:600;">Space group: \t
            </span>%s<span style=" font-weight:600;"><br/>Unit cell: \t</span>%s</p></body></html> '''
            self.SpaceGrpUCellText.setText(htmlstr %(spacegroup_info, unitcell_info) )
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
            #self.MillerComboBox.addItems( [""] + self.millerarraylabels )
            self.MillerComboBox.addItem("", userData=-1)
            for k,lbl in enumerate(self.millerarraylabels):
              self.MillerComboBox.addItem(lbl, userData=k)

            self.MillerComboBox.setCurrentIndex(0) # select the first item which is no miller array
            self.comboviewwidth = 0
            for e in self.millerarraylabels:
              self.comboviewwidth = max(self.comboviewwidth, self.MillerComboBox.fontMetrics().width( e) )
            self.MillerComboBox.view().setMinimumWidth(self.comboviewwidth)

            self.millertable.clearContents()
            self.millertable.setRowCount(len(self.array_infotpls))
            for row,millarr in enumerate(self.array_infotpls):
              for col,elm in enumerate(millarr):
                self.millertable.setItem(row, col, QTableWidgetItem(str(elm)))
            #self.functionTabWidget.setDisabled(True)
            self.NewFileLoaded = False

          if self.NewHKLscenes:
            self.NewHKLscenes = False
      except Exception as e:
        errmsg = str(e)
        if "Resource temporarily unavailable" not in errmsg:
          print( errmsg  +  traceback.format_exc(limit=10) )
        #print(errmsg)


  def UpdateGUI(self):
    self.unfeedback = True
    #if not self.isfirsttime:
    self.ManualPowerScalecheckbox.setChecked( self.currentphilstringdict['viewer.nth_power_scale_radii'] >= 0.0 )
    self.power_scale_spinBox.setEnabled( self.currentphilstringdict['viewer.nth_power_scale_radii'] >= 0.0 )
    self.radii_scale_spinBox.setValue( self.currentphilstringdict['viewer.scale'])
    self.showsliceGroupCheckbox.setChecked( self.currentphilstringdict['viewer.slice_mode'])
    self.expandP1checkbox.setChecked( self.currentphilstringdict['viewer.expand_to_p1'])
    self.expandAnomalouscheckbox.setChecked( self.currentphilstringdict['viewer.expand_anomalous'])
    self.sysabsentcheckbox.setChecked( self.currentphilstringdict['viewer.show_systematic_absences'])
    self.ttipalpha_spinBox.setValue( self.currentphilstringdict['NGL.tooltip_alpha'])
    self.mousemoveslider.setValue( 2000*self.currentphilstringdict['NGL.mouse_sensitivity'])
    vecnr, dgr = self.currentphilstringdict['clip_plane.angle_around_vector']
    self.rotavecangle_labeltxt.setText("Reflections rotated around Vector with Angle: %2.fº" %dgr)


    self.sliceindexspinBox.setValue( self.currentphilstringdict['viewer.slice_index'])
    self.Nbins_spinBox.setValue( self.currentphilstringdict['nbins'])
    if self.currentphilstringdict['spacegroup_choice'] is not None:
      self.SpaceGroupComboBox.setCurrentIndex(  self.currentphilstringdict['spacegroup_choice'] )
    #self.clipParallelBtn.setChecked( self.currentphilstringdict['clip_plane.is_parallel'])
    self.missingcheckbox.setChecked( self.currentphilstringdict['viewer.show_missing'])
    self.onlymissingcheckbox.setEnabled( self.currentphilstringdict['viewer.show_missing'] )
    axidx = -1
    for axidx,c in enumerate(self.sliceaxis.values()):
      if c in self.currentphilstringdict['viewer.slice_axis']:
        break

    self.SliceLabelComboBox.setCurrentIndex( axidx )
    self.cameraPerspectCheckBox.setChecked( "perspective" in self.currentphilstringdict['NGL.camera_type'])
    #self.ClipPlaneChkGroupBox.setChecked( self.currentphilstringdict['clip_plane.clipwidth'] != None and not self.showsliceGroupCheckbox.isChecked() )
    if self.currentphilstringdict['clip_plane.clipwidth']:
      self.clipwidth_spinBox.setValue( self.currentphilstringdict['clip_plane.clipwidth'])
    self.hkldist_spinBox.setValue( self.currentphilstringdict['clip_plane.hkldist'])
    self.PlaneParallelCheckbox.setChecked( self.currentphilstringdict['viewer.fixorientation']== "reflection_slice")
    self.AlignVectorGroupBox.setChecked( self.currentphilstringdict['viewer.fixorientation']== "vector")
    self.onlymissingcheckbox.setChecked( self.currentphilstringdict['viewer.show_only_missing'])
    if self.currentphilstringdict['real_space_unit_cell_scale_fraction'] is not None:
      self.DrawRealUnitCellBox.setChecked(True)
      self.unitcellslider.setValue( self.currentphilstringdict['real_space_unit_cell_scale_fraction'] * self.unitcellslider.maximum())
    else:
      self.DrawRealUnitCellBox.setChecked(False)
    if self.currentphilstringdict['reciprocal_unit_cell_scale_fraction'] is not None:
      self.DrawReciprocUnitCellBox.setChecked(True)
      self.reciprocunitcellslider.setValue( self.currentphilstringdict['reciprocal_unit_cell_scale_fraction'] * self.reciprocunitcellslider.maximum())
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
      self.millerarraytable.resizeColumnsToContents()
      return
    idx = i-1
    if type(self.millerarraytablemodel._data[0][idx]) is str:
      print("Cannot sort this column.")
      return
    if self.millerarraytableform.sortChkbox.checkState() == Qt.Unchecked:
      self.millerarraytable.sortByColumn(idx, Qt.SortOrder.DescendingOrder)
    else:
      self.millerarraytable.sortByColumn(idx, Qt.SortOrder.AscendingOrder)
    self.millerarraytable.resizeColumnsToContents()

  def onSortChkbox(self):
    self.onSortComboBoxSelchange(self.millerarraytableform.SortComboBox.currentIndex() )


  def onPrecisionChanged(self, val):
    self.millerarraytablemodel.precision = val
    self.millerarraytable.resizeColumnsToContents()


  def onFinalMouseSensitivity(self):
    val = self.mousemoveslider.value()/2000.0
    self.PhilToJsRender('NGL.mouse_sensitivity = %f' %val)


  def onMouseSensitivity(self):
    val = self.mousemoveslider.value()/2000.0
    self.mousesensitxtbox.setText("%2.2f" %val )


  def onTooltipAlphaChanged(self, val):
    if self.unfeedback:
      return
    self.ttipalpha = val
    self.PhilToJsRender('NGL.tooltip_alpha = %f' %val)


  def onShowTooltips(self, val):
    if self.ttipClickradio.isChecked() or val=="click":
      self.PhilToJsRender("NGL.show_tooltips = click")
      self.ttip_click_invoke = "click"
    if self.ttipHoverradio.isChecked() or val=="hover":
      self.PhilToJsRender("NGL.show_tooltips = hover")
      self.ttip_click_invoke = "hover"


  def onFontsizeChanged(self, val):
    font = self.app.font()
    font.setPointSize(val);
    self.fontsize = val
    self.app.setFont(font);
    self.settingsform.setFixedSize( self.settingsform.sizeHint() )
    self.aboutform.setFixedSize( self.aboutform.sizeHint() )
    self.ColourMapSelectDlg.setFixedSize( self.ColourMapSelectDlg.sizeHint() )


  def onBrowserFontsizeChanged(self, val):
    self.browserfontsize = val
    self.PhilToJsRender("NGL.fontsize = %d" %val)


  def onClearTextBuffer(self):
    self.textInfo.clear()
    self.infostr = ""


  def onCameraPerspect(self,val):
    if self.cameraPerspectCheckBox.isChecked():
      self.PhilToJsRender("NGL.camera_type = perspective")
    else:
      self.PhilToJsRender("NGL.camera_type = orthographic")


  def ExpandRefls(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      if self.ExpandReflsGroupBox.isChecked():
        self.PhilToJsRender('''
        viewer.expand_to_p1 = True
        viewer.expand_anomalous = True
        viewer.inbrowser = False
                        ''' )
      else:
        self.PhilToJsRender('''
        viewer.expand_to_p1 = False
        viewer.expand_anomalous = False
        viewer.inbrowser = False
                        ''' )
    else:
      if self.ExpandReflsGroupBox.isChecked():
        self.PhilToJsRender('''
        viewer.expand_to_p1 = True
        viewer.expand_anomalous = True
        viewer.inbrowser = True
                        ''' )
      else:
        self.PhilToJsRender('''
        viewer.expand_to_p1 = False
        viewer.expand_anomalous = False
        viewer.inbrowser = True
                        ''' )


  def ExpandToP1(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      if self.expandP1checkbox.isChecked():
        self.PhilToJsRender('''
        viewer.expand_to_p1 = True
        viewer.inbrowser = False
                        ''' )
      else:
        self.PhilToJsRender('''
        viewer.expand_to_p1 = False
        viewer.inbrowser = False
                        ''' )
    else:
      if self.expandP1checkbox.isChecked():
        self.PhilToJsRender('''
        viewer.expand_to_p1 = True
        viewer.inbrowser = True
                        ''' )
      else:
        self.PhilToJsRender('''
        viewer.expand_to_p1 = False
        viewer.inbrowser = True
                        ''' )


  def ExpandAnomalous(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      if self.expandAnomalouscheckbox.isChecked():
        self.PhilToJsRender('''
        viewer.expand_anomalous = True
        viewer.inbrowser = False
                        ''' )
      else:
        self.PhilToJsRender('''
        viewer.expand_anomalous = False
        viewer.inbrowser = False
                        ''' )
    else:
      if self.expandAnomalouscheckbox.isChecked():
        self.PhilToJsRender('''
        viewer.expand_anomalous = True
        viewer.inbrowser = True
                        ''' )
      else:
        self.PhilToJsRender('''
        viewer.expand_anomalous = False
        viewer.inbrowser = True
                        ''' )

  def showSysAbsent(self):
    if self.unfeedback:
      return
    if self.sysabsentcheckbox.isChecked():
      self.PhilToJsRender('viewer.show_systematic_absences = True')
    else:
      self.PhilToJsRender('viewer.show_systematic_absences = False')


  def showMissing(self):
    if self.unfeedback:
      return
    if self.missingcheckbox.isChecked():
      self.PhilToJsRender('viewer.show_missing = True')
      self.onlymissingcheckbox.setEnabled(True)
    else:
      self.PhilToJsRender("""viewer.show_missing = False
                             viewer.show_only_missing = False
                          """)
      self.onlymissingcheckbox.setEnabled(False)


  def showOnlyMissing(self):
    if self.unfeedback:
      return
    if self.onlymissingcheckbox.isChecked():
      self.PhilToJsRender('viewer.show_only_missing = True')
    else:
      self.PhilToJsRender('viewer.show_only_missing = False')


  def showSlice(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      #self.ClipPlaneChkGroupBox.setChecked(False)
      self.PhilToJsRender("""viewer {
  slice_mode = True
  inbrowser = False
  fixorientation = "reflection_slice"
  is_parallel = False
}
    """)
      i = self.SliceLabelComboBox.currentIndex()
      rmin = self.array_infotpls[self.currentmillarray_idx][3][0][i]
      rmax = self.array_infotpls[self.currentmillarray_idx][3][1][i]
      self.sliceindexspinBox.setRange(rmin, rmax)
    else:
      #self.ClipPlaneChkGroupBox.setChecked(True)
      self.PhilToJsRender("""viewer {
  slice_mode = False
  inbrowser = True
  fixorientation = "None"
}
       """)


  def onSliceComboSelchange(self,i):
    if self.unfeedback:
      return
    # 3th element in each table row is the min-max span of hkls.
    rmin = self.array_infotpls[self.currentmillarray_idx][3][0][i]
    rmax = self.array_infotpls[self.currentmillarray_idx][3][1][i]
    self.sliceindexspinBox.setRange(rmin, rmax)
    val = "None"
    if self.PlaneParallelCheckbox.isChecked():
      val = "reflection_slice"
    self.PhilToJsRender("""viewer {
    slice_axis = %s
    is_parallel = False
    fixorientation = "%s"
}""" %(self.sliceaxis[i], val))


  def onSliceIndexChanged(self, val):
    if self.unfeedback:
      return
    self.sliceindex = val
    self.PhilToJsRender("viewer.slice_index = %d" %self.sliceindex)


  def onBindataComboSelchange(self, i):
    if self.BinDataComboBox.currentText():
      binner_idx = self.BinDataComboBox.currentIndex()
      self.PhilToJsRender('binner_idx = %d' % binner_idx )
      bin_opacitieslst = []
      for j in range(self.nbins):
        bin_opacitieslst.append((1.0, j))
      self.bin_opacities = str(bin_opacitieslst)
      self.OpaqueAllCheckbox.setCheckState(Qt.Checked)


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


  def onBinsTableItemChanged(self, item):
    #print( "in itemChanged %s,  %s" %(item.text(), str( item.checkState())) )
    row = item.row()
    col = item.column()
    try:
      if not self.bin_opacities:
        return
      bin_opacitieslst = eval(self.bin_opacities)
      alpha = max(0.0, min(1.0, float(item.text()) ) ) # between 0 and 1 only
      try:
        (oldalpha, row) = bin_opacitieslst[row]
        if oldalpha == float(item.text()):
          if item.checkState()==Qt.Unchecked:
            alpha = 0.0
          else:
            alpha = 1.0
      except Exception as e:
        pass

      if col==3 and self.binstable_isready: # changing opacity
        bin_opacitieslst[row] = (alpha, row)
        self.bin_opacities = str(bin_opacitieslst)
        self.SetAllOpaqueCheckboxes()
        self.PhilToJsRender('NGL.bin_opacities = "%s"' %self.bin_opacities )

      if col==1 and self.binstable_isready: # changing scene_bin_thresholds
        aboveitem = self.binstable.item(row-1, 1)
        belowitem = self.binstable.item(row+1, 1)
        if aboveitem is None:
          aboveval = -9e99
        else:
          aboveval = float(aboveitem.text())
        if belowitem is None:
          belowval = 9e99
        else:
          belowval = float(belowitem.text())
        # the new value must be between above and below values only
        newval = min(belowval, max(aboveval, float(item.text()) ) )
        # but the other way round if binning against resolution
        if self.BinDataComboBox.currentIndex() == 0:
          newval = min(aboveval, max(belowval, float(item.text()) ) )
        self.binstable.item(row,col).setText(str(newval))
        #nbins = len(self.bin_infotpls)
        self.lowerbinvals[row] = newval
        allbinvals = self.lowerbinvals + [ self.upperbinvals[-1] ]
        nbins = len(allbinvals)
        self.PhilToJsRender('''
        scene_bin_thresholds = \"%s\"
        nbins = %d
        ''' %(allbinvals, nbins) )

    except Exception as e:
      print( str(e)  +  traceback.format_exc(limit=10) )


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
    self.PhilToJsRender('NGL.bin_opacities = "%s"' %self.bin_opacities)
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
      self.PhilToJsRender("nbins = %d" %self.nbins)


  def onRadiiScaleChanged(self, val):
    if self.unfeedback:
      return
    self.PhilToJsRender("viewer.scale = %f" %self.radii_scale_spinBox.value() )


  def onPowerScaleChanged(self, val):
    if self.unfeedback:
      return
    self.PhilToJsRender("viewer.nth_power_scale_radii = %f" %self.power_scale_spinBox.value() )


  def onManualPowerScale(self, val=None):
    """
    if val is not None:
      self.PhilToJsRender('viewer.nth_power_scale_radii = %f' %val)
      return
    """
    if self.unfeedback:
      return
    if self.ManualPowerScalecheckbox.isChecked():
      self.PhilToJsRender('viewer.nth_power_scale_radii = %f' %self.power_scale_spinBox.value())
      #self.power_scale_spinBox.setEnabled(True)
    else:
      self.PhilToJsRender('viewer.nth_power_scale_radii = -1.0')
      #self.power_scale_spinBox.setEnabled(False)


  def onShowAllVectors(self):
    #self.unfeedback = True
    if self.ShowAllVectorsBtn.checkState()==Qt.Checked:
      for rvrow in range(self.vectortable2.rowCount()):
        self.vectortable2.item(rvrow, 0).setCheckState(Qt.Checked)
    if self.ShowAllVectorsBtn.checkState()==Qt.Unchecked:
      for rvrow in range(self.vectortable2.rowCount()):
        self.vectortable2.item(rvrow, 0).setCheckState(Qt.Unchecked)
    #self.unfeedback = False


  def onVectorTableItemChanged(self, item):
    if self.unfeedback:
      return
    row = item.row()
    col = item.column()
    try:
      rc = self.vectortable2.rowCount()
      label = None
      if self.vectortable2.item(row, 0) is not None:
        label = self.vectortable2.item(row, 0).text()

      if row < len(self.all_vectors):
        if label is None:
          return

        if col==0:
          if item.checkState()==Qt.Checked:
            self.PhilToJsRender(" viewer.show_vector = '[%d, %s]'" %(row, True ))
          else:
            self.PhilToJsRender(" viewer.show_vector = '[%d, %s]'" %(row, False ))

          if self.rotvec is not None: # reset any imposed angle to 0 whenever checking or unchecking a vector
              self.PhilToJsRender("""clip_plane {
                angle_around_vector = '[%d, 0]'
                bequiet = False
              }""" %self.rotvec)
              self.rotavecangle_slider.setValue(0)

          self.rotvec = None
          sum = 0
          for rvrow in range(self.vectortable2.rowCount()):
            if self.vectortable2.item(rvrow, 0) is not None:
              if self.vectortable2.item(rvrow, 0).checkState()==Qt.Checked:
                self.rotvec = rvrow
                sum +=1
          if sum > 1 or sum == 0: # can only use one vector to rotate around. so if more are selected then deselect them altogether
            self.PhilToJsRender("""clip_plane {
  animate_rotation_around_vector = '[%d, %f]'
}""" %(0, -1.0)) #
            self.PhilToJsRender('viewer.fixorientation = "None"')
            self.AnimaRotCheckBox.setCheckState(Qt.Unchecked)
            self.rotvec = None

          if self.rotvec is not None:
            self.RotateAroundframe.setEnabled(True)
            # notify cctbx which is the curently selected vector
            self.PhilToJsRender("clip_plane.angle_around_vector = '[%d, 0.0]'" %self.rotvec)
          else:
            self.RotateAroundframe.setDisabled(True)
          if not self.unfeedback:
            if sum >= rc:
              self.ShowAllVectorsBtn.setCheckState(Qt.Checked)
            if sum == 0:
              self.ShowAllVectorsBtn.setCheckState(Qt.Unchecked)
            if sum >0.0 and sum < rc:
              self.ShowAllVectorsBtn.setCheckState(Qt.PartiallyChecked)

      if row==(rc-1) and label !="" and label != "new vector": # a user defined vector
        if col==1:
          hklop = self.vectortable2.item(row, 1).text()
          self.PhilToJsRender("""viewer.add_user_vector_hkl_op = '%s'
          viewer.user_label = %s """ %(hklop, label))
        if col==2:
          hklvec = self.vectortable2.item(row, 2).text()
          self.PhilToJsRender("""viewer.add_user_vector_hkl = '(%s)'
          viewer.user_label = %s """ %(hklvec, label))
        if col==3:
          abcvec = self.vectortable2.item(row, 3).text()
          self.PhilToJsRender("""viewer.add_user_vector_abc = '(%s)'
          viewer.user_label = %s """ %(abcvec, label))

    except Exception as e:
      print( str(e)  +  traceback.format_exc(limit=10) )


  def createExpansionBox(self):
    self.SpaceGroupComboBox.activated.connect(self.SpacegroupSelchange)
    self.expandP1checkbox.clicked.connect(self.ExpandToP1)
    self.expandAnomalouscheckbox.clicked.connect(self.ExpandAnomalous)
    self.ExpandReflsGroupBox.clicked.connect(self.ExpandRefls)
    self.sysabsentcheckbox.clicked.connect(self.showSysAbsent)
    self.missingcheckbox.clicked.connect(self.showMissing)
    self.onlymissingcheckbox.clicked.connect(self.showOnlyMissing)


  def CreateSliceTabs(self):
    #self.SliceReflectionsBox.clicked.connect(self.onSliceReflectionsBoxclicked)
    self.showsliceGroupCheckbox.clicked.connect(self.showSlice)

    self.sliceindex = 0
    self.sliceindexspinBox.setValue(self.sliceindex)
    self.sliceindexspinBox.setSingleStep(1)
    self.sliceindexspinBox.setRange(-100, 100)
    self.sliceindexspinBox.valueChanged.connect(self.onSliceIndexChanged)
    self.SliceLabelComboBox.activated.connect(self.onSliceComboSelchange)
    self.sliceaxis = { 0:"h", 1:"k", 2:"l" }
    self.SliceLabelComboBox.addItems( list( self.sliceaxis.values()) )

    self.PlaneParallelCheckbox.clicked.connect(self.onPlaneParallelCheckbox)
    vprec = 2
    self.hkldistval = 0.0
    self.hkldist_spinBox.setValue(self.hkldistval)
    self.hkldist_spinBox.setDecimals(vprec)
    self.hkldist_spinBox.setSingleStep(0.5)
    self.hkldist_spinBox.setRange(-100.0, 100.0)
    self.hkldist_spinBox.valueChanged.connect(self.onHKLdistChanged)
    self.clipwidth_spinBox.setValue(2 )
    self.clipwidth_spinBox.setDecimals(vprec)
    self.clipwidth_spinBox.setSingleStep(0.1)
    self.clipwidth_spinBox.setRange(0.0, 100.0)
    self.clipwidth_spinBox.valueChanged.connect(self.onClipwidthChanged)
    self.rotavecangle_labeltxt.setText("Reflections rotated around Vector with Angle: 0º")
    #self.rotavecangle_slider.setValue(0)
    #self.rotavecangle_slider.setSingleStep(1)
    self.rotavecangle_slider.sliderReleased.connect(self.onFinalRotaVecAngle)
    self.rotavecangle_slider.valueChanged.connect(self.onRotaVecAngleChanged)
    self.rotavecangle_slider.setTracking(False)
    self.ClipPlaneChkGroupBox.clicked.connect(self.onClipPlaneChkBox)


  def onAlignedVector(self):
    if self.unfeedback:
      return
    val = "None"
    if self.AlignVectorGroupBox.isChecked():
      val = "vector"
    philstr = """viewer {
        is_parallel = %s
        fixorientation = "%s"
      } """ %(str(self.AlignParallelBtn.isChecked()), val )
    self.PhilToJsRender(philstr)


  def onClipPlaneChkBox(self):
    if self.unfeedback:
      return
    hkldist, clipwidth = 0.0, None
    if self.ClipPlaneChkGroupBox.isChecked():
      if self.showsliceGroupCheckbox.isChecked() == False:
        self.showsliceGroupCheckbox.setChecked(False)
        self.showsliceGroupCheckbox.setChecked(False)
        self.PhilToJsRender("""viewer {
                                                      slice_mode = False
                                                      inbrowser = True
                                                  }
                            """)
      hkldist, clipwidth = self.hkldistval, self.clipwidth_spinBox.value()
    philstr = """clip_plane {
  clipwidth = %s
}
        """ %clipwidth

    self.PhilToJsRender(philstr)


  def onRotaVecAngleChanged(self, val):
    if self.unfeedback or self.rotvec is None:
      return
    #self.rotavecangle_labeltxt.setText("Reflections rotated around Vector with Angle: %dº" %val)
    self.PhilToJsRender("""clip_plane {
    angle_around_vector = '[%d, %f]'
    bequiet = False
}""" %(self.rotvec, val))


  def onFinalRotaVecAngle(self):
    if self.unfeedback or self.rotvec is None:
      return
    val = self.rotavecangle_slider.value()
    self.PhilToJsRender("""clip_plane {
    angle_around_vector = '[%d, %f]'
    bequiet = False
}""" %(self.rotvec, val))


  def onAnimateRotation(self):
    if self.AnimaRotCheckBox.isChecked() == True:
      self.AnimateSpeedSlider.setEnabled(True)
      self.rotavecangle_slider.setDisabled(True)
      speed = self.AnimateSpeedSlider.value()
      self.PhilToJsRender("""clip_plane {
      animate_rotation_around_vector = '[%d, %f]'
}""" %(self.rotvec, speed))
    else:
      self.rotavecangle_slider.setEnabled(True)
      self.AnimateSpeedSlider.setDisabled(True)
      self.PhilToJsRender("""clip_plane {
      animate_rotation_around_vector = '[%d, %f]'
}""" %(self.rotvec, -1.0))


  def onClipwidthChanged(self, val):
    if not self.unfeedback:
      self.PhilToJsRender("clip_plane.clipwidth = %f" %self.clipwidth_spinBox.value())


  def onHKLdistChanged(self, val):
    if not self.unfeedback:
      self.hkldistval = val
      self.PhilToJsRender("clip_plane.hkldist = %f" %self.hkldistval)


  def onHvecChanged(self, val):
    if not self.unfeedback:
      self.PhilToJsRender("clip_plane.h = %f" %self.hvec_spinBox.value())


  def onKvecChanged(self, val):
    if not self.unfeedback:
      self.PhilToJsRender("clip_plane.k = %f" %self.kvec_spinBox.value())


  def onLvecChanged(self, val):
    if not self.unfeedback:
      self.PhilToJsRender("clip_plane.l = %f" %self.lvec_spinBox.value())


  def onPlaneParallelCheckbox(self):
    if not self.unfeedback:
      val = "None"
      if self.PlaneParallelCheckbox.isChecked():
        val = "reflection_slice"
      self.PhilToJsRender('viewer.fixorientation = "%s"' %val)


  def onMillerTableCellPressed(self, row, col):
    #print( "in millertable CellPressed " + self.millertable.currentItem().text() )
    if self.millertable.mousebutton == Qt.RightButton:
      self.MillerTableContextMenuHandler(QCursor.pos(), row)
    if self.millertable.mousebutton == QEvent.MouseButtonDblClick:
      # quickly display data with a double click
      #for sceneid,(scenelabel,labeltype,arrayid,sceneid) in enumerate(self.scenearraylabeltypes):
      for scenelabel,labeltype,arrayid,sceneid in self.scenearraylabeltypes:
        if row == arrayid:
          self.DisplayData(sceneid, row)
          break


  def onMillerTableitemSelectionChanged(self):
    self.millertable.selectedrows = list(set([ e.row() for e in self.millertable.selectedItems() ]))


  def MillerTableContextMenuHandler(self, pos, row):
    self.millertablemenu.clear()
    # Tag menu items with data being int or a (string, int) tuple.
    # These are being checked for in onMillerTableMenuAction() and appropriate
    # action taken
    for i,(scenelabel,labeltype,arrayid,sceneid) in enumerate(self.scenearraylabeltypes): # loop over scenes
      scenelabelstr = ",".join(scenelabel)
      if self.millerarraylabels[row] == scenelabelstr or self.millerarraylabels[row] + " + " in scenelabelstr:
        if labeltype == "hassigmas":
          myqa = QAction("Display data of %s" %scenelabelstr, self.window, triggered=self.testaction)
          myqa.setData((i, row))
          self.millertablemenu.addAction(myqa)
          myqa = QAction("Display sigmas of %s" %scenelabelstr, self.window, triggered=self.testaction)
          myqa.setData((i + 1000, row)) # want to show the sigmas rather than the data if we add 1000
          self.millertablemenu.addAction(myqa)
        else:
          myqa = QAction("Display %s" %scenelabelstr, self.window, triggered=self.testaction)
          myqa.setData((i, row))
          self.millertablemenu.addAction(myqa)
    myqa = QAction("Make new data as a function of this data and another data set...", self.window, triggered=self.testaction)
    myqa.setData( ("newdata", row ))
    self.millertablemenu.addAction(myqa)

    if len(self.millertable.selectedrows) > 0:
      arraystr = ""
      labels = []
      for i,r in enumerate(self.millertable.selectedrows):
        labels.extend( self.millerarraylabels[r].split(",") ) # to cope with I,SigI or other multiple labels
      myqa = QAction("Show a table of %s data ..." %  " and ".join(labels), self.window, triggered=self.testaction)
      myqa.setData( ("tabulate_data", labels ))
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
        self.operate_arrayidx2 = -1 # i.e. no second miller array selected yet
        if strval=="newdata":
          #self.operationlabeltxt.setText("Enter a python expression of " + self.millerarraylabels[idx] + " 'data1' and or 'sigmas1' variable")
          self.operationlabeltxt.setText("Enter a python expression of " + self.millerarraylabels[idx]
           + " 'data1' and/or 'sigmas1' variable and optionally 'data2' and/or 'sigmas2'"
           +" variable from the miller array below:")
          self.MillerComboBox.setEnabled(True)
          self.MillerLabel3.setText("Example: 'newdata = data1 - flex.pow(data2); newsigmas= sigmas1 - data2 / sigmas1' ")
          self.makenewdataform.show()
        if strval=="tabulate_data":
          self.PhilToJsRender('tabulate_miller_array_ids = "%s"' %str(idx))


  def DisplayData(self, idx, row):
    # want to show the sigmas rather than the data if idx we add 1000
    if (idx - 1000) >= 0:
      idx = idx - 1000
      philstr = """
      viewer
      {
        sigma_radius = True
        sigma_color = True
        scene_id = %d
      }
      """ %idx
    else:
      philstr = """
      viewer
      {
        sigma_radius = False
        sigma_color = False
        scene_id = %d
      }
      """ %idx
    self.PhilToJsRender(philstr)
    self.currentmillarray_idx = row
    if self.fileisvalid:
      self.functionTabWidget.setEnabled(True)
      self.expandAnomalouscheckbox.setEnabled(True)
      self.expandP1checkbox.setEnabled(True)
      # don't allow anomalous expansion for data that's already anomalous
      arrayinfo = self.array_infotpls[self.currentmillarray_idx]
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
    self.PhilToJsRender('miller_array_operations = "[ %s ]"' %str(mtpl) )
    self.makenewdataform.accept()


  def onMillerComboSelchange(self, i):
    self.operate_arrayidx2 = self.MillerComboBox.itemData(i)
    # -1 if first item. Otherwise itemdata is index of the list of miller arrays


  def testaction(self):
    pass


  def createFileInfoBox(self):
    labels = ["Column label", "Type", "# HKLs", "Span of HKLs",
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
    labels = ["# HKLs", "lower bin value", "upper bin value", "opacity"]
    self.binstable.setHorizontalHeaderLabels(labels)
    self.binstable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    self.Nbins_spinBox.valueChanged.connect(self.onNbinsChanged)
    self.OpaqueAllCheckbox.clicked.connect(self.onOpaqueAll)
    self.binstable.itemChanged.connect(self.onBinsTableItemChanged  )
    self.binstable.itemPressed.connect(self.onBinsTableitemPressed  )
    self.BinDataComboBox.activated.connect(self.onBindataComboSelchange)


  def CreateVectorsBox(self):
    self.DrawRealUnitCellBox.clicked.connect(self.onDrawUnitCellBoxClick)
    self.DrawReciprocUnitCellBox.clicked.connect(self.onDrawReciprocUnitCellBoxClick)
    self.unitcellslider.sliderReleased.connect(self.onUnitcellScale)
    self.reciprocunitcellslider.sliderReleased.connect(self.onReciprocUnitcellScale)
    labels = ["draw", "rotation", "as hkl", "as abc"]
    self.vectortable2.setHorizontalHeaderLabels(labels)
    self.vectortable2.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    self.vectortable2.itemChanged.connect(self.onVectorTableItemChanged  )
    #self.vectortable2.itemPressed.connect(self.onVectorTableItemChanged  )
    self.RotateAroundframe.setDisabled(True)
    self.ShowAllVectorsBtn.clicked.connect(self.onShowAllVectors)
    self.AlignParallelBtn.clicked.connect(self.onAlignedVector)
    self.AlignNormalBtn.clicked.connect(self.onAlignedVector)
    self.AlignVectorGroupBox.clicked.connect(self.onAlignedVector)
    self.AnimaRotCheckBox.clicked.connect(self.onAnimateRotation)
    self.AnimateSpeedSlider.sliderReleased.connect(self.onAnimateRotation)
    self.AnimateSpeedSlider.setDisabled(True)


  #def onResetView(self):
  #  self.PhilToJsRender('action = reset_view')


  def onSaveImage(self):
    if self.cctbxpythonversion == 'cctbx.python.version: 3': # streaming image over websockets
      options = QFileDialog.Options()
      fileName, filtr = QFileDialog.getSaveFileName(self.window,
              "Save screenshot to file", "",
              "PNG Files (*.png);;All Files (*)", "", options)
      if fileName:
        self.PhilToJsRender('save_image_name = "%s" '%fileName)
    else:
      self.PhilToJsRender('save_image_name = "dummy.png" ')
      # eventual file name prompted to us by Browser_download_requested(


  def onDrawReciprocUnitCellBoxClick(self):
    if not self.unfeedback:
      if self.DrawReciprocUnitCellBox.isChecked():
        val = self.reciprocunitcellslider.value()/self.reciprocunitcellslider.maximum()
        self.PhilToJsRender("reciprocal_unit_cell_scale_fraction = %f" %val)
      else:
        self.PhilToJsRender("reciprocal_unit_cell_scale_fraction = None")


  def onDrawUnitCellBoxClick(self):
    if not self.unfeedback:
      if self.DrawRealUnitCellBox.isChecked():
        val = self.unitcellslider.value()/self.unitcellslider.maximum()
        self.PhilToJsRender("real_space_unit_cell_scale_fraction = %f" %val)
      else:
        self.PhilToJsRender("real_space_unit_cell_scale_fraction = None")


  def onUnitcellScale(self):
    if self.unfeedback:
      return
    val = self.unitcellslider.value()/self.unitcellslider.maximum()
    self.PhilToJsRender("real_space_unit_cell_scale_fraction = %f" %val)


  def onReciprocUnitcellScale(self):
    if self.unfeedback:
      return
    val = self.reciprocunitcellslider.value()/self.reciprocunitcellslider.maximum()
    self.PhilToJsRender("reciprocal_unit_cell_scale_fraction = %f" %val)


  def HighlightReflection(self, hkl):
    self.PhilToJsRender("viewer.show_hkl = '%s'" %str(hkl))


  def DebugInteractively(self):
    import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )


  def SpacegroupSelchange(self,i):
    self.PhilToJsRender("spacegroup_choice = %d" %i)


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

    guiargs = [ 'useGuiSocket=' + str(self.sockport),
               'high_quality=True',
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
  import time
  #time.sleep(10)
  try:
    import PySide2.QtCore
    Qtversion = str(PySide2.QtCore.qVersion())
    debugtrue = False
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = " "
    for e in sys.argv:
      if "devmode" in e or "debug" in e and not "UseOSBrowser" in e:
        debugtrue = True
        print("Qt version " + Qtversion)
        # some useful flags as per https://doc.qt.io/qt-5/qtwebengine-debugging.html
        if "debug" in e:
          os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--remote-debugging-port=9741 --single-process --js-flags='--expose_gc'"
        if "devmode" in e:  # --single-process will freeze the WebEngineDebugForm at breakpoints
          os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--js-flags='--expose_gc'"

    settings = QSettings("CCTBX", "HKLviewer" )
    # In case of more than one PySide2 installation tag the settings by version number of PySide2
    # as different versions occasionally use slightly different metrics for font and window sizes
    settings.beginGroup("PySide2_" + Qtversion)
    QWebEngineViewFlags = settings.value("QWebEngineViewFlags", None)
    fontsize = settings.value("FontSize", None)
    browserfontsize = settings.value("BrowserFontSize", None)
    power_scale_value = settings.value("PowerScaleValue", None)
    radii_scale_value = settings.value("RadiiScaleValue", None)
    ttip_click_invoke = settings.value("ttip_click_invoke", None)
    windowsize = settings.value("windowsize", None)
    splitter1sizes = settings.value("splitter1Sizes", None)
    splitter2sizes = settings.value("splitter2Sizes", None)
    settings.endGroup()
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )

    if QWebEngineViewFlags is None: # avoid doing this test over and over again on the same PC
      QWebEngineViewFlags = " --disable-web-security" # for chromium
      print("testing if WebGL works in QWebEngineView....")
      cmdargs = [ sys.executable, QtChromiumCheck.__file__ ]
      webglproc = subprocess.Popen( cmdargs, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      procout, procerr = webglproc.communicate()
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      if not "WebGL works" in procout.decode():
        QWebEngineViewFlags = " --enable-webgl-software-rendering --ignore-gpu-blacklist "
    if "verbose" in sys.argv[1:]:
      print("using flags for QWebEngineView: " + QWebEngineViewFlags)
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] += QWebEngineViewFlags

    from PySide2.QtWidgets import QApplication
    # ensure QWebEngineView scales correctly on a screen with high DPI
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    guiobj = NGL_HKLViewer(app)
    time.sleep(1) # make time for zmq_listen loop to start in cctbx subprocess

    def MyAppClosing():
      settings.beginGroup("PySide2_" + Qtversion )
      settings.setValue("QWebEngineViewFlags", QWebEngineViewFlags)
      settings.setValue("FontSize", guiobj.fontsize )
      settings.setValue("BrowserFontSize", guiobj.browserfontsize )
      settings.setValue("PowerScaleValue", guiobj.currentphilstringdict["viewer.nth_power_scale_radii"])
      settings.setValue("RadiiScaleValue", guiobj.radii_scale_spinBox.value())
      settings.setValue("ttip_click_invoke", guiobj.ttip_click_invoke)
      settings.setValue("windowsize", guiobj.window.size())
      settings.setValue("splitter1Sizes", guiobj.splitter.saveState())
      settings.setValue("splitter2Sizes", guiobj.splitter_2.saveState())
      settings.endGroup()

    app.lastWindowClosed.connect(MyAppClosing)

    timer = QTimer()
    timer.setInterval(20)
    timer.timeout.connect(guiobj.ProcessMessages)
    timer.start()

    if power_scale_value is not None:
      psv = eval(power_scale_value)
      if psv >= 0.0:
        guiobj.power_scale_spinBox.setValue(psv)
      guiobj.ManualPowerScalecheckbox.setChecked(psv >= 0.0)
      guiobj.onManualPowerScale() # disables power_scale_spinBox if psv is less than 0.0
    if radii_scale_value is not None:
      guiobj.radii_scale_spinBox.setValue(eval(radii_scale_value))
      guiobj.onRadiiScaleChanged(None)
    if fontsize is not None:
      guiobj.onFontsizeChanged(int(fontsize))
      guiobj.fontspinBox.setValue(int(fontsize))
    if browserfontsize is not None:
      guiobj.onBrowserFontsizeChanged(int(browserfontsize))
      guiobj.browserfontspinBox.setValue(int(browserfontsize))
    if ttip_click_invoke is not None:
      guiobj.onShowTooltips(ttip_click_invoke)
      guiobj.ttipClickradio.setChecked(ttip_click_invoke == "click")
      guiobj.ttipHoverradio.setChecked(ttip_click_invoke == "hover")
    if splitter1sizes is not None and splitter2sizes is not None and windowsize is not None:
      guiobj.window.resize(windowsize)
      if guiobj.webpagedebugform and guiobj.devmode:
        guiobj.webpagedebugform.resize( guiobj.window.size())
      guiobj.splitter.restoreState(splitter1sizes)
      guiobj.splitter_2.restoreState(splitter2sizes)

    ret = app.exec_()

  except Exception as e:
    print( str(e)  +  traceback.format_exc(limit=10) )


if (__name__ == "__main__") :
  run()
