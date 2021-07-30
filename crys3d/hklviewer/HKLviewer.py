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

import sys, zmq, subprocess, time, traceback, zlib, io, os, math, os.path
if sys.version_info[0] < 3:
  print("HKLviewer GUI must be run from Python 3")
  sys.exit(-42)

from .qt import Qt, QtCore, QCoreApplication, QEvent, QItemSelectionModel, QSize, QSettings, QTimer, QUrl
from .qt import (  QAction, QCheckBox, QComboBox, QDialog, QDoubleSpinBox,
    QFileDialog, QFrame, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
    QMainWindow, QMenu, QMenuBar, QProgressBar, QPushButton, QRadioButton, QRect,
    QScrollBar, QSizePolicy, QSlider, QSpinBox, QStyleFactory, QStatusBar, QTableView, QTableWidget,
    QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget )

from .qt import QColor, QFont, QCursor, QDesktopServices
from .qt import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )
from . import HKLviewerGui
try: # if invoked by cctbx.python or some such
  from crys3d.hklview import HKLviewerGui
  from crys3d.hklview.helpers import ( MillerArrayTableView, MillerArrayTableForm,
                                     MillerArrayTableModel, MPLColourSchemes )
except Exception as e: # if invoked by a generic python that doesn't know cctbx modules
  from . import HKLviewerGui
  from .helpers import MillerArrayTableView, MillerArrayTableForm, MillerArrayTableModel, MPLColourSchemes

class MakeNewDataForm(QDialog):
  def __init__(self, parent=None):
    super(MakeNewDataForm, self).__init__(parent.window)
    self.setWindowFlag(Qt.WindowContextHelpButtonHint,False);
    self.setWindowTitle("Create a new reflection dataset")
    self.setSizeGripEnabled(True)
    layout = QGridLayout()
    sp = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    sp.setVerticalStretch(10)
    parent.operationtxtbox.setSizePolicy(sp)
    parent.operationtxtbox.setMinimumSize(QSize(0, 30))

    sp2 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    sp2.setVerticalStretch(0)
    parent.operationlabeltxt.setSizePolicy(sp2)

    layout.addWidget(parent.operationlabeltxt,     0, 0, 1, 2)
    layout.addWidget(parent.MillerComboBox,        1, 0, 1, 1)
    layout.addWidget(parent.MillerLabel2,          1, 1, 1, 1)
    layout.addWidget(parent.MillerLabel3,          2, 0, 1, 2)
    layout.addWidget(parent.operationtxtbox,       3, 0, 1, 2)
    layout.addWidget(parent.newlabelLabel,         4, 0, 1, 1)
    layout.addWidget(parent.newlabeltxtbox,        4, 1, 1, 1)
    layout.addWidget(parent.operationbutton,       5, 0, 1, 2)
    layout.setRowStretch (0, 1)
    layout.setRowStretch (1 ,0)
    self.setLayout(layout)
    m = self.fontMetrics().width( "asdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdfasdf")


class AboutForm(QDialog):
  def __init__(self, parent=None):
    super(AboutForm, self).__init__(parent.window)
    self.setWindowTitle("About HKLviewer")
    self.setWindowFlags(Qt.Tool)
    mainLayout = QGridLayout()
    self.aboutlabel = QLabel()
    self.aboutlabel.setWordWrap(True)
    self.aboutlabel.setTextInteractionFlags(Qt.TextBrowserInteraction);
    self.aboutlabel.setOpenExternalLinks(True);
    self.writeAboutstr("")
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
  def writeAboutstr(self, versionstr):
    aboutstr = """<html><head/><body><p>
    <span style=" font-weight:600;">HKLviewer, </span>
    CCTBX version: """ + versionstr + \
    """
    <br/>A reflection data viewer for crystallography
    <br/>Developers: Dr. Robert D. Oeffner<br/>
    Cambridge Institute for Medical Research, University of Cambridge.<br/>
    HKLviewer is part of the <a href="http://cci.lbl.gov/docs/cctbx/"> CCTBX library</a>
    as well as derived software thereof.<br/>
    HKLviewer uses functionality provided by the
    <a href="https://github.com/nglviewer/ngl">NGL Viewer</a> project and
    the <a href="https://github.com/niklasvh/html2canvas">html2canvas</a> project.
    Refer to rdo20@cam.ac.uk or cctbx@cci.lbl.gov for queries or bug reports.
    </p></body></html>"""
    self.aboutlabel.setText(aboutstr)
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
    self.setWindowTitle("Chromium QWebEngineDebug")
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


class MyQMainDialog(QDialog):
  def __init__(self, parent):
    super(MyQMainDialog, self).__init__()
    self.parent = parent

  def closeEvent(self, event):
    self.parent.closeEvent(event)
    event.accept()


class NGL_HKLViewer(HKLviewerGui.Ui_MainWindow):
  settings = QSettings("CCTBX", "HKLviewer" )
  # qversion() comes out like '5.12.5'. We just want '5.12'
  Qtversion = "Qt" + ".".join( QtCore.qVersion().split(".")[0:2])

  def __init__(self, thisapp, isembedded=False): #, cctbxpython=None):
    self.datatypedict = { }
    self.browserfontsize = None
    self.isembedded = isembedded
    print("version " + self.Qtversion)

    self.ReadPersistedQsettings()

    if isembedded:
      self.window = MyQMainDialog(self)
      self.window.hide()
    else:
      self.window = MyQMainWindow(self)
    self.setupUi(self.window)
    if isembedded:
      mainLayout = QGridLayout()
      mainLayout.addWidget(self.widget, 0, 0)
      self.window.setLayout(mainLayout)
    else:
      self.window.setCentralWidget(self.centralwidget)
      self.menubar = QMenuBar(self.window)
      self.menubar.setObjectName(u"menubar")
      self.menubar.setGeometry(QRect(0, 0, 1093, 22))
      self.menuFile = QMenu(self.menubar)
      self.menuFile.setObjectName(u"menuFile")
      self.menuHelp = QMenu(self.menubar)
      self.menuHelp.setObjectName(u"menuHelp")
      self.window.setMenuBar(self.menubar)
      self.statusBar = QStatusBar(self.window)
      self.statusBar.setObjectName(u"statusBar")
      self.window.setStatusBar(self.statusBar)
      self.menubar.addAction(self.menuFile.menuAction())
      self.menubar.addAction(self.menuHelp.menuAction())
      self.menuFile.addAction(self.actionOpen_reflection_file)
      self.menuFile.addAction(self.actionSave_reflection_file)
      self.menuFile.addAction(self.actionSettings)
      self.menuFile.addAction(self.actiondebug)
      self.menuFile.addAction(self.actionColour_Gradient)
      self.menuFile.addAction(self.actionSave_Current_Image)
      self.menuFile.addAction(self.actionExit)
      self.menuHelp.addAction(self.actionLocal_Help)
      self.menuHelp.addAction(self.actionHKLviewer_Tutorial)
      self.menuHelp.addAction(self.actionCCTBXwebsite)
      self.menuHelp.addAction(self.actionAbout)
      self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
      self.menuHelp.setTitle(QCoreApplication.translate("MainWindow", u"Help", None))

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
    #self.cctbxpython = cctbxpython
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

    self.browserfontspinBox = QDoubleSpinBox()
    self.browserfontspinBox.setSingleStep(1)
    self.browserfontspinBox.setRange(4, 50)
    self.browserfontspinBox.setValue(self.font.pointSize())
    self.browserfontspinBox.valueChanged.connect(self.onBrowserFontsizeChanged)
    self.BrowserFontsize_labeltxt = QLabel()
    self.BrowserFontsize_labeltxt.setText("Browser font size:")

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
    #self.ttip_click_invoke = "hover"

    self.ColourMapSelectDlg = MPLColourSchemes(self)
    self.ColourMapSelectDlg.setWindowTitle("HKLviewer Colour Gradient Maps")
    # colour schemes and radii mapping for types of datasets stored in jsview_3d.py but persisted here:
    # colourmap=brg, colourpower=1, powerscale=1, radiiscale=1
    self.settingsform = SettingsForm(self)
    self.aboutform = AboutForm(self)
    self.webpagedebugform = None

    self.MillerComboBox = QComboBox()
    self.MillerComboBox.activated.connect(self.onMillerComboSelchange)
    self.operationlabeltxt = QLabel()
    self.operationlabeltxt.setWordWrap(True)
    self.MillerLabel2 = QLabel()
    self.MillerLabel2.setText("is represented by array2")
    self.MillerLabel3 = QLabel()
    self.MillerLabel3.setText("""<html><head/><body><p>
    For examples on creating a dataset from existing ones see
    <a href="http://cci.lbl.gov/docs/cctbx/doc_hklviewer/#making-a-new-dataset">Making a new dataset</a>.
    <br>
    For details on python scripting cctbx.miller.array see
    <a href="https://cctbx.github.io/cctbx/cctbx.miller.html#the-miller-array">cctbx.miller arrays</a>.
    """)
    self.MillerLabel3.setTextInteractionFlags(Qt.TextBrowserInteraction);
    self.MillerLabel3.setOpenExternalLinks(True);

    self.newlabelLabel = QLabel()
    self.newlabelLabel.setText("Column label for new reflection dataset:")
    self.newlabeltxtbox = QLineEdit('')
    self.operationtxtbox = QTextEdit('')
    self.operationtxtbox.setPlaceholderText("""Example:
dat = array1.data()/array1.sigmas() * array2.normalize().data()
sigs = 2*flex.exp(1/array2.sigmas())
newarray._data = dat
newarray._sigmas = sigs
    """)
    self.operationtxtbox.setToolTip(u"<html><head/><body><p>Rather than using math.exp(), math.pow(), etc. " +
                                  "use flex.exp(), flex.pow() etc when operating on flex.array variables.</p></body></html>")
    self.operationbutton = QPushButton("Compute new reflection dataset")
    self.operationbutton.clicked.connect(self.onMakeNewData)
    self.makenewdataform = MakeNewDataForm(self)
    self.makenewdataform.setModal(True)

    self.millerarraytable = MillerArrayTableView(self.window)
    self.millerarraytable.setSortingEnabled(False)
    self.millerarraytable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    self.millerarraytableform = MillerArrayTableForm(self)
    self.millerarraytablemodel = None
    self.millerarraytable.installEventFilter(self.millerarraytableform) # for keyboard copying to clipboard
    self.millerarraytable.setToolTip("Double-click the column and row of a reflection for the "+
                                     "viewer to zoom in on that.\nRight-click a reflection in " +
                                     "the viewer to locate its entry in this table.")

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
    self.functionTabWidget.setDisabled(True)
    self.Statusbartxtbox = None
    self.chimeraxprocmsghandler = None
    if not isembedded:
      self.window.statusBar().showMessage("")
      self.hklLabel = QLabel()
      self.hklLabel.setText("HKL vectors along X-axis, Y-axis, Z-axis:")
      self.Statusbartxtbox = QLineEdit('')
      self.Statusbartxtbox.setReadOnly(True)
      self.Statusbartxtbox.setAlignment(Qt.AlignRight)
      self.window.statusBar().addPermanentWidget(self.hklLabel)
      self.window.statusBar().addPermanentWidget(self.Statusbartxtbox, 1)
      self.actionOpen_reflection_file.triggered.connect(self.onOpenReflectionFile)
      self.actionLocal_Help.triggered.connect(self.onOpenHelpFile)
      self.actionHKLviewer_Tutorial.triggered.connect(self.onOpenTutorialHelpFile)
      self.actionCCTBXwebsite.triggered.connect(self.onOpenCCTBXwebsite)
      self.actiondebug.triggered.connect(self.DebugInteractively)
      self.actionSave_Current_Image.triggered.connect(self.onSaveImage)
      self.actionSettings.triggered.connect(self.SettingsDialog)
      self.actionAbout.triggered.connect(self.aboutform.show)
      self.actionExit.triggered.connect(self.window.close)
      self.actionSave_reflection_file.triggered.connect(self.onSaveReflectionFile)
      self.actionColour_Gradient.triggered.connect(self.ColourMapSelectDlg.show)
    else:
      self.textInfo.setVisible(False) # stdout sent to chimeraX's console instead

    self.functionTabWidget.setCurrentIndex(0) # if accidentally set to a different tab in the Qtdesigner
    self.window.show()


  def onOpenHelpFile(self):
    QDesktopServices.openUrl("http://cci.lbl.gov/docs/cctbx/doc_hklviewer/")


  def onOpenTutorialHelpFile(self):
    QDesktopServices.openUrl("http://cci.lbl.gov/docs/cctbx/tuto_hklviewer/")


  def onOpenCCTBXwebsite(self):
    QDesktopServices.openUrl("http://cci.lbl.gov/docs/cctbx/")


  def closeEvent(self, event):
    self.send_message('action = is_terminating')
    self.closing = True
    #self.window.setVisible(False)
    if self.UseOSBrowser == False:
      self.webpage.deleteLater() # avoid "Release of profile requested but WebEnginePage still not deleted. Expect troubles !"
    print("HKLviewer closing down...")
    nc = 0
    sleeptime = 0.2
    maxtime = 3
    while not self.canexit and nc < maxtime: # until cctbx.python has finished or after maxtime sec
      time.sleep(sleeptime)
      self.ProcessMessages()
      nc += sleeptime
    try:
      #self.cctbxproc.terminate()
      self.out, self.err = self.cctbxproc.communicate(input="exit()", timeout=maxtime)
      print(str(self.out) + "\n" + str(self.err))
    except Exception as e:
      print("Terminating hanging cctbx.python process...")
      import psutil
      parent_pid = self.cctbxproc.pid   # my example
      parent = psutil.Process(parent_pid)
      for child in parent.children(recursive=True):  # or parent.children() for recursive=False
        child.kill()
      parent.kill()
    if self.UseOSBrowser == False:
      if self.webpagedebugform and self.devmode:
        self.webpagedebugform.close()
        self.webpagedebugform.deleteLater()
      self.BrowserBox.close()
      self.BrowserBox.deleteLater()
    self.PersistQsettings()
    if not self.isembedded:
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
      #self.webpage.setUrl("https://cctbx.github.io/")
      self.webpage.setUrl(QUrl("http://cci.lbl.gov/docs/cctbx/doc_hklviewer/"))
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
            "MTZ Files (*.mtz);;CIF Files (*.cif);;HKL Files (*.hkl);;SCA Files (*.sca);;All Files (*)", "", options)
    if fileName:
      #self.HKLnameedit.setText(fileName)
      self.window.setWindowTitle("HKLviewer: " + fileName)
      #self.infostr = ""
      self.textInfo.setPlainText("")
      self.fileisvalid = False
      self.send_message('openfilename = "%s"' %fileName )
      self.MillerComboBox.clear()
      self.BinDataComboBox.clear()
      self.millertable.clearContents()
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
            "Save datasets to a new reflection file", "",
            "MTZ Files (*.mtz);;CIF Files (*.cif);;All Files (*)", "", options)
    if fileName:
      self.send_message('savefilename = "%s"' %fileName )


  def SettingsDialog(self):
    self.settingsform.show()
    # don't know why valueChanged.connect() method only takes effect from here on
    self.fontspinBox.valueChanged.connect(self.onFontsizeChanged)


  def onColourChartSelect(self, selcolmap, colourpowscale):
    # called when user clicks OK in helpers.MPLColourSchemes.EnactColourMapSelection()
    if selcolmap != "":
      self.send_message("""viewer.color_scheme = %s
viewer.color_powscale = %s""" %(selcolmap, colourpowscale) )


  def ProcessMessages(self):
    """
    Deal with the messages posted to this GUI by cmdlineframes.py
    """
    if self.webpagedebugform is not None:
      self.webpagedebugform.update()
    if self.zmq_context:
      try:
        binmsg = self.socket.recv(flags=zmq.NOBLOCK) #To empty the socket from previous messages
        msg = zlib.decompress(binmsg)
        nan = float("nan") # workaround for "evaluating" any NaN values in the messages received
        msgstr = msg.decode()

        if "cctbx.python.version:" in msgstr:
          self.cctbxpythonversion = msgstr
          self.send_message("""NGL {
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
                txts = txts + " "*20 + copyrighttitle + ":\n\n"
                txts = txts + f.read()
                txts = txts + "\n" + "#" * 50  + "\n"
            self.aboutform.copyrightstxt.setText(txts)
            self.aboutform.setFixedSize( self.aboutform.sizeHint() )
          if self.infodict.get("cctbxversion"):
            self.aboutform.writeAboutstr( self.infodict["cctbxversion"])

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
                  item.setToolTip("Change a bin threshold by entering a preferred value in the " +
                                  "\"lower bin value\" column for a particular bin.")
                if col == 3:
                  item = QTableWidgetItem()
                  item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
                  item.setCheckState(Qt.Checked)
                  item.setToolTip("Change visibility of bins either by ticking or unticking the check " +
                                  "boxes or by entering opacity values between 0 and 1. Reflections with "+
                                  "values less than 0.3 will not respond to mouse clicks.")
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
              self.BrowserBox.setUrl(QUrl(self.html_url))
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
          if self.infodict.get("all_vectors") is not None:
            self.rotvec = None
            self.all_vectors = self.infodict.get("all_vectors",[])

            self.vectortable2.clearContents()
            self.vectortable2.setRowCount(len(self.all_vectors)+1)
            for row, (opnr, label, order, v, hklop, hkls, abcs) in enumerate(self.all_vectors):
              for col,elm in enumerate((label, hklop, hkls, abcs)):
                item = QTableWidgetItem(str(elm))
                if col == 0:
                  item.setFlags((Qt.ItemIsUserCheckable | Qt.ItemIsEnabled) ^ Qt.ItemIsEditable)
                  item.setCheckState(Qt.Unchecked)
                item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                self.vectortable2.setItem(row, col, item)

            rc = self.vectortable2.rowCount()-1 # last row is for user defined vector
            item = QTableWidgetItem("new vector")
            item.setFlags((Qt.ItemIsUserCheckable | Qt.ItemIsEnabled) ^ Qt.ItemIsEditable)
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
            self.currentmillarray_idx = None
            self.vectortable2.clearContents()
            self.ShowAllVectorsBtn.setCheckState(Qt.Unchecked)
            self.functionTabWidget.setDisabled(True)
            self.AlignVectorGroupBox.setChecked( False)

          if self.infodict.get("NewHKLscenes"):
            self.NewHKLscenes = self.infodict.get("NewHKLscenes",False)

          if self.infodict.get("NewMillerArray"):
            self.NewMillerArray = self.infodict.get("NewMillerArray",False)

          if self.infodict.get("StatusBar") and self.Statusbartxtbox is not None:
            self.Statusbartxtbox.setText(self.infodict.get("StatusBar", "") )

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
            self.ColourMapSelectDlg.selcolmap = self.infodict.get("ColourChart", "brg")
            self.ColourMapSelectDlg.setPowerScaleSliderVal( self.infodict.get("ColourPowerScale", 1.0))
            if self.infodict.get("ShowColourMapDialog"):
              arrayinfo = self.array_infotpls[self.currentmillarray_idx]
              self.ColourMapSelectDlg.setDataType(arrayinfo[1])
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


          if self.infodict.get("datatype_dict"):
            self.datatypedict = self.infodict.get("datatype_dict", {} )

          if self.infodict.get("spacegroup_info"):
            spacegroup_info = self.infodict.get("spacegroup_info",False)
            unitcell_info = self.infodict.get("unitcell_info",False)
            htmlstr = '''<html><head/><body><p><span style=" font-weight:600;">Space group: \t
            </span>%s<span style=" font-weight:600;"><br/>Unit cell: \t</span>%s</p></body></html> '''
            self.SpaceGrpUCellText.setText(htmlstr %(spacegroup_info, unitcell_info) )
          self.fileisvalid = True
          #print("ngl_hkl_infodict: " + str(ngl_hkl_infodict))

          if currentinfostr:
            if self.isembedded:
              print(currentinfostr)
            else:
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
        if "Resource temporarily unavailable" not in errmsg: # ignore errors from no connection to ZMQ socket
          print( errmsg  +  traceback.format_exc(limit=10) )


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
    self.rotavecangle_labeltxt.setText("Reflections rotated around Vector with Angle: %3.1fº" %dgr)

    self.ColourMapSelectDlg.selcolmap = self.currentphilstringdict["viewer.color_scheme"]
    if self.currentmillarray_idx is not None:
      arrayinfo = self.array_infotpls[self.currentmillarray_idx]
      self.ColourMapSelectDlg.setDataType(arrayinfo[1])
    self.ColourMapSelectDlg.setPowerScaleSliderVal( self.currentphilstringdict["viewer.color_powscale"] )

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
    self.send_message('NGL.mouse_sensitivity = %f' %val)


  def onMouseSensitivity(self):
    val = self.mousemoveslider.value()/2000.0
    self.mousesensitxtbox.setText("%2.2f" %val )


  def onTooltipAlphaChanged(self, val):
    if self.unfeedback:
      return
    self.ttipalpha = val
    self.send_message('NGL.tooltip_alpha = %f' %val)


  def onShowTooltips(self, val):
    if self.ttipClickradio.isChecked() or val=="click":
      self.send_message("NGL.show_tooltips = click")
      self.ttip_click_invoke = "click"
    if self.ttipHoverradio.isChecked() or val=="hover":
      self.send_message("NGL.show_tooltips = hover")
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
    self.send_message("NGL.fontsize = %d" %val)


  def onClearTextBuffer(self):
    self.textInfo.clear()
    self.infostr = ""


  def onCameraPerspect(self,val):
    if self.cameraPerspectCheckBox.isChecked():
      self.send_message("NGL.camera_type = perspective")
    else:
      self.send_message("NGL.camera_type = orthographic")


  def ExpandRefls(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      if self.ExpandReflsGroupBox.isChecked():
        self.send_message('''
        viewer.expand_to_p1 = True
        viewer.expand_anomalous = True
        viewer.inbrowser = False
                        ''' )
      else:
        self.send_message('''
        viewer.expand_to_p1 = False
        viewer.expand_anomalous = False
        viewer.inbrowser = False
                        ''' )
    else:
      if self.ExpandReflsGroupBox.isChecked():
        self.send_message('''
        viewer.expand_to_p1 = True
        viewer.expand_anomalous = True
        viewer.inbrowser = True
                        ''' )
      else:
        self.send_message('''
        viewer.expand_to_p1 = False
        viewer.expand_anomalous = False
        viewer.inbrowser = True
                        ''' )


  def ExpandToP1(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      if self.expandP1checkbox.isChecked():
        self.send_message('''
        viewer.expand_to_p1 = True
        viewer.inbrowser = False
                        ''' )
      else:
        self.send_message('''
        viewer.expand_to_p1 = False
        viewer.inbrowser = False
                        ''' )
    else:
      if self.expandP1checkbox.isChecked():
        self.send_message('''
        viewer.expand_to_p1 = True
        viewer.inbrowser = True
                        ''' )
      else:
        self.send_message('''
        viewer.expand_to_p1 = False
        viewer.inbrowser = True
                        ''' )


  def ExpandAnomalous(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      if self.expandAnomalouscheckbox.isChecked():
        self.send_message('''
        viewer.expand_anomalous = True
        viewer.inbrowser = False
                        ''' )
      else:
        self.send_message('''
        viewer.expand_anomalous = False
        viewer.inbrowser = False
                        ''' )
    else:
      if self.expandAnomalouscheckbox.isChecked():
        self.send_message('''
        viewer.expand_anomalous = True
        viewer.inbrowser = True
                        ''' )
      else:
        self.send_message('''
        viewer.expand_anomalous = False
        viewer.inbrowser = True
                        ''' )

  def showSysAbsent(self):
    if self.unfeedback:
      return
    if self.sysabsentcheckbox.isChecked():
      self.send_message('viewer.show_systematic_absences = True')
    else:
      self.send_message('viewer.show_systematic_absences = False')


  def showMissing(self):
    if self.unfeedback:
      return
    if self.missingcheckbox.isChecked():
      self.send_message('viewer.show_missing = True')
      self.onlymissingcheckbox.setEnabled(True)
    else:
      self.send_message("""viewer.show_missing = False
                             viewer.show_only_missing = False
                          """)
      self.onlymissingcheckbox.setEnabled(False)


  def showOnlyMissing(self):
    if self.unfeedback:
      return
    if self.onlymissingcheckbox.isChecked():
      self.send_message('viewer.show_only_missing = True')
    else:
      self.send_message('viewer.show_only_missing = False')


  def showSlice(self):
    if self.unfeedback:
      return
    if self.showsliceGroupCheckbox.isChecked():
      #self.ClipPlaneChkGroupBox.setChecked(False)
      self.send_message("""viewer {
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
      self.send_message("""viewer {
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
    self.send_message("""viewer {
    slice_axis = %s
    is_parallel = False
    fixorientation = "%s"
}""" %(self.sliceaxis[i], val))


  def onSliceIndexChanged(self, val):
    if self.unfeedback:
      return
    self.sliceindex = val
    self.send_message("viewer.slice_index = %d" %self.sliceindex)


  def onBindataComboSelchange(self, i):
    if self.BinDataComboBox.currentText():
      binner_idx = self.BinDataComboBox.currentIndex()
      self.send_message('binner_idx = %d' % binner_idx )
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
        self.send_message('NGL.bin_opacities = "%s"' %self.bin_opacities )
      if col==1 and self.binstable_isready: # changing scene_bin_thresholds
        aboveitem = self.binstable.item(row-1, 1)
        belowitem = self.binstable.item(row+1, 1)
        rightitem = self.binstable.item(row, 2)
        if aboveitem is None:
          aboveval = -9e99
        else:
          aboveval = float(aboveitem.text())
        if belowitem is None:
          belowval = 9e99
        else:
          belowval = float(belowitem.text())
          if math.isnan(belowval):
            belowval = float(rightitem.text())
        # the new value must be between above and below values only
        newval = min(belowval, max(aboveval, float(item.text()) ) )
        # but the other way round if binning against resolution
        if self.BinDataComboBox.currentIndex() == 0:
          newval = min(aboveval, max(belowval, float(item.text()) ) )
        self.binstable.item(row,col).setText(str(newval))
        self.lowerbinvals[row] = newval
        allbinvals = self.lowerbinvals + [ self.upperbinvals[-1] ]
        nbins = len(allbinvals)
        self.send_message('''
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
    self.send_message('NGL.bin_opacities = "%s"' %self.bin_opacities)
    self.binstableitemchanges = False
    self.binstable_isready = True


  def onNbinsChanged(self, val):
    if self.unfeedback:
      return
    self.nbins = val
    self.send_message("nbins = %d" %self.nbins)


  def onRadiiScaleChanged(self, val):
    if self.unfeedback:
      return
    self.send_message("viewer.scale = %f" %self.radii_scale_spinBox.value() )


  def onPowerScaleChanged(self, val):
    if self.unfeedback:
      return
    self.send_message("viewer.nth_power_scale_radii = %f" %self.power_scale_spinBox.value() )


  def onManualPowerScale(self, val=None):
    if self.unfeedback:
      return
    if self.ManualPowerScalecheckbox.isChecked():
      self.send_message('viewer.nth_power_scale_radii = %f' %self.power_scale_spinBox.value())
    else:
      self.send_message('viewer.nth_power_scale_radii = -1.0')


  def onShowAllVectors(self):
    if self.ShowAllVectorsBtn.checkState()==Qt.Checked:
      for rvrow in range(self.vectortable2.rowCount()):
        if self.vectortable2.item(rvrow, 0).text() !=  "new vector":
          self.vectortable2.item(rvrow, 0).setCheckState(Qt.Checked)
    if self.ShowAllVectorsBtn.checkState()==Qt.Unchecked:
      for rvrow in range(self.vectortable2.rowCount()):
        self.vectortable2.item(rvrow, 0).setCheckState(Qt.Unchecked)


  def onVectorTableItemChanged(self, item):
    if self.unfeedback:
      return
    row = item.row()
    col = item.column()
    try:
      rc = len(self.all_vectors)
      label = None
      if self.vectortable2.item(row, 0) is not None:
        label = self.vectortable2.item(row, 0).text()
      if row < rc:
        if label is None:
          return
        if col==0:
          if item.checkState()==Qt.Checked:
            self.send_message(" viewer.show_vector = '[%d, %s]'" %(row, True ))
          else:
            self.send_message(" viewer.show_vector = '[%d, %s]'" %(row, False ))
          if self.rotvec is not None: # reset any imposed angle to 0 whenever checking or unchecking a vector
              self.send_message("""clip_plane {
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
            self.send_message("""clip_plane {
  animate_rotation_around_vector = '[%d, %f]'
}""" %(0, -1.0)) #
            self.send_message('viewer.fixorientation = "None"')
            self.AnimaRotCheckBox.setCheckState(Qt.Unchecked)
            self.rotvec = None

          if self.rotvec is not None:
            self.RotateAroundframe.setEnabled(True)
            # notify cctbx which is the curently selected vector
            self.send_message("clip_plane.angle_around_vector = '[%d, 0.0]'" %self.rotvec)
          else:
            self.RotateAroundframe.setDisabled(True)
          if not self.unfeedback:
            if sum >= rc:
              self.ShowAllVectorsBtn.setCheckState(Qt.Checked)
            if sum == 0:
              self.ShowAllVectorsBtn.setCheckState(Qt.Unchecked)
            if sum >0.0 and sum < rc:
              self.ShowAllVectorsBtn.setCheckState(Qt.PartiallyChecked)
      if row==rc and label !="" and label != "new vector": # a user defined vector
        if col==1:
          hklop = self.vectortable2.item(row, 1).text()
          self.send_message("""viewer.add_user_vector_hkl_op = '%s'
          viewer.user_label = %s """ %(hklop, label))
        if col==2:
          hklvec = self.vectortable2.item(row, 2).text()
          self.send_message("""viewer.add_user_vector_hkl = '(%s)'
          viewer.user_label = %s """ %(hklvec, label))
        if col==3:
          abcvec = self.vectortable2.item(row, 3).text()
          self.send_message("""viewer.add_user_vector_abc = '(%s)'
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
    self.yHKLrotBtn.clicked.connect(self.onYangleHKLrotate)
    self.xHKLrotBtn.clicked.connect(self.onXangleHKLrotate)
    self.zHKLrotBtn.clicked.connect(self.onZangleHKLrotate)
    self.yHKLbackrotBtn.clicked.connect(self.onYangleHKLrotateback)
    self.xHKLbackrotBtn.clicked.connect(self.onXangleHKLrotateback)
    self.zHKLbackrotBtn.clicked.connect(self.onZangleHKLrotateback)


  def onXangleHKLrotate(self):
    self.send_message("viewer.angle_around_XHKL_vector = %f" %self.angleStepHKLrotSpinBox.value() )


  def onYangleHKLrotate(self):
    self.send_message("viewer.angle_around_YHKL_vector = %f" %self.angleStepHKLrotSpinBox.value() )


  def onZangleHKLrotate(self):
    self.send_message("viewer.angle_around_ZHKL_vector = %f" %self.angleStepHKLrotSpinBox.value() )


  def onXangleHKLrotateback(self):
    self.send_message("viewer.angle_around_XHKL_vector = %f" %(-1*self.angleStepHKLrotSpinBox.value()) )


  def onYangleHKLrotateback(self):
    self.send_message("viewer.angle_around_YHKL_vector = %f" %(-1*self.angleStepHKLrotSpinBox.value()) )


  def onZangleHKLrotateback(self):
    self.send_message("viewer.angle_around_ZHKL_vector = %f" %(-1*self.angleStepHKLrotSpinBox.value()) )


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
    self.send_message(philstr)


  def onClipPlaneChkBox(self):
    if self.unfeedback:
      return
    hkldist, clipwidth = 0.0, None
    if self.ClipPlaneChkGroupBox.isChecked():
      if self.showsliceGroupCheckbox.isChecked() == False:
        self.showsliceGroupCheckbox.setChecked(False)
        self.showsliceGroupCheckbox.setChecked(False)
        self.send_message("""viewer {
                                                      slice_mode = False
                                                      inbrowser = True
                                                  }
                            """)
      hkldist, clipwidth = self.hkldistval, self.clipwidth_spinBox.value()
    philstr = """clip_plane {
  clipwidth = %s
}
        """ %clipwidth

    self.send_message(philstr)


  def onRotaVecAngleChanged(self, val):
    if self.unfeedback or self.rotvec is None:
      return
    self.send_message("""clip_plane {
    angle_around_vector = '[%d, %f]'
}""" %(self.rotvec, val*0.5))


  def onFinalRotaVecAngle(self):
    if self.unfeedback or self.rotvec is None:
      return
    val = self.rotavecangle_slider.value()*0.5
    self.send_message("""clip_plane {
    angle_around_vector = '[%d, %f]'
}""" %(self.rotvec, val))


  def onAnimateRotation(self):
    if self.AnimaRotCheckBox.isChecked() == True:
      self.AnimateSpeedSlider.setEnabled(True)
      self.rotavecangle_slider.setDisabled(True)
      speed = self.AnimateSpeedSlider.value()
      self.send_message("""clip_plane {
      animate_rotation_around_vector = '[%d, %f]'
}""" %(self.rotvec, speed))
    else:
      self.rotavecangle_slider.setEnabled(True)
      self.AnimateSpeedSlider.setDisabled(True)
      self.send_message("""clip_plane {
      animate_rotation_around_vector = '[%d, %f]'
}""" %(self.rotvec, -1.0))


  def onClipwidthChanged(self, val):
    if not self.unfeedback:
      self.send_message("clip_plane.clipwidth = %f" %self.clipwidth_spinBox.value())


  def onHKLdistChanged(self, val):
    if not self.unfeedback:
      self.hkldistval = val
      self.send_message("clip_plane.hkldist = %f" %self.hkldistval)


  def onHvecChanged(self, val):
    if not self.unfeedback:
      self.send_message("clip_plane.h = %f" %self.hvec_spinBox.value())


  def onKvecChanged(self, val):
    if not self.unfeedback:
      self.send_message("clip_plane.k = %f" %self.kvec_spinBox.value())


  def onLvecChanged(self, val):
    if not self.unfeedback:
      self.send_message("clip_plane.l = %f" %self.lvec_spinBox.value())


  def onPlaneParallelCheckbox(self):
    if not self.unfeedback:
      val = "None"
      if self.PlaneParallelCheckbox.isChecked():
        val = "reflection_slice"
      self.send_message('viewer.fixorientation = "%s"' %val)


  def onMillerTableCellPressed(self, row, col):
    #print( "in millertable CellPressed " + self.millertable.currentItem().text() )
    if self.millertable.mousebutton == Qt.RightButton:
      self.MillerTableContextMenuHandler(QCursor.pos(), row)
    if self.millertable.mousebutton == QEvent.MouseButtonDblClick:
      # quickly display data with a double click
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
          myqa.setData((sceneid, row))
          self.millertablemenu.addAction(myqa)
          myqa = QAction("Display sigmas of %s" %scenelabelstr, self.window, triggered=self.testaction)
          myqa.setData((sceneid + 1000, row)) # want to show the sigmas rather than the data if we add 1000
          self.millertablemenu.addAction(myqa)
        else:
          myqa = QAction("Display %s" %scenelabelstr, self.window, triggered=self.testaction)
          myqa.setData((sceneid, row))
          self.millertablemenu.addAction(myqa)
    myqa = QAction("Make a new dataset from this dataset and another dataset...",
                   self.window, triggered=self.testaction)
    myqa.setData( ("newdata", row ))
    self.millertablemenu.addAction(myqa)

    if len(self.millertable.selectedrows) > 0:
      arraystr = ""
      labels = []
      for i,r in enumerate(self.millertable.selectedrows):
        labels.extend( self.millerarraylabels[r].split(",") ) # to cope with I,SigI or other multiple labels
      myqa = QAction("Show a table of the %s dataset ..." %  " and ".join(labels),
                     self.window, triggered=self.testaction)
      lbls =[] # group any crystal_id=1, wavelength_id, scale_group_code with labels in lists
      for i,r in enumerate(self.millertable.selectedrows):
        lbls.extend( [ self.millerarraylabels[r].split(",") ] ) # to cope with I,SigI or other multiple labels
      myqa.setData( ("tabulate_data", lbls ))
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
          self.operationlabeltxt.setText("Define a new cctbx.miller.array object, \"newarray\", "
           + "by entering a python expression for \"newarray\" or by assigning \"newarray._data\" "
           + "(and optionally \"newarray._sigmas\") to a function of the cctbx.miller.array object, "
           + "\"array1\", representing the " + self.millerarraylabels[idx] + " dataset. Optionally "
           + "also include the cctbx.miller.array object, \"array2\" representing a dataset "
           + "selected from the dropdown list below."
           )
          self.makenewdataform.show()
        if strval=="tabulate_data":
          self.send_message('tabulate_miller_array_ids = "%s"' %str(idx))


  def DisplayData(self, idx, row):
    # want to show the sigmas rather than the data if idx we add 1000
    self.currentmillarray_idx = row
    arrayinfo = self.array_infotpls[self.currentmillarray_idx]
    datatype = arrayinfo[1]
    if (idx - 1000) >= 0:
      idx = idx - 1000
      philstr = """
      viewer.sigma_color_radius = True
      viewer.scene_id = %d
      """ %idx
    else:
      philstr = """
      viewer.sigma_color_radius = False
      viewer.scene_id = %d
      """ %idx
    self.send_message(philstr)
    if self.fileisvalid:
      self.functionTabWidget.setEnabled(True)
      self.expandAnomalouscheckbox.setEnabled(True)
      self.expandP1checkbox.setEnabled(True)
      # don't allow anomalous expansion for data that's already anomalous
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
    mtpl = (self.operationtxtbox.toPlainText(), self.newlabeltxtbox.text() ,
              self.operate_arrayidx1, self.operate_arrayidx2 )
    self.send_message('miller_array_operations = "[ %s ]"' %str(mtpl) )
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
    self.all_vectors = []
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
    self.rotavecangle_labeltxt.setText("Reflections rotated around Vector with Angle: 0º")
    self.rotavecangle_slider.sliderReleased.connect(self.onFinalRotaVecAngle)
    self.rotavecangle_slider.valueChanged.connect(self.onRotaVecAngleChanged)
    self.rotavecangle_slider.setTracking(False)
    self.ClipPlaneChkGroupBox.clicked.connect(self.onClipPlaneChkBox)


  def onSaveImage(self):
    if self.cctbxpythonversion == 'cctbx.python.version: 3': # streaming image over websockets
      options = QFileDialog.Options()
      fileName, filtr = QFileDialog.getSaveFileName(self.window,
              "Save screenshot to file", "",
              "PNG Files (*.png);;All Files (*)", "", options)
      if fileName:
        self.send_message('save_image_name = "%s" '%fileName)
    else:
      self.send_message('save_image_name = "dummy.png" ')
      # eventual file name prompted to us by Browser_download_requested(


  def onDrawReciprocUnitCellBoxClick(self):
    if not self.unfeedback:
      if self.DrawReciprocUnitCellBox.isChecked():
        val = self.reciprocunitcellslider.value()/self.reciprocunitcellslider.maximum()
        self.send_message("reciprocal_unit_cell_scale_fraction = %f" %val)
      else:
        self.send_message("reciprocal_unit_cell_scale_fraction = None")


  def onDrawUnitCellBoxClick(self):
    if not self.unfeedback:
      if self.DrawRealUnitCellBox.isChecked():
        val = self.unitcellslider.value()/self.unitcellslider.maximum()
        self.send_message("real_space_unit_cell_scale_fraction = %f" %val)
      else:
        self.send_message("real_space_unit_cell_scale_fraction = None")


  def onUnitcellScale(self):
    if self.unfeedback:
      return
    val = self.unitcellslider.value()/self.unitcellslider.maximum()
    self.send_message("real_space_unit_cell_scale_fraction = %f" %val)


  def onReciprocUnitcellScale(self):
    if self.unfeedback:
      return
    val = self.reciprocunitcellslider.value()/self.reciprocunitcellslider.maximum()
    self.send_message("reciprocal_unit_cell_scale_fraction = %f" %val)


  def HighlightReflection(self, hkl):
    self.send_message("viewer.show_hkl = '%s'" %str(hkl))


  def DebugInteractively(self):
    import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )


  def SpacegroupSelchange(self,i):
    self.send_message("spacegroup_choice = %d" %i)


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
    assert self.cctbxpython is not None
    # subprocess will not create interactive programs ( using popen.communicate() will simply terminate
    # the subprocess after execution). Since we need cmdlineframes.run() to be interactive
    # we start it with shell=True and flags -i -c for cmdlineframes.run() to remain running.
    # Care must be taken when closing HKLviewer to ensure the shell and its child process are both closed.
    args = ' '.join( guiargs + sys.argv[1:])
    cmdargs =  self.cctbxpython + ' -i -c "from crys3d.hklviewer import cmdlineframes;' \
     + ' cmdlineframes.run()" ' + args
    self.cctbxproc = subprocess.Popen( cmdargs, shell=True,
                                      universal_newlines=True,
                                      stdin=subprocess.PIPE,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    # Wait for connection from the zmq socket in CCTBX by testing if we can send an empty string
    t=0.0; dt = 0.3; timeout = 5
    err = zmq.EAGAIN
    while err == zmq.EAGAIN:
      try:
        err = 0
        self.socket.send(bytes("","utf-8"), zmq.NOBLOCK)
      except Exception as e:
        err = e.errno # if EAGAIN then message cannot be sent at the moment.
      time.sleep(dt)
      t += dt
      if  t > timeout:
        raise Exception("HKLviewer GUI failed making the ZMQ socket connection to the CCTBX.")
        break



  def send_message(self, cmdstr, msgtype="philstr"):
    try:
      msg = str([msgtype, cmdstr])
      if sys.version_info.major==3:
        self.socket.send(bytes(msg,"utf-8"), zmq.NOBLOCK)
      else:
        self.socket.send(bytes(msg), zmq.NOBLOCK)
      return True
    except Exception as e:
      print( str(e) + "\nFailed sending message to the CCTBX\n" + traceback.format_exc(limit=10))
      return False


  def setDatatypedict(self, datatypedict):
    self.datatypedict = datatypedict
    # send persisted colour schemes and raddi mappings to jsview_3d.py
    return self.send_message(str(self.datatypedict), msgtype="dict")


  def PersistQsettings(self):
    self.settings.setValue("PythonPath", self.cctbxpython )
    self.settings.beginGroup(self.Qtversion )
    self.settings.setValue("QWebEngineViewFlags", self.QWebEngineViewFlags)
    self.settings.setValue("FontSize", self.fontsize )
    self.settings.setValue("BrowserFontSize", self.browserfontsize )
    self.settings.setValue("ttip_click_invoke", self.ttip_click_invoke)
    self.settings.setValue("windowsize", self.window.size())
    self.settings.setValue("splitter1Sizes", self.splitter.saveState())
    self.settings.setValue("splitter2Sizes", self.splitter_2.saveState())

    self.settings.beginGroup("DataTypesGroups")
    datatypesgroups = self.settings.childGroups()
    for datatype in list(self.datatypedict.keys()):
      self.settings.setValue(datatype + "/ColourChart", self.datatypedict[ datatype ][0] )
      self.settings.setValue(datatype + "/ColourPowerScale", self.datatypedict[ datatype ][1] )
      self.settings.setValue(datatype + "/PowerScale", self.datatypedict[ datatype ][2])
      self.settings.setValue(datatype + "/RadiiScale", self.datatypedict[ datatype ][3])
    self.settings.endGroup() # DataTypesGroups
    self.settings.endGroup() # PySide2_ + Qtversion


  def ReadPersistedQsettings(self):
    # read the users persisted settings from disc
    # Locate cctbx.python. If not in the Qsettings then try if in the executable path environment
    self.cctbxpython = self.settings.value("PythonPath", "")
    if not os.path.isfile(self.cctbxpython):
      wherecmd = "which"
      if sys.platform == 'win32':
        wherecmd = "where"
      proc = subprocess.Popen([wherecmd, "cctbx.python"],
                              universal_newlines=True, # avoid them annoying byte strings
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE).communicate()
      if proc[0] != "":
        self.cctbxpython = proc[0].strip()
    if not os.path.isfile(self.cctbxpython):
      from .qt import QInputDialog
      self.cctbxpython, ok = QInputDialog.getText(None, "cctbx.python needs specifying",
              'The HKLviewer GUI needs to know where the cctbx.python is located.\n' +
              'Enter the full path for the executable cctbx.python dispatcher file.\n' +
              'Tip: Use the "which" or "where" command from a shell with an active CCTBX environment.')
    if not os.path.isfile(self.cctbxpython):
      raise Exception("The file, %s, does not exists!\n" %self.cctbxpython)
    print("HKLviewer using cctbx.python from: %s" %self.cctbxpython)
    # In case of more than one PySide2 installation tag the settings by version number of PySide2
    # as different versions may use different metrics for font and window sizes
    self.settings.beginGroup(self.Qtversion)
    self.settings.beginGroup("DataTypesGroups")
    datatypes = self.settings.childGroups()
    #datatypedict = { }
    if datatypes:
      for datatype in datatypes:
        self.datatypedict[ datatype ] = [ self.settings.value(datatype + "/ColourChart", "brg"),
                                      float(self.settings.value(datatype + "/ColourPowerScale", 1.0)),
                                      float(self.settings.value(datatype + "/PowerScale", 1.0)),
                                      float(self.settings.value(datatype + "/RadiiScale", 1.0)),
                                    ]
    self.settings.endGroup()
    self.QWebEngineViewFlags = self.settings.value("QWebEngineViewFlags", None)
    self.fontsize = self.settings.value("FontSize", 10)
    self.browserfontsize = self.settings.value("BrowserFontSize", 9)
    self.ttip_click_invoke = self.settings.value("ttip_click_invoke", None)
    self.windowsize = self.settings.value("windowsize", None)
    self.splitter1sizes = self.settings.value("splitter1Sizes", None)
    self.splitter2sizes = self.settings.value("splitter2Sizes", None)
    self.settings.endGroup()
    # test for any necessary flags for WebGL to work on this platform
    if self.QWebEngineViewFlags is None: # avoid doing this test over and over again on the same PC
      self.QWebEngineViewFlags = " --disable-web-security" # for chromium
      if not self.isembedded:
        print("testing if WebGL works in QWebEngineView....")
        QtChromiumCheck_fpath = os.path.join(os.path.split(HKLviewerGui.__file__)[0], "QtChromiumCheck.py")
        cmdargs = [ sys.executable, QtChromiumCheck_fpath ]
        webglproc = subprocess.Popen( cmdargs, stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        procout, procerr = webglproc.communicate()
        if not "WebGL works" in procout.decode():
          self.QWebEngineViewFlags = " --enable-webgl-software-rendering --ignore-gpu-blacklist "
    if "verbose" in sys.argv[1:]:
      print("using flags for QWebEngineView: " + self.QWebEngineViewFlags)
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] += self.QWebEngineViewFlags


  def UsePersistedQsettings(self):
    # Now assign the users persisted settings to the GUI
    if self.fontsize is not None:
      self.onFontsizeChanged(int(self.fontsize))
      self.fontspinBox.setValue(int(self.fontsize))
    if self.browserfontsize is not None:
      self.onBrowserFontsizeChanged(int(self.browserfontsize))
      self.browserfontspinBox.setValue(int(self.browserfontsize))
    if self.ttip_click_invoke is not None:
      self.onShowTooltips(self.ttip_click_invoke)
      self.ttipClickradio.setChecked(self.ttip_click_invoke == "click")
      self.ttipHoverradio.setChecked(self.ttip_click_invoke == "hover")
    if self.splitter1sizes is not None and self.splitter2sizes is not None and self.windowsize is not None:
      self.window.resize(self.windowsize)
      if self.webpagedebugform and self.devmode:
        self.webpagedebugform.resize( self.window.size())
      self.splitter.restoreState(self.splitter1sizes)
      self.splitter_2.restoreState(self.splitter2sizes)
    self.setDatatypedict(self.datatypedict)


  @staticmethod
  def RemoveQsettings(all=False):
    mstr = NGL_HKLViewer.Qtversion
    if all:
      mstr = ""
    NGL_HKLViewer.settings.remove(mstr)


def run(isembedded=False, chimeraxsession=None):
  import time
  #time.sleep(10) # enough time for attaching debugger
  try:
    debugtrue = False
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = " "
    for e in sys.argv:
      if "devmode" in e or "debug" in e and not "UseOSBrowser" in e:
        debugtrue = True
        # some useful flags as per https://doc.qt.io/qt-5/qtwebengine-debugging.html
        if "debug" in e:
          os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--remote-debugging-port=9741 --single-process --js-flags='--expose_gc'"
        if "devmode" in e:  # --single-process will freeze the WebEngineDebugForm at breakpoints
          os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--js-flags='--expose_gc'"
      if "remove_settings" in e:
        NGL_HKLViewer.RemoveQsettings()
        print("HKLViewer settings for %s removed." %NGL_HKLViewer.Qtversion)
        exit()
      if "remove_all_settings" in e:
        NGL_HKLViewer.RemoveQsettings(all=True)
        print("All HKLViewer settings removed.")
        exit()

    from .qt import QApplication
    # ensure QWebEngineView scales correctly on a screen with high DPI
    if not isembedded:
      QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)

    HKLguiobj = NGL_HKLViewer(app, isembedded)

    if not isembedded:
      timer = QTimer()
      timer.setInterval(20)
      timer.timeout.connect(HKLguiobj.ProcessMessages)
      timer.start()
    else:
      start_time = [time.time()]

      def ChXTimer(trigger, trigger_data):
        elapsed_time = time.time()-start_time[0]
        if elapsed_time > 0.02:
          start_time[0] = time.time()
          HKLguiobj.ProcessMessages()

      HKLguiobj.chimeraxprocmsghandler = chimeraxsession.triggers.add_handler('new frame', ChXTimer)

    HKLguiobj.UsePersistedQsettings()

    if isembedded:
      return HKLguiobj

    ret = app.exec_()

  except Exception as e:
    print( str(e)  +  traceback.format_exc(limit=10) )
