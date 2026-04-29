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

import sys, zmq, subprocess, time, traceback, zlib, io, os, math, os.path, re
if sys.version_info[0] < 3:
  print("HKLviewer GUI must be run from Python 3")
  sys.exit(-42)

os.environ['QT_MAC_WANTS_LAYER'] = '1'

from .qt import Qt, QtCore, QCoreApplication, QEvent, QItemSelectionModel, QSize, QSettings, QTimer, QUrl
from .qt import (  QAction, QAbstractScrollArea, QCheckBox, QColorDialog, QComboBox, QDialog, QDoubleSpinBox,
    QFileDialog, QFrame, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QInputDialog, QLabel, QLineEdit, QCloseEvent,
    QMainWindow, QMenu, QMenuBar, QMessageBox, QPalette, QPlainTextEdit, QProgressBar, QPushButton, QRadioButton, QRect,
    QScrollBar, QSizePolicy, QSlider, QSpinBox, QSplitter, QStyleFactory, QStatusBar, QTableView, QTableWidget,
    QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget )

from .qt import QColor, QFont, QIcon, QCursor, QDesktopServices
from .qt import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )


from . import hklviewer_gui
from .helpers import ( MillerArrayTableView, MillerArrayTableForm, MyhorizontalHeader, MyQPlainTextEdit,
    MillerArrayTableModel, MPLColourSchemes, MillerTableColumnHeaderDialog, MyQDoubleSpinBox, FindDialog )


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
    self.aboutlabel.setMaximumSize(QSize(16777215, 16777215))
    self.writeAboutstr("")
    self.copyrightstxt = QTextEdit()
    self.copyrightstxt.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    self.copyrightstxt.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    self.copyrightstxt.setReadOnly(True)
    self.copyrightstxt.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
    self.copyrightstxt.setMaximumSize(QSize(16777215, 16777215))
    sp = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
    self.copyrightstxt.setMinimumSize(QSize(350, 150))
    self.copyrightstxt.setSizePolicy(sp)
    self.aboutlabel.setSizePolicy(sp)
    self.OKbtn = QPushButton()
    self.OKbtn.setText("OK")
    self.OKbtn.clicked.connect(self.onOK)
    mainLayout.addWidget(self.aboutlabel,  0, 0, 1, 3)
    mainLayout.addWidget(self.copyrightstxt,  1, 0, 1, 3)
    mainLayout.addWidget(self.OKbtn,  2, 1, 1, 1)
    self.setLayout(mainLayout)
    self.setMinimumSize(QSize(400, 300))
    self.setSizePolicy(sp)
    self.setFixedSize( self.sizeHint() )
  def writeAboutstr(self, versionstr):
    aboutstr = """<html><head/><body><p>
    <span style=" font-weight:600;">HKLviewer, </span>
    CCTBX version: """ + versionstr + \
    """
    <br/>A reflection data viewer for crystallography
    <br/>Developers: Robert D. Oeffner<br/>
    Cambridge Institute for Medical Research, University of Cambridge.<br/>
    HKLviewer is part of the <a href="http://cci.lbl.gov/docs/cctbx/"> CCTBX library</a>
    as well as derived software thereof.<br/>
    HKLviewer uses functionality provided by the
    <a href="https://github.com/nglviewer/ngl">NGL Viewer</a> project and
    the <a href="https://github.com/niklasvh/html2canvas">html2canvas</a> project.
    Refer to robert@oeffner.net or cctbx@cci.lbl.gov for queries or bug reports.
    </p></body></html>"""
    self.aboutlabel.setText(aboutstr)
  def onOK(self):
    self.hide()


class SettingsForm(QDialog):
  def __init__(self, parent=None):
    super(SettingsForm, self).__init__(parent.window)
    self.setWindowTitle("HKLviewer Settings")
    self.setWindowFlags(Qt.Tool)
    myGroupBox = QGroupBox()
    layout = QGridLayout()
    layout.addWidget(parent.mousespeed_labeltxt,     0, 0, 1, 1)
    layout.addWidget(parent.mousemoveslider,         0, 1, 1, 3)
    layout.addWidget(parent.mousesensitxtbox,        0, 4, 1, 1)
    layout.addWidget(parent.Fontsize_labeltxt,       1, 0, 1, 1)
    layout.addWidget(parent.fontspinBox,             1, 4, 1, 1)
    layout.addWidget(parent.BrowserFontsize_labeltxt, 2, 0, 1, 1)
    layout.addWidget(parent.browserfontspinBox,      2, 4, 1, 1)
    layout.addWidget(parent.vectorWidth_labeltxt,    3, 0, 1, 1)
    layout.addWidget(parent.vectorWidthspinBox,      3, 4, 1, 1)
    layout.addWidget(parent.cameraPerspectCheckBox,  4, 0, 1, 1)
    layout.addWidget(parent.bufsize_labeltxt,        5, 0, 1, 1)
    layout.addWidget(parent.bufsizespinBox,          5, 4, 1, 1)
    layout.addWidget(parent.clearbufbtn,             6, 0, 1, 1)
    layout.addWidget(parent.wraptextbtn,             6, 2, 1, 2)
    layout.addWidget(parent.ttiplabeltxt,            7, 0, 1, 1)
    layout.addWidget(parent.ttipClickradio,          7, 1, 1, 1)
    layout.addWidget(parent.ttipHoverradio,          7, 2, 1, 1)
    layout.addWidget(parent.ttipalphalabeltxt,       7, 3, 1, 1)
    layout.addWidget(parent.ttipalpha_spinBox,       7, 4, 1, 1)

    layout.setRowStretch (0, 1)
    layout.setRowStretch (1 ,0)
    myGroupBox.setLayout(layout)
    mainLayout = QGridLayout()
    mainLayout.addWidget(myGroupBox,                       0, 0, 1, 3)
    mainLayout.addWidget(parent.resetlabeltxt,             1, 0, 1, 2)
    mainLayout.addWidget(parent.resetFactoryDefaultbtn,    1, 2, 1, 1)
    mainLayout.addWidget(parent.showphillabeltxt,          2, 0, 1, 2)
    mainLayout.addWidget(parent.showphilbtn,               2, 2, 1, 1)
    self.setLayout(mainLayout)
    self.setFixedSize( self.sizeHint() )


class WebEngineDebugForm(QDialog):
  def __init__(self, parent=None):
    super(WebEngineDebugForm, self).__init__(None,
                        Qt.WindowMinimizeButtonHint | # Want minimise and maximise buttons on window.
                        Qt.WindowMaximizeButtonHint | # As they are not the default for QDialog we
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

"""
# hklviewer_gui.py should be generated from the Qtdesigner output file, HKLviewer4.ui, with the commandline:

"C:\\Program Files\\Python38\\Lib\\site-packages\\PySide2\\uic.exe" HKLviewer6.ui -o hklviewer_gui.py -g python

# Replace the autogenerated imports in hklviewer_gui.py with:

from .qt import QWebEngineView
try: # if invoked by cctbx.python or some such
  from crys3d.hklviewer.helpers import HeaderDataTableWidget, MyQDoubleSpinBox # implicit import
except Exception as e: # if invoked by a generic python that doesn't know cctbx modules
  from .helpers import HeaderDataTableWidget, MyQDoubleSpinBox # implicit import

from .qt import ( QCoreApplication, QMetaObject, QRect, QSize, Qt,  # implicit import
 QFont, QAbstractItemView, QAction, QCheckBox, QComboBox, QLineEdit, QDockWidget,
 QDoubleSpinBox, QFrame, QGridLayout, QGroupBox, QLabel, QPlainTextEdit,
 QPushButton, QRadioButton, QScrollArea, QSlider, QSplitter, QSizePolicy, QSpinBox,
 QTableWidget, QTabWidget, QTextEdit, QWidget, QIcon, QAbstractScrollArea, )

# For allowing embedding in ChimeraX comment out the line:

MainWindow.setCentralWidget(self.centralwidget)

from  the function hklviewer_gui.Ui_MainWindow.setupUi()

"""
timer = QTimer()

class NGL_HKLViewer(hklviewer_gui.Ui_MainWindow):
  # qversion() comes out like '5.12.5'. We just want '5.12'
  Qtversion = "Qt" + ".".join( QtCore.qVersion().split(".")[0:2])
  settings = QSettings("CCTBX", "HKLviewer" )
  reset_to_factorydefaults = False

  def __init__(self, thisapp, isembedded=False): #, cctbxpython=None):
    self.datatypedict = { }
    self.browserfontsize = None
    self.vectorwidth = None
    self.mousespeedscale = 2000
    self.isembedded = isembedded
    self.philfname = "" # for regression tests
    self.image_fname = "testimage.png"
    print("version " + self.Qtversion)
    self.colnames_select_dict = {}
    self.lasttime = time.monotonic()
    self.lastwaittime = time.monotonic()
    self.backgroundcolour = QColor(127,127,127)

    if isembedded:
      self.window = MyQMainDialog(self)
      self.window.hide()
    else:
      self.window = MyQMainWindow(self)
    self.setupUi(self.window)
    if isembedded: # in ChimeraX
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
      self.menuFile.addAction(self.actionBackground_Colour)
      self.menuFile.addAction(self.actionSave_Current_Image)
      self.menuFile.addAction(self.actionExit)
      self.menuHelp.addAction(self.actionLocal_Help)
      self.menuHelp.addAction(self.actionHKLviewer_Tutorial)
      self.menuHelp.addAction(self.actionCCTBXwebsite)
      self.menuHelp.addAction(self.actionNGLmousebindings)
      self.menuHelp.addAction(self.actionAbout)
      self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
      self.menuHelp.setTitle(QCoreApplication.translate("MainWindow", u"Help", None))
    self.nsplitters = 0
    for a in dir(self):
      if isinstance( self.__getattribute__(a), QSplitter):
        self.nsplitters += 1
    self.ntabs = 0
    for a in dir(self):
      if isinstance( self.__getattribute__(a), QTabWidget):
        self.ntabs += 1
    self.factorydefaultfname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "HKLviewerDefaults.ini")
    if not self.ReadPersistedQsettings():
      sys.exit()
    self.buttonsdeflist =[]
    self.app = thisapp
    self.cctbxversion = "unversioned"
    self.actiondebug.setVisible(False)
    self.UseOSBrowser = False
    self.devmode = False
    self.make_new_factory_default_settings = False
    self.XtricorderBtn.setVisible(False)
    self.XtriageBtn.setVisible(True)
    for e in sys.argv:
      if "UseOSBrowser" in e:
        self.UseOSBrowser = True
      if "devmode" in e:
        self.devmode = True
      if  "debug" in e or "devmode" in e:
        self.actiondebug.setVisible(True)
      if "new_factory_defaults" in e:
        self.make_new_factory_default_settings= True

    self.zmq_context = None
    self.unfeedback = False
    self.cctbxpythonversion = None
    self.mousespeed_labeltxt = QLabel()
    self.mousespeed_labeltxt.setText("Mouse speed:")
    self.mousemoveslider = QSlider(Qt.Horizontal)
    self.mousemoveslider.setMinimum(1)
    self.mousemoveslider.setMaximum(200)
    self.mousemoveslider.setValue(float(self.mousespeed)*self.mousespeedscale)
    self.mousemoveslider.sliderReleased.connect(self.onFinalMouseSensitivity)
    self.mousemoveslider.valueChanged.connect(self.onMouseSensitivity)
    self.mousesensitxtbox = QLineEdit('')
    self.mousesensitxtbox.setReadOnly(True)

    self.fontspinBox = MyQDoubleSpinBox()
    self.fontspinBox.setSingleStep(1)
    self.fontspinBox.setRange(4, 50)
    self.font = QFont()
    self.font.setFamily(self.font.defaultFamily())
    self.fontspinBox.setValue(self.font.pointSize())
    self.fontspinBox.editingFinished.connect(self.onFontsizeChanged)
    # valueChanged signal is buggy and gets triggered twice when onFontsizeChanged takes too long
    # This may cause font step=2 rather than step=1
    # Workaround is to use onMouseRelease invoked by MyQDoubleSpinBox.mouseReleaseEvent
    self.fontspinBox.onMouseRelease = self.onFontsizeChanged
    self.texttabfont = QFont("Courier New")
    self.texttabfont.setPointSize(self.font.pointSize())
    self.texttabfont.setBold(True)
    self.finddlg = FindDialog(self.tabText) # show near tabText when making finddlg visible
    self.textAlerts.finddlg = self.finddlg
    self.textInfo.finddlg = self.finddlg
    self.Fontsize_labeltxt = QLabel()
    self.Fontsize_labeltxt.setText("Font size:")
    self.waiting = False

    self.browserfontspinBox = QDoubleSpinBox()
    self.browserfontspinBox.setSingleStep(1)
    self.browserfontspinBox.setRange(4, 50)
    self.browserfontspinBox.setValue(self.font.pointSize())
    self.browserfontspinBox.valueChanged.connect(self.onBrowserFontsizeChanged)
    self.BrowserFontsize_labeltxt = QLabel()
    self.BrowserFontsize_labeltxt.setText("Browser font size:")

    self.vectorWidthspinBox = QSpinBox()
    self.vectorWidthspinBox.setSingleStep(1)
    self.vectorWidthspinBox.setRange(1, 30)
    self.vectorWidthspinBox.setValue(self.vectorwidth)
    self.vectorWidthspinBox.valueChanged.connect(self.onVectorWidthChanged)
    self.vectorWidth_labeltxt = QLabel()
    self.vectorWidth_labeltxt.setText("Thickness of vectors and axes:")

    self.cameraPerspectCheckBox = QCheckBox()
    self.cameraPerspectCheckBox.setText("Perspective camera")
    self.cameraPerspectCheckBox.clicked.connect(self.onCameraPerspect)
    self.cameraPerspectCheckBox.setCheckState(Qt.Unchecked)

    self.bufsizespinBox = QSpinBox()
    self.bufsizespinBox.setSingleStep(1)
    self.bufsizespinBox.setRange(1, 100)
    self.bufsizespinBox.setValue(int(self.textinfosize))
    self.bufsizespinBox.valueChanged.connect(self.onTextbufferSizeChanged)
    self.bufsize_labeltxt = QLabel()
    self.bufsize_labeltxt.setText("Text buffer size (Kbytes):")
    self.clearbufbtn = QPushButton()
    self.clearbufbtn.setText("Clear all text")
    self.clearbufbtn.clicked.connect(self.onClearTextBuffer)
    self.wraptextbtn = QCheckBox()
    self.wraptextbtn.setText("Word wrap text")
    self.wraptextbtn.clicked.connect(self.onWrapTextBuffer)
    self.wraptextbtn.setCheckState(Qt.Checked if self.wraptextinfo else Qt.Unchecked)
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
    self.resetlabeltxt = QLabel()
    self.resetlabeltxt.setWordWrap(True)
    self.resetlabeltxt.setText("Delete user settings and revert to factory defaults for GUI, colour and radii scheme")
    self.resetFactoryDefaultbtn = QPushButton()
    self.resetFactoryDefaultbtn.setText("Reset Settings")
    self.resetFactoryDefaultbtn.clicked.connect(self.onResetFactoryDefault)
    self.showphillabeltxt = QLabel()
    self.showphillabeltxt.setWordWrap(True)
    self.showphillabeltxt.setText("Show current non-default phil parameters")
    self.showphilbtn = QPushButton()
    self.showphilbtn.setText("Show phil")
    self.showphilbtn.clicked.connect(self.onDebugShowPhil)
    # Set the rich text of the SpaceGrpUCellText here rather than in QtDesigner which on windows
    # include MS Font in it. MS Font are not understood by MacOS
    htmlstr = '''<html><head/><body><p><span style=" font-weight:600;">Space group: \t
    </span>P21 21 21<span style=" font-weight:600;"><br/>Unit cell(s): \t</span>(98.371, 98.371, 263.131, 90, 90, 120)</p></body></html> '''
    self.SpaceGrpUCellText.setText(htmlstr )
    self.ColourMapSelectDlg = MPLColourSchemes(self)
    self.select_millertable_column_dlg = MillerTableColumnHeaderDialog(self)
    self.ColourMapSelectDlg.setWindowTitle("HKLviewer Colour Gradient Maps")
    self.ColourMapSelectDlg.hide()
    # colour schemes and radii mapping for types of datasets stored in jsview_3d.py but persisted here:
    # colourmap=brg, colourpower=1, powerscale=1, radiiscale=1
    self.settingsform = SettingsForm(self)
    self.aboutform = AboutForm(self)
    self.webpagedebugform = None

    self.BgrndColourDlg = QColorDialog(self.window)
    self.BgrndColourDlg.setOption(QColorDialog.NoButtons)
    self.BgrndColourDlg.setWindowTitle("HKLviewer Background Colour")
    self.BgrndColourDlg.currentColorChanged.connect(self.onBackgroundColourChanged)

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
    self.operationtxtbox.setAcceptRichText(False)
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
    self.ano_spg_tpls = []
    self.colnames_select_lst = []
    self.currentmillarray_idx = None
    self.matching_arrays = []
    self.bin_infotpls = None
    self.bin_opacities= None
    self.nbins = 1
    self.lowerbinvals = []
    self.upperbinvals = []
    self.html_url = None
    self.spacegroups = {}
    self.info = []
    self.current_labels = None
    self.infostr = ""
    self.alertstr = ""
    self.fileisvalid = False
    self.currentfileName = None
    self.NewFileLoaded = False
    self.NewMillerArray = False
    self.NewHKLscenes = False
    self.binstableitemchanges = False
    self.canexit = False
    self.filterlst = ["MTZ Files (*.mtz)", "CIF Files (*.cif)", "HKL Files (*.hkl)",
                      "SCA Files (*.sca)", "All Files (*)"]
    self.ipresetbtn = -1
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
    self.chimeraxsession = None
    self.XtricorderBtn.clicked.connect(self.onXtricorderRun)
    self.XtriageBtn.clicked.connect(self.onXtriageRun)
    self.tabText.setCurrentIndex(0)
    self.tabText.currentChanged.connect( self.onTabtextChanged )
    if not isembedded:
      self.window.statusBar().showMessage("")
      self.hklLabel = QLabel()
      self.hklLabel.setText("X,Y,Z axes as HKL coordinates:")
      self.Statusbartxtbox = QLineEdit('')
      self.Statusbartxtbox.setReadOnly(True)
      self.Statusbartxtbox.setAlignment(Qt.AlignRight)
      self.window.statusBar().addPermanentWidget(self.hklLabel)
      self.window.statusBar().addPermanentWidget(self.Statusbartxtbox, 1)
      self.actionOpen_reflection_file.triggered.connect(self.onOpenReflectionFile)
      self.actionLocal_Help.triggered.connect(self.onOpenHelpFile)
      self.actionHKLviewer_Tutorial.triggered.connect(self.onOpenTutorialHelpFile)
      self.actionCCTBXwebsite.triggered.connect(self.onOpenCCTBXwebsite)
      self.actionNGLmousebindings.triggered.connect(self.onNGLmouseControlswebsite)
      self.actiondebug.triggered.connect(self.DebugInteractively)
      self.actionSave_Current_Image.triggered.connect(self.onSaveImage)
      self.actionSave_Current_Image.setDisabled(True)
      self.actionSettings.triggered.connect(self.SettingsDialog)
      self.actionAbout.triggered.connect(self.aboutform.show)
      self.actionExit.triggered.connect(self.window.close)
      self.actionSave_reflection_file.triggered.connect(self.onSaveReflectionFile)
      self.actionSave_reflection_file.setDisabled(True)
      self.actionColour_Gradient.triggered.connect(self.ColourMapSelectDlg.show)
      self.actionBackground_Colour.triggered.connect(self.BgrndColourDlg.show)
    else:
      self.tabText.setVisible(False) # stdout sent to chimeraX's console instead
    self.functionTabWidget.setCurrentIndex(0) # if accidentally set to a different tab in the Qtdesigner
    self.tabWidget.setCurrentIndex(0) # if accidentally set to a different tab in the Qtdesigner
    self.window.show()


  def onOpenHelpFile(self):
    QDesktopServices.openUrl(QUrl("https://cci.lbl.gov/docs/cctbx/doc_hklviewer/"))


  def onOpenTutorialHelpFile(self):
    QDesktopServices.openUrl(QUrl("https://cci.lbl.gov/docs/cctbx/tuto_hklviewer/"))


  def onOpenCCTBXwebsite(self):
    QDesktopServices.openUrl(QUrl("https://cci.lbl.gov/docs/cctbx/"))


  def onNGLmouseControlswebsite(self):
    QDesktopServices.openUrl(QUrl("https://nglviewer.org/ngl/api/manual/interaction-controls.html#controls"))


  def onPresetbtn_click(self):
    for i,((btnname, label, philstr), isenabled, datalabel, moniker_veclabels) in enumerate(self.buttonsdeflist):
      vectorscombobox = None
      #if moniker_veclabels:
      if moniker_veclabels is not None and moniker_veclabels[0] != "" and len(moniker_veclabels[1]):
        vectorscombobox = self.__getattribute__(btnname + "_vectors")
        vectorscombobox.setEnabled(False)
      if self.__getattribute__(btnname).isChecked():
        if vectorscombobox is not None:
          vectorscombobox.setEnabled(True)
          if vectorscombobox.currentIndex() > -1:
            currentvecname = vectorscombobox.currentText()
            philstr = re.sub(moniker_veclabels[0], currentvecname, philstr)
        self.send_message(philstr, msgtype = "preset_philstr")
        self.ipresetbtn = i
      else:
        if vectorscombobox is not None:
          vectorscombobox.setEnabled(False)


  def onVectorsComboSelchange(self, i):
    self.onPresetbtn_click()


  def closeEvent(self, event=QCloseEvent()): # provide default value when called explicitly below
    self.send_message('action = is_terminating')
    self.closing = True
    self.finddlg.setVisible(False)
    self.settingsform.setVisible(False)
    self.aboutform.setVisible(False)
    self.millerarraytableform.setVisible(False)
    self.ColourMapSelectDlg.setVisible(False)
    self.select_millertable_column_dlg.setVisible(False)
    self.ControlsWidget.setVisible(False)
    self.InfoWidget.setVisible(False)
    self.window.setVisible(False)
    self.BgrndColourDlg.close()

    if self.UseOSBrowser == False:
      self.webpage.deleteLater() # avoid "Release of profile requested but WebEnginePage still not deleted. Expect troubles !"
    print("HKLviewer closing down")
    nc = 0
    sleeptime = 0.2
    maxtime = 3
    while not self.canexit and nc < maxtime: # until cctbx.python has finished or after maxtime sec
      time.sleep(sleeptime)
      print(".", end='', flush=True)
      self.ProcessMessages()
      nc += sleeptime
    try:
      #self.cctbxproc.terminate()
      self.out, self.err = self.cctbxproc.communicate(input="exit()", timeout=maxtime)
      print(str(self.out) + "\n" + str(self.err))
    except Exception as e:
      print("\nExterminating unresponsive cctbx.python process, unconditionally, at will, with impunity, effective immediately!")
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
    if not self.reset_to_factorydefaults:
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


  def setWindowFilenameTitles(self, fname):
    self.window.setWindowTitle("HKLviewer: " + fname)
    self.ControlsWidget.setWindowTitle("HKLviewer Controls: " + fname)
    self.InfoWidget.setWindowTitle("HKLviewer Info: " + fname)
    self.current_labels = None
    self.textInfo.setPlainText("")
    self.textAlerts.setPlainText("")
    self.fileisvalid = False
    self.MillerComboBox.clear()
    self.millertable.clearContents()
    self.expandP1checkbox.setChecked(False)
    self.expandAnomalouscheckbox.setChecked(False)
    self.sysabsentcheckbox.setChecked(False)
    self.missingcheckbox.setChecked(False)
    self.onlymissingcheckbox.setChecked(False)


  def onOpenReflectionFile(self):
    options = QFileDialog.Options()
    self.currentfileName, filtr = QFileDialog.getOpenFileName(self.window,
            "Open a reflection file", "",
            ";;".join(self.filterlst), "", options)
    if self.openReflectionFile():
      # Rearrange filters to use the current filter as default for next time we open a file
      # The default filter is the first one. So promote the one chosen to be first in the list for next time
      idx = self.filterlst.index(filtr) # find its place in the list of filters
      self.filterlst.pop(idx ) # remove it from the list
      # make a new list with it as the first element
      newfilterlst = [filtr]
      newfilterlst.extend(self.filterlst)
      self.filterlst = newfilterlst


  def openReflectionFile(self):
    if self.currentfileName:
      self.setWindowFilenameTitles( self.currentfileName)
      self.send_message('openfilename = "%s"' %self.currentfileName )
      return True
    return False


  def onSaveReflectionFile(self):
    if len(self.millertable.selectedrows) ==0:
      QMessageBox.warning(self.window, "HKLviewer",
        "First highlight one or more datasets in the table of datasets on the details tab to save them as a new datafile", buttons=QMessageBox.Ok)
      return
    datasets2savestr = "', '".join([ self.millertable.item(i, 0).text() for i in self.millertable.selectedrows ] )
    options = QFileDialog.Options()
    fileName, filtr = QFileDialog.getSaveFileName(self.window,
            "Save datasets '%s' to a reflection file" %datasets2savestr, "",
            "MTZ Files (*.mtz);;CIF Files (*.cif);;All Files (*)", "", options)
    if fileName:
      self.send_message('''
savefilename = "%s"
datasets_to_save = %s'''
                        %(fileName, " ".join([ str(e) for e in self.millertable.selectedrows] ))
                        )


  def SettingsDialog(self):
    self.settingsform.show()
    # don't know why valueChanged.connect() method only takes effect from here on
    #self.fontspinBox.valueChanged.connect(self.onFontsizeChanged)
    self.settingsform.activateWindow()


  def onColourChartSelect(self, selcolmap, colourpowscale):
    # called when user clicks OK in helpers.MPLColourSchemes.EnactColourMapSelection()
    if selcolmap != "":
      self.send_message("""hkls.color_scheme = %s
hkls.color_powscale = %s""" %(selcolmap, colourpowscale) )


  def onBackgroundColourChanged(self, color):
    self.backgroundcolour = color
    pal = QPalette()
    pal.setColor(self.BrowserBox.backgroundRole(), self.backgroundcolour)
    self.BrowserBox.setPalette(pal)
    self.send_message("NGL.background_colour = 'rgb(%d, %d, %d)'" %(color.red(), color.green(), color.blue()) )


  def onSelect_millertable_column_dlg(self):
    """
    Dialog for choosing what columns of the miller table should be displayed
    Invoked by helpers.MyhorizontalHeader()
    """
    self.select_millertable_column_dlg.show()
    self.select_millertable_column_dlg.activateWindow()


  def SaveImage(self):
    self.send_message('save_image_name = "%s"' %self.image_fname)


  def SetFirstScene(self):
    self.send_message("viewer.scene_id = 0")


  def SetStateFromPHILfile(self):
    with open(self.philfname, "r") as f:
      self.AddAlertsText("Processing PHIL file: %s\n" %self.philfname)
      PHILstr = f.read()
      self.send_message(PHILstr, msgtype = "preset_philstr")


  def onXtricorderRun(self):
    if not self.currentfileName:
      QMessageBox.warning(self.window, "HKLviewer",
        "To run Xtricorder you must first open a datafile in the HKLviewer", buttons=QMessageBox.Ok)
      return
    if "_xtricorder.mtz" in self.currentfileName:
      self.AddInfoText("File looks like it has already been processed by Xtricorder. Try loading another reflection file.\n")
      return
    if "expanded" in self.currentfileName:
      self.AddInfoText("Datasets appear to have been expanded to a subgroup. Save these to a new file and load that to run Xtricorder\n")
      return
    self.XtricorderBtn.setEnabled(False)
    self.XtriageBtn.setEnabled(False)
    self.send_message("external_cmd = 'runXtricorder'" )
    self.AddInfoText("Running Xtricorder")
    self.AddAlertsText("Running Xtricorder")
    self.waiting = True


  def onXtriageRun(self):
    if not self.current_labels:
      QMessageBox.warning(self.window, "HKLviewer",
        "To run Xtriage you must first display an observation dataset", buttons=QMessageBox.Ok)
      return
    if "expanded" in self.currentfileName:
      self.AddInfoText("Datasets appear to have been expanded to a subgroup. Save these to a new file and load that to run Xtriage\n")
      return
    self.XtricorderBtn.setEnabled(False)
    self.XtriageBtn.setEnabled(False)
    self.send_message("external_cmd = 'runXtriage'" )
    self.AddInfoText("Running Xtriage")
    self.AddAlertsText("Running Xtriage")
    self.waiting = True


  def ProcessMessages(self):
    """
    Deal with the messages posted to this GUI by hklview_frame.HKLViewFrame.SendInfoToGUI()
    """
    if self.webpagedebugform is not None:
      try: # During shutdown this may fail if webpagedebugform exits before message handler terminates
        self.webpagedebugform.update()
      except Exception as e:
        pass
    if self.zmq_context:
      if self.cctbxproc.poll() is not None and not self.closing:
        print("Critical Error: CCTBX process has terminated. Please restart HKLviewer.", flush=True)
        timer.stop()
        self.closeEvent()

      if (time.monotonic() - 5) > self.lasttime: # send Isoldes clipper data every 5 sec
        self.lasttime = time.monotonic()
        if self.chimeraxsession is not None and self.chimeraxsession.HKLviewer is not None \
         and hasattr(self.chimeraxsession, "isolde"):
          self.chimeraxsession.HKLviewer.isolde_clipper_data_to_dict()
          self.send_message(str(self.chimeraxsession.HKLviewer.clipper_crystdict),
                            msgtype="clipper_crystdict")

      if self.waiting and (time.monotonic() - 0.5) > self.lastwaittime:
        self.lastwaittime = time.monotonic() # reassure the user we have not crashed
        self.AddInfoText(".")
        self.AddAlertsText(".")

      try:
        binmsg = self.socket.recv(flags=zmq.NOBLOCK) #To empty the socket from previous messages
        msg = zlib.decompress(binmsg)
        nan = float("nan") # workaround for "evaluating" any NaN or inf values in the messages received
        inf = math.inf

        msgstr = msg.decode()

        if "cctbx.python.version:" in msgstr:
          self.cctbxpythonversion = msgstr
          self.send_message("""NGL {
  fontsize = %s
  vector_width = %s
  show_tooltips = %s
}
""" %(self.browserfontsize, self.vectorwidth, self.ttip_click_invoke) )

          if self.cctbxpythonversion == 'cctbx.python.version: 2':
            # use NGL's download feature for images since websocket_server fails to handle large streams
            self.webpage.profile().downloadRequested.connect(self.Browser_download_requested)
          return
        self.infodict = eval(msgstr)
        if self.infodict:

          if self.infodict.get("AddXtricorderButton", False):
            self.XtricorderBtn.setVisible(True)

          if self.infodict.get("closing_time"): # notified by cctbx in regression tests
            QTimer.singleShot(10000, self.closeEvent )

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
            self.cctbxversion = self.infodict["cctbxversion"]
            self.aboutform.writeAboutstr( self.cctbxversion )

          if self.infodict.get("scene_array_label_types"):
            self.scenearraylabeltypes = self.infodict.get("scene_array_label_types", [])

          if self.infodict.get("array_infotpls", None) is not None:
            self.array_infotpls = self.infodict["array_infotpls"]
            self.millerarraylabels =  [e[1][0] for e in self.array_infotpls]
            self.make_new_millertable()

            self.MillerComboBox.clear()
            self.MillerComboBox.addItem("", userData=[-1, ""])
            for k,lbl in enumerate(self.millerarraylabels):
              self.MillerComboBox.addItem(lbl, userData= [k, self.scenearraylabeltypes[k][1]] )

            self.MillerComboBox.setCurrentIndex(0) # select the first item which is no miller array
            self.comboviewwidth = 0
            for e in self.millerarraylabels:
              self.comboviewwidth = max(self.comboviewwidth, self.MillerComboBox.fontMetrics().width( e) )
            self.MillerComboBox.view().setMinimumWidth(self.comboviewwidth)

          if self.infodict.get("show_log_file_from_external_cmd"):
            if self.infodict.get("show_log_file_from_external_cmd") != -42:
              tabname, fname = self.infodict.get("show_log_file_from_external_cmd")
              mstr = ""
              with open(fname, 'r') as f:
                mstr += f.read() + '\\n'
              self.add_another_text_tab(fname, mstr, os.path.abspath(fname))
            self.waiting = False
            self.XtricorderBtn.setEnabled(True)
            self.XtriageBtn.setEnabled(True)
            self.AddInfoText("\n")
            self.AddAlertsText("\n")

          if self.infodict.get("ano_spg_tpls"):
            # needed for determining if expansion checkbox for P1 and friedel are enabled or disabled
            self.ano_spg_tpls = self.infodict.get("ano_spg_tpls",[])

          if self.infodict.get("current_labels", None) is not None:
            # needed in case we run xtriage on this miller array
            self.current_labels = self.infodict.get("current_labels")

          if self.infodict.get("colnames_select_lst", None) is not None:
            self.colnames_select_lst = self.infodict["colnames_select_lst"]
            self.select_millertable_column_dlg.make_new_selection_table()

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
                  item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                if col==1:
                  self.lowerbinvals.append(elm)
                if col==2:
                  item = QTableWidgetItem(str(elm))
                  item.setFlags(item.flags() ^ Qt.ItemIsEditable)
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

          if self.infodict.get("bin_opacity"):
            self.bin_opacities = self.infodict["bin_opacity"]
            if self.binstable.rowCount() > 0:
              self.update_table_opacities()

          if self.infodict.get("html_url") and self.html_url is None:
            self.html_url = self.infodict["html_url"]
            if self.UseOSBrowser==False:
              self.BrowserBox.setUrl(QUrl(self.html_url))
              # workaround for background colour bug in chromium
              # https://bugreports.qt.io/browse/QTBUG-41960
              self.BrowserBox.page().setBackgroundColor(QColor(127, 127, 127, 1) )
              pal = QPalette()
              self.BrowserBox.setAutoFillBackground(True)
              pal.setColor(self.BrowserBox.backgroundRole(), self.backgroundcolour)
              self.BrowserBox.setPalette(pal)
              self.send_message("Loading %s in QwebEngine" %self.html_url, msgtype="debug_info")

          if self.infodict.get("spacegroups"):
            spgs = self.infodict.get("spacegroups",[])
            self.spacegroups = { i : e for i,e in enumerate(spgs) }
            self.SpaceGroupComboBox.clear()
            self.SpaceGroupComboBox.addItems( list(self.spacegroups.values()) )

          currentinfostr = ""
          if self.infodict.get("info"):
            currentinfostr = self.infodict.get("info","")
            if "Destroying HKLViewFrame" in currentinfostr:
              self.canexit = True

          currentalertstr = ""
          if self.infodict.get("alert") or self.infodict.get("info"):
            currentalertstr = self.infodict.get("alert","") + self.infodict.get("info","")
            if self.closing: # Qt GUI has gone so print oending strings to the shell
              print(currentalertstr)

          if self.infodict.get("tabulate_miller_array"):
            currentinfostr = "Received table data\n"
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
            self.millerarraytable.resizeColumnsToContents()

            self.millerarraytableform.SortComboBox.clear()
            self.millerarraytableform.SortComboBox.addItems(["unsorted"] + labels )
            self.millerarraytableform.SortComboBox.view().setMinimumWidth(self.comboviewwidth)
            self.millerarraytableform.resize(tablewidth, self.millerarraytable.rowHeight(0)*15)
            self.millerarraytableform.show()
            self.millerarraytableform.activateWindow()

          self.unfeedback = True
          if self.infodict.get("all_vectors") is not None:
            self.rotvec = None
            self.all_vectors = self.infodict.get("all_vectors",[])

            self.clipplane_normal_vector_combo.clear()
            self.vectortable2.clearContents()
            self.vectortable2.setRowCount(len(self.all_vectors)+1)
            cw = 0
            for row, (opnr, label, order, cartvec, hklop, hkls, abcs, length) in enumerate(self.all_vectors):
              for col,elm in enumerate((label, hklop, hkls, abcs)):
                item = QTableWidgetItem(str(elm))
                if col == 0:
                  item.setFlags((Qt.ItemIsUserCheckable | Qt.ItemIsEnabled) ^ Qt.ItemIsEditable)
                  item.setCheckState(Qt.Unchecked)
                  item.setData(Qt.UserRole, opnr)
                item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                self.vectortable2.setItem(row, col, item)
              self.clipplane_normal_vector_combo.addItem(label, userData=length)
              cw = max(cw, self.clipplane_normal_vector_combo.fontMetrics().width( label) )
            self.clipplane_normal_vector_combo.view().setMinimumWidth(cw)

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
            self.currentfileName = self.infodict.get("file_name", "")
            self.setWindowFilenameTitles( self.currentfileName)
            for tabname in ["xtricorder.log",  "xtriage.log"]:
              self.removeNamedTab(tabname)

          if self.infodict.get("NewFileLoaded"):
            self.NewFileLoaded = self.infodict.get("NewFileLoaded",False)
            self.currentmillarray_idx = None
            self.vectortable2.clearContents()
            self.ShowAllVectorsBtn.setCheckState(Qt.Unchecked)
            self.functionTabWidget.setDisabled(True)
            self.AlignVectorGroupBox.setChecked( False)
            for tabname in ["xtricorder.log",  "xtriage.log"]:
              self.removeNamedTab(tabname)
            # display only those miller table columns that have been persisted if any
            stored_colnames_select_lst = []
            current_philstr = "selected_info {\n"
            for philname, caption, value in self.colnames_select_lst:
              is_selected = value
              if self.colnames_select_dict.get(philname, -42) != -42:
                is_selected = bool(self.colnames_select_dict[philname])
              stored_colnames_select_lst.append( (philname, caption, is_selected) )
              current_philstr += "  %s = %s\n" %(philname, is_selected)
            current_philstr += "}\n"
            self.send_message(current_philstr)
            self.colnames_select_lst = stored_colnames_select_lst
            self.select_millertable_column_dlg.make_new_selection_table()
    # Notify CCTBX that GUI has been initiated and it can now process messages.
    # This is critical as it releases a waiting semaphore in CCTBX
            self.send_message("", msgtype="initiated_gui")

          if self.infodict.get("NewHKLscenes"):
            self.NewHKLscenes = self.infodict.get("NewHKLscenes",False)

          if self.infodict.get("NewMillerArray"):
            self.NewMillerArray = self.infodict.get("NewMillerArray",False)

          if self.infodict.get("StatusBar") and self.Statusbartxtbox is not None:
            self.Statusbartxtbox.setText(self.infodict.get("StatusBar", "") )

          if self.infodict.get("include_tooltip_lst"):
            self.include_tooltip_lst = self.infodict.get("include_tooltip_lst", [])

          if self.infodict.get("clicked_HKL"):
            (h,k,l) = self.infodict.get("clicked_HKL", ( ) )
          if self.infodict.get("orig_hkl_ids"):
            if self.millerarraytablemodel is not None and self.millerarraytable.isVisible() :
              if self.millerarraytableform.SortComboBox.currentIndex() == 0:
                # can only match hkls in the unsorted table
                orig_hkl_ids = self.infodict.get("orig_hkl_ids", [])
                mode = QItemSelectionModel.Select | QItemSelectionModel.Rows
                self.millerarraytable.activateWindow()
                self.millerarraytable.setFocus()
                for ids in orig_hkl_ids:
                  self.millerarraytable.selectRow(ids)
              else:
                self.AddInfoText("Un-sort the table in order to highlight data corresponding " \
                  "to reflection that was clicked by the mouse.\n")

          if self.infodict.get("ColourChart") and self.infodict.get("ColourPowerScale"):
            self.ColourMapSelectDlg.selcolmap = self.infodict.get("ColourChart", "brg")
            self.ColourMapSelectDlg.setPowerScaleSliderVal( self.infodict.get("ColourPowerScale", 1.0))
            if self.infodict.get("ShowColourMapDialog"):
              self.ColourMapSelectDlg.show()
              self.ColourMapSelectDlg.activateWindow()

          if self.infodict.get("CurrentDatatype"):
            # sent by jsview_3d.HKLview_3d.DrawNGLJavaScript() once reflection have been rendered
            self.ColourMapSelectDlg.setDataType(self.infodict.get("CurrentDatatype", ""))
            self.actionSave_Current_Image.setDisabled(False)

          if self.infodict.get("bin_labels_type_idxs"):
            bin_labels_type_idxs = self.infodict.get("bin_labels_type_idxs",False)
            self.BinDataComboBox.clear()
            # fill combobox with labels of data that can be used for binning
            for label,labeltype,idx in bin_labels_type_idxs:
              self.BinDataComboBox.addItem(label, (labeltype, idx) )
            self.BinDataComboBox.view().setMinimumWidth(self.comboviewwidth)

          if self.infodict.get("binner_idx"):
            self.BinDataComboBox.setCurrentIndex(self.infodict.get("binner_idx", 0))

          if self.infodict.get("used_nth_power_scale_radii", None) is not None:
            self.unfeedback = True
            self.power_scale_spinBox.setValue( self.infodict.get("used_nth_power_scale_radii", 0.0))
            self.unfeedback = False

          if self.infodict.get("datatype_dict"):
            self.datatypedict = self.infodict.get("datatype_dict", {} )

          if self.infodict.get("enable_disable_preset_buttons"):
            newlist = eval(self.infodict.get("enable_disable_preset_buttons", "[]" ))
            if newlist != self.buttonsdeflist:
              self.buttonsdeflist = newlist

              for i in reversed(range(self.gridLayout_24.count())):
                # first delete any previous widgets from last time a file was loaded
                widgetToRemove = self.gridLayout_24.itemAt(i).widget()
                self.gridLayout_24.removeWidget(widgetToRemove)
                widgetToRemove.setParent(None)
              # programmatically create preset buttons on the self.gridLayout_24 from the QtDesigner
              for i,((btnname, label, _), datalabel, tooltip, moniker_veclabels) in enumerate(self.buttonsdeflist):
                moniker = ""
                veclabels = []
                if moniker_veclabels:
                  (moniker, veclabels) = moniker_veclabels
                self.__dict__[btnname] = QRadioButton(self.PresetButtonsFrame)
                pbutton = self.__getattribute__(btnname)
                pbutton.setObjectName(btnname)
                # since QRadioButton cannot wrap text the text for pbutton goes in
                # btnlabel below which is a QLabel where text can be wrapped if needed
                pbutton.setText("")
                pbutton.setToolTip(tooltip)
                pbutton.clicked.connect(self.onPresetbtn_click)
                sizePolicy1 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
                sizePolicy1.setHorizontalStretch(0)
                sizePolicy1.setVerticalStretch(0)
                sizePolicy1.setHeightForWidth(pbutton.sizePolicy().hasHeightForWidth())
                pbutton.setSizePolicy(sizePolicy1)
                self.gridLayout_24.addWidget(pbutton, i, 0, 1, 1)
                btnlabel = QLabel(self.widget)
                sizePolicy2 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
                sizePolicy2.setHorizontalStretch(1)
                sizePolicy2.setVerticalStretch(1)
                sizePolicy2.setHeightForWidth(btnlabel.sizePolicy().hasHeightForWidth())
                btnlabel.setSizePolicy(sizePolicy2)
                btnlabel.setWordWrap(True)
                btnlabel.setToolTip(tooltip)
                btnlabel.setText(label + " (using " + datalabel + ")")
                self.gridLayout_24.addWidget(btnlabel, i, 1, 1, 1)

                if moniker != "" and len(veclabels):
                  comboboxname = btnname + "_vectors"
                  self.__dict__[comboboxname] = QComboBox(self.PresetButtonsFrame)
                  vectorscombobox = self.__getattribute__(comboboxname)
                  vectorscombobox.addItems( veclabels )
                  vectorscombobox.activated.connect(self.onVectorsComboSelchange)
                  self.gridLayout_24.addWidget(vectorscombobox, i, 2, 1, 1)

                if self.ipresetbtn == i:
                  pbutton.setChecked(True)

          if self.infodict.get("spacegroup_info"):
            spacegroup_info = self.infodict.get("spacegroup_info",False)
            unitcell_info = self.infodict.get("unitcell_info",False)
            htmlstr = '''<html><head/><body><p><span style=" font-weight:600;">Space group: \t
            </span>%s<span style=" font-weight:600;"><br/>Unit cell(s): \t</span>%s</p></body></html> '''
            self.SpaceGrpUCellText.setText(htmlstr %(spacegroup_info, "<br/>".join(unitcell_info)) )
          self.fileisvalid = True

          if currentinfostr:
            self.AddInfoText(currentinfostr)

          if currentalertstr:
            self.AddAlertsText(currentalertstr)

          if (self.NewFileLoaded or self.NewMillerArray) or self.NewHKLscenes:
            self.NewMillerArray = False
            if self.millerarraytablemodel:
              self.millerarraytablemodel.clear()
              self.millerarraytablemodel = MillerArrayTableModel([[]], [], self)

            self.make_new_millertable()
            #self.UsePersistedQsettings
            self.NewFileLoaded = False
            self.actionSave_reflection_file.setDisabled(False)

          if self.NewHKLscenes:
            self.NewHKLscenes = False

      except Exception as e:
        errmsg = str(e)
        if "Resource temporarily unavailable" not in errmsg: # ignore errors from no connection to ZMQ socket
          print( errmsg  +  traceback.format_exc(limit=10) )


  def tabTextScrollDownShow(self, textbox):
    # If tab is not the current one shown scroling down to the bottom of the textbox is exceeded by one pagestep
    # So subtract a pagestep from the scroll value (and then some) to display the latest addition to the text output
    textbox.verticalScrollBar().setValue( textbox.verticalScrollBar().maximum() - (textbox.verticalScrollBar().pageStep()-2) )
    ctextbox = self.tabText.widget( self.tabText.currentIndex() ).children()[1]
    ctextbox.verticalScrollBar().setValue( ctextbox.verticalScrollBar().maximum()  )


  def onTabtextChanged(self, i):
    self.finddlg.hide() # if it is present from a recent string search


  def AddInfoText(self, currentinfostr):
    if self.isembedded:
      print(currentinfostr)
    else:
      self.infostr += currentinfostr
      # display no more than self.bufsize Kbytes of text
      self.infostr = self.infostr[-1000*self.textinfosize:]
      self.textInfo.setPlainText(self.infostr)
      self.tabTextScrollDownShow(self.textInfo)


  def AddAlertsText(self, currentalertstr):
    if self.isembedded:
      print(currentalertstr)
    else:
      self.alertstr += currentalertstr
      # display no more than self.bufsize Kbytes of text
      self.alertstr = self.alertstr[-1000*self.textinfosize:]
      self.textAlerts.setPlainText(self.alertstr)
      self.tabTextScrollDownShow(self.textAlerts)


  def removeNamedTab(self, tabsubstr):
    for n in range(self.tabText.count()):
      tabname = self.tabText.tabText(n)
      if tabsubstr in tabname:
        self.tabText.removeTab( self.tabText.indexOf(self.__dict__[tabname]) )
        self.__dict__[tabname].setParent(None)


  def add_another_text_tab(self, tabname, mstr, ttip):
    self.removeNamedTab(tabname)
    self.__dict__[tabname] = QWidget()
    gridLayout = QGridLayout(self.__dict__[tabname])
    gridLayout.setSpacing(4)
    gridLayout.setContentsMargins(3, 3, 3, 3)
    gridLayout.setContentsMargins(0, 0, 0, 0)
    newtabedit = MyQPlainTextEdit(self.__dict__[tabname])
    sp = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    sp.setHeightForWidth(newtabedit.sizePolicy().hasHeightForWidth())
    newtabedit.setSizePolicy(sp)
    newtabedit.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    newtabedit.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    newtabedit.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
    newtabedit.setLineWrapMode(QPlainTextEdit.NoWrap)
    newtabedit.setReadOnly(True)
    newtabedit.setPlainText(mstr)
    newtabedit.setFont(self.texttabfont)
    newtabedit.finddlg = self.finddlg
    gridLayout.addWidget(newtabedit, 0, 0, 1, 1)
    self.tabText.addTab(self.__dict__[tabname], tabname)
    idx = self.tabText.indexOf(self.__dict__[tabname])
    self.tabText.setCurrentIndex( idx )
    self.tabText.setTabToolTip(idx, ttip)


  def make_new_millertable(self):
    self.millertable.clearContents()
    if len(self.array_infotpls) == 0:
      return
    self.millertable.setRowCount(len(self.array_infotpls))
    labels = [ e.strip("|").strip() for e in self.array_infotpls[0][0] ]
    self.millertable.setColumnCount(len(labels))
    self.millertable.setHorizontalHeaderLabels(labels)
    self.millertable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    for row,(headerlst,infotpls,fmtstrtpls,fmtstr2tpls) in enumerate(self.array_infotpls):
      for col,elm in enumerate(infotpls):
        self.millertable.setItem(row, col, QTableWidgetItem(fmtstrtpls[col].format(elm)))


  def UpdateGUI(self):
    self.unfeedback = True
    self.ManualPowerScalecheckbox.setChecked( math.isnan( self.currentphilstringdict['hkls.nth_power_scale_radii'] )==False )
    self.power_scale_spinBox.setEnabled( self.ManualPowerScalecheckbox.isChecked() )
    self.radii_scale_spinBox.setValue( self.currentphilstringdict['hkls.scale'])
    self.expandP1checkbox.setChecked( self.currentphilstringdict['hkls.expand_to_p1'])
    self.expandAnomalouscheckbox.setChecked( self.currentphilstringdict['hkls.expand_anomalous'])
    self.ExpandReflsGroupBox.setChecked(self.expandP1checkbox.isChecked() or self.expandAnomalouscheckbox.isChecked())
    self.sysabsentcheckbox.setChecked( self.currentphilstringdict['hkls.show_systematic_absences'])
    self.ttipalpha_spinBox.setValue( self.currentphilstringdict['NGL.tooltip_alpha'])
    self.mousemoveslider.setValue( self.mousespeedscale*self.currentphilstringdict['NGL.mouse_sensitivity'])
    if self.currentphilstringdict['viewer.angle_around_vector'] is not None:
      self.rotvec, dgr = self.currentphilstringdict['viewer.angle_around_vector']
      self.rotavecangle_labeltxt.setText("Reflections rotated around Vector with Angle: %3.1fº" %dgr)

    self.ColourMapSelectDlg.selcolmap = self.currentphilstringdict["hkls.color_scheme"]
    self.ColourMapSelectDlg.setPowerScaleSliderVal( self.currentphilstringdict["hkls.color_powscale"] )

    self.Nbins_spinBox.setValue( self.currentphilstringdict['binning.nbins'])
    if self.currentphilstringdict['spacegroup_choice'] is not None:
      self.SpaceGroupComboBox.setCurrentIndex(  self.currentphilstringdict['spacegroup_choice'] )
    #self.clipParallelBtn.setChecked( self.currentphilstringdict['clip_plane.is_parallel'])
    self.missingcheckbox.setChecked( self.currentphilstringdict['hkls.show_missing'])
    self.onlymissingcheckbox.setEnabled( self.currentphilstringdict['hkls.show_missing'] )
    if self.currentphilstringdict['viewer.scene_id'] is not None:
      self.functionTabWidget.setEnabled(True)
    self.cameraPerspectCheckBox.setChecked( "perspective" in self.currentphilstringdict['NGL.camera_type'])
    if self.currentphilstringdict['clip_plane.clip_width']:
      self.clipwidth_spinBox.setValue( self.currentphilstringdict['clip_plane.clip_width'])
    self.hkldist_spinBox.setValue( self.currentphilstringdict['clip_plane.hkldist'])
    if self.currentphilstringdict['viewer.fixorientation'] == "vector":
      self.AlignVectorGroupBox.setChecked( True)
      self.RotateAroundframe.setEnabled(True)
    else:
      self.AlignVectorGroupBox.setChecked( False)
    self.onlymissingcheckbox.setChecked( self.currentphilstringdict['hkls.show_only_missing'])

    self.DrawRealUnitCellBox.setChecked(self.currentphilstringdict['draw_real_space_unit_cell'])
    self.unitcellslider.setValue( self.currentphilstringdict['real_space_unit_cell_scale_fraction'] * self.unitcellslider.maximum())
    self.DrawReciprocUnitCellBox.setChecked(self.currentphilstringdict['draw_reciprocal_unit_cell'])
    self.reciprocunitcellslider.setValue( self.currentphilstringdict['reciprocal_unit_cell_scale_fraction'] * self.reciprocunitcellslider.maximum())

    if self.currentphilstringdict['viewer.animate_rotation_around_vector'] is not None:
      self.rotvec,speed = self.currentphilstringdict['viewer.animate_rotation_around_vector']
      self.AnimaRotCheckBox.setChecked( speed > 0 )

    self.ClipPlaneChkGroupBox.setChecked(self.currentphilstringdict['clip_plane.clip_width'] != None)
    self.AutoClipWidthCheckBox.setChecked(self.currentphilstringdict['clip_plane.auto_clip_width'] == True)
    self.clipwidth_spinBox.setDisabled(self.AutoClipWidthCheckBox.isChecked() )

    lbl = self.currentphilstringdict['clip_plane.normal_vector']
    idx = -1
    if len(self.all_vectors) > 0:
      for i,(opnr,label,order,cartvec,hklop,hkls,abcs,length) in enumerate(self.all_vectors):
        if lbl == label:
          idx = i
          break
      opnr,label,order,cartvec,hklop,hkls,abcs,length = self.all_vectors[idx]
      if hkls == "" or not self.ClipPlaneChkGroupBox.isChecked():
        self.normal_realspace_vec_btn.setEnabled(False)
        self.normal_realspace_vec_label.setEnabled(False)
      else:
        self.normal_realspace_vec_btn.setEnabled(True)
        self.normal_realspace_vec_label.setEnabled(True)

    if self.currentphilstringdict['viewer.fixorientation'] is not None:
      self.parallel_current_orientation_btn.setChecked( "None" in self.currentphilstringdict['viewer.fixorientation'] \
         or self.currentphilstringdict['viewer.is_parallel'] )
      self.normal_vec_btn.setChecked( "vector" in self.currentphilstringdict['viewer.fixorientation'] and \
        not self.currentphilstringdict['viewer.is_parallel'] and \
        not self.currentphilstringdict['clip_plane.is_assoc_real_space_vector'])
      self.normal_realspace_vec_btn.setChecked( "vector" in self.currentphilstringdict['viewer.fixorientation'] and \
        not self.currentphilstringdict['viewer.is_parallel'] and \
        self.currentphilstringdict['clip_plane.is_assoc_real_space_vector'])

      self.clipplane_normal_vector_combo.setCurrentIndex(idx )
      if isinstance(self.clipplane_normal_vector_combo.currentData(), float) or isinstance(self.clipplane_normal_vector_combo.currentData(), int):
        self.clipplane_normal_vector_length.setText("{:.6g}".format(self.clipplane_normal_vector_combo.currentData()))

    for i in range(self.vectortable2.rowCount()):
      try: # checkboxes are not present first time gui is started
        self.vectortable2.item(i, 0).setCheckState(Qt.Unchecked)
      except Exception as e:
        pass

    if self.currentphilstringdict['viewer.show_vector'] is not None:
      ivecs = self.currentphilstringdict['viewer.show_vector']
      for ivec in ivecs:
        try: # checkboxes are not present first time gui is started
          [i,b] = eval(ivec)
          if i < self.vectortable2.rowCount():
            if b:
              self.vectortable2.item(i, 0).setCheckState(Qt.Checked)
            else:
              self.vectortable2.item(i, 0).setCheckState(Qt.Unchecked)
        except Exception as e:
          pass

    if self.currentphilstringdict['viewer.show_all_vectors'] is not None:
      show_all_vectors = self.currentphilstringdict['viewer.show_all_vectors']
      for rvrow in range(self.vectortable2.rowCount()):
        if self.vectortable2.item(rvrow, 0) is not None:
          if show_all_vectors == 1:
            self.vectortable2.item(rvrow, 0).setCheckState(Qt.Checked)
          if show_all_vectors == -1:
            self.vectortable2.item(rvrow, 0).setCheckState(Qt.Unchecked)

    if self.currentphilstringdict.get('viewer.user_vector', None) is not None:
      self.user_vectors = self.currentphilstringdict['viewer.user_vector']

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


  def onTextbufferSizeChanged(self,val):
    if self.unfeedback:
      return
    self.textinfosize = val


  def onFinalMouseSensitivity(self):
    self.mousespeed = self.mousemoveslider.value()/self.mousespeedscale
    self.send_message('NGL.mouse_sensitivity = "%2.3f"' %self.mousespeed)


  def onMouseSensitivity(self):
    self.mousespeed = self.mousemoveslider.value()/self.mousespeedscale
    self.mousesensitxtbox.setText("%2.1f" %(self.mousemoveslider.value()*10.0/self.mousemoveslider.maximum()) )


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


  def onFontsizeChanged(self):
    val = self.fontspinBox.value()
    font = self.app.font()
    font.setPointSize(val)
    self.fontsize = val
    self.app.setFont(font)
    self.settingsform.setFixedSize( self.settingsform.sizeHint() )
    self.aboutform.setFixedSize( self.aboutform.sizeHint() )
    self.BgrndColourDlg.setFixedSize( self.BgrndColourDlg.sizeHint() )
    self.ColourMapSelectDlg.setFixedSize( self.ColourMapSelectDlg.sizeHint() )
    self.select_millertable_column_dlg.resize()
    self.SpaceGrpUCellText.setFont(font)
    self.texttabfont = QFont("Courier New")
    self.texttabfont.setBold(True)
    self.texttabfont.setPointSize(val)
    for i in range(self.tabText.count()):
      self.tabText.widget(i).children()[1].setFont(self.texttabfont)


  def onBrowserFontsizeChanged(self, val):
    self.browserfontsize = val
    self.send_message("NGL.fontsize = %d" %val)


  def onVectorWidthChanged(self, val):
    self.vectorwidth = val
    self.send_message("NGL.vector_width = %d" %val)


  def onClearTextBuffer(self):
    self.textInfo.clear()
    self.textAlerts.clear()
    self.infostr = ""
    self.alertstr = ""


  def onDebugShowPhil(self):
    return self.send_message("", msgtype="debug_show_phil")


  def onResetFactoryDefault(self):
    ret = QMessageBox.warning(self.window, "Reset to factory defaults next time HKLviewer starts",
                              "Are you sure?",
                              buttons=QMessageBox.Yes|QMessageBox.No, defaultButton=QMessageBox.No)
    if ret == QMessageBox.Yes:
      self.RemoveQsettings()
      msg = "User settings for %s have been removed. Factory defaults will be used after restart." %self.Qtversion
      self.AddInfoText(msg)
      self.resetFactoryDefaultbtn.setEnabled(False)


  def onWrapTextBuffer(self):
    if self.wraptextbtn.isChecked():
      self.wraptextinfo = True
      self.textInfo.setLineWrapMode(QPlainTextEdit.WidgetWidth)
      self.textAlerts.setLineWrapMode(QPlainTextEdit.WidgetWidth)
    else:
      self.wraptextinfo = False
      self.textInfo.setLineWrapMode(QPlainTextEdit.NoWrap)
      self.textAlerts.setLineWrapMode(QPlainTextEdit.NoWrap)


  def onCameraPerspect(self,val):
    if self.cameraPerspectCheckBox.isChecked():
      self.send_message("NGL.camera_type = perspective")
    else:
      self.send_message("NGL.camera_type = orthographic")


  def ExpandRefls(self):
    if self.unfeedback:
      return
    if self.ExpandReflsGroupBox.isChecked():
      self.send_message('''
      hkls.expand_to_p1 = True
      hkls.expand_anomalous = True
                      ''' )
    else:
      self.send_message('''
      hkls.expand_to_p1 = False
      hkls.expand_anomalous = False
                      ''' )


  def ExpandToP1(self):
    if self.unfeedback:
      return
    if self.expandP1checkbox.isChecked():
      self.send_message('hkls.expand_to_p1 = True' )
    else:
      self.send_message('hkls.expand_to_p1 = False' )


  def ExpandAnomalous(self):
    if self.unfeedback:
      return
    if self.expandAnomalouscheckbox.isChecked():
      self.send_message('hkls.expand_anomalous = True' )
    else:
      self.send_message('hkls.expand_anomalous = False' )


  def showSysAbsent(self):
    if self.unfeedback:
      return
    if self.sysabsentcheckbox.isChecked():
      self.send_message('hkls.show_systematic_absences = True')
    else:
      self.send_message('hkls.show_systematic_absences = False')


  def showMissing(self):
    if self.unfeedback:
      return
    if self.missingcheckbox.isChecked():
      self.send_message('hkls.show_missing = True')
      self.onlymissingcheckbox.setEnabled(True)
    else:
      self.send_message("""hkls.show_missing = False
                             hkls.show_only_missing = False
                          """)
      self.onlymissingcheckbox.setEnabled(False)


  def showOnlyMissing(self):
    if self.unfeedback:
      return
    if self.onlymissingcheckbox.isChecked():
      self.send_message('hkls.show_only_missing = True')
    else:
      self.send_message('hkls.show_only_missing = False')


  def onBindataComboSelchange(self, i):
    if self.BinDataComboBox.currentText():
      binner_idx = self.BinDataComboBox.currentIndex()
      self.send_message('binning.binner_idx = %d' % binner_idx )
      bin_opacitieslst = []
      for j in range(self.nbins):
        bin_opacitieslst.append((1.0, j))
      self.bin_opacities = bin_opacitieslst
      self.OpaqueAllCheckbox.setCheckState(Qt.Checked)


  def update_table_opacities(self, allalpha=None):
    bin_opacitieslst = self.bin_opacities
    self.binstable_isready = False
    for binopacity in bin_opacitieslst:
      if not allalpha:
        alpha = binopacity[0]
      else:
        alpha = allalpha
      bin = binopacity[1]
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
    bin_opacitieslst = self.bin_opacities
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
    self.binTableCheckState = item.checkState()
    self.bintableAlpha = float(item.text())


  def onBinsTableItemChanged(self, item):
    row = item.row()
    col = item.column()
    try:
      if not self.bin_opacities:
        return
      bin_opacitieslst = self.bin_opacities
      alpha = max(0.0, min(1.0, float(item.text()) ) ) # between 0 and 1 only
      try:
        (oldalpha, row) = bin_opacitieslst[row]
        row = int(row)
        if oldalpha == float(item.text()):
          if item.checkState()==Qt.Unchecked:
            alpha = 0.0
          else:
            alpha = 1.0
      except Exception as e:
        pass
      if col==3 and self.binstable_isready: # changing opacity
        bin_opacitieslst[row] = (alpha, row)
        self.bin_opacities = bin_opacitieslst
        self.SetAllOpaqueCheckboxes()
        philstr = ""
        for opa,bin in self.bin_opacities:
          philstr += 'binning.bin_opacity = %s %s\n' %(opa, bin)
        self.send_message(philstr)
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
        newval = float(item.text())
        self.binstable.item(row,col).setText(str(newval))
        self.lowerbinvals[row] = newval
        upperbinvalsnonan = [e for e in self.upperbinvals if not math.isnan(e)]
        allbinvals = self.lowerbinvals + [ upperbinvalsnonan[-1] ]
        self.send_message("binning.scene_bin_thresholds = %s" %" ".join([ str(e) for e in allbinvals]))
    except Exception as e:
      print( str(e)  +  traceback.format_exc(limit=10) )


  def onOpaqueAll(self):
    self.binstableitemchanges = True
    bin_opacitieslst = self.bin_opacities
    nbins = len(bin_opacitieslst)
    bin_opacitieslst = []
    self.binstable_isready = False
    if self.OpaqueAllCheckbox.isChecked():
      for i in range(nbins):
        bin_opacitieslst.append((1.0, i))  #  ("1.0, %d" %i)
    else:
      for i in range(nbins):
        bin_opacitieslst.append((0.0, i))  #   ("0.0, %d" %i)
    self.bin_opacities = bin_opacitieslst
    philstr = ""
    for opa,bin in self.bin_opacities:
      philstr += 'binning.bin_opacity = %s %s\n' %(opa, bin)
    self.send_message(philstr)
    self.binstableitemchanges = False
    self.binstable_isready = True


  def onNbinsChanged(self, val):
    if self.unfeedback:
      return
    self.nbins = val
    self.send_message("binning.nbins = %d" %self.nbins)


  def onRadiiScaleEditFinished(self):
    if self.unfeedback:
      return
    self.send_message("hkls.scale = %f" %self.radii_scale_spinBox.value() )


  def onPowerScaleEditFinished(self):
    if self.unfeedback:
      return
    self.send_message("hkls.nth_power_scale_radii = %f" %self.power_scale_spinBox.value() )


  def onManualPowerScale(self, val=None):
    if self.unfeedback:
      return
    if self.ManualPowerScalecheckbox.isChecked():
      self.send_message('hkls.nth_power_scale_radii = %f' %self.power_scale_spinBox.value())
    else:
      self.send_message('hkls.nth_power_scale_radii = %s' %float("nan"))


  def onShowAllVectors(self):
    if self.unfeedback:
      return
    if self.ShowAllVectorsBtn.checkState()==Qt.Checked:
      self.send_message("viewer.show_all_vectors = 1")
    if self.ShowAllVectorsBtn.checkState()==Qt.Unchecked:
      self.send_message("viewer.show_all_vectors = -1")
    self.unfeedback = False


  def onVectorTableItemChanged(self, item):
    if self.unfeedback:
      return
    row = item.row()
    col = item.column()
    try:
      rc = len(self.all_vectors)
      label = None
      opnr = -1
      if self.vectortable2.item(row, 0) is not None:
        label = self.vectortable2.item(row, 0).text()
      if row < rc:
        if label is None:
          return
        if col==0:
          philstr = ""
          if self.rotvec is not None: # reset any imposed angle to 0 whenever checking or unchecking a vector
              philstr += "viewer.angle_around_vector = '[%d, 0]'\n" %self.rotvec
              self.rotavecangle_slider.setValue(0)
          self.rotvec = None
          sum = 0
          ivec= []
          for rvrow in range(self.vectortable2.rowCount()):
            if self.vectortable2.item(rvrow, 0) is not None:
              opnr = self.vectortable2.item(rvrow, 0).data(Qt.UserRole)
              if self.vectortable2.item(rvrow, 0).checkState()==Qt.Checked:
                self.rotvec = rvrow
                sum +=1
                ivec = [opnr, True]
              else:
                ivec = [opnr, False]
              philstr += "viewer.show_vector = " + str(ivec) + "\n"
          if sum > 1 or sum == 0: # can only use one vector to rotate around. so if more are selected then deselect them altogether
            philstr += "viewer.animate_rotation_around_vector = '[%d, %f]'\n" %(0, -1.0)
            philstr += 'viewer.fixorientation = *None\n'
            self.AnimaRotCheckBox.setCheckState(Qt.Unchecked)
            self.rotvec = None
          if self.rotvec is not None:
            self.RotateAroundframe.setEnabled(True)
            # notify cctbx which is the curently selected vector
            philstr += "viewer.angle_around_vector = '[%d, 0.0]'\n" %self.rotvec
          else:
            self.RotateAroundframe.setDisabled(True)
          if sum >= rc:
            self.ShowAllVectorsBtn.setCheckState(Qt.Checked)
          if sum == 0:
            self.ShowAllVectorsBtn.setCheckState(Qt.Unchecked)
          if sum >0.0 and sum < rc:
            self.ShowAllVectorsBtn.setCheckState(Qt.PartiallyChecked)
          self.send_message(philstr)
      if row==rc and label !="" and label != "new vector": # a user defined vector
        if col==1:
          hklop = self.vectortable2.item(row, 1).text()
          self.send_message("""
viewer.user_vector {
  hkl_op = '%s'
  label = %s
}""" %(hklop, label))
        if col==2:
          hklvec = self.vectortable2.item(row, 2).text()
          self.send_message("""
viewer.user_vector {
  hkl = '%s'
  label = %s
}""" %(hklvec, label))
        if col==3:
          abcvec = self.vectortable2.item(row, 3).text()
          self.send_message("""
viewer.user_vector {
  abc = '%s'
  label = %s
}""" %(abcvec, label))
    except Exception as e:
      print( str(e)  +  traceback.format_exc(limit=10) )


  def createExpansionBox(self):
    self.SpaceGroupComboBox.activated.connect(self.onSpacegroupSelchange)
    self.expandP1checkbox.clicked.connect(self.ExpandToP1)
    self.expandAnomalouscheckbox.clicked.connect(self.ExpandAnomalous)
    self.ExpandReflsGroupBox.clicked.connect(self.ExpandRefls)
    self.commitSubgroupExpansionBtn.clicked.connect(self.onCommitSubgroupExpansion)
    self.commitSubgroupExpansionBtn.setEnabled(False)
    self.sysabsentcheckbox.clicked.connect(self.showSysAbsent)
    self.missingcheckbox.clicked.connect(self.showMissing)
    self.onlymissingcheckbox.clicked.connect(self.showOnlyMissing)


  def CreateSliceTabs(self):
    vprec = 2
    self.hkldistval = 0.0
    self.hkldist_spinBox.setValue(self.hkldistval)
    self.hkldist_spinBox.setDecimals(1)
    self.hkldist_spinBox.setSingleStep(1)
    self.hkldist_spinBox.setRange(-1000.0, 1000.0)
    self.hkldist_spinBox.editingFinished.connect(self.onHKLdistEditFinished)
    self.hkldist_spinBox.onStepBy = self.onHKLdistEditFinished

    self.clipwidth_spinBox.setValue(0.35 )
    self.clipwidth_spinBox.setDecimals(3)
    self.clipwidth_spinBox.setSingleStep(0.05)
    self.clipwidth_spinBox.setRange(0.0, 100.0)
    self.clipwidth_spinBox.editingFinished.connect(self.onClipwidthEditFinished)
    self.clipwidth_spinBox.onStepBy = self.onClipwidthEditFinished

    self.yHKLrotBtn.clicked.connect(self.onYangleHKLrotate)
    self.xHKLrotBtn.clicked.connect(self.onXangleHKLrotate)
    self.zHKLrotBtn.clicked.connect(self.onZangleHKLrotate)
    self.yHKLbackrotBtn.clicked.connect(self.onYangleHKLrotateback)
    self.xHKLbackrotBtn.clicked.connect(self.onXangleHKLrotateback)
    self.zHKLbackrotBtn.clicked.connect(self.onZangleHKLrotateback)
    self.clipplane_normal_vector_length.editingFinished.connect(self.onClipwidthNormalVecLengthEditFinished)
    self.clipplane_normal_vector_combo.activated.connect(self.onClipPlaneNormalVecSelchange)

    self.parallel_current_orientation_btn.clicked.connect(self.onParallel_current_orientation_btn_click)
    self.normal_realspace_vec_btn.clicked.connect(self.onNormal_realspace_vec_btn_click)
    self.normal_vec_btn.clicked.connect(self.onNormal_vec_btn_click)


  def onClipwidthNormalVecLengthEditFinished(self):
    if self.unfeedback:
      return
    try:
      val = eval(self.clipplane_normal_vector_length.text())
      philstr = "clip_plane.normal_vector_length_scale = %s"  %val
      self.send_message(philstr)
    except Exception as e:
      print( str(e) )


  def onClipPlaneNormalVecSelchange(self):
    if self.unfeedback:
      return
    self.clipplane_normal_vector_length.setText("{:.6g}".format(self.clipplane_normal_vector_combo.currentData()))
    philstr = """viewer.fixorientation = *vector
clip_plane.clip_width = %f
clip_plane.normal_vector = "%s"
clip_plane.normal_vector_length_scale = -1
""" %(self.clipwidth_spinBox.value(), self.clipplane_normal_vector_combo.currentText())
    self.send_message(philstr)


  def onParallel_current_orientation_btn_click(self):
    if self.unfeedback:
      return
    self.clipplane_normal_vector_combo.setEnabled(False)
    self.RotateGroupBox.setEnabled(True)
    philstr = """viewer.fixorientation = *None
clip_plane.clip_width = %f
clip_plane.normal_vector = ""
""" %self.clipwidth_spinBox.value()
    self.send_message(philstr)


  def onNormal_vec_btn_click(self):
    if self.unfeedback:
      return
    self.clipplane_normal_vector_combo.setEnabled(True)
    self.RotateGroupBox.setEnabled(False)
    philstr = """viewer.fixorientation = *vector
viewer.is_parallel = False
clip_plane.clip_width = %f
clip_plane.normal_vector = "%s"
clip_plane.is_assoc_real_space_vector = False
clip_plane.normal_vector_length_scale = -1
""" %( self.clipwidth_spinBox.value(), self.clipplane_normal_vector_combo.currentText())
    self.send_message(philstr)


  def onNormal_realspace_vec_btn_click(self):
    if self.unfeedback:
      return
    self.clipplane_normal_vector_combo.setEnabled(True)
    self.RotateGroupBox.setEnabled(False)
    philstr = """viewer.fixorientation = *vector
viewer.is_parallel = False
clip_plane.clip_width = %f
clip_plane.normal_vector = "%s"
clip_plane.is_assoc_real_space_vector = True
clip_plane.normal_vector_length_scale = -1
""" %( self.clipwidth_spinBox.value(), self.clipplane_normal_vector_combo.currentText())
    self.send_message(philstr)


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
    val = "*None"
    if self.AlignVectorGroupBox.isChecked():
      val = "*vector"
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
      if self.normal_realspace_vec_btn.isChecked():
        self.clipplane_normal_vector_combo.setEnabled(True)
        self.RotateGroupBox.setEnabled(False)
        philstr = """hkls {
  slice_mode = False
}
viewer.is_parallel = False
viewer.fixorientation = *vector
clip_plane.clip_width = %f
clip_plane.normal_vector = "%s"
clip_plane.is_assoc_real_space_vector = True
clip_plane.normal_vector_length_scale = -1
""" %( self.clipwidth_spinBox.value(), self.clipplane_normal_vector_combo.currentText())
      elif self.normal_vec_btn.isChecked():
        self.clipplane_normal_vector_combo.setEnabled(True)
        self.RotateGroupBox.setEnabled(False)
        philstr = """hkls {
  slice_mode = False
}
viewer.is_parallel = False
viewer.fixorientation = *vector
clip_plane.clip_width = %f
clip_plane.normal_vector = "%s"
clip_plane.is_assoc_real_space_vector = False
clip_plane.normal_vector_length_scale = -1
""" %( self.clipwidth_spinBox.value(), self.clipplane_normal_vector_combo.currentText())
      else:
        self.clipplane_normal_vector_combo.setEnabled(False)
        self.RotateGroupBox.setEnabled(True)
        philstr = """hkls {
  slice_mode = False
}
viewer.is_parallel = True
viewer.fixorientation = *None
clip_plane.clip_width = %f
clip_plane.normal_vector = ""
          """ %self.clipwidth_spinBox.value()
    else:
      philstr = """hkls {
  slice_mode = False
}
viewer.fixorientation = *None
clip_plane {
  normal_vector = ""
  clip_width = None
}
       """
    self.send_message(philstr)


  def onAutoClipWidthChkBox(self):
    if self.unfeedback:
      return
    if self.AutoClipWidthCheckBox.isChecked():
      philstr = "clip_plane.auto_clip_width = True"
    else:
      philstr = "clip_plane.auto_clip_width = False"
    self.send_message(philstr)


  def onRotaVecAngleChanged(self, val):
    if self.unfeedback or self.rotvec is None:
      return
    self.send_message("viewer.angle_around_vector = '[%d, %f]'" %(self.rotvec, val*0.5))


  def onFinalRotaVecAngle(self):
    if self.unfeedback or self.rotvec is None:
      return
    val = self.rotavecangle_slider.value()*0.5
    self.send_message("viewer.angle_around_vector = '[%d, %f]'" %(self.rotvec, val))


  def onAnimateRotation(self):
    if not self.unfeedback:
      if self.AnimaRotCheckBox.isChecked() == True:
        self.AnimateSpeedSlider.setEnabled(True)
        self.rotavecangle_slider.setDisabled(True)
        speed = self.AnimateSpeedSlider.value()
        for rvrow in range(self.vectortable2.rowCount()):
          if self.vectortable2.item(rvrow, 0) is not None:
            if self.vectortable2.item(rvrow, 0).checkState()==Qt.Checked:
              self.rotvec = rvrow
              break
        self.send_message("viewer.animate_rotation_around_vector = '[%d, %f]'" %(self.rotvec, speed))
      else:
        self.rotavecangle_slider.setEnabled(True)
        self.AnimateSpeedSlider.setDisabled(True)
        self.send_message("viewer.animate_rotation_around_vector = '[%d, %f]'" %(self.rotvec, -1.0))


  def onClipwidthEditFinished(self):
    if not self.unfeedback:
      self.send_message("clip_plane.clip_width = %f" %self.clipwidth_spinBox.value())


  def onHKLdistEditFinished(self):
    if not self.unfeedback:
      self.hkldistval = self.hkldist_spinBox.value()
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


  def onMillerTableCellPressed(self, row, col):
    #print( "in millertable CellPressed " + self.millertable.currentItem().text() )
    if self.millertable.mousebutton == Qt.RightButton:
      self.MillerTableContextMenuHandler(QCursor.pos(), row)
    if self.millertable.mousebutton == QEvent.MouseButtonDblClick:
      # quickly display data with a double click
      for scenelabel,labeltype,arrayid,hassigmas,sceneid in self.scenearraylabeltypes:
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
    for i,(scenelabel,labeltype,arrayid,hassigmas,sceneid) in enumerate(self.scenearraylabeltypes): # loop over scenes
      scenelabelstr = scenelabel
      if self.millerarraylabels[row] == scenelabelstr or self.millerarraylabels[row] + " + " in scenelabelstr:
        if hassigmas: # sigmas are present for this array
          myqa = QAction("Display data of %s" %scenelabelstr, self.window, triggered=self.testaction)
          myqa.setData(("display_data",[sceneid, row]))
          self.millertablemenu.addAction(myqa)
          myqa = QAction("Display sigmas of %s" %scenelabelstr, self.window, triggered=self.testaction)
          myqa.setData(("display_data",[sceneid + 1000, row])) # want to show the sigmas rather than the data if we add 1000
          self.millertablemenu.addAction(myqa)
        else: # no sigmas present for this array
          myqa = QAction("Display %s" %scenelabelstr, self.window, triggered=self.testaction)
          myqa.setData(("display_data",[sceneid, row]))
          self.millertablemenu.addAction(myqa)
    myqa = QAction("Make a new dataset from this dataset and another dataset...",
                   self.window, triggered=self.testaction)
    myqa.setData(("make_newdata", row ))
    self.millertablemenu.addAction(myqa)
    if len(self.millertable.selectedrows) ==1:
      myqa = QAction("Display tooltips for %s" %self.millerarraylabels[row], self.window, triggered=self.testaction)
      myqa.setCheckable(True)
      myqa.setChecked(self.include_tooltip_lst[row] )
      myqa.setData(("tooltip_data", row))
      self.millertablemenu.addAction(myqa)

    if len(self.millertable.selectedrows) > 1:
      arraystr = ""
      labels = []
      sum = 0
      for i,row in enumerate(self.millertable.selectedrows):
        labels.append( self.millerarraylabels[row] )
        sum += int(self.include_tooltip_lst[row])
      myqa = QAction("Display tooltips for %s" %  " and ".join(labels),
                     self.window, triggered=self.testaction)
      myqa.setCheckable(True)
      if len(self.millertable.selectedrows) == sum:
        myqa.setChecked(True )
      if sum == 0:
        myqa.setChecked(False )
      if len(self.millertable.selectedrows) > sum and sum > 0:
        myqa = QAction(QIcon( os.path.join(os.path.dirname(os.path.abspath(__file__)), "partiallychecked.png")),
                       "Display tooltips for %s" %  " and ".join(labels),
                       self.window, triggered=self.testaction)
        myqa.setCheckable(True)

      myqa.setData(("tooltip_data", self.millertable.selectedrows))
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
      strval, val = data
      if strval=="display_data":
        [idx, row] = val
        self.DisplayData(idx,row)
      if strval=="make_newdata":
        self.operate_arrayidx1 = val
        self.operate_arraytype1 = self.scenearraylabeltypes[val][1]
        self.operate_arrayidx2 = -1 # i.e. no second miller array selected yet in MillerComboBox
        self.operate_arraytype2 = ""
        self.operationlabeltxt.setText("Define a new cctbx.miller.array object, \"newarray\", "
          + "by entering a python expression for \"newarray\" or by assigning \"newarray._data\" "
          + "(and optionally \"newarray._sigmas\") to a function of the cctbx.miller.array object, "
          + "\"array1\", representing the " + self.millerarraylabels[val] + " dataset. Optionally "
          + "also include the cctbx.miller.array object, \"array2\" representing a dataset "
          + "selected from the dropdown list below."
          )
        self.makenewdataform.show()
      if strval=="tabulate_data":
        self.send_message('tabulate_miller_array_ids = "%s"' %str(val))
      if strval=="tooltip_data":
        if isinstance(val, list): # if we highlighted more rows then toggle tooltips on or off by
          sum = 0   # taking the opposite of the rounded integer average of the tooltip visibility
          for row in self.millertable.selectedrows:
            sum += int(self.include_tooltip_lst[row])
          bval = bool(round(float(sum)/len(self.millertable.selectedrows)))
          for row in self.millertable.selectedrows:
            self.include_tooltip_lst[row] = not bval
            self.send_message('tooltip_data = "%s"' %str([row, self.include_tooltip_lst[row]]))
        else: # only rightclicked one row so toggle its tooltip visibility
          self.include_tooltip_lst[val] = not self.include_tooltip_lst[val]
          self.send_message('tooltip_data = "%s"' %str([val, self.include_tooltip_lst[val]]))


  def DisplayData(self, idx, row):
    # want to show the sigmas rather than the data if we have added 1000 to idx
    self.currentmillarray_idx = row
    arrayinfo = self.array_infotpls[self.currentmillarray_idx]
    if (idx - 1000) >= 0:
      idx = idx - 1000
      philstr = """
      hkls.sigma_color_radius = True
      viewer.scene_id = %d
      """ %idx
    else:
      philstr = """
      hkls.sigma_color_radius = False
      viewer.scene_id = %d
      """ %idx
    self.send_message(philstr)
    if self.fileisvalid:
      self.functionTabWidget.setEnabled(True)
      self.expandAnomalouscheckbox.setEnabled(True)
      self.expandP1checkbox.setEnabled(True)
      # don't allow anomalous expansion for data that's already anomalous
      isanomalous, spacegroup = self.ano_spg_tpls[self.currentmillarray_idx]
      label = arrayinfo[0][0]
      if isanomalous:
        self.expandAnomalouscheckbox.setDisabled(True)
      if spacegroup=='P 1 (No. 1)':
        self.expandP1checkbox.setDisabled(True)
    else:
      self.functionTabWidget.setDisabled(True)
    self.SpaceGroupComboBox.clear()
    self.SpaceGroupComboBox.addItems( list(self.spacegroups.values() ))


  def onMakeNewData(self):
    lbltype2 = ["",""]
    if self.operate_arrayidx2 >= 0: # a dataset was selected in the millercombobox. Get label and type
      lbltype2 = [self.millerarraylabels[self.operate_arrayidx2],
                  self.operate_arraytype2]
    mtpl = (self.operationtxtbox.toPlainText(),
              self.newlabeltxtbox.text() ,
              [self.millerarraylabels[self.operate_arrayidx1],
               self.operate_arraytype1],
               lbltype2
            )
    self.send_message('miller_array_operation = "%s"' %str(mtpl) )
    self.makenewdataform.accept()


  def onMillerComboSelchange(self, i):
    self.operate_arrayidx2, self.operate_arraytype2 = self.MillerComboBox.itemData(i)
    # -1 if first item. Otherwise itemdata is index of the list of miller arrays


  def testaction(self):
    pass


  def onCommitSubgroupExpansion(self):
    ret = QMessageBox.warning(self.window, "Expand datasets to selected subgroup",
                              "This will irreversibly expand all datasets to the selected subgroup\nAre you sure?",
                              buttons=QMessageBox.Yes|QMessageBox.No, defaultButton=QMessageBox.No)
    if ret == QMessageBox.Yes:
      self.send_message('commit_subgroup_datasets = True')
      self.commitSubgroupExpansionBtn.setEnabled(False)


  def onAddDataset(self):
    label, ok = QInputDialog().getText(self.window, "Enter a unique label for the new dataset",
                                         "Create new dataset of visible reflections with this label:")
    if ok and label:
       self.send_message('visible_dataset_label = "%s"' %label)


  def createFileInfoBox(self):
    labels = ["Column label", "Type", "λ(Å)", "# HKLs", "Span of HKLs",
       "Min Max data", "Min Max sigmas", "d_min, d_max (Å)", "Anomalous", "Symmetry unique"]
    self.millertable.setColumnCount(len(labels))
    self.millertable.setHorizontalHeader( MyhorizontalHeader(self.window) )
    self.millertable.setHorizontalHeaderLabels(labels)
    self.millertable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    # don't allow editing this table
    self.millertable.setEditTriggers(QTableWidget.NoEditTriggers)
    self.millertable.cellPressed.connect(self.onMillerTableCellPressed)
    self.millertable.cellDoubleClicked.connect(self.onMillerTableCellPressed)
    self.millertable.itemSelectionChanged.connect(self.onMillerTableitemSelectionChanged)


  def createRadiiScaleGroupBox(self):
    self.ManualPowerScalecheckbox.clicked.connect(self.onManualPowerScale)
    self.power_scale_spinBox.editingFinished.connect(self.onPowerScaleEditFinished)
    self.power_scale_spinBox.onStepBy = self.onPowerScaleEditFinished
    self.radii_scale_spinBox.editingFinished.connect(self.onRadiiScaleEditFinished)
    self.radii_scale_spinBox.onStepBy = self.onRadiiScaleEditFinished


  def createBinsBox(self):
    self.binstable_isready = False
    labels = ["# HKLs", "lower bin value", "upper bin value", "opacity"]
    self.binstable.setHorizontalHeaderLabels(labels)
    self.binstable.horizontalHeader().setDefaultAlignment(Qt.AlignLeft)
    self.Nbins_spinBox.valueChanged.connect(self.onNbinsChanged)
    self.OpaqueAllCheckbox.clicked.connect(self.onOpaqueAll)
    self.addDatasetBtn.clicked.connect(self.onAddDataset)
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
    self.AutoClipWidthCheckBox.clicked.connect(self.onAutoClipWidthChkBox)


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
        self.send_message("draw_reciprocal_unit_cell = True")
      else:
        self.send_message("draw_reciprocal_unit_cell = False")


  def onDrawUnitCellBoxClick(self):
    if not self.unfeedback:
      if self.DrawRealUnitCellBox.isChecked():
        self.send_message("draw_real_space_unit_cell = True")
      else:
        self.send_message("draw_real_space_unit_cell = False")


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


  def onSpacegroupSelchange(self,i):
    self.send_message("spacegroup_choice = %d" %i)
    if i==0:
      self.commitSubgroupExpansionBtn.setEnabled(False)
    else:
      self.commitSubgroupExpansionBtn.setEnabled(True)


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
                                      stdin=subprocess.PIPE)
    # Wait for connection from the zmq socket in CCTBX by testing if we can send an empty string
    print("Establishing ZMQ socket to CCTBX process", end='')
    t=0.0; dt = 0.5; timeout = 60 # more than enough time for connecting
    err = zmq.EAGAIN
    while err == zmq.EAGAIN:
      try:
        err = 0
        self.socket.send(bytes("","utf-8"), zmq.NOBLOCK)
      except Exception as e:
        err = e.errno # if EAGAIN then message cannot be sent at the moment.
      time.sleep(dt)
      t += dt
      print(".", end='', flush=True)
      if  t > timeout:
        raise Exception("\nHKLviewer GUI failed making ZMQ socket connection to CCTBX process.")
        break
    print("\nDone.")


  def send_message(self, cmdstr, msgtype="philstr"):
    try:
      msg = str([msgtype, cmdstr])
      if sys.version_info.major==3:
        self.socket.send(bytes(msg,"utf-8"), zmq.NOBLOCK)
      else:
        self.socket.send(bytes(msg), zmq.NOBLOCK)
      if msgtype == "philstr":
        self.ipresetbtn = -1 # so we don't tick any of the preset radio buttons after remaking them
      return True
    except Exception as e:
      print( str(e) + "\nFailed sending message to the CCTBX\n" + traceback.format_exc(limit=10))
      return False


  def setDatatypedict(self, datatypedict):
    self.datatypedict = datatypedict
    # send persisted colour schemes and raddi mappings to jsview_3d.py
    return self.send_message(str(self.datatypedict), msgtype="datatypedict")


  def PersistQsettings(self, write_factory_default_settings = False):
    Qtversion = self.Qtversion
    if write_factory_default_settings:
      # For developers only:
      # When supplying the new_factory_defaults argument on the command line all the settings
      # usually stored in the users Qsettings database will instead be stored in the file
      # crys3d\HKLviewer\HKLviewerDefaults.ini which can then be committed to git repo if desired and
      # work as new factory default settings. Adjust colours and sizes for datasets as well
      # as other settings as desired to write a HKLviewerDefaults.ini with those new settings
      # when exiting HKLviewer. The file should then be committed to the cctbx_project repo.
      self.settings = QSettings(self.factorydefaultfname,  QSettings.IniFormat)
      self.AddInfoText("Writing factory defaults to %s\n" %self.factorydefaultfname)
      self.AddAlertsText("Writing factory defaults to %s\n" %self.factorydefaultfname)
      Qtversion = "Qt"
    if not write_factory_default_settings:  # don't store system specific value as a default
      self.settings.setValue("PythonPath", self.cctbxpython )
    self.settings.beginGroup("MillerTableColumnHeader")
    ###### Code below should be done elsewhere such as in conda installation scripts for CCTBX
    # This presumes Qt is present in the CCTBX build and will then provide a path for GUI dispatchers
    cctbxbindir = os.path.split(self.cctbxpython)[0]
    cctbxversionsettings = QSettings("CCTBX", self.cctbxversion ) # version number sent from cctbx process
    cctbxversionsettings.setValue("DispatcherPath", cctbxbindir )
    #####################
    for philname, dummy, value in self.colnames_select_lst:
      self.settings.setValue(philname, int(value) )
    if len(self.colnames_select_lst) == 0:
      # No hkl file was opened so just save whatever MillerTableColumnHeader was already on disc
      # ReadPersistedQsettings() stored this in self.colnames_select_dict initially
      for philname in list(self.colnames_select_dict.keys()):
        self.settings.setValue(philname, int(self.colnames_select_dict[philname]) )
    self.settings.endGroup() # MillerTableColumnHeader

    self.settings.beginGroup(Qtversion )
    if not write_factory_default_settings: # don't store system specific value as a default
      self.settings.setValue("QWebEngineViewFlags", self.QWebEngineViewFlags)
    self.settings.setValue("QSplitter_number", self.nsplitters )
    self.settings.setValue("QTabWidget_number", self.ntabs )
    self.settings.setValue("FontSize", self.fontsize )
    self.settings.setValue("BackgroundColour", str(self.backgroundcolour.getRgb()) )
    self.settings.setValue("WordWrapTextInfo", int(self.wraptextinfo ))
    self.settings.setValue("MouseSpeed", self.mousespeed )
    self.settings.setValue("TextBufferSize", self.textinfosize )
    self.settings.setValue("BrowserFontSize", self.browserfontsize )
    self.settings.setValue("VectorWidth", self.vectorwidth )
    self.settings.setValue("ttip_click_invoke", self.ttip_click_invoke)
    self.settings.setValue("geometry", self.window.saveGeometry())
    self.settings.setValue("windowstate", self.window.saveState())
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
    if write_factory_default_settings: # reset path for when closing app will save settings again
      self.settings = QSettings("CCTBX", "HKLviewer" )


  def ReadPersistedQsettings(self):
    # Read the user's persisted settings from disc
    self.settings.beginGroup(self.Qtversion)
    use_factory_default_settings = False
    # First see if there are any. If not then use factory defaults stored in .ini file
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    if len(self.settings.allKeys()) == 0:
       # no settings for this Qt version. Use defaults then
      use_factory_default_settings = True
      # Numbers of splitters and tabs in the GUI are a very crude indication of
      # GUI complexity. If numbers differs from what is stored in the settings on disk the
      # settings are likely from a newer or older GUI version and should be ignored to prevent
      # messing up GUI layout. Use the defaults instead
    if self.nsplitters !=  int(self.settings.value("QSplitter_number", 0)):
      use_factory_default_settings = True
    if self.ntabs != int(self.settings.value("QTabWidget_number", 0)):
      use_factory_default_settings = True
    self.settings.endGroup()

    Qtversion = self.Qtversion
    if use_factory_default_settings:
      print("Reading factory defaults from " + self.factorydefaultfname)
      self.settings = QSettings(self.factorydefaultfname, QSettings.IniFormat)
      Qtversion = "Qt"

    # Locate cctbx.python. If not in the Qsettings then try if in the executable path environment
    cctbxpython_from_settings = self.settings.value("PythonPath", "")
    cctbxpython_from_env = ""

    wherecmd = "which"
    if sys.platform == 'win32':
      wherecmd = "where"
    proc = subprocess.Popen([wherecmd, "cctbx.python"],
                            universal_newlines=True, # avoid them annoying byte strings
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE).communicate()
    if proc[0] != "":
      cctbxpython_from_env = proc[0].strip()

    self.cctbxpython = cctbxpython_from_env
    if not os.path.isfile(self.cctbxpython):
      self.cctbxpython = cctbxpython_from_settings
    if not os.path.isfile(self.cctbxpython):
      from .qt import QInputDialog
      self.cctbxpython, ok = QInputDialog.getText(None, "cctbx.python needs specifying",
              'The HKLviewer GUI needs to know where the cctbx.python is located.\n' +
              'Enter the full path for the executable cctbx.python dispatcher file.\n' +
              'Tip: Use the "which" or "where" command from a shell with an active CCTBX environment.')
    if not os.path.isfile(self.cctbxpython):
      raise Exception("The file, %s, does not exists!\n" %self.cctbxpython)
    print("HKLviewer using cctbx.python from: %s" %self.cctbxpython)

    self.settings.beginGroup("MillerTableColumnHeader")
    keys = self.settings.childKeys()
    for philname in keys:
      self.colnames_select_dict[philname] = int(self.settings.value(philname, 1))
    self.settings.endGroup() # MillerTableColumnHeader
    # In case of more than one PySide2 installation tag the settings by version number of PySide2
    # as different versions seem too use different scaling for font and window sizes
    # But in case of using factory defaults then just look for the "Qt" group
    self.settings.beginGroup(Qtversion)
    self.settings.beginGroup("DataTypesGroups")
    datatypes = self.settings.childGroups()
    if datatypes:
      for datatype in datatypes:
        self.datatypedict[ datatype ] = [ self.settings.value(datatype + "/ColourChart", "brg"),
                                      float(self.settings.value(datatype + "/ColourPowerScale", 1.0)),
                                      float(self.settings.value(datatype + "/PowerScale", 1.0)),
                                      float(self.settings.value(datatype + "/RadiiScale", 1.0)),
                                    ]
    self.settings.endGroup()
    self.QWebEngineViewFlags = self.settings.value("QWebEngineViewFlags", None)
    self.mousespeed = float(self.settings.value("MouseSpeed", 0.3))
    self.textinfosize = int(self.settings.value("TextBufferSize", 30))
    # Get a boolean from the stored value which on windows is an int but a string on linux.
    # Do this by casting it into a string and compare with "1"
    self.wraptextinfo = (str(self.settings.value("WordWrapTextInfo", "0")) == "1")
    self.fontsize = float(self.settings.value("FontSize", 10))
    bcol = eval(self.settings.value("BackgroundColour", "(127,127,127,255)"))
    self.backgroundcolour = QColor(*bcol)
    self.browserfontsize = float(self.settings.value("BrowserFontSize", 9))
    self.vectorwidth = float(self.settings.value("VectorWidth", 5))
    self.ttip_click_invoke = self.settings.value("ttip_click_invoke", None)
    self.geometry = self.settings.value("geometry", None)
    self.windowstate = self.settings.value("windowstate", None)
    self.splitter1sizes = self.settings.value("splitter1Sizes", None)
    self.splitter2sizes = self.settings.value("splitter2Sizes", None)
    self.settings.endGroup()
    if use_factory_default_settings:
      # Revert to storing settings in the default Qsettings location such as
      # Windows: HKEY_CURRENT_USER\Software\ , Linux: $HOME/.config/ or MacOS: $HOME/Library/Preferences/
      # Do this by constructing a Qsettings object with our program scope
      self.settings = QSettings("CCTBX", "HKLviewer" )
    if self.QWebEngineViewFlags is None: # avoid doing this test over and over again on the same PC
      flgs = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "")
      self.QWebEngineViewFlags = flgs + " --disable-web-security" #\
      #  + " --enable-webgl-software-rendering --disable-gpu-compositing" \
      #  + " --disable_chromium_framebuffer_multisample --use-gl=swiftshader" \
      #  + " --swiftshader --swiftshader-webgl --ignore-gpu-blacklist"

      os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = self.QWebEngineViewFlags
      # QTWEBENGINE_CHROMIUM_FLAGS environment is now set for TestWebGL()
      if not self.isembedded:
        print("Testing if WebGL works in QWebEngineView....")
        #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
        if not TestWebGL(): # try software rendering
          self.QWebEngineViewFlags = flgs + " --disable-web-security" \
            + " --enable-webgl-software-rendering --disable-gpu-compositing" \
            + " --disable_chromium_framebuffer_multisample --use-gl=swiftshader" \
            + " --swiftshader --swiftshader-webgl --ignore-gpu-blocklist"

          os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = self.QWebEngineViewFlags
          if not TestWebGL():
            print("FATAL ERROR: WebGL does not work in QWebEngineView on this platform!")
            return False
        print(" It does!")
    for arg in sys.argv[1:]:
      if "verbose" in arg:
         print("using flags for QWebEngineView: " + os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] )
    return True


  def UsePersistedQsettings(self):
    # Now assign the users persisted settings to the GUI
    if self.backgroundcolour is not None:
      self.BgrndColourDlg.setCurrentColor(self.backgroundcolour)
      self.onBackgroundColourChanged(self.backgroundcolour)
    if self.mousespeed is not None:
      self.mousemoveslider.setValue(float(self.mousespeed)*self.mousespeedscale)
      self.mousesensitxtbox.setText("%2.1f" %(self.mousemoveslider.value()*10.0/self.mousemoveslider.maximum()) )
      self.send_message('NGL.mouse_sensitivity = %s' %self.mousespeed)
    if self.wraptextinfo is not None:
      self.textInfo.setLineWrapMode(QPlainTextEdit.WidgetWidth if self.wraptextinfo else QPlainTextEdit.NoWrap)
    if self.textinfosize is not None:
      self.onTextbufferSizeChanged(int(self.textinfosize))
      self.bufsizespinBox.setValue(int(self.textinfosize))
    if self.fontsize is not None:
      #self.onFontsizeChanged(int(self.fontsize))
      self.fontspinBox.setValue(int(self.fontsize))
      self.onFontsizeChanged()
    if self.browserfontsize is not None:
      self.onBrowserFontsizeChanged(int(self.browserfontsize))
      self.browserfontspinBox.setValue(int(self.browserfontsize))
    if self.vectorwidth is not None:
      self.onVectorWidthChanged(int(self.vectorwidth))
      self.vectorWidthspinBox.setValue(int(self.vectorwidth))
    if self.ttip_click_invoke is not None:
      self.onShowTooltips(self.ttip_click_invoke)
      self.ttipClickradio.setChecked(self.ttip_click_invoke == "click")
      self.ttipHoverradio.setChecked(self.ttip_click_invoke == "hover")
    if self.splitter1sizes is not None and self.splitter2sizes is not None and \
     self.windowstate is not None:
      self.window.restoreGeometry(self.geometry) # this restores layout of main window
      self.window.restoreState(self.windowstate) # as well as ControlsWidget
      if self.webpagedebugform and self.devmode:
        self.webpagedebugform.resize( self.window.size())
      self.splitter.restoreState(self.splitter1sizes)
      self.splitter_2.restoreState(self.splitter2sizes)
      # must make ControlsWidget visible explicitly for restoring it
      self.ControlsWidget.setVisible(True)
      self.InfoWidget.setVisible(True)
    self.setDatatypedict(self.datatypedict)
    if self.make_new_factory_default_settings:
      # Create a new Factory default settings .ini file to be stored alongside this source file.
      self.PersistQsettings(True)


  @staticmethod
  def RemoveQsettings(all=False):
    mstr = NGL_HKLViewer.Qtversion
    if all:
      mstr = ""
    NGL_HKLViewer.settings.remove(mstr)
    print("HKLviewer settings removed. Using factory defaults next time.")
    NGL_HKLViewer.reset_to_factorydefaults = True


# test for any necessary flags for WebGL to work on this platform
def TestWebGL():
  #print("QTWEBENGINE_CHROMIUM_FLAGS = %s" %os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] )
  QtChromiumCheck_fpath = os.path.join(os.path.split(hklviewer_gui.__file__)[0], "qt_chromium_check.py")
  cmdargs = [ sys.executable, QtChromiumCheck_fpath ]
  webglproc = subprocess.Popen( cmdargs, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  procout, procerr = webglproc.communicate()
  if "WebGL works" in procout.decode():
    return True
  print(procout.decode())
  return False


def run(isembedded=False, chimeraxsession=None):
  import time
  #time.sleep(15) # enough time for attaching debugger

  # TODO: rewrite this function at some point to use python's argparse module
  try:
    kwargs = dict(arg.split('=') for arg in sys.argv if '=' in arg)
    # if an argument is a filename then have it as a keyword argument and assume it's a reflection file
    for arg in sys.argv[1:]:
      if '=' not in arg:
      # if so add it as a keyword argument
        if os.path.isfile(arg):
          if not kwargs.get('hklin', False):
            kwargs['hklin'] = arg

    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "")
    argstr = " ".join(sys.argv[1:])

    if "remove_settings" in argstr:
      NGL_HKLViewer.RemoveQsettings()
      sys.exit()
    if "devmode" in argstr or "debug" in argstr:
      os.environ["PYTHONASYNCIODEBUG"] = "1" # print debug output from asyncio used in webbrowser_messenger_py3
    if "devmode" in argstr or "debug" in argstr and not "UseOSBrowser" in argstr:
      # some useful flags as per https://doc.qt.io/qt-5/qtwebengine-debugging.html
      if "debug" in argstr:
        os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--remote-debugging-port=9741 --single-process --js-flags='--expose_gc'"
      if "devmode" in argstr: # Also start our WebEngineDebugForm
# Don't use --single-process as it will freeze the WebEngineDebugForm when reaching user defined JavaScript breakpoints
        os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--js-flags='--expose_gc'"
    closingtime = None
    if kwargs.get('closing_time', False): # close when time is up during regression tests
      closingtime = int(kwargs['closing_time']) * 2000 # miliseconds

    from .qt import QApplication
    # ensure QWebEngineView scales correctly on a screen with high DPI
    if not isembedded:
      QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    HKLguiobj = NGL_HKLViewer(app, isembedded)
    if not isembedded:
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
      HKLguiobj.chimeraxsession = chimeraxsession
    # Call HKLguiobj.UsePersistedQsettings() but through QTimer so it happens after
    # the QApplication eventloop has started as to ensure resizing according to persisted
    # font size is done properly
    QTimer.singleShot(1000, HKLguiobj.UsePersistedQsettings)
    # For regression tests close us after a specified time
    if closingtime:
      QTimer.singleShot(closingtime, HKLguiobj.closeEvent)

    if isembedded:
      return HKLguiobj

    ret = app.exec_()

  except Exception as e:
    print( str(e)  +  traceback.format_exc(limit=10) )
