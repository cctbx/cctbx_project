from __future__ import absolute_import, division
import sys

if 'PyQt5' in sys.modules:
  qt_mod = PyQT5
else:
  qt_mod = PySide2

from qt_mod import ( QtCore, QtWidgets )
from qt_mod.QtCore import  QAbstractTableModel, QCoreApplication, QMetaObject, QModelIndex,  # special import
      QRect, Qt, QEvent, QItemSelectionModel, QUrl, QSize, QSettings, QTimer
from qt_mod.QtWidgets import ( QAction, QAbstractItemView, QApplication, QCheckBox, QComboBox,  # special import
      QDialog, QDoubleSpinBox, QFileDialog, QFrame, QGridLayout, QGroupBox, QHeaderView,
      QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMenu, QMenuBar, QPlainTextEdit,
      QProgressBar, QPushButton, QRadioButton, QScrollArea, QScrollBar, QSizePolicy,
      QSlider, QSplitter, QSpinBox, QStatusBar, QStyleFactory, QTableView, QTableWidget,
      QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget, QVBoxLayout )
from qt_mod.QtGui import ( QColor, QFont, QCursor, QDesktopServices, QIcon, QKeySequence )   # implicit import
from qt_mod.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )  # special import
