from __future__ import absolute_import, division

try:
  import PyQt5
  # PyQt5
  from PyQt5 import ( QtCore, QtWidgets )   # special import
  from PyQt5.QtCore import  QAbstractTableModel, QCoreApplication, QMetaObject, QModelIndex, \
        QRect, Qt, QEvent, QItemSelectionModel, QUrl, QSize, QSettings, QTimer  # special import
  from PyQt5.QtWidgets import ( QAction, QAbstractItemView, QApplication, QCheckBox, QColorDialog, QComboBox,   # special import
        QDialog, QDockWidget, QDoubleSpinBox, QFileDialog, QFrame, QGridLayout, QGroupBox, QHeaderView,
        QHBoxLayout, QInputDialog, QLabel, QLineEdit, QMainWindow, QMenu, QMenuBar, QMessageBox,
        QPlainTextEdit, QAbstractScrollArea,
        QProgressBar, QPushButton, QRadioButton, QScrollArea, QScrollBar, QSizePolicy,
        QSlider, QSplitter, QSpinBox, QStatusBar, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget, QVBoxLayout )
  from PyQt5.QtGui import ( QBrush, QCloseEvent, QColor, QFont, QCursor, QDesktopServices, QIcon,    # special import
                           QKeySequence, QPalette, QTextDocument )
  from PyQt5.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )   # special import

except Exception:
  # PySide2
  from PySide2 import ( QtCore, QtWidgets )   # special import
  from PySide2.QtCore import Qt, QEvent, QAbstractTableModel, QCoreApplication, QMetaObject, \
        QModelIndex, QUrl, QRect, QItemSelectionModel, QSize, QSettings, QTimer   # special import
  from PySide2.QtWidgets import ( QAction, QAbstractItemView, QApplication, QCheckBox, QColorDialog, QComboBox,   # special import
        QDialog, QDockWidget, QDoubleSpinBox, QFileDialog, QFrame, QGridLayout, QGroupBox, QHeaderView,
        QHBoxLayout, QInputDialog, QLabel, QLineEdit, QMainWindow, QMenu, QMenuBar,  QMessageBox,
        QPlainTextEdit, QAbstractScrollArea,
        QProgressBar, QPushButton, QRadioButton, QScrollArea, QScrollBar, QSizePolicy,
        QSlider, QSplitter, QSpinBox, QStatusBar, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget, QVBoxLayout )
  from PySide2.QtGui import ( QBrush, QCloseEvent, QColor, QFont, QCursor, QDesktopServices, QIcon,    # special import
                             QKeySequence, QPalette, QTextDocument )
  from PySide2.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )   # special import
