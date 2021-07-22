from __future__ import division
import sys

if 'PyQt5' in sys.modules:
  # PyQt5
  from PyQt5.QtWebEngineWidgets import QWebEngineView # import dependency
  from PyQt5.QtCore import ( QCoreApplication, QMetaObject, QModelIndex, QRect, QSize, Qt # import dependency
  from PyQt5.QtWidgets import QAbstractItemView, QAction, QApplication, QCheckBox, QComboBox, \
    QDoubleSpinBox, QFrame, QGridLayout, QGroupBox, QLabel, QMenu, QMenuBar, \
    QPlainTextEdit, QPushButton, QRadioButton, QSlider, QSplitter, QSizePolicy, QSpinBox, \
    QStatusBar, QTableWidget, QTabWidget, QTextEdit, QWidget ) # import dependency
  from PyQt5.QtGui import QIcon, QCursor, QKeySequence # import dependency
  from PyQt5 import QtCore, QtWidgets # import dependency
  from PyQt5.QtCore import QAbstractTableModel, Qt, QEvent, QItemSelectionModel, QUrl, QSize, QSettings, QTimer # import dependency
  from PyQt5.QtWidgets import (  QAction, QCheckBox, QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMainWindow, QMenu, QProgressBar, QPushButton, QRadioButton, QScrollArea, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget, QVBoxLayout ) # import dependency
  from PyQt5.QtGui import QColor, QFont, QCursor, QDesktopServices # import dependency
  from PyQt5.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage ) # import dependency


else:
  # PySide2
  from PySide2.QtWebEngineWidgets import QWebEngineView # import dependency
  from PySide2.QtCore import QAbstractTableModel, QCoreApplication, QMetaObject, QModelIndex, QUrl, QRect, QSize, Qt # import dependency
  from PySide2.QtWidgets import ( QAbstractItemView, QAction, QApplication, QCheckBox, QComboBox, \
    QDoubleSpinBox, QFrame, QGridLayout, QGroupBox, QLabel, QMenu, QMenuBar, \
    QPlainTextEdit, QPushButton, QRadioButton, QSlider, QSplitter, QSizePolicy, QSpinBox, \
    QStatusBar, QTableWidget, QTabWidget, QTextEdit, QWidget )# import dependency
  from PySide2.QtGui import QIcon, QCursor, QKeySequence # import dependency
  from PySide2 import QtCore, QtWidgets # import dependency
  from PySide2.QtCore import Qt, QEvent, QItemSelectionModel, QSize, QSettings, QTimer # import dependency
  from PySide2.QtWidgets import (  QAction, QCheckBox, QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMainWindow, QMenu, QProgressBar, QPushButton, QRadioButton, QScrollArea, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget, QVBoxLayout ) # import dependency
  from PySide2.QtGui import QColor, QFont, QCursor, QDesktopServices # import dependency
  from PySide2.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage ) # import dependency

