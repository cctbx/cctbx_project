
import sys

if 'PyQt5' in sys.modules:
  # PyQt5
  from PyQt5.QtWebEngineWidgets import QWebEngineView
  from PyQt5.QtCore import QCoreApplication, QMetaObject, QModelIndex, QRect, QSize, Qt
  from PyQt5.QtWidgets import QAbstractItemView, QAction, QApplication, QCheckBox, QComboBox, \
    QDoubleSpinBox, QFrame, QGridLayout, QGroupBox, QLabel, QMenu, QMenuBar, \
    QPlainTextEdit, QPushButton, QRadioButton, QSlider, QSplitter, QSizePolicy, QSpinBox, \
    QStatusBar, QTableWidget, QTabWidget, QTextEdit, QWidget
  from PyQt5.QtGui import QIcon, QCursor, QKeySequence
  from PyQt5 import QtCore, QtWidgets
  from PyQt5.QtCore import QAbstractTableModel, Qt, QEvent, QItemSelectionModel, QUrl, QSize, QSettings, QTimer
  from PyQt5.QtWidgets import (  QAction, QCheckBox, QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMainWindow, QMenu, QProgressBar, QPushButton, QRadioButton, QScrollArea, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget, QVBoxLayout )
  from PyQt5.QtGui import QColor, QFont, QCursor, QDesktopServices
  from PyQt5.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )


else:
  # PySide2
  from PySide2.QtWebEngineWidgets import QWebEngineView
  from PySide2.QtCore import QAbstractTableModel, QCoreApplication, QMetaObject, QModelIndex, QUrl, QRect, QSize, Qt
  from PySide2.QtWidgets import QAbstractItemView, QAction, QApplication, QCheckBox, QComboBox, \
    QDoubleSpinBox, QFrame, QGridLayout, QGroupBox, QLabel, QMenu, QMenuBar, \
    QPlainTextEdit, QPushButton, QRadioButton, QSlider, QSplitter, QSizePolicy, QSpinBox, \
    QStatusBar, QTableWidget, QTabWidget, QTextEdit, QWidget
  from PySide2.QtGui import QIcon, QCursor, QKeySequence
  from PySide2 import QtCore, QtWidgets
  from PySide2.QtCore import Qt, QEvent, QItemSelectionModel, QSize, QSettings, QTimer
  from PySide2.QtWidgets import (  QAction, QCheckBox, QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMainWindow, QMenu, QProgressBar, QPushButton, QRadioButton, QScrollArea, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget, QVBoxLayout )
  from PySide2.QtGui import QColor, QFont, QCursor, QDesktopServices
  from PySide2.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )

