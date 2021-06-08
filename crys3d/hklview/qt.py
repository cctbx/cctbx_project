
import sys

if 'PyQt5' in sys.modules:
  # PyQt5
  from PyQt5.QtWebEngineWidgets import QWebEngineView
  from PyQt5.QtCore import QCoreApplication, QMetaObject, QRect, QSize, Qt
  from PyQt5.QtWidgets import QAbstractItemView, QAction, QCheckBox, QComboBox, \
    QDoubleSpinBox, QFrame, QGridLayout, QGroupBox, QLabel, QMenu, QMenuBar, \
    QPlainTextEdit, QPushButton, QRadioButton, QSlider, QSplitter, QSizePolicy, QSpinBox, \
    QStatusBar, QTableWidget, QTabWidget, QTextEdit, QWidget
  from PyQt5.QtGui import QIcon
  from PyQt5.QtCore import Qt, QEvent, QItemSelectionModel, QSize, QSettings, QTimer
  from PyQt5.QtWidgets import (  QAction, QCheckBox, QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMainWindow, QMenu, QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget )
  from PyQt5.QtGui import QColor, QFont, QCursor, QDesktopServices
  from PyQt5.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )
else:
  # PySide2
  from PySide2.QtWebEngineWidgets import QWebEngineView
  from PySide2.QtCore import QCoreApplication, QMetaObject, QRect, QSize, Qt
  from PySide2.QtWidgets import QAbstractItemView, QAction, QCheckBox, QComboBox, \
    QDoubleSpinBox, QFrame, QGridLayout, QGroupBox, QLabel, QMenu, QMenuBar, \
    QPlainTextEdit, QPushButton, QRadioButton, QSlider, QSplitter, QSizePolicy, QSpinBox, \
    QStatusBar, QTableWidget, QTabWidget, QTextEdit, QWidget
  from PySide2.QtGui import QIcon

  from PySide2.QtCore import Qt, QEvent, QItemSelectionModel, QSize, QSettings, QTimer
  from PySide2.QtWidgets import (  QAction, QCheckBox, QComboBox, QDialog,
        QFileDialog, QGridLayout, QGroupBox, QHeaderView, QHBoxLayout, QLabel, QLineEdit,
        QMainWindow, QMenu, QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableView, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QTextBrowser, QWidget )
  from PySide2.QtGui import QColor, QFont, QCursor, QDesktopServices
  from PySide2.QtWebEngineWidgets import ( QWebEngineView, QWebEngineProfile, QWebEnginePage )

