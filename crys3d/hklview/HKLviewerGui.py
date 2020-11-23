# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'HKLviewer.ui'
##
## Created by: Qt User Interface Compiler version 5.14.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from __future__ import absolute_import, division, print_function

from PySide2.QtWebEngineWidgets import QWebEngineView
try: # if invoked by cctbx.python or some such
  from crys3d.hklview.helpers import HeaderDataTableWidget
except Exception as e: # if invoked by a generic python that doesn't know cctbx modules
  from helpers import HeaderDataTableWidge

from PySide2.QtCore import (QCoreApplication, QMetaObject, QObject, QPoint,
    QRect, QSize, QUrl, Qt)
from PySide2.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont,
    QFontDatabase, QIcon, QLinearGradient, QPalette, QPainter, QPixmap,
    QRadialGradient)
from PySide2.QtWidgets import *


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1066, 840)
        self.actionOpen_reflection_file = QAction(MainWindow)
        self.actionOpen_reflection_file.setObjectName(u"actionOpen_reflection_file")
        self.actionSettings = QAction(MainWindow)
        self.actionSettings.setObjectName(u"actionSettings")
        self.actiondebug = QAction(MainWindow)
        self.actiondebug.setObjectName(u"actiondebug")
        self.actionExit = QAction(MainWindow)
        self.actionExit.setObjectName(u"actionExit")
        self.actionSave_reflection_file = QAction(MainWindow)
        self.actionSave_reflection_file.setObjectName(u"actionSave_reflection_file")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.gridLayout_2 = QGridLayout(self.centralwidget)
        self.gridLayout_2.setSpacing(4)
        self.gridLayout_2.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.widget = QWidget(self.centralwidget)
        self.widget.setObjectName(u"widget")
        self.gridLayout = QGridLayout(self.widget)
        self.gridLayout.setSpacing(4)
        self.gridLayout.setContentsMargins(3, 3, 3, 3)
        self.gridLayout.setObjectName(u"gridLayout")
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.splitter = QSplitter(self.widget)
        self.splitter.setObjectName(u"splitter")
        self.splitter.setOrientation(Qt.Horizontal)
        self.widget_2 = QWidget(self.splitter)
        self.widget_2.setObjectName(u"widget_2")
        self.gridLayout_4 = QGridLayout(self.widget_2)
        self.gridLayout_4.setSpacing(4)
        self.gridLayout_4.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_4.setObjectName(u"gridLayout_4")
        self.splitter_2 = QSplitter(self.widget_2)
        self.splitter_2.setObjectName(u"splitter_2")
        self.splitter_2.setFrameShape(QFrame.NoFrame)
        self.splitter_2.setFrameShadow(QFrame.Plain)
        self.splitter_2.setMidLineWidth(0)
        self.splitter_2.setOrientation(Qt.Vertical)
        self.widget_4 = QWidget(self.splitter_2)
        self.widget_4.setObjectName(u"widget_4")
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(3)
        sizePolicy.setHeightForWidth(self.widget_4.sizePolicy().hasHeightForWidth())
        self.widget_4.setSizePolicy(sizePolicy)
        self.gridLayout_5 = QGridLayout(self.widget_4)
        self.gridLayout_5.setSpacing(4)
        self.gridLayout_5.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_5.setObjectName(u"gridLayout_5")
        self.gridLayout_5.setContentsMargins(0, 6, 0, -1)
        self.label = QLabel(self.widget_4)
        self.label.setObjectName(u"label")

        self.gridLayout_5.addWidget(self.label, 0, 0, 1, 1)

        self.millertable = HeaderDataTableWidget(self.widget_4)
        if (self.millertable.columnCount() < 10):
            self.millertable.setColumnCount(10)
        if (self.millertable.rowCount() < 1):
            self.millertable.setRowCount(1)
        self.millertable.setObjectName(u"millertable")
        sizePolicy1 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.millertable.sizePolicy().hasHeightForWidth())
        self.millertable.setSizePolicy(sizePolicy1)
        self.millertable.setVerticalScrollMode(QAbstractItemView.ScrollPerPixel)
        self.millertable.setHorizontalScrollMode(QAbstractItemView.ScrollPerPixel)
        self.millertable.setRowCount(1)
        self.millertable.setColumnCount(10)
        self.millertable.horizontalHeader().setMinimumSectionSize(5)

        self.gridLayout_5.addWidget(self.millertable, 1, 0, 1, 1)

        self.splitter_2.addWidget(self.widget_4)
        self.functionTabWidget = QTabWidget(self.splitter_2)
        self.functionTabWidget.setObjectName(u"functionTabWidget")
        sizePolicy2 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(1)
        sizePolicy2.setHeightForWidth(self.functionTabWidget.sizePolicy().hasHeightForWidth())
        self.functionTabWidget.setSizePolicy(sizePolicy2)
        self.tab = QWidget()
        self.tab.setObjectName(u"tab")
        self.gridLayout_28 = QGridLayout(self.tab)
        self.gridLayout_28.setSpacing(4)
        self.gridLayout_28.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_28.setObjectName(u"gridLayout_28")
        self.ExpandReflsGroupBox = QGroupBox(self.tab)
        self.ExpandReflsGroupBox.setObjectName(u"ExpandReflsGroupBox")
        sizePolicy3 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.ExpandReflsGroupBox.sizePolicy().hasHeightForWidth())
        self.ExpandReflsGroupBox.setSizePolicy(sizePolicy3)
        self.ExpandReflsGroupBox.setCheckable(True)
        self.ExpandReflsGroupBox.setChecked(False)
        self.gridLayout_3 = QGridLayout(self.ExpandReflsGroupBox)
        self.gridLayout_3.setSpacing(4)
        self.gridLayout_3.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.gridLayout_3.setContentsMargins(6, 6, 6, 6)
        self.expandP1checkbox = QCheckBox(self.ExpandReflsGroupBox)
        self.expandP1checkbox.setObjectName(u"expandP1checkbox")

        self.gridLayout_3.addWidget(self.expandP1checkbox, 0, 0, 1, 1)

        self.expandAnomalouscheckbox = QCheckBox(self.ExpandReflsGroupBox)
        self.expandAnomalouscheckbox.setObjectName(u"expandAnomalouscheckbox")

        self.gridLayout_3.addWidget(self.expandAnomalouscheckbox, 0, 1, 1, 1)


        self.gridLayout_28.addWidget(self.ExpandReflsGroupBox, 0, 0, 1, 3)

        self.SpacegroupLabel = QLabel(self.tab)
        self.SpacegroupLabel.setObjectName(u"SpacegroupLabel")
        sizePolicy4 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(2)
        sizePolicy4.setHeightForWidth(self.SpacegroupLabel.sizePolicy().hasHeightForWidth())
        self.SpacegroupLabel.setSizePolicy(sizePolicy4)

        self.gridLayout_28.addWidget(self.SpacegroupLabel, 2, 0, 1, 1)

        self.sysabsentcheckbox = QCheckBox(self.tab)
        self.sysabsentcheckbox.setObjectName(u"sysabsentcheckbox")

        self.gridLayout_28.addWidget(self.sysabsentcheckbox, 1, 0, 1, 1)

        self.widget_3 = QWidget(self.tab)
        self.widget_3.setObjectName(u"widget_3")
        self.gridLayout_20 = QGridLayout(self.widget_3)
        self.gridLayout_20.setSpacing(4)
        self.gridLayout_20.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_20.setObjectName(u"gridLayout_20")
        self.gridLayout_20.setContentsMargins(0, 0, 0, 0)
        self.missingcheckbox = QCheckBox(self.widget_3)
        self.missingcheckbox.setObjectName(u"missingcheckbox")

        self.gridLayout_20.addWidget(self.missingcheckbox, 0, 0, 2, 2)

        self.onlymissingcheckbox = QCheckBox(self.widget_3)
        self.onlymissingcheckbox.setObjectName(u"onlymissingcheckbox")

        self.gridLayout_20.addWidget(self.onlymissingcheckbox, 2, 0, 1, 1)


        self.gridLayout_28.addWidget(self.widget_3, 1, 1, 1, 2)

        self.SpaceGroupComboBox = QComboBox(self.tab)
        self.SpaceGroupComboBox.setObjectName(u"SpaceGroupComboBox")
        sizePolicy5 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(2)
        sizePolicy5.setHeightForWidth(self.SpaceGroupComboBox.sizePolicy().hasHeightForWidth())
        self.SpaceGroupComboBox.setSizePolicy(sizePolicy5)

        self.gridLayout_28.addWidget(self.SpaceGroupComboBox, 2, 1, 1, 2)

        self.functionTabWidget.addTab(self.tab, "")
        self.tab_2 = QWidget()
        self.tab_2.setObjectName(u"tab_2")
        self.gridLayout_9 = QGridLayout(self.tab_2)
        self.gridLayout_9.setSpacing(4)
        self.gridLayout_9.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_9.setObjectName(u"gridLayout_9")
        self.SliceReflectionsBox = QCheckBox(self.tab_2)
        self.SliceReflectionsBox.setObjectName(u"SliceReflectionsBox")
        sizePolicy6 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        sizePolicy6.setHorizontalStretch(1)
        sizePolicy6.setVerticalStretch(10)
        sizePolicy6.setHeightForWidth(self.SliceReflectionsBox.sizePolicy().hasHeightForWidth())
        self.SliceReflectionsBox.setSizePolicy(sizePolicy6)

        self.gridLayout_9.addWidget(self.SliceReflectionsBox, 0, 0, 1, 1)

        self.fixedorientcheckbox = QCheckBox(self.tab_2)
        self.fixedorientcheckbox.setObjectName(u"fixedorientcheckbox")
        sizePolicy7 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        sizePolicy7.setHorizontalStretch(3)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.fixedorientcheckbox.sizePolicy().hasHeightForWidth())
        self.fixedorientcheckbox.setSizePolicy(sizePolicy7)

        self.gridLayout_9.addWidget(self.fixedorientcheckbox, 0, 1, 1, 1)

        self.ClipPlaneChkGroupBox = QGroupBox(self.tab_2)
        self.ClipPlaneChkGroupBox.setObjectName(u"ClipPlaneChkGroupBox")
        sizePolicy8 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy8.setHorizontalStretch(0)
        sizePolicy8.setVerticalStretch(1)
        sizePolicy8.setHeightForWidth(self.ClipPlaneChkGroupBox.sizePolicy().hasHeightForWidth())
        self.ClipPlaneChkGroupBox.setSizePolicy(sizePolicy8)
        self.ClipPlaneChkGroupBox.setCheckable(True)
        self.gridLayout_7 = QGridLayout(self.ClipPlaneChkGroupBox)
        self.gridLayout_7.setSpacing(4)
        self.gridLayout_7.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.gridLayout_7.setContentsMargins(3, 3, 3, 6)
        self.clipNormalBtn = QRadioButton(self.ClipPlaneChkGroupBox)
        self.clipNormalBtn.setObjectName(u"clipNormalBtn")
        sizePolicy9 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy9.setHorizontalStretch(0)
        sizePolicy9.setVerticalStretch(0)
        sizePolicy9.setHeightForWidth(self.clipNormalBtn.sizePolicy().hasHeightForWidth())
        self.clipNormalBtn.setSizePolicy(sizePolicy9)

        self.gridLayout_7.addWidget(self.clipNormalBtn, 0, 1, 1, 1)

        self.clipParallelBtn = QRadioButton(self.ClipPlaneChkGroupBox)
        self.clipParallelBtn.setObjectName(u"clipParallelBtn")
        sizePolicy9.setHeightForWidth(self.clipParallelBtn.sizePolicy().hasHeightForWidth())
        self.clipParallelBtn.setSizePolicy(sizePolicy9)

        self.gridLayout_7.addWidget(self.clipParallelBtn, 0, 0, 1, 1)

        self.groupBox_2 = QGroupBox(self.ClipPlaneChkGroupBox)
        self.groupBox_2.setObjectName(u"groupBox_2")
        self.gridLayout_19 = QGridLayout(self.groupBox_2)
        self.gridLayout_19.setSpacing(4)
        self.gridLayout_19.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_19.setObjectName(u"gridLayout_19")
        self.frame_2 = QFrame(self.groupBox_2)
        self.frame_2.setObjectName(u"frame_2")
        sizePolicy3.setHeightForWidth(self.frame_2.sizePolicy().hasHeightForWidth())
        self.frame_2.setSizePolicy(sizePolicy3)
        self.frame_2.setFrameShape(QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QFrame.Raised)
        self.gridLayout_15 = QGridLayout(self.frame_2)
        self.gridLayout_15.setSpacing(4)
        self.gridLayout_15.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_15.setObjectName(u"gridLayout_15")
        self.gridLayout_15.setContentsMargins(0, 0, 0, 0)
        self.recipvecBtn = QRadioButton(self.frame_2)
        self.recipvecBtn.setObjectName(u"recipvecBtn")

        self.gridLayout_15.addWidget(self.recipvecBtn, 0, 0, 1, 1)

        self.clipTNCSBtn = QRadioButton(self.frame_2)
        self.clipTNCSBtn.setObjectName(u"clipTNCSBtn")

        self.gridLayout_15.addWidget(self.clipTNCSBtn, 0, 1, 1, 1)

        self.realspacevecBtn = QRadioButton(self.frame_2)
        self.realspacevecBtn.setObjectName(u"realspacevecBtn")

        self.gridLayout_15.addWidget(self.realspacevecBtn, 1, 0, 1, 1)


        self.gridLayout_19.addWidget(self.frame_2, 0, 0, 1, 1)

        self.frame_6 = QFrame(self.groupBox_2)
        self.frame_6.setObjectName(u"frame_6")
        self.frame_6.setFrameShape(QFrame.StyledPanel)
        self.frame_6.setFrameShadow(QFrame.Raised)
        self.gridLayout_18 = QGridLayout(self.frame_6)
        self.gridLayout_18.setSpacing(4)
        self.gridLayout_18.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_18.setObjectName(u"gridLayout_18")
        self.gridLayout_18.setContentsMargins(0, 0, 0, 0)
        self.rotavecangle_labeltxt = QLabel(self.frame_6)
        self.rotavecangle_labeltxt.setObjectName(u"rotavecangle_labeltxt")
        sizePolicy10 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy10.setHorizontalStretch(0)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(self.rotavecangle_labeltxt.sizePolicy().hasHeightForWidth())
        self.rotavecangle_labeltxt.setSizePolicy(sizePolicy10)

        self.gridLayout_18.addWidget(self.rotavecangle_labeltxt, 0, 0, 1, 1)

        self.rotavecangle_slider = QSlider(self.frame_6)
        self.rotavecangle_slider.setObjectName(u"rotavecangle_slider")
        sizePolicy11 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy11.setHorizontalStretch(1)
        sizePolicy11.setVerticalStretch(0)
        sizePolicy11.setHeightForWidth(self.rotavecangle_slider.sizePolicy().hasHeightForWidth())
        self.rotavecangle_slider.setSizePolicy(sizePolicy11)
        self.rotavecangle_slider.setMaximum(359)
        self.rotavecangle_slider.setOrientation(Qt.Horizontal)
        self.rotavecangle_slider.setTickPosition(QSlider.TicksAbove)
        self.rotavecangle_slider.setTickInterval(18)

        self.gridLayout_18.addWidget(self.rotavecangle_slider, 0, 1, 1, 1)


        self.gridLayout_19.addWidget(self.frame_6, 3, 0, 1, 1)

        self.frame_3 = QFrame(self.groupBox_2)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setFrameShape(QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QFrame.Raised)
        self.gridLayout_16 = QGridLayout(self.frame_3)
        self.gridLayout_16.setSpacing(4)
        self.gridLayout_16.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_16.setObjectName(u"gridLayout_16")
        self.gridLayout_16.setContentsMargins(0, 0, 0, 0)
        self.frame_9 = QFrame(self.frame_3)
        self.frame_9.setObjectName(u"frame_9")
        sizePolicy12 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        sizePolicy12.setHorizontalStretch(0)
        sizePolicy12.setVerticalStretch(0)
        sizePolicy12.setHeightForWidth(self.frame_9.sizePolicy().hasHeightForWidth())
        self.frame_9.setSizePolicy(sizePolicy12)
        self.frame_9.setFrameShape(QFrame.StyledPanel)
        self.frame_9.setFrameShadow(QFrame.Raised)
        self.gridLayout_22 = QGridLayout(self.frame_9)
        self.gridLayout_22.setSpacing(4)
        self.gridLayout_22.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_22.setObjectName(u"gridLayout_22")
        self.gridLayout_22.setContentsMargins(0, 0, 0, 0)
        self.label_21 = QLabel(self.frame_9)
        self.label_21.setObjectName(u"label_21")
        sizePolicy13 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Preferred)
        sizePolicy13.setHorizontalStretch(0)
        sizePolicy13.setVerticalStretch(0)
        sizePolicy13.setHeightForWidth(self.label_21.sizePolicy().hasHeightForWidth())
        self.label_21.setSizePolicy(sizePolicy13)

        self.gridLayout_22.addWidget(self.label_21, 0, 0, 1, 1)

        self.hvec_spinBox = QDoubleSpinBox(self.frame_9)
        self.hvec_spinBox.setObjectName(u"hvec_spinBox")
        sizePolicy13.setHeightForWidth(self.hvec_spinBox.sizePolicy().hasHeightForWidth())
        self.hvec_spinBox.setSizePolicy(sizePolicy13)

        self.gridLayout_22.addWidget(self.hvec_spinBox, 0, 1, 1, 1)


        self.gridLayout_16.addWidget(self.frame_9, 0, 0, 1, 1)

        self.frame_10 = QFrame(self.frame_3)
        self.frame_10.setObjectName(u"frame_10")
        sizePolicy12.setHeightForWidth(self.frame_10.sizePolicy().hasHeightForWidth())
        self.frame_10.setSizePolicy(sizePolicy12)
        self.frame_10.setFrameShape(QFrame.StyledPanel)
        self.frame_10.setFrameShadow(QFrame.Raised)
        self.gridLayout_23 = QGridLayout(self.frame_10)
        self.gridLayout_23.setSpacing(4)
        self.gridLayout_23.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_23.setObjectName(u"gridLayout_23")
        self.gridLayout_23.setContentsMargins(0, 0, 0, 0)
        self.label_22 = QLabel(self.frame_10)
        self.label_22.setObjectName(u"label_22")
        sizePolicy13.setHeightForWidth(self.label_22.sizePolicy().hasHeightForWidth())
        self.label_22.setSizePolicy(sizePolicy13)

        self.gridLayout_23.addWidget(self.label_22, 0, 0, 1, 1)

        self.kvec_spinBox = QDoubleSpinBox(self.frame_10)
        self.kvec_spinBox.setObjectName(u"kvec_spinBox")
        sizePolicy13.setHeightForWidth(self.kvec_spinBox.sizePolicy().hasHeightForWidth())
        self.kvec_spinBox.setSizePolicy(sizePolicy13)

        self.gridLayout_23.addWidget(self.kvec_spinBox, 0, 1, 1, 1)


        self.gridLayout_16.addWidget(self.frame_10, 0, 1, 1, 1)

        self.frame_11 = QFrame(self.frame_3)
        self.frame_11.setObjectName(u"frame_11")
        sizePolicy12.setHeightForWidth(self.frame_11.sizePolicy().hasHeightForWidth())
        self.frame_11.setSizePolicy(sizePolicy12)
        self.frame_11.setFrameShape(QFrame.StyledPanel)
        self.frame_11.setFrameShadow(QFrame.Raised)
        self.gridLayout_24 = QGridLayout(self.frame_11)
        self.gridLayout_24.setSpacing(4)
        self.gridLayout_24.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_24.setObjectName(u"gridLayout_24")
        self.gridLayout_24.setContentsMargins(0, 0, 0, 0)
        self.label_23 = QLabel(self.frame_11)
        self.label_23.setObjectName(u"label_23")
        sizePolicy13.setHeightForWidth(self.label_23.sizePolicy().hasHeightForWidth())
        self.label_23.setSizePolicy(sizePolicy13)

        self.gridLayout_24.addWidget(self.label_23, 0, 0, 1, 1)

        self.lvec_spinBox = QDoubleSpinBox(self.frame_11)
        self.lvec_spinBox.setObjectName(u"lvec_spinBox")
        sizePolicy13.setHeightForWidth(self.lvec_spinBox.sizePolicy().hasHeightForWidth())
        self.lvec_spinBox.setSizePolicy(sizePolicy13)

        self.gridLayout_24.addWidget(self.lvec_spinBox, 0, 1, 1, 1)


        self.gridLayout_16.addWidget(self.frame_11, 0, 2, 1, 1)


        self.gridLayout_19.addWidget(self.frame_3, 1, 0, 1, 1)

        self.frame = QFrame(self.groupBox_2)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.gridLayout_29 = QGridLayout(self.frame)
        self.gridLayout_29.setSpacing(4)
        self.gridLayout_29.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_29.setObjectName(u"gridLayout_29")
        self.gridLayout_29.setContentsMargins(0, 0, 0, 0)
        self.label_15 = QLabel(self.frame)
        self.label_15.setObjectName(u"label_15")
        sizePolicy12.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy12)
        self.label_15.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.gridLayout_29.addWidget(self.label_15, 0, 3, 1, 1)

        self.hkldist_spinBox = QDoubleSpinBox(self.frame)
        self.hkldist_spinBox.setObjectName(u"hkldist_spinBox")
        sizePolicy12.setHeightForWidth(self.hkldist_spinBox.sizePolicy().hasHeightForWidth())
        self.hkldist_spinBox.setSizePolicy(sizePolicy12)

        self.gridLayout_29.addWidget(self.hkldist_spinBox, 0, 1, 1, 1)

        self.clipwidth_spinBox = QDoubleSpinBox(self.frame)
        self.clipwidth_spinBox.setObjectName(u"clipwidth_spinBox")
        sizePolicy12.setHeightForWidth(self.clipwidth_spinBox.sizePolicy().hasHeightForWidth())
        self.clipwidth_spinBox.setSizePolicy(sizePolicy12)

        self.gridLayout_29.addWidget(self.clipwidth_spinBox, 0, 4, 1, 1)

        self.label_14 = QLabel(self.frame)
        self.label_14.setObjectName(u"label_14")

        self.gridLayout_29.addWidget(self.label_14, 0, 0, 1, 1)


        self.gridLayout_19.addWidget(self.frame, 2, 0, 1, 1)


        self.gridLayout_7.addWidget(self.groupBox_2, 1, 0, 1, 3)


        self.gridLayout_9.addWidget(self.ClipPlaneChkGroupBox, 1, 0, 1, 2)

        self.showsliceGroupCheckbox = QGroupBox(self.tab_2)
        self.showsliceGroupCheckbox.setObjectName(u"showsliceGroupCheckbox")
        sizePolicy14 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        sizePolicy14.setHorizontalStretch(0)
        sizePolicy14.setVerticalStretch(0)
        sizePolicy14.setHeightForWidth(self.showsliceGroupCheckbox.sizePolicy().hasHeightForWidth())
        self.showsliceGroupCheckbox.setSizePolicy(sizePolicy14)
        self.showsliceGroupCheckbox.setCheckable(True)
        self.showsliceGroupCheckbox.setChecked(True)
        self.gridLayout_6 = QGridLayout(self.showsliceGroupCheckbox)
        self.gridLayout_6.setSpacing(4)
        self.gridLayout_6.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_6.setObjectName(u"gridLayout_6")
        self.gridLayout_6.setContentsMargins(3, 3, 3, 3)
        self.SliceLabelComboBox = QComboBox(self.showsliceGroupCheckbox)
        self.SliceLabelComboBox.setObjectName(u"SliceLabelComboBox")
        sizePolicy3.setHeightForWidth(self.SliceLabelComboBox.sizePolicy().hasHeightForWidth())
        self.SliceLabelComboBox.setSizePolicy(sizePolicy3)

        self.gridLayout_6.addWidget(self.SliceLabelComboBox, 0, 0, 1, 1)

        self.label_3 = QLabel(self.showsliceGroupCheckbox)
        self.label_3.setObjectName(u"label_3")
        sizePolicy10.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy10)
        self.label_3.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.label_3.setIndent(1)

        self.gridLayout_6.addWidget(self.label_3, 0, 1, 1, 1)

        self.sliceindexspinBox = QSpinBox(self.showsliceGroupCheckbox)
        self.sliceindexspinBox.setObjectName(u"sliceindexspinBox")
        sizePolicy10.setHeightForWidth(self.sliceindexspinBox.sizePolicy().hasHeightForWidth())
        self.sliceindexspinBox.setSizePolicy(sizePolicy10)

        self.gridLayout_6.addWidget(self.sliceindexspinBox, 0, 2, 1, 1)


        self.gridLayout_9.addWidget(self.showsliceGroupCheckbox, 2, 0, 1, 2)

        self.functionTabWidget.addTab(self.tab_2, "")
        self.tab_3 = QWidget()
        self.tab_3.setObjectName(u"tab_3")
        self.gridLayout_27 = QGridLayout(self.tab_3)
        self.gridLayout_27.setSpacing(4)
        self.gridLayout_27.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_27.setObjectName(u"gridLayout_27")
        self.groupBox_3 = QGroupBox(self.tab_3)
        self.groupBox_3.setObjectName(u"groupBox_3")
        sizePolicy15 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy15.setHorizontalStretch(1)
        sizePolicy15.setVerticalStretch(0)
        sizePolicy15.setHeightForWidth(self.groupBox_3.sizePolicy().hasHeightForWidth())
        self.groupBox_3.setSizePolicy(sizePolicy15)
        self.groupBox_3.setCheckable(False)
        self.gridLayout_10 = QGridLayout(self.groupBox_3)
        self.gridLayout_10.setSpacing(4)
        self.gridLayout_10.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.label_7 = QLabel(self.groupBox_3)
        self.label_7.setObjectName(u"label_7")
        sizePolicy10.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy10)
        self.label_7.setTextFormat(Qt.AutoText)
        self.label_7.setAlignment(Qt.AlignLeading|Qt.AlignLeft|Qt.AlignVCenter)
        self.label_7.setWordWrap(True)

        self.gridLayout_10.addWidget(self.label_7, 0, 0, 1, 4)

        self.ManualPowerScalecheckbox = QCheckBox(self.groupBox_3)
        self.ManualPowerScalecheckbox.setObjectName(u"ManualPowerScalecheckbox")
        sizePolicy16 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.MinimumExpanding)
        sizePolicy16.setHorizontalStretch(0)
        sizePolicy16.setVerticalStretch(0)
        sizePolicy16.setHeightForWidth(self.ManualPowerScalecheckbox.sizePolicy().hasHeightForWidth())
        self.ManualPowerScalecheckbox.setSizePolicy(sizePolicy16)

        self.gridLayout_10.addWidget(self.ManualPowerScalecheckbox, 1, 0, 1, 2)

        self.label_10 = QLabel(self.groupBox_3)
        self.label_10.setObjectName(u"label_10")
        sizePolicy16.setHeightForWidth(self.label_10.sizePolicy().hasHeightForWidth())
        self.label_10.setSizePolicy(sizePolicy16)
        self.label_10.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.label_10.setIndent(1)

        self.gridLayout_10.addWidget(self.label_10, 1, 2, 1, 1)

        self.power_scale_spinBox = QDoubleSpinBox(self.groupBox_3)
        self.power_scale_spinBox.setObjectName(u"power_scale_spinBox")
        sizePolicy17 = QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        sizePolicy17.setHorizontalStretch(0)
        sizePolicy17.setVerticalStretch(0)
        sizePolicy17.setHeightForWidth(self.power_scale_spinBox.sizePolicy().hasHeightForWidth())
        self.power_scale_spinBox.setSizePolicy(sizePolicy17)
        self.power_scale_spinBox.setMaximum(1.000000000000000)
        self.power_scale_spinBox.setSingleStep(0.050000000000000)
        self.power_scale_spinBox.setValue(0.350000000000000)

        self.gridLayout_10.addWidget(self.power_scale_spinBox, 1, 3, 1, 1)

        self.label_11 = QLabel(self.groupBox_3)
        self.label_11.setObjectName(u"label_11")
        sizePolicy16.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy16)

        self.gridLayout_10.addWidget(self.label_11, 2, 0, 1, 1)

        self.radii_scale_spinBox = QDoubleSpinBox(self.groupBox_3)
        self.radii_scale_spinBox.setObjectName(u"radii_scale_spinBox")
        sizePolicy17.setHeightForWidth(self.radii_scale_spinBox.sizePolicy().hasHeightForWidth())
        self.radii_scale_spinBox.setSizePolicy(sizePolicy17)
        self.radii_scale_spinBox.setMaximum(5.000000000000000)
        self.radii_scale_spinBox.setSingleStep(0.100000000000000)

        self.gridLayout_10.addWidget(self.radii_scale_spinBox, 2, 1, 1, 1)


        self.gridLayout_27.addWidget(self.groupBox_3, 0, 0, 1, 1)

        self.functionTabWidget.addTab(self.tab_3, "")
        self.tab_4 = QWidget()
        self.tab_4.setObjectName(u"tab_4")
        self.gridLayout_13 = QGridLayout(self.tab_4)
        self.gridLayout_13.setSpacing(4)
        self.gridLayout_13.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_13.setObjectName(u"gridLayout_13")
        self.groupBox_5 = QGroupBox(self.tab_4)
        self.groupBox_5.setObjectName(u"groupBox_5")
        sizePolicy3.setHeightForWidth(self.groupBox_5.sizePolicy().hasHeightForWidth())
        self.groupBox_5.setSizePolicy(sizePolicy3)
        self.gridLayout_11 = QGridLayout(self.groupBox_5)
        self.gridLayout_11.setSpacing(4)
        self.gridLayout_11.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_11.setObjectName(u"gridLayout_11")
        self.BinDataComboBox = QComboBox(self.groupBox_5)
        self.BinDataComboBox.setObjectName(u"BinDataComboBox")
        sizePolicy18 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy18.setHorizontalStretch(2)
        sizePolicy18.setVerticalStretch(0)
        sizePolicy18.setHeightForWidth(self.BinDataComboBox.sizePolicy().hasHeightForWidth())
        self.BinDataComboBox.setSizePolicy(sizePolicy18)

        self.gridLayout_11.addWidget(self.BinDataComboBox, 0, 0, 1, 1)

        self.label_12 = QLabel(self.groupBox_5)
        self.label_12.setObjectName(u"label_12")

        self.gridLayout_11.addWidget(self.label_12, 0, 1, 1, 1)

        self.Nbins_spinBox = QSpinBox(self.groupBox_5)
        self.Nbins_spinBox.setObjectName(u"Nbins_spinBox")
        sizePolicy13.setHeightForWidth(self.Nbins_spinBox.sizePolicy().hasHeightForWidth())
        self.Nbins_spinBox.setSizePolicy(sizePolicy13)
        self.Nbins_spinBox.setMinimum(1)
        self.Nbins_spinBox.setMaximum(20)

        self.gridLayout_11.addWidget(self.Nbins_spinBox, 0, 2, 1, 1)


        self.gridLayout_13.addWidget(self.groupBox_5, 0, 0, 1, 1)

        self.widget_6 = QWidget(self.tab_4)
        self.widget_6.setObjectName(u"widget_6")
        sizePolicy8.setHeightForWidth(self.widget_6.sizePolicy().hasHeightForWidth())
        self.widget_6.setSizePolicy(sizePolicy8)
        self.gridLayout_12 = QGridLayout(self.widget_6)
        self.gridLayout_12.setSpacing(4)
        self.gridLayout_12.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_12.setObjectName(u"gridLayout_12")
        self.gridLayout_12.setContentsMargins(3, 3, 3, 3)
        self.OpaqueAllCheckbox = QCheckBox(self.widget_6)
        self.OpaqueAllCheckbox.setObjectName(u"OpaqueAllCheckbox")

        self.gridLayout_12.addWidget(self.OpaqueAllCheckbox, 0, 0, 1, 1)

        self.binstable = QTableWidget(self.widget_6)
        if (self.binstable.columnCount() < 4):
            self.binstable.setColumnCount(4)
        if (self.binstable.rowCount() < 5):
            self.binstable.setRowCount(5)
        self.binstable.setObjectName(u"binstable")
        sizePolicy19 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy19.setHorizontalStretch(0)
        sizePolicy19.setVerticalStretch(1)
        sizePolicy19.setHeightForWidth(self.binstable.sizePolicy().hasHeightForWidth())
        self.binstable.setSizePolicy(sizePolicy19)
        self.binstable.setRowCount(5)
        self.binstable.setColumnCount(4)

        self.gridLayout_12.addWidget(self.binstable, 1, 0, 1, 1)


        self.gridLayout_13.addWidget(self.widget_6, 1, 0, 1, 1)

        self.functionTabWidget.addTab(self.tab_4, "")
        self.tab_5 = QWidget()
        self.tab_5.setObjectName(u"tab_5")
        self.gridLayout_8 = QGridLayout(self.tab_5)
        self.gridLayout_8.setSpacing(4)
        self.gridLayout_8.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_8.setObjectName(u"gridLayout_8")
        self.groupBox_6 = QGroupBox(self.tab_5)
        self.groupBox_6.setObjectName(u"groupBox_6")
        sizePolicy10.setHeightForWidth(self.groupBox_6.sizePolicy().hasHeightForWidth())
        self.groupBox_6.setSizePolicy(sizePolicy10)
        self.groupBox_6.setAlignment(Qt.AlignLeading|Qt.AlignLeft|Qt.AlignVCenter)
        self.gridLayout_21 = QGridLayout(self.groupBox_6)
        self.gridLayout_21.setSpacing(4)
        self.gridLayout_21.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_21.setObjectName(u"gridLayout_21")
        self.DrawReciprocUnitCellBox = QGroupBox(self.groupBox_6)
        self.DrawReciprocUnitCellBox.setObjectName(u"DrawReciprocUnitCellBox")
        self.DrawReciprocUnitCellBox.setCheckable(True)
        self.DrawReciprocUnitCellBox.setChecked(False)
        self.gridLayout_25 = QGridLayout(self.DrawReciprocUnitCellBox)
        self.gridLayout_25.setSpacing(4)
        self.gridLayout_25.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_25.setObjectName(u"gridLayout_25")
        self.label_2 = QLabel(self.DrawReciprocUnitCellBox)
        self.label_2.setObjectName(u"label_2")
        sizePolicy10.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy10)
        self.label_2.setWordWrap(True)

        self.gridLayout_25.addWidget(self.label_2, 0, 0, 1, 1)

        self.reciprocunitcellslider = QSlider(self.DrawReciprocUnitCellBox)
        self.reciprocunitcellslider.setObjectName(u"reciprocunitcellslider")
        sizePolicy20 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy20.setHorizontalStretch(0)
        sizePolicy20.setVerticalStretch(0)
        sizePolicy20.setHeightForWidth(self.reciprocunitcellslider.sizePolicy().hasHeightForWidth())
        self.reciprocunitcellslider.setSizePolicy(sizePolicy20)
        self.reciprocunitcellslider.setMinimum(0)
        self.reciprocunitcellslider.setMaximum(100)
        self.reciprocunitcellslider.setOrientation(Qt.Horizontal)

        self.gridLayout_25.addWidget(self.reciprocunitcellslider, 0, 1, 1, 1)

        self.label_4 = QLabel(self.DrawReciprocUnitCellBox)
        self.label_4.setObjectName(u"label_4")
        self.label_4.setWordWrap(True)

        self.gridLayout_25.addWidget(self.label_4, 0, 2, 1, 1)


        self.gridLayout_21.addWidget(self.DrawReciprocUnitCellBox, 1, 0, 1, 1)

        self.DrawRealUnitCellBox = QGroupBox(self.groupBox_6)
        self.DrawRealUnitCellBox.setObjectName(u"DrawRealUnitCellBox")
        self.DrawRealUnitCellBox.setCheckable(True)
        self.DrawRealUnitCellBox.setChecked(False)
        self.gridLayout_26 = QGridLayout(self.DrawRealUnitCellBox)
        self.gridLayout_26.setSpacing(4)
        self.gridLayout_26.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_26.setObjectName(u"gridLayout_26")
        self.label_6 = QLabel(self.DrawRealUnitCellBox)
        self.label_6.setObjectName(u"label_6")
        self.label_6.setWordWrap(True)

        self.gridLayout_26.addWidget(self.label_6, 0, 0, 1, 1)

        self.unitcellslider = QSlider(self.DrawRealUnitCellBox)
        self.unitcellslider.setObjectName(u"unitcellslider")
        self.unitcellslider.setMinimum(0)
        self.unitcellslider.setMaximum(100)
        self.unitcellslider.setOrientation(Qt.Horizontal)

        self.gridLayout_26.addWidget(self.unitcellslider, 0, 1, 1, 1)

        self.label_5 = QLabel(self.DrawRealUnitCellBox)
        self.label_5.setObjectName(u"label_5")
        self.label_5.setWordWrap(True)

        self.gridLayout_26.addWidget(self.label_5, 0, 2, 1, 1)


        self.gridLayout_21.addWidget(self.DrawRealUnitCellBox, 2, 0, 1, 1)

        self.frame_4 = QFrame(self.groupBox_6)
        self.frame_4.setObjectName(u"frame_4")
        sizePolicy3.setHeightForWidth(self.frame_4.sizePolicy().hasHeightForWidth())
        self.frame_4.setSizePolicy(sizePolicy3)
        self.frame_4.setFrameShape(QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QFrame.Raised)
        self.gridLayout_14 = QGridLayout(self.frame_4)
        self.gridLayout_14.setSpacing(4)
        self.gridLayout_14.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_14.setObjectName(u"gridLayout_14")
        self.ResetViewBtn = QPushButton(self.frame_4)
        self.ResetViewBtn.setObjectName(u"ResetViewBtn")
        sizePolicy9.setHeightForWidth(self.ResetViewBtn.sizePolicy().hasHeightForWidth())
        self.ResetViewBtn.setSizePolicy(sizePolicy9)

        self.gridLayout_14.addWidget(self.ResetViewBtn, 0, 0, 1, 1)

        self.SaveImageBtn = QPushButton(self.frame_4)
        self.SaveImageBtn.setObjectName(u"SaveImageBtn")
        sizePolicy9.setHeightForWidth(self.SaveImageBtn.sizePolicy().hasHeightForWidth())
        self.SaveImageBtn.setSizePolicy(sizePolicy9)

        self.gridLayout_14.addWidget(self.SaveImageBtn, 0, 1, 1, 1)


        self.gridLayout_21.addWidget(self.frame_4, 0, 0, 1, 1)


        self.gridLayout_8.addWidget(self.groupBox_6, 0, 0, 1, 1)

        self.functionTabWidget.addTab(self.tab_5, "")
        self.splitter_2.addWidget(self.functionTabWidget)
        self.textInfo = QPlainTextEdit(self.splitter_2)
        self.textInfo.setObjectName(u"textInfo")
        sizePolicy21 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy21.setHorizontalStretch(0)
        sizePolicy21.setVerticalStretch(2)
        sizePolicy21.setHeightForWidth(self.textInfo.sizePolicy().hasHeightForWidth())
        self.textInfo.setSizePolicy(sizePolicy21)
        self.textInfo.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.textInfo.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.textInfo.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.splitter_2.addWidget(self.textInfo)

        self.gridLayout_4.addWidget(self.splitter_2, 0, 0, 1, 1)

        self.splitter.addWidget(self.widget_2)
        self.BrowserBox = QWebEngineView(self.splitter)
        self.BrowserBox.setObjectName(u"BrowserBox")
        sizePolicy15.setHeightForWidth(self.BrowserBox.sizePolicy().hasHeightForWidth())
        self.BrowserBox.setSizePolicy(sizePolicy15)
        self.splitter.addWidget(self.BrowserBox)

        self.gridLayout.addWidget(self.splitter, 0, 0, 1, 1)


        self.gridLayout_2.addWidget(self.widget, 0, 0, 1, 1)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1066, 22))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusBar = QStatusBar(MainWindow)
        self.statusBar.setObjectName(u"statusBar")
        MainWindow.setStatusBar(self.statusBar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menuFile.addAction(self.actionOpen_reflection_file)
        self.menuFile.addAction(self.actionSave_reflection_file)
        self.menuFile.addAction(self.actionSettings)
        self.menuFile.addAction(self.actiondebug)
        self.menuFile.addAction(self.actionExit)

        self.retranslateUi(MainWindow)

        self.functionTabWidget.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"MainWindow", None))
        self.actionOpen_reflection_file.setText(QCoreApplication.translate("MainWindow", u"Open reflection file...", None))
        self.actionSettings.setText(QCoreApplication.translate("MainWindow", u"Settings", None))
        self.actiondebug.setText(QCoreApplication.translate("MainWindow", u"Debug", None))
        self.actionExit.setText(QCoreApplication.translate("MainWindow", u"Exit", None))
        self.actionSave_reflection_file.setText(QCoreApplication.translate("MainWindow", u"Save reflection file", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"Display a data set with a double-click. Right-click table for more options.", None))
        self.ExpandReflsGroupBox.setTitle(QCoreApplication.translate("MainWindow", u"Expand reflections", None))
        self.expandP1checkbox.setText(QCoreApplication.translate("MainWindow", u"Expand to P1", None))
        self.expandAnomalouscheckbox.setText(QCoreApplication.translate("MainWindow", u"Show Friedel pairs", None))
        self.SpacegroupLabel.setText(QCoreApplication.translate("MainWindow", u"Space Subgroups", None))
        self.sysabsentcheckbox.setText(QCoreApplication.translate("MainWindow", u"Show Systematic Absences", None))
        self.missingcheckbox.setText(QCoreApplication.translate("MainWindow", u"Show Missing Reflections", None))
        self.onlymissingcheckbox.setText(QCoreApplication.translate("MainWindow", u"Only", None))
        self.functionTabWidget.setTabText(self.functionTabWidget.indexOf(self.tab), QCoreApplication.translate("MainWindow", u"Expansion", None))
        self.SliceReflectionsBox.setText(QCoreApplication.translate("MainWindow", u"Slice Reflections", None))
        self.fixedorientcheckbox.setText(QCoreApplication.translate("MainWindow", u"Fix orientation but allow zoom and translation", None))
        self.ClipPlaneChkGroupBox.setTitle(QCoreApplication.translate("MainWindow", u"with a clip plane oriented", None))
        self.clipNormalBtn.setText(QCoreApplication.translate("MainWindow", u"Perpendicular to vector below", None))
        self.clipParallelBtn.setText(QCoreApplication.translate("MainWindow", u"Parallel to vector below", None))
        self.groupBox_2.setTitle(QCoreApplication.translate("MainWindow", u"Vector Components (r1, r2, r3) are specified", None))
        self.recipvecBtn.setText(QCoreApplication.translate("MainWindow", u"in reciprocal unit cell coordinates", None))
        self.clipTNCSBtn.setText(QCoreApplication.translate("MainWindow", u"as tNCS vector", None))
        self.realspacevecBtn.setText(QCoreApplication.translate("MainWindow", u"in real space unit cell coordinates", None))
        self.rotavecangle_labeltxt.setText(QCoreApplication.translate("MainWindow", u"Angle rotated around vector ", None))
        self.label_21.setText(QCoreApplication.translate("MainWindow", u"r1:", None))
        self.label_22.setText(QCoreApplication.translate("MainWindow", u"r2:", None))
        self.label_23.setText(QCoreApplication.translate("MainWindow", u"r3:", None))
        self.label_15.setText(QCoreApplication.translate("MainWindow", u"Clip Plane Width:", None))
        self.label_14.setText(QCoreApplication.translate("MainWindow", u"Distance from origin:", None))
        self.showsliceGroupCheckbox.setTitle(QCoreApplication.translate("MainWindow", u"Explicitly with a plane at", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"equal to", None))
        self.functionTabWidget.setTabText(self.functionTabWidget.indexOf(self.tab_2), QCoreApplication.translate("MainWindow", u"Slicing", None))
        self.groupBox_3.setTitle(QCoreApplication.translate("MainWindow", u"Radii size of HKL spheres", None))
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Tick the &quot;User defined power scaling&quot; box for chosing a power factor. Untick this box for an automatically computed power factor. Automatic power factor will render reflections so that the reflection with the smallest data value is rendered about 10 times smaller than the reflection with the largest data value. <br/>Entering a power factor of 1 means the size of a sphere for each reflection is proportional to its data value. For a typical data set of intensities with values differing by some orders of magnitude only the strongest reflections would then be visible. Entering a power factor of 0 means all spheres of reflections will have the same size irrespective of their data values. This means even very weak reflections would then be visible. </p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.ManualPowerScalecheckbox.setToolTip("")
#endif // QT_CONFIG(tooltip)
        self.ManualPowerScalecheckbox.setText(QCoreApplication.translate("MainWindow", u"User defined power scaling", None))
        self.label_10.setText(QCoreApplication.translate("MainWindow", u"Power scale factor:", None))
#if QT_CONFIG(tooltip)
        self.power_scale_spinBox.setToolTip("")
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(tooltip)
        self.label_11.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>A linear scale factor of 1 means that if two reflectiions with the largest data values happen to be next to each other their spheres will touch but not overlap.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_11.setText(QCoreApplication.translate("MainWindow", u"Linear scale factor", None))
#if QT_CONFIG(tooltip)
        self.radii_scale_spinBox.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>A linear scale factor of 1 means that if two reflectiions with the largest data values happen to be next to each other their spheres will touch but not overlap.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.functionTabWidget.setTabText(self.functionTabWidget.indexOf(self.tab_3), QCoreApplication.translate("MainWindow", u"Sizing", None))
        self.groupBox_5.setTitle(QCoreApplication.translate("MainWindow", u"Bin according to", None))
#if QT_CONFIG(tooltip)
        self.BinDataComboBox.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Order the curently displayed data set in bins according to the values of the selected data set. Reflections in the currently displayed data set not matching any reflections in the selected data set will be put in a bin labelled with &quot;NaN&quot;.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_12.setText(QCoreApplication.translate("MainWindow", u"Number of bins:", None))
#if QT_CONFIG(tooltip)
        self.Nbins_spinBox.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Select the number of bins for ordering the currently displayed data set. More than about 10-15 bins will slow down the HKLviewer when working with large data sets.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.OpaqueAllCheckbox.setText(QCoreApplication.translate("MainWindow", u"Show/hide all data", None))
#if QT_CONFIG(tooltip)
        self.binstable.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Change the visibility of selected bins either by ticking or unticking the check boxes or by entering an opacity value between 0 and 1</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.functionTabWidget.setTabText(self.functionTabWidget.indexOf(self.tab_4), QCoreApplication.translate("MainWindow", u"Binning", None))
        self.groupBox_6.setTitle(QCoreApplication.translate("MainWindow", u"Geometry", None))
#if QT_CONFIG(tooltip)
        self.DrawReciprocUnitCellBox.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The reciprocal unit cell spanning a*,b*,c* from (0,0,0) to (1,1,1) is located in the centre of the displayed reflections. and is therefore quite small. Adjust the slider to scale its outline to become compatible with the sphere of displayed reflections.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.DrawReciprocUnitCellBox.setTitle(QCoreApplication.translate("MainWindow", u"Show Reciprocal Unit Cell", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"True to Scale of HKL grid", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"Close to radius of sphere of reflections", None))
#if QT_CONFIG(tooltip)
        self.DrawRealUnitCellBox.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The real space unit cell is often much larger than the sphere of displayed reflections because the reciprocal unit cell size is inversely proportional to the size of the real space unit cell. Adjust the slider to scale its outline to become compatible with the sphere of displayed reflections.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.DrawRealUnitCellBox.setTitle(QCoreApplication.translate("MainWindow", u"Show Real Space Unit Cell", None))
        self.label_6.setText(QCoreApplication.translate("MainWindow", u"True to Scale of HKL grid", None))
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"Close to radius of sphere of reflections", None))
#if QT_CONFIG(tooltip)
        self.ResetViewBtn.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>Place the displayed reflections in the centre of the view.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.ResetViewBtn.setText(QCoreApplication.translate("MainWindow", u"Reset View", None))
        self.SaveImageBtn.setText(QCoreApplication.translate("MainWindow", u"Save Image", None))
        self.functionTabWidget.setTabText(self.functionTabWidget.indexOf(self.tab_5), QCoreApplication.translate("MainWindow", u"Extras", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
    # retranslateUi

