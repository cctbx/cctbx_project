# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'HKLviewer2.ui'
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
        MainWindow.resize(790, 904)
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
        self.actionReset_View = QAction(MainWindow)
        self.actionReset_View.setObjectName(u"actionReset_View")
        self.actionSave_Current_Image = QAction(MainWindow)
        self.actionSave_Current_Image.setObjectName(u"actionSave_Current_Image")
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
        sizePolicy1 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy1)

        self.gridLayout_5.addWidget(self.label, 0, 0, 1, 1)

        self.millertable = HeaderDataTableWidget(self.widget_4)
        if (self.millertable.columnCount() < 10):
            self.millertable.setColumnCount(10)
        if (self.millertable.rowCount() < 1):
            self.millertable.setRowCount(1)
        self.millertable.setObjectName(u"millertable")
        sizePolicy2 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.millertable.sizePolicy().hasHeightForWidth())
        self.millertable.setSizePolicy(sizePolicy2)
        self.millertable.setVerticalScrollMode(QAbstractItemView.ScrollPerPixel)
        self.millertable.setHorizontalScrollMode(QAbstractItemView.ScrollPerPixel)
        self.millertable.setRowCount(1)
        self.millertable.setColumnCount(10)
        self.millertable.horizontalHeader().setMinimumSectionSize(5)

        self.gridLayout_5.addWidget(self.millertable, 1, 0, 1, 1)

        self.splitter_2.addWidget(self.widget_4)
        self.functionTabWidget = QTabWidget(self.splitter_2)
        self.functionTabWidget.setObjectName(u"functionTabWidget")
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
        self.showsliceGroupCheckbox = QGroupBox(self.tab_2)
        self.showsliceGroupCheckbox.setObjectName(u"showsliceGroupCheckbox")
        sizePolicy6 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(1)
        sizePolicy6.setHeightForWidth(self.showsliceGroupCheckbox.sizePolicy().hasHeightForWidth())
        self.showsliceGroupCheckbox.setSizePolicy(sizePolicy6)
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
        sizePolicy7 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy7.setHorizontalStretch(0)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy7)
        self.label_3.setAlignment(Qt.AlignLeading|Qt.AlignLeft|Qt.AlignTop)
        self.label_3.setIndent(1)

        self.gridLayout_6.addWidget(self.label_3, 0, 1, 1, 1)

        self.sliceindexspinBox = QSpinBox(self.showsliceGroupCheckbox)
        self.sliceindexspinBox.setObjectName(u"sliceindexspinBox")
        sizePolicy7.setHeightForWidth(self.sliceindexspinBox.sizePolicy().hasHeightForWidth())
        self.sliceindexspinBox.setSizePolicy(sizePolicy7)
        self.sliceindexspinBox.setAlignment(Qt.AlignLeading|Qt.AlignLeft|Qt.AlignTop)

        self.gridLayout_6.addWidget(self.sliceindexspinBox, 0, 2, 1, 1)


        self.gridLayout_9.addWidget(self.showsliceGroupCheckbox, 1, 0, 1, 2)

        self.ClipPlaneChkGroupBox = QGroupBox(self.tab_2)
        self.ClipPlaneChkGroupBox.setObjectName(u"ClipPlaneChkGroupBox")
        sizePolicy7.setHeightForWidth(self.ClipPlaneChkGroupBox.sizePolicy().hasHeightForWidth())
        self.ClipPlaneChkGroupBox.setSizePolicy(sizePolicy7)
        self.ClipPlaneChkGroupBox.setCheckable(True)
        self.gridLayout_7 = QGridLayout(self.ClipPlaneChkGroupBox)
        self.gridLayout_7.setSpacing(4)
        self.gridLayout_7.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.frame_5 = QFrame(self.ClipPlaneChkGroupBox)
        self.frame_5.setObjectName(u"frame_5")
        self.frame_5.setFrameShape(QFrame.StyledPanel)
        self.frame_5.setFrameShadow(QFrame.Raised)
        self.gridLayout_17 = QGridLayout(self.frame_5)
        self.gridLayout_17.setSpacing(4)
        self.gridLayout_17.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_17.setObjectName(u"gridLayout_17")
        self.gridLayout_17.setContentsMargins(0, 0, 0, 0)
        self.frame_12 = QFrame(self.frame_5)
        self.frame_12.setObjectName(u"frame_12")
        sizePolicy8 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        sizePolicy8.setHorizontalStretch(0)
        sizePolicy8.setVerticalStretch(0)
        sizePolicy8.setHeightForWidth(self.frame_12.sizePolicy().hasHeightForWidth())
        self.frame_12.setSizePolicy(sizePolicy8)
        self.frame_12.setFrameShape(QFrame.StyledPanel)
        self.frame_12.setFrameShadow(QFrame.Raised)
        self.gridLayout_30 = QGridLayout(self.frame_12)
        self.gridLayout_30.setSpacing(4)
        self.gridLayout_30.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_30.setObjectName(u"gridLayout_30")
        self.gridLayout_30.setContentsMargins(0, 0, 0, 0)
        self.label_24 = QLabel(self.frame_12)
        self.label_24.setObjectName(u"label_24")
        sizePolicy9 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Preferred)
        sizePolicy9.setHorizontalStretch(0)
        sizePolicy9.setVerticalStretch(0)
        sizePolicy9.setHeightForWidth(self.label_24.sizePolicy().hasHeightForWidth())
        self.label_24.setSizePolicy(sizePolicy9)

        self.gridLayout_30.addWidget(self.label_24, 0, 0, 1, 1)

        self.hvec_spinBox = QDoubleSpinBox(self.frame_12)
        self.hvec_spinBox.setObjectName(u"hvec_spinBox")
        sizePolicy9.setHeightForWidth(self.hvec_spinBox.sizePolicy().hasHeightForWidth())
        self.hvec_spinBox.setSizePolicy(sizePolicy9)

        self.gridLayout_30.addWidget(self.hvec_spinBox, 0, 1, 1, 1)


        self.gridLayout_17.addWidget(self.frame_12, 0, 0, 1, 1)

        self.frame_13 = QFrame(self.frame_5)
        self.frame_13.setObjectName(u"frame_13")
        sizePolicy8.setHeightForWidth(self.frame_13.sizePolicy().hasHeightForWidth())
        self.frame_13.setSizePolicy(sizePolicy8)
        self.frame_13.setFrameShape(QFrame.StyledPanel)
        self.frame_13.setFrameShadow(QFrame.Raised)
        self.gridLayout_31 = QGridLayout(self.frame_13)
        self.gridLayout_31.setSpacing(4)
        self.gridLayout_31.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_31.setObjectName(u"gridLayout_31")
        self.gridLayout_31.setContentsMargins(0, 0, 0, 0)
        self.label_25 = QLabel(self.frame_13)
        self.label_25.setObjectName(u"label_25")
        sizePolicy9.setHeightForWidth(self.label_25.sizePolicy().hasHeightForWidth())
        self.label_25.setSizePolicy(sizePolicy9)

        self.gridLayout_31.addWidget(self.label_25, 0, 0, 1, 1)

        self.kvec_spinBox = QDoubleSpinBox(self.frame_13)
        self.kvec_spinBox.setObjectName(u"kvec_spinBox")
        sizePolicy9.setHeightForWidth(self.kvec_spinBox.sizePolicy().hasHeightForWidth())
        self.kvec_spinBox.setSizePolicy(sizePolicy9)

        self.gridLayout_31.addWidget(self.kvec_spinBox, 0, 1, 1, 1)


        self.gridLayout_17.addWidget(self.frame_13, 0, 1, 1, 1)

        self.frame_14 = QFrame(self.frame_5)
        self.frame_14.setObjectName(u"frame_14")
        sizePolicy8.setHeightForWidth(self.frame_14.sizePolicy().hasHeightForWidth())
        self.frame_14.setSizePolicy(sizePolicy8)
        self.frame_14.setFrameShape(QFrame.StyledPanel)
        self.frame_14.setFrameShadow(QFrame.Raised)
        self.gridLayout_32 = QGridLayout(self.frame_14)
        self.gridLayout_32.setSpacing(4)
        self.gridLayout_32.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_32.setObjectName(u"gridLayout_32")
        self.gridLayout_32.setContentsMargins(0, 0, 0, 0)
        self.label_26 = QLabel(self.frame_14)
        self.label_26.setObjectName(u"label_26")
        sizePolicy9.setHeightForWidth(self.label_26.sizePolicy().hasHeightForWidth())
        self.label_26.setSizePolicy(sizePolicy9)

        self.gridLayout_32.addWidget(self.label_26, 0, 0, 1, 1)

        self.lvec_spinBox = QDoubleSpinBox(self.frame_14)
        self.lvec_spinBox.setObjectName(u"lvec_spinBox")
        sizePolicy9.setHeightForWidth(self.lvec_spinBox.sizePolicy().hasHeightForWidth())
        self.lvec_spinBox.setSizePolicy(sizePolicy9)

        self.gridLayout_32.addWidget(self.lvec_spinBox, 0, 1, 1, 1)


        self.gridLayout_17.addWidget(self.frame_14, 0, 2, 1, 1)


        self.gridLayout_7.addWidget(self.frame_5, 0, 0, 1, 1)

        self.fixedorientcheckbox = QCheckBox(self.ClipPlaneChkGroupBox)
        self.fixedorientcheckbox.setObjectName(u"fixedorientcheckbox")
        self.fixedorientcheckbox.setChecked(True)

        self.gridLayout_7.addWidget(self.fixedorientcheckbox, 0, 1, 1, 1)

        self.frame_2 = QFrame(self.ClipPlaneChkGroupBox)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QFrame.Raised)
        self.gridLayout_33 = QGridLayout(self.frame_2)
        self.gridLayout_33.setSpacing(4)
        self.gridLayout_33.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_33.setObjectName(u"gridLayout_33")
        self.gridLayout_33.setContentsMargins(0, 0, 0, 0)
        self.label_16 = QLabel(self.frame_2)
        self.label_16.setObjectName(u"label_16")
        sizePolicy8.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy8)
        self.label_16.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.gridLayout_33.addWidget(self.label_16, 0, 3, 1, 1)

        self.hkldist_spinBox = QDoubleSpinBox(self.frame_2)
        self.hkldist_spinBox.setObjectName(u"hkldist_spinBox")
        sizePolicy8.setHeightForWidth(self.hkldist_spinBox.sizePolicy().hasHeightForWidth())
        self.hkldist_spinBox.setSizePolicy(sizePolicy8)

        self.gridLayout_33.addWidget(self.hkldist_spinBox, 0, 1, 1, 1)

        self.clipwidth_spinBox = QDoubleSpinBox(self.frame_2)
        self.clipwidth_spinBox.setObjectName(u"clipwidth_spinBox")
        sizePolicy8.setHeightForWidth(self.clipwidth_spinBox.sizePolicy().hasHeightForWidth())
        self.clipwidth_spinBox.setSizePolicy(sizePolicy8)

        self.gridLayout_33.addWidget(self.clipwidth_spinBox, 0, 4, 1, 1)

        self.label_17 = QLabel(self.frame_2)
        self.label_17.setObjectName(u"label_17")

        self.gridLayout_33.addWidget(self.label_17, 0, 0, 1, 1)


        self.gridLayout_7.addWidget(self.frame_2, 1, 0, 1, 2)


        self.gridLayout_9.addWidget(self.ClipPlaneChkGroupBox, 0, 0, 1, 2)

        self.functionTabWidget.addTab(self.tab_2, "")
        self.tab_3 = QWidget()
        self.tab_3.setObjectName(u"tab_3")
        self.gridLayout_27 = QGridLayout(self.tab_3)
        self.gridLayout_27.setSpacing(4)
        self.gridLayout_27.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_27.setObjectName(u"gridLayout_27")
        self.groupBox_3 = QGroupBox(self.tab_3)
        self.groupBox_3.setObjectName(u"groupBox_3")
        sizePolicy10 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy10.setHorizontalStretch(1)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(self.groupBox_3.sizePolicy().hasHeightForWidth())
        self.groupBox_3.setSizePolicy(sizePolicy10)
        self.groupBox_3.setCheckable(False)
        self.gridLayout_10 = QGridLayout(self.groupBox_3)
        self.gridLayout_10.setSpacing(4)
        self.gridLayout_10.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.label_7 = QLabel(self.groupBox_3)
        self.label_7.setObjectName(u"label_7")
        sizePolicy7.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy7)
        self.label_7.setTextFormat(Qt.RichText)
        self.label_7.setScaledContents(False)
        self.label_7.setAlignment(Qt.AlignLeading|Qt.AlignLeft|Qt.AlignVCenter)
        self.label_7.setWordWrap(True)

        self.gridLayout_10.addWidget(self.label_7, 0, 0, 1, 4)

        self.ManualPowerScalecheckbox = QCheckBox(self.groupBox_3)
        self.ManualPowerScalecheckbox.setObjectName(u"ManualPowerScalecheckbox")
        sizePolicy11 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.MinimumExpanding)
        sizePolicy11.setHorizontalStretch(0)
        sizePolicy11.setVerticalStretch(0)
        sizePolicy11.setHeightForWidth(self.ManualPowerScalecheckbox.sizePolicy().hasHeightForWidth())
        self.ManualPowerScalecheckbox.setSizePolicy(sizePolicy11)

        self.gridLayout_10.addWidget(self.ManualPowerScalecheckbox, 1, 0, 1, 2)

        self.label_10 = QLabel(self.groupBox_3)
        self.label_10.setObjectName(u"label_10")
        sizePolicy11.setHeightForWidth(self.label_10.sizePolicy().hasHeightForWidth())
        self.label_10.setSizePolicy(sizePolicy11)
        self.label_10.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.label_10.setIndent(1)

        self.gridLayout_10.addWidget(self.label_10, 1, 2, 1, 1)

        self.power_scale_spinBox = QDoubleSpinBox(self.groupBox_3)
        self.power_scale_spinBox.setObjectName(u"power_scale_spinBox")
        sizePolicy12 = QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        sizePolicy12.setHorizontalStretch(0)
        sizePolicy12.setVerticalStretch(0)
        sizePolicy12.setHeightForWidth(self.power_scale_spinBox.sizePolicy().hasHeightForWidth())
        self.power_scale_spinBox.setSizePolicy(sizePolicy12)
        self.power_scale_spinBox.setMaximum(1.000000000000000)
        self.power_scale_spinBox.setSingleStep(0.050000000000000)
        self.power_scale_spinBox.setValue(0.350000000000000)

        self.gridLayout_10.addWidget(self.power_scale_spinBox, 1, 3, 1, 1)

        self.label_11 = QLabel(self.groupBox_3)
        self.label_11.setObjectName(u"label_11")
        sizePolicy11.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy11)

        self.gridLayout_10.addWidget(self.label_11, 2, 0, 1, 1)

        self.radii_scale_spinBox = QDoubleSpinBox(self.groupBox_3)
        self.radii_scale_spinBox.setObjectName(u"radii_scale_spinBox")
        sizePolicy12.setHeightForWidth(self.radii_scale_spinBox.sizePolicy().hasHeightForWidth())
        self.radii_scale_spinBox.setSizePolicy(sizePolicy12)
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
        sizePolicy13 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy13.setHorizontalStretch(2)
        sizePolicy13.setVerticalStretch(0)
        sizePolicy13.setHeightForWidth(self.BinDataComboBox.sizePolicy().hasHeightForWidth())
        self.BinDataComboBox.setSizePolicy(sizePolicy13)

        self.gridLayout_11.addWidget(self.BinDataComboBox, 0, 0, 1, 1)

        self.label_12 = QLabel(self.groupBox_5)
        self.label_12.setObjectName(u"label_12")

        self.gridLayout_11.addWidget(self.label_12, 0, 1, 1, 1)

        self.Nbins_spinBox = QSpinBox(self.groupBox_5)
        self.Nbins_spinBox.setObjectName(u"Nbins_spinBox")
        sizePolicy9.setHeightForWidth(self.Nbins_spinBox.sizePolicy().hasHeightForWidth())
        self.Nbins_spinBox.setSizePolicy(sizePolicy9)
        self.Nbins_spinBox.setMinimum(1)
        self.Nbins_spinBox.setMaximum(20)

        self.gridLayout_11.addWidget(self.Nbins_spinBox, 0, 2, 1, 1)


        self.gridLayout_13.addWidget(self.groupBox_5, 0, 0, 1, 1)

        self.widget_6 = QWidget(self.tab_4)
        self.widget_6.setObjectName(u"widget_6")
        sizePolicy14 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy14.setHorizontalStretch(0)
        sizePolicy14.setVerticalStretch(1)
        sizePolicy14.setHeightForWidth(self.widget_6.sizePolicy().hasHeightForWidth())
        self.widget_6.setSizePolicy(sizePolicy14)
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
        sizePolicy15 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy15.setHorizontalStretch(0)
        sizePolicy15.setVerticalStretch(1)
        sizePolicy15.setHeightForWidth(self.binstable.sizePolicy().hasHeightForWidth())
        self.binstable.setSizePolicy(sizePolicy15)
        self.binstable.setRowCount(5)
        self.binstable.setColumnCount(4)

        self.gridLayout_12.addWidget(self.binstable, 1, 0, 1, 1)


        self.gridLayout_13.addWidget(self.widget_6, 1, 0, 1, 1)

        self.functionTabWidget.addTab(self.tab_4, "")
        self.tab_5 = QWidget()
        self.tab_5.setObjectName(u"tab_5")
        self.gridLayout_14 = QGridLayout(self.tab_5)
        self.gridLayout_14.setSpacing(4)
        self.gridLayout_14.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_14.setObjectName(u"gridLayout_14")
        self.groupBox_2 = QGroupBox(self.tab_5)
        self.groupBox_2.setObjectName(u"groupBox_2")
        self.groupBox_2.setMaximumSize(QSize(16777215, 16777215))
        self.groupBox_2.setCheckable(True)
        self.gridLayout_16 = QGridLayout(self.groupBox_2)
        self.gridLayout_16.setSpacing(4)
        self.gridLayout_16.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_16.setObjectName(u"gridLayout_16")
        self.gridLayout_16.setContentsMargins(3, -1, 3, -1)
        self.frame_7 = QFrame(self.groupBox_2)
        self.frame_7.setObjectName(u"frame_7")
        sizePolicy7.setHeightForWidth(self.frame_7.sizePolicy().hasHeightForWidth())
        self.frame_7.setSizePolicy(sizePolicy7)
        self.frame_7.setFrameShape(QFrame.StyledPanel)
        self.frame_7.setFrameShadow(QFrame.Raised)
        self.gridLayout_25 = QGridLayout(self.frame_7)
        self.gridLayout_25.setSpacing(4)
        self.gridLayout_25.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_25.setObjectName(u"gridLayout_25")
        self.gridLayout_25.setContentsMargins(0, 0, 0, 0)
        self.frame_15 = QFrame(self.frame_7)
        self.frame_15.setObjectName(u"frame_15")
        sizePolicy8.setHeightForWidth(self.frame_15.sizePolicy().hasHeightForWidth())
        self.frame_15.setSizePolicy(sizePolicy8)
        self.frame_15.setFrameShape(QFrame.StyledPanel)
        self.frame_15.setFrameShadow(QFrame.Raised)
        self.gridLayout_34 = QGridLayout(self.frame_15)
        self.gridLayout_34.setSpacing(4)
        self.gridLayout_34.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_34.setObjectName(u"gridLayout_34")
        self.gridLayout_34.setContentsMargins(0, 0, 0, 0)
        self.label_27 = QLabel(self.frame_15)
        self.label_27.setObjectName(u"label_27")
        sizePolicy9.setHeightForWidth(self.label_27.sizePolicy().hasHeightForWidth())
        self.label_27.setSizePolicy(sizePolicy9)

        self.gridLayout_34.addWidget(self.label_27, 0, 0, 1, 1)

        self.R1_spinBox = QDoubleSpinBox(self.frame_15)
        self.R1_spinBox.setObjectName(u"R1_spinBox")
        sizePolicy9.setHeightForWidth(self.R1_spinBox.sizePolicy().hasHeightForWidth())
        self.R1_spinBox.setSizePolicy(sizePolicy9)

        self.gridLayout_34.addWidget(self.R1_spinBox, 0, 1, 1, 1)


        self.gridLayout_25.addWidget(self.frame_15, 0, 0, 1, 1)

        self.frame_16 = QFrame(self.frame_7)
        self.frame_16.setObjectName(u"frame_16")
        sizePolicy8.setHeightForWidth(self.frame_16.sizePolicy().hasHeightForWidth())
        self.frame_16.setSizePolicy(sizePolicy8)
        self.frame_16.setFrameShape(QFrame.StyledPanel)
        self.frame_16.setFrameShadow(QFrame.Raised)
        self.gridLayout_35 = QGridLayout(self.frame_16)
        self.gridLayout_35.setSpacing(4)
        self.gridLayout_35.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_35.setObjectName(u"gridLayout_35")
        self.gridLayout_35.setContentsMargins(0, 0, 0, 0)
        self.label_28 = QLabel(self.frame_16)
        self.label_28.setObjectName(u"label_28")
        sizePolicy9.setHeightForWidth(self.label_28.sizePolicy().hasHeightForWidth())
        self.label_28.setSizePolicy(sizePolicy9)

        self.gridLayout_35.addWidget(self.label_28, 0, 0, 1, 1)

        self.R2_spinBox = QDoubleSpinBox(self.frame_16)
        self.R2_spinBox.setObjectName(u"R2_spinBox")
        sizePolicy9.setHeightForWidth(self.R2_spinBox.sizePolicy().hasHeightForWidth())
        self.R2_spinBox.setSizePolicy(sizePolicy9)

        self.gridLayout_35.addWidget(self.R2_spinBox, 0, 1, 1, 1)


        self.gridLayout_25.addWidget(self.frame_16, 0, 1, 1, 1)

        self.frame_17 = QFrame(self.frame_7)
        self.frame_17.setObjectName(u"frame_17")
        sizePolicy8.setHeightForWidth(self.frame_17.sizePolicy().hasHeightForWidth())
        self.frame_17.setSizePolicy(sizePolicy8)
        self.frame_17.setFrameShape(QFrame.StyledPanel)
        self.frame_17.setFrameShadow(QFrame.Raised)
        self.gridLayout_36 = QGridLayout(self.frame_17)
        self.gridLayout_36.setSpacing(4)
        self.gridLayout_36.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_36.setObjectName(u"gridLayout_36")
        self.gridLayout_36.setContentsMargins(0, 0, 0, 0)
        self.R3_spinBox = QDoubleSpinBox(self.frame_17)
        self.R3_spinBox.setObjectName(u"R3_spinBox")
        sizePolicy9.setHeightForWidth(self.R3_spinBox.sizePolicy().hasHeightForWidth())
        self.R3_spinBox.setSizePolicy(sizePolicy9)

        self.gridLayout_36.addWidget(self.R3_spinBox, 0, 1, 1, 1)

        self.label_29 = QLabel(self.frame_17)
        self.label_29.setObjectName(u"label_29")
        sizePolicy9.setHeightForWidth(self.label_29.sizePolicy().hasHeightForWidth())
        self.label_29.setSizePolicy(sizePolicy9)

        self.gridLayout_36.addWidget(self.label_29, 0, 0, 1, 1)


        self.gridLayout_25.addWidget(self.frame_17, 0, 2, 1, 1)


        self.gridLayout_16.addWidget(self.frame_7, 0, 0, 1, 1)

        self.TwinAxisComboBox = QComboBox(self.groupBox_2)
        self.TwinAxisComboBox.setObjectName(u"TwinAxisComboBox")

        self.gridLayout_16.addWidget(self.TwinAxisComboBox, 0, 1, 1, 1)

        self.frame_4 = QFrame(self.groupBox_2)
        self.frame_4.setObjectName(u"frame_4")
        sizePolicy7.setHeightForWidth(self.frame_4.sizePolicy().hasHeightForWidth())
        self.frame_4.setSizePolicy(sizePolicy7)
        self.frame_4.setFrameShape(QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QFrame.Raised)
        self.gridLayout_8 = QGridLayout(self.frame_4)
        self.gridLayout_8.setSpacing(4)
        self.gridLayout_8.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_8.setObjectName(u"gridLayout_8")
        self.gridLayout_8.setContentsMargins(3, 3, 3, 3)
        self.recipvecBtn = QRadioButton(self.frame_4)
        self.recipvecBtn.setObjectName(u"recipvecBtn")
        self.recipvecBtn.setChecked(True)

        self.gridLayout_8.addWidget(self.recipvecBtn, 0, 0, 1, 1)

        self.groupBox = QGroupBox(self.frame_4)
        self.groupBox.setObjectName(u"groupBox")
        self.groupBox.setCheckable(True)
        self.verticalLayout = QVBoxLayout(self.groupBox)
        self.verticalLayout.setSpacing(4)
        self.verticalLayout.setContentsMargins(3, 3, 3, 3)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.verticalLayout.setContentsMargins(-1, 3, -1, 3)
        self.clipParallelBtn = QRadioButton(self.groupBox)
        self.clipParallelBtn.setObjectName(u"clipParallelBtn")

        self.verticalLayout.addWidget(self.clipParallelBtn)

        self.clipNormalBtn = QRadioButton(self.groupBox)
        self.clipNormalBtn.setObjectName(u"clipNormalBtn")

        self.verticalLayout.addWidget(self.clipNormalBtn)


        self.gridLayout_8.addWidget(self.groupBox, 0, 1, 3, 1)

        self.realspacevecBtn = QRadioButton(self.frame_4)
        self.realspacevecBtn.setObjectName(u"realspacevecBtn")

        self.gridLayout_8.addWidget(self.realspacevecBtn, 1, 0, 1, 1)

        self.clipTNCSBtn = QRadioButton(self.frame_4)
        self.clipTNCSBtn.setObjectName(u"clipTNCSBtn")

        self.gridLayout_8.addWidget(self.clipTNCSBtn, 2, 0, 1, 1)


        self.gridLayout_16.addWidget(self.frame_4, 1, 0, 1, 2)

        self.frame_6 = QFrame(self.groupBox_2)
        self.frame_6.setObjectName(u"frame_6")
        sizePolicy14.setHeightForWidth(self.frame_6.sizePolicy().hasHeightForWidth())
        self.frame_6.setSizePolicy(sizePolicy14)
        self.frame_6.setFrameShape(QFrame.StyledPanel)
        self.frame_6.setFrameShadow(QFrame.Raised)
        self.gridLayout_18 = QGridLayout(self.frame_6)
        self.gridLayout_18.setSpacing(4)
        self.gridLayout_18.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_18.setObjectName(u"gridLayout_18")
        self.gridLayout_18.setContentsMargins(0, 0, 0, 0)
        self.rotavecangle_slider = QSlider(self.frame_6)
        self.rotavecangle_slider.setObjectName(u"rotavecangle_slider")
        sizePolicy16 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy16.setHorizontalStretch(1)
        sizePolicy16.setVerticalStretch(0)
        sizePolicy16.setHeightForWidth(self.rotavecangle_slider.sizePolicy().hasHeightForWidth())
        self.rotavecangle_slider.setSizePolicy(sizePolicy16)
        self.rotavecangle_slider.setMaximum(359)
        self.rotavecangle_slider.setOrientation(Qt.Horizontal)
        self.rotavecangle_slider.setTickPosition(QSlider.TicksAbove)
        self.rotavecangle_slider.setTickInterval(18)

        self.gridLayout_18.addWidget(self.rotavecangle_slider, 0, 1, 1, 1)

        self.rotavecangle_labeltxt = QLabel(self.frame_6)
        self.rotavecangle_labeltxt.setObjectName(u"rotavecangle_labeltxt")
        sizePolicy17 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        sizePolicy17.setHorizontalStretch(1)
        sizePolicy17.setVerticalStretch(0)
        sizePolicy17.setHeightForWidth(self.rotavecangle_labeltxt.sizePolicy().hasHeightForWidth())
        self.rotavecangle_labeltxt.setSizePolicy(sizePolicy17)
        self.rotavecangle_labeltxt.setWordWrap(True)

        self.gridLayout_18.addWidget(self.rotavecangle_labeltxt, 0, 0, 1, 1)


        self.gridLayout_16.addWidget(self.frame_6, 2, 0, 1, 2)

        self.AnimaRotCheckBox = QCheckBox(self.groupBox_2)
        self.AnimaRotCheckBox.setObjectName(u"AnimaRotCheckBox")
        sizePolicy18 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        sizePolicy18.setHorizontalStretch(0)
        sizePolicy18.setVerticalStretch(1)
        sizePolicy18.setHeightForWidth(self.AnimaRotCheckBox.sizePolicy().hasHeightForWidth())
        self.AnimaRotCheckBox.setSizePolicy(sizePolicy18)

        self.gridLayout_16.addWidget(self.AnimaRotCheckBox, 3, 0, 1, 1)


        self.gridLayout_14.addWidget(self.groupBox_2, 0, 0, 1, 1)

        self.DrawReciprocUnitCellBox = QGroupBox(self.tab_5)
        self.DrawReciprocUnitCellBox.setObjectName(u"DrawReciprocUnitCellBox")
        self.DrawReciprocUnitCellBox.setCheckable(True)
        self.DrawReciprocUnitCellBox.setChecked(False)
        self.gridLayout_15 = QGridLayout(self.DrawReciprocUnitCellBox)
        self.gridLayout_15.setSpacing(4)
        self.gridLayout_15.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_15.setObjectName(u"gridLayout_15")
        self.gridLayout_15.setContentsMargins(-1, 3, -1, 3)
        self.label_2 = QLabel(self.DrawReciprocUnitCellBox)
        self.label_2.setObjectName(u"label_2")
        sizePolicy7.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy7)
        self.label_2.setWordWrap(True)

        self.gridLayout_15.addWidget(self.label_2, 0, 0, 1, 1)

        self.reciprocunitcellslider = QSlider(self.DrawReciprocUnitCellBox)
        self.reciprocunitcellslider.setObjectName(u"reciprocunitcellslider")
        sizePolicy19 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy19.setHorizontalStretch(0)
        sizePolicy19.setVerticalStretch(0)
        sizePolicy19.setHeightForWidth(self.reciprocunitcellslider.sizePolicy().hasHeightForWidth())
        self.reciprocunitcellslider.setSizePolicy(sizePolicy19)
        self.reciprocunitcellslider.setMinimum(0)
        self.reciprocunitcellslider.setMaximum(100)
        self.reciprocunitcellslider.setOrientation(Qt.Horizontal)

        self.gridLayout_15.addWidget(self.reciprocunitcellslider, 0, 1, 1, 1)

        self.label_4 = QLabel(self.DrawReciprocUnitCellBox)
        self.label_4.setObjectName(u"label_4")
        self.label_4.setWordWrap(True)

        self.gridLayout_15.addWidget(self.label_4, 0, 2, 1, 1)


        self.gridLayout_14.addWidget(self.DrawReciprocUnitCellBox, 1, 0, 1, 1)

        self.DrawRealUnitCellBox = QGroupBox(self.tab_5)
        self.DrawRealUnitCellBox.setObjectName(u"DrawRealUnitCellBox")
        self.DrawRealUnitCellBox.setCheckable(True)
        self.DrawRealUnitCellBox.setChecked(False)
        self.gridLayout_26 = QGridLayout(self.DrawRealUnitCellBox)
        self.gridLayout_26.setSpacing(4)
        self.gridLayout_26.setContentsMargins(3, 3, 3, 3)
        self.gridLayout_26.setObjectName(u"gridLayout_26")
        self.gridLayout_26.setContentsMargins(-1, 3, -1, 3)
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


        self.gridLayout_14.addWidget(self.DrawRealUnitCellBox, 2, 0, 1, 1)

        self.functionTabWidget.addTab(self.tab_5, "")
        self.splitter_2.addWidget(self.functionTabWidget)
        self.textInfo = QPlainTextEdit(self.splitter_2)
        self.textInfo.setObjectName(u"textInfo")
        sizePolicy15.setHeightForWidth(self.textInfo.sizePolicy().hasHeightForWidth())
        self.textInfo.setSizePolicy(sizePolicy15)
        self.textInfo.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.textInfo.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.textInfo.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.splitter_2.addWidget(self.textInfo)

        self.gridLayout_4.addWidget(self.splitter_2, 0, 0, 1, 1)

        self.splitter.addWidget(self.widget_2)
        self.BrowserBox = QWebEngineView(self.splitter)
        self.BrowserBox.setObjectName(u"BrowserBox")
        sizePolicy10.setHeightForWidth(self.BrowserBox.sizePolicy().hasHeightForWidth())
        self.BrowserBox.setSizePolicy(sizePolicy10)
        self.splitter.addWidget(self.BrowserBox)

        self.gridLayout.addWidget(self.splitter, 0, 0, 1, 1)


        self.gridLayout_2.addWidget(self.widget, 0, 0, 1, 1)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 790, 22))
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
        self.menuFile.addAction(self.actionReset_View)
        self.menuFile.addAction(self.actionSave_Current_Image)
        self.menuFile.addAction(self.actionExit)

        self.retranslateUi(MainWindow)

        self.functionTabWidget.setCurrentIndex(4)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"MainWindow", None))
        self.actionOpen_reflection_file.setText(QCoreApplication.translate("MainWindow", u"Open reflection file...", None))
        self.actionSettings.setText(QCoreApplication.translate("MainWindow", u"Settings", None))
        self.actiondebug.setText(QCoreApplication.translate("MainWindow", u"Debug", None))
        self.actionExit.setText(QCoreApplication.translate("MainWindow", u"Exit", None))
        self.actionSave_reflection_file.setText(QCoreApplication.translate("MainWindow", u"Save reflection file", None))
        self.actionReset_View.setText(QCoreApplication.translate("MainWindow", u"Reset View", None))
        self.actionSave_Current_Image.setText(QCoreApplication.translate("MainWindow", u"Save Current Image", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"Display a data set with a double-click. Right-click table for more options.", None))
        self.ExpandReflsGroupBox.setTitle(QCoreApplication.translate("MainWindow", u"Expand reflections", None))
        self.expandP1checkbox.setText(QCoreApplication.translate("MainWindow", u"Expand to P1", None))
        self.expandAnomalouscheckbox.setText(QCoreApplication.translate("MainWindow", u"Show Friedel pairs", None))
        self.SpacegroupLabel.setText(QCoreApplication.translate("MainWindow", u"Space Subgroups", None))
        self.sysabsentcheckbox.setText(QCoreApplication.translate("MainWindow", u"Show Systematic Absences", None))
        self.missingcheckbox.setText(QCoreApplication.translate("MainWindow", u"Show Missing Reflections", None))
        self.onlymissingcheckbox.setText(QCoreApplication.translate("MainWindow", u"Only", None))
        self.functionTabWidget.setTabText(self.functionTabWidget.indexOf(self.tab), QCoreApplication.translate("MainWindow", u"Expansion", None))
        self.showsliceGroupCheckbox.setTitle(QCoreApplication.translate("MainWindow", u"Explicitly with a plane at", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"equal to", None))
        self.ClipPlaneChkGroupBox.setTitle(QCoreApplication.translate("MainWindow", u"Slice with a clip plane normal to vector", None))
        self.label_24.setText(QCoreApplication.translate("MainWindow", u"h:", None))
        self.label_25.setText(QCoreApplication.translate("MainWindow", u"k:", None))
        self.label_26.setText(QCoreApplication.translate("MainWindow", u"l:", None))
        self.fixedorientcheckbox.setText(QCoreApplication.translate("MainWindow", u"Fix orientation", None))
        self.label_16.setText(QCoreApplication.translate("MainWindow", u"Clip Plane Width:", None))
        self.label_17.setText(QCoreApplication.translate("MainWindow", u"Distance from origin:", None))
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
        self.groupBox_2.setTitle(QCoreApplication.translate("MainWindow", u"Add Vector from Origin", None))
        self.label_27.setText(QCoreApplication.translate("MainWindow", u"r1:", None))
        self.label_28.setText(QCoreApplication.translate("MainWindow", u"r2:", None))
        self.label_29.setText(QCoreApplication.translate("MainWindow", u"r3:", None))
        self.TwinAxisComboBox.setCurrentText("")
        self.recipvecBtn.setText(QCoreApplication.translate("MainWindow", u"in Reciprocal Space Coordinates", None))
        self.groupBox.setTitle(QCoreApplication.translate("MainWindow", u"Align Vector to Screen", None))
        self.clipParallelBtn.setText(QCoreApplication.translate("MainWindow", u"Parallel", None))
        self.clipNormalBtn.setText(QCoreApplication.translate("MainWindow", u"Perpendicular", None))
        self.realspacevecBtn.setText(QCoreApplication.translate("MainWindow", u"in Real Space Coordinates", None))
        self.clipTNCSBtn.setText(QCoreApplication.translate("MainWindow", u"TNCS Vector in File Header", None))
        self.rotavecangle_labeltxt.setText(QCoreApplication.translate("MainWindow", u"Rotate Reflections around Vector with Angle", None))
        self.AnimaRotCheckBox.setText(QCoreApplication.translate("MainWindow", u"Rotate continuously", None))
#if QT_CONFIG(tooltip)
        self.DrawReciprocUnitCellBox.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The reciprocal unit cell spanning a*,b*,c* from (0,0,0) to (1,1,1) is located in the centre of the displayed reflections. and is therefore quite small. Adjust the slider to scale its outline to become compatible with the sphere of displayed reflections.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.DrawReciprocUnitCellBox.setTitle(QCoreApplication.translate("MainWindow", u"Show Reciprocal Unit Cell", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"True to Scale", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"Close to radius of sphere of reflections", None))
#if QT_CONFIG(tooltip)
        self.DrawRealUnitCellBox.setToolTip(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>The real space unit cell is often much larger than the sphere of displayed reflections because the reciprocal unit cell size is inversely proportional to the size of the real space unit cell. Adjust the slider to scale its outline to become compatible with the sphere of displayed reflections.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.DrawRealUnitCellBox.setTitle(QCoreApplication.translate("MainWindow", u"Show Real Space Unit Cell", None))
        self.label_6.setText(QCoreApplication.translate("MainWindow", u"True to Scale", None))
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"Close to radius of sphere of reflections", None))
        self.functionTabWidget.setTabText(self.functionTabWidget.indexOf(self.tab_5), QCoreApplication.translate("MainWindow", u"Vector", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
    # retranslateUi

