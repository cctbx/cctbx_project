from __future__ import absolute_import, division, print_function

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

from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import ( QApplication, QCheckBox, QComboBox,
        QDial, QDialog, QFileDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QDoubleSpinBox, QSpinBox, QStyleFactory, QTableWidget,
        QTableWidgetItem, QTabWidget, QTextEdit, QVBoxLayout, QWidget )

import sys, zmq, subprocess, time

class NGL_HKLViewer(QDialog):
  def __init__(self, parent=None):
    super(NGL_HKLViewer, self).__init__(parent)
    self.context = None

    self.originalPalette = QApplication.palette()

    self.openFileNameButton = QPushButton("Load reflection file")
    self.openFileNameButton.setDefault(True)
    self.openFileNameButton.clicked.connect(self.OpenReflectionsFile)

    self.flatPushButton = QPushButton("Flat Push Button")
    self.flatPushButton.setFlat(True)
    self.flatPushButton.clicked.connect(self.DoSomething)

    self.MillerComboBox = QComboBox()
    self.MillerComboBox.activated.connect(self.MillerComboSelchange)

    self.MillerLabel = QLabel()
    self.MillerLabel.setText("Selected Miller Array")

    self.FOMComboBox = QComboBox()
    self.FOMComboBox.activated.connect(self.FOMComboSelchange)

    self.FOMLabel = QLabel()
    self.FOMLabel.setText("Use Figure of Merits")

    self.SpaceGroupComboBox = QComboBox()
    self.SpaceGroupComboBox.activated.connect(self.SpacegroupSelchange)

    self.SpacegroupLabel = QLabel()
    self.SpacegroupLabel.setText("Space Subgroups")

    self.mergecheckbox = QCheckBox()
    self.mergecheckbox.setText("Merge data")
    #self.mergecheckbox.setTristate (True)
    self.mergecheckbox.clicked.connect(self.MergeData)

    self.expandP1checkbox = QCheckBox()
    self.expandP1checkbox.setText("Expand to P1")
    self.expandP1checkbox.clicked.connect(self.ExpandToP1)

    self.expandAnomalouscheckbox = QCheckBox()
    self.expandAnomalouscheckbox.setText("Show Friedel pairs")
    self.expandAnomalouscheckbox.clicked.connect(self.ExpandAnomalous)

    self.sysabsentcheckbox = QCheckBox()
    self.sysabsentcheckbox.setText("Show Systematic Absences")
    self.sysabsentcheckbox.clicked.connect(self.showSysAbsent)

    self.missingcheckbox = QCheckBox()
    self.missingcheckbox.setText("Show Missing")
    self.missingcheckbox.clicked.connect(self.showMissing)

    self.onlymissingcheckbox = QCheckBox()
    self.onlymissingcheckbox.setText("Only Show Missing")
    self.onlymissingcheckbox.clicked.connect(self.showOnlyMissing)

    self.showslicecheckbox = QCheckBox()
    self.showslicecheckbox.setText("Show Slice")
    self.showslicecheckbox.clicked.connect(self.showSlice)

    self.sliceindexspinBox = QDoubleSpinBox()
    self.sliceindex = 0
    self.sliceindexspinBox.setValue(self.sliceindex)
    self.sliceindexspinBox.setDecimals(0)
    self.sliceindexspinBox.setSingleStep(1)
    self.sliceindexspinBox.setRange(0, 20)
    self.sliceindexspinBox.valueChanged.connect(self.onSliceIndexChanged)

    self.SliceLabelComboBox = QComboBox()
    self.SliceLabelComboBox.activated.connect(self.onSliceComboSelchange)
    self.sliceaxis = [ "h", "k", "l" ]
    self.SliceLabelComboBox.addItems( self.sliceaxis )

    self.HKLnameedit = QLineEdit('')
    self.HKLnameedit.setReadOnly(True)
    self.textInfo = QTextEdit()
    self.textInfo.setReadOnly(True)

    self.RadiiScaleGroupBox = QGroupBox("Radii Size of HKL Spheres")

    #self.PowerScaleGroupBox = QGroupBox("Manual Power Scaling of Sphere Radii")

    self.ManualPowerScalecheckbox = QCheckBox()
    self.ManualPowerScalecheckbox.setText("Manual Power Scaling of Sphere Radii")
    self.ManualPowerScalecheckbox.clicked.connect(self.onManualPowerScale)

    self.power_scale_spinBox = QDoubleSpinBox(self.RadiiScaleGroupBox)
    self.nth_power_scale = 0.5
    self.power_scale_spinBox.setValue(self.nth_power_scale)
    self.power_scale_spinBox.setDecimals(2)
    self.power_scale_spinBox.setSingleStep(0.05)
    self.power_scale_spinBox.setRange(0.0, 1.0)
    self.power_scale_spinBox.valueChanged.connect(self.onPowerScaleChanged)
    self.powerscaleLabel = QLabel()
    self.powerscaleLabel.setText("Power scale Factor")

    self.radii_scale_spinBox = QDoubleSpinBox(self.RadiiScaleGroupBox)
    self.radii_scale = 1.0
    self.radii_scale_spinBox.setValue(self.radii_scale)
    self.radii_scale_spinBox.setDecimals(1)
    self.radii_scale_spinBox.setSingleStep(0.1)
    self.radii_scale_spinBox.setRange(0.2, 2.0)
    self.radii_scale_spinBox.valueChanged.connect(self.onRadiiScaleChanged)
    self.radiiscaleLabel = QLabel()
    self.radiiscaleLabel.setText("Linear Scale Factor")

    self.millertable = QTableWidget(0, 8)
    labels = ["label", "type", "no. of HKLs", "span of HKLs",
       "min max data", "min max sigmas", "d_min, d_max", "symmetry unique"]
    self.millertable.setHorizontalHeaderLabels(labels)
    # don't allow editing the miller array info
    self.millertable.setEditTriggers(QTableWidget.NoEditTriggers)

    self.createTopLeftGroupBox()
    #self.createTopRightGroupBox()
    self.createBottomLeftTabWidget()
    self.createRadiiScaleGroupBox()

    #topLayout = QHBoxLayout()
    #topLayout.addWidget(self.openFileNameButton)
    #topLayout.addStretch(1)

    mainLayout = QGridLayout()
    mainLayout.addWidget(self.openFileNameButton,  0, 0, 1, 1)
    mainLayout.addWidget(self.HKLnameedit,         1, 0, 1, 1)

    mainLayout.addWidget(self.topLeftGroupBox,     2, 0)
    #mainLayout.addWidget(self.topRightGroupBox,    1, 1)
    mainLayout.addWidget(self.RadiiScaleGroupBox,  3, 0)
    mainLayout.addWidget(self.bottomLeftGroupBox,  4, 0)
    mainLayout.setRowStretch(0, 0)
    mainLayout.setRowStretch(1, 0)
    mainLayout.setRowStretch(2, 0)
    mainLayout.setRowStretch(3, 0)
    mainLayout.setRowStretch(4, 1)
    #mainLayout.setColumnStretch(0, 1)
    #mainLayout.setColumnStretch(1, 0)
    self.setLayout(mainLayout)

    self.setWindowTitle("NGL-HKL-viewer")
    self.cctbxproc = None
    self.LaunchCCTBXPython()
    self.out = None
    self.err = None
    self.miller_arrays = None
    self.matching_arrays = None
    self.bin_info = None
    self.html_url = None
    self.spacegroups = None
    self.info = None
    self.infostr = ""
    self.fileisvalid = False
    self.NewFileLoaded = False

    #self.msgqueuethrd = threading.Thread(target = self.update )
    #self.msgqueuethrd.daemon = True
    #self.msgqueuethrd.start()

    self.show()


  def update(self):
    #while 1:
    #  time.sleep(1)
    if self.cctbxproc:
      if self.cctbxproc.stdout:
        print(self.cctbxproc.stdout.read().decode("utf-8"))
      if self.cctbxproc.stderr:
        print(self.cctbxproc.stderr.read().decode("utf-8"))
    if self.out:
      print(self.out.decode("utf-8"))
    if self.err:
      print(self.err.decode("utf-8"))
    #print("in update\n")
    if self.context:
      try:
        msg = self.socket.recv(flags=zmq.NOBLOCK) #To empty the socket from previous messages
        #msg = self.socket.recv()
        self.info = eval(msg.decode())

        ngl_hkl_infodict = self.info
        if ngl_hkl_infodict:
          self.miller_arrays = ngl_hkl_infodict["miller_arrays"]
          self.matching_arrays = ngl_hkl_infodict["matching_arrays"]
          self.bin_info = ngl_hkl_infodict["bin_info"]
          self.html_url = ngl_hkl_infodict["html_url"]
          self.spacegroups = ngl_hkl_infodict["spacegroups"]
          self.mergedata = ngl_hkl_infodict["mergedata"]
          self.infostr = ngl_hkl_infodict["info"]
          self.NewFileLoaded = ngl_hkl_infodict["NewFileLoaded"]
          self.fileisvalid = True

          if self.infostr:
            print(self.infostr)
            self.textInfo.setPlainText(self.infostr)
            #self.mergecheckbox.setEnabled(True)
          else:
            self.textInfo.setPlainText("")
            #self.mergecheckbox.setEnabled(False)
          if self.NewFileLoaded:
            if self.mergedata == True : val = 2
            if self.mergedata == None : val = 1
            if self.mergedata == False : val = 0
            self.mergecheckbox.setCheckState(val )

            self.MillerComboBox.clear()
            self.MillerComboBox.addItems( [ (str(e[0]) + " (" + str(e[1]) +")" )
                                             for e in self.miller_arrays ] )
            self.FOMComboBox.clear()
            self.FOMComboBox.addItems( [ (str(e[0]) + " (" + str(e[1]) +")" )
                                           for e in self.miller_arrays ] )
            self.SpaceGroupComboBox.clear()
            self.SpaceGroupComboBox.addItems( self.spacegroups )

            self.millertable.setRowCount(len(self.miller_arrays))
            #self.millertable.setColumnCount(8)
            for n,millarr in enumerate(self.miller_arrays):
              for m,elm in enumerate(millarr):
                self.millertable.setItem(n, m, QTableWidgetItem(str(elm)))

      except Exception as e:
        #print( str(e) )
        pass




  def MergeData(self):
    if self.mergecheckbox.checkState()== 2:
      self.NGL_HKL_command('NGL_HKLviewer.mergedata = True')
    if self.mergecheckbox.checkState()== 1:
      self.NGL_HKL_command('NGL_HKLviewer.mergedata = None')
    if self.mergecheckbox.checkState()== 0:
      self.NGL_HKL_command('NGL_HKLviewer.mergedata = False')


  def ExpandToP1(self):
    if self.expandP1checkbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.expand_to_p1 = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.expand_to_p1 = False')


  def ExpandAnomalous(self):
    if self.expandAnomalouscheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.expand_anomalous = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.expand_anomalous = False')


  def showSysAbsent(self):
    if self.sysabsentcheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_systematic_absences = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_systematic_absences = False')


  def showMissing(self):
    if self.missingcheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_missing = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_missing = False')


  def showOnlyMissing(self):
    if self.onlymissingcheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_only_missing = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.show_only_missing = False')


  def showSlice(self):
    if self.showslicecheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.slice_mode = True')
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.slice_mode = False')


  def onSliceComboSelchange(self,i):
    rmin = self.miller_arrays[self.MillerComboBox.currentIndex()][3][0][i]
    rmax = self.miller_arrays[self.MillerComboBox.currentIndex()][3][1][i]
    self.sliceindexspinBox.setRange(rmin, rmax)
    self.NGL_HKL_command("NGL_HKLviewer.viewer.slice_axis = %s" % self.sliceaxis[i] )


  def onSliceIndexChanged(self, val):
    self.sliceindex = val
    self.NGL_HKL_command("NGL_HKLviewer.viewer.slice_index = %d" %self.sliceindex)


  def onRadiiScaleChanged(self, val):
    self.radii_scale = val
    self.NGL_HKL_command("""
      NGL_HKLviewer.viewer {
        nth_power_scale_radii = %f
        scale = %f
      }
      """ %(self.nth_power_scale, self.radii_scale)
    )


  def onPowerScaleChanged(self, val):
    self.nth_power_scale = val
    self.NGL_HKL_command("""
      NGL_HKLviewer.viewer {
        nth_power_scale_radii = %f
        scale = %f
      }
      """ %(self.nth_power_scale, self.radii_scale)
    )


  def onManualPowerScale(self):
    if self.ManualPowerScalecheckbox.isChecked():
      self.NGL_HKL_command('NGL_HKLviewer.viewer.nth_power_scale_radii = %f' %self.nth_power_scale)
      self.power_scale_spinBox.setEnabled(True)
    else:
      self.NGL_HKL_command('NGL_HKLviewer.viewer.nth_power_scale_radii = -1.0')
      self.power_scale_spinBox.setEnabled(False)


  def OpenReflectionsFile(self):
    options = QFileDialog.Options()
    fileName, filtr = QFileDialog.getOpenFileName(self,
            "Load reflections file",
            "",
            "All Files (*);;MTZ Files (*.mtz);;CIF (*.cif)", "", options)
    if fileName:
      self.HKLnameedit.setText(fileName)
      self.fileisvalid = False
      self.NGL_HKL_command('NGL_HKLviewer.filename = "%s"' %fileName )
      self.MillerComboBox.clear()
      self.FOMComboBox.clear()
      #while not self.fileisvalid:
      #  time.sleep(1)
        #print("file not valid")


  def createTopLeftGroupBox(self):
    self.topLeftGroupBox = QGroupBox("Group 1")

    layout = QGridLayout()
    layout.addWidget(self.MillerComboBox,            1, 1, 1, 1)
    layout.addWidget(self.MillerLabel,               1, 0, 1, 1)
    layout.addWidget(self.FOMComboBox,               2, 1, 1, 1)
    layout.addWidget(self.FOMLabel,                  2, 0, 1, 1)
    layout.addWidget(self.SpaceGroupComboBox,        3, 1, 1, 1)
    layout.addWidget(self.SpacegroupLabel,           3, 0, 1, 1)
    layout.addWidget(self.mergecheckbox,             4, 0, 1, 1)
    layout.addWidget(self.expandP1checkbox,          4, 1, 1, 1)
    layout.addWidget(self.expandAnomalouscheckbox,   5, 0, 1, 1)
    layout.addWidget(self.sysabsentcheckbox,         5, 1, 1, 1)
    layout.addWidget(self.missingcheckbox,           6, 0, 1, 1)
    layout.addWidget(self.onlymissingcheckbox,       6, 1, 1, 1)

    layout.addWidget(self.showslicecheckbox,         7, 0, 1, 1)
    layout.addWidget(self.SliceLabelComboBox,        7, 1, 1, 1)
    layout.addWidget(self.sliceindexspinBox,         7, 2, 1, 1)
    #layout.addStretch(1)
    self.topLeftGroupBox.setLayout(layout)


  def createTopRightGroupBox(self):
    self.topRightGroupBox = QGroupBox("Group 2")
    togglePushButton = QPushButton("Toggle Push Button")
    togglePushButton.setCheckable(True)
    togglePushButton.setChecked(True)

    slider = QSlider(Qt.Horizontal, self.RadiiScaleGroupBox)
    slider.setValue(40)

    scrollBar = QScrollBar(Qt.Horizontal, self.RadiiScaleGroupBox)
    scrollBar.setValue(60)

    dial = QDial(self.RadiiScaleGroupBox)
    dial.setValue(30)
    dial.setNotchesVisible(True)

    layout = QVBoxLayout()
    layout.addWidget(self.openFileNameButton)
    layout.addWidget(togglePushButton)
    layout.addWidget(self.flatPushButton)

    layout.addWidget(slider)
    layout.addWidget(scrollBar)
    layout.addWidget(dial)

    layout.addStretch(1)
    self.topRightGroupBox.setLayout(layout)


  def DoSomething(self):
    print( self.miller_arrays )
    print( self.matching_arrays )
    print( self.bin_info )
    print( self.html_url )
    print( self.spacegroups )
    print(self.info)
    import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )



  def createBottomLeftTabWidget(self):
    self.bottomLeftGroupBox = QGroupBox("Group 3")
    layout = QGridLayout()

    """
    self.bottomLeftTabWidget = QTabWidget()
    self.bottomLeftTabWidget.setSizePolicy(QSizePolicy.Preferred,  QSizePolicy.Ignored)
    tab1 = QWidget()
    tab1hbox = QHBoxLayout()
    tab1hbox.setContentsMargins(5, 5, 5, 5)
    tab1hbox.addWidget(self.millertable)
    tab1.setLayout(tab1hbox)
    tab2 = QWidget()

    tab2hbox = QHBoxLayout()
    tab2hbox.setContentsMargins(5, 5, 5, 5)
    tab2hbox.addWidget(self.textInfo)
    tab2.setLayout(tab2hbox)
    self.bottomLeftTabWidget.addTab(tab1, "&Miller Arrays")
    self.bottomLeftTabWidget.addTab(tab2, "Information")
    """

    layout.addWidget(self.millertable,            0, 0, 1, 1)
    layout.addWidget(self.textInfo,               1, 0, 1, 1)
    layout.setRowStretch (0, 1)
    layout.setRowStretch (1 ,0)
    self.bottomLeftGroupBox.setLayout(layout)


  def MillerComboSelchange(self,i):
    self.NGL_HKL_command("NGL_HKLviewer.column = %d" %i)
    if self.miller_arrays[ i ][1] == 'Map coeffs':
      self.FOMComboBox.setEnabled(True)
    else:
      self.FOMComboBox.setEnabled(False)

    self.SpaceGroupComboBox.clear()
    self.SpaceGroupComboBox.addItems( self.spacegroups )

    if self.miller_arrays[ i ][7] == False:
      self.mergecheckbox.setEnabled(True)
    else:
      self.mergecheckbox.setEnabled(False)


  def FOMComboSelchange(self,i):
    self.NGL_HKL_command("NGL_HKLviewer.fomcolumn = %d" %i)


  def SpacegroupSelchange(self,i):
    self.NGL_HKL_command("NGL_HKLviewer.spacegroupchoice = %d" %i)


  def createRadiiScaleGroupBox(self):
    layout = QGridLayout()
    layout.addWidget(self.ManualPowerScalecheckbox, 1, 0, 1, 2)
    layout.addWidget(self.powerscaleLabel,          2, 0, 1, 2)
    layout.addWidget(self.power_scale_spinBox,      2, 1, 1, 2)
    layout.addWidget(self.radiiscaleLabel,          3, 0, 1, 2)
    layout.addWidget(self.radii_scale_spinBox,      3, 1, 1, 2)
    layout.setColumnStretch (0, 1)
    layout.setColumnStretch (1 ,0)
    self.RadiiScaleGroupBox.setLayout(layout)



  def LaunchCCTBXPython(self):
    self.context = zmq.Context()
    self.socket = self.context.socket(zmq.PAIR)
    self.socket.bind("tcp://127.0.0.1:7895")
    try: msg = self.socket.recv(flags=zmq.NOBLOCK) #To empty the socket from previous messages
    except Exception: pass

    cmdargs = 'cctbx.python.bat -i -c "from crys3d.hklview import cmdlineframes; myHKLview = cmdlineframes.HKLViewFrame(useSocket=True, verbose=False)"\n'
    #self.cctbxproc = subprocess.Popen( cmdargs, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    self.cctbxproc = subprocess.Popen( cmdargs, shell=True, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr)
    #self.cctbxproc = subprocess.Popen( cmdargs, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    time.sleep(1)
    """
    self.NGL_HKL_command('''
    NGL_HKLviewer {
      filename = "C:\\Users\\oeffner\\Buser\\Tests\LLGperResidue\\map_1six.1.mtz"
      column=0
      viewer.expand_anomalous=True
    }
    ''')
    """

  def NGL_HKL_command(self, cmdstr):
    #stdinstr = "myHKLview.ExecutePhilString("+  cmdstr + ")\n"
    #self.cctbxproc.stdin.write( stdinstr.encode() )
    print("sending:\n" + cmdstr)
    self.socket.send(bytes(cmdstr,"utf-8"))
    print("stuff sent")

    #( self.out, self.err ) = self.cctbxproc.communicate(input = bytes(cmdstr, 'utf-8') )

    #print(self.err)


if __name__ == '__main__':
  app = QApplication(sys.argv)
  guiobj = NGL_HKLViewer()

  timer = QTimer()
  timer.setInterval(0.1)
  timer.timeout.connect(guiobj.update)
  timer.start(500)

  if guiobj.cctbxproc:
    guiobj.cctbxproc.terminate()
  sys.exit(app.exec_())
