import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PySide2.QtCore import QObject, QPoint, QTimer, Qt, Signal
from PySide2.QtGui import QMouseEvent, QPainter
from PySide2.QtGui import QCursor
from PySide2.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListView,
    QPushButton,
    QSlider,
    QStyle,
    QStyleOptionComboBox,
    QTextEdit,
    QVBoxLayout,
    QWidget
)

_active_toasts = []


class EditsEditDialog(QDialog):
  input_names = ["ideal", "sigma"]

  def __init__(self, parent=None, input_names=None, defaults_dict=None,action="add"):
    super().__init__(parent)
    
    if action == "add":
      title = "Add"
    elif action == "mod":
      title = "Modify"
    self.action = action
    self.setWindowTitle(f'Geometry Edit: {title}')
    mainLayout = QVBoxLayout()
    self.inputsLayout = QHBoxLayout()

    if input_names is None:
      input_names = self.input_names
    if defaults_dict is None:
      defaults_dict = {name: "" for name in input_names}

    # Dictionary to store QLineEdit references
    self.inputFields = {}

    for name in input_names:
      inputLayout = QVBoxLayout()
      inputLabel = QLabel(name.capitalize())
      inputField = QLineEdit(self)

      default_value = str(defaults_dict.get(name, ""))
      inputField.setText(default_value)

      inputLayout.addWidget(inputLabel)
      inputLayout.addWidget(inputField)

      self.inputsLayout.addLayout(inputLayout)

      # Store QLineEdit reference with its name as key
      self.inputFields[name] = inputField

    mainLayout.addLayout(self.inputsLayout)
    self.setLayout(mainLayout)

    self.buttonsLayout = QHBoxLayout()
    self.cancelButton = QPushButton('Cancel', self)
    self.cancelButton.clicked.connect(self.reject)
    self.buttonsLayout.addWidget(self.cancelButton)

    self.acceptButton = QPushButton('Accept', self)
    self.acceptButton.clicked.connect(self.accept)
    self.buttonsLayout.addWidget(self.acceptButton)

    mainLayout.addLayout(self.buttonsLayout)

  def collectInputValues(self):
    # Collect values from QLineEdit widgets
    values = {name: field.text() for name, field in self.inputFields.items()}
    values["action"] = self.action
    return values

  def exec_(self):
    # Get the cursor's current position
    cursorPos = QCursor.pos()

    # Move the dialog to the cursor's position initially
    self.move(cursorPos.x(), cursorPos.y())

    # Adjust position to ensure the dialog is fully visible on the screen
    screen = QApplication.desktop().screenNumber(cursorPos)
    screen_geom = QApplication.desktop().screenGeometry(screen)

    # Calculate the dialog's geometry after the initial move
    dialog_geom = self.geometry()
    dialog_x = cursorPos.x()
    dialog_y = cursorPos.y()
    dialog_width = dialog_geom.width()
    dialog_height = dialog_geom.height()

    # Check right boundary
    if dialog_x + dialog_width > screen_geom.right():
        dialog_x = screen_geom.right() - dialog_width

    # Check bottom boundary
    if dialog_y + dialog_height > screen_geom.bottom():
        dialog_y = screen_geom.bottom() - dialog_height

    # Check left boundary
    if dialog_x < screen_geom.left():
        dialog_x = screen_geom.left()

    # Check top boundary
    if dialog_y < screen_geom.top():
        dialog_y = screen_geom.top()

    # Move dialog to the adjusted position
    self.move(dialog_x, dialog_y)
    
    result = super().exec_()
    # Collect values if dialog was accepted
    if result == QDialog.Accepted:
      return self.collectInputValues()
    else:
      return None

class BondEditDialog(EditsEditDialog):
   input_names = ["ideal","sigma"]
class AngleEditDialog(EditsEditDialog):
   input_names = ["ideal","sigma"]
class DihedralEditDialog(EditsEditDialog):
   input_names = ["ideal","sigma","weight","harmonic"]
class ChiralEditDialog(EditsEditDialog):
   input_names = ["ideal","sigma"]
class PlaneEditDialog(EditsEditDialog):
   input_names = []

# End edits dialog
class InfoDialog(QDialog):
  def __init__(self, text, title="Info",parent=None):
    super().__init__(parent)
    self.setWindowTitle(title)

    self.text_edit = QTextEdit(self)
    self.text_edit.setReadOnly(True)
    self.text_edit.setText(text)

    layout = QVBoxLayout(self)
    layout.addWidget(self.text_edit)


class Toast(QWidget):
  def __init__(self, message, duration=2000, active_toasts=None, parent_widget=None):
    super(Toast, self).__init__(parent_widget)
    self.setAutoFillBackground(True)
    self.setWindowFlags(Qt.SplashScreen | Qt.FramelessWindowHint)
    self.setAttribute(Qt.WA_TranslucentBackground)

    layout = QVBoxLayout()
    label = QLabel(message)
    layout.addWidget(label)
    self.setLayout(layout)

    # Set maximum size
    self.setMaximumSize(200, 100)

    if active_toasts is not None:
      active_toasts.append(self)
      self.destroyed.connect(lambda: active_toasts.remove(self))

    # Set position relative to parent_widget
    if parent_widget:
      parent_pos = parent_widget.mapToGlobal(QPoint(0, 0))
      parent_geom = parent_widget.frameGeometry()
      self.show()  # Need to show first to get frameGeometry
      toast_geom = self.frameGeometry()
      x = parent_pos.x() + (parent_geom.width() - toast_geom.width()) // 2
      y = parent_pos.y() + (parent_geom.height() - toast_geom.height()) // 2
      self.move(QPoint(x, y))

    QTimer.singleShot(duration, self.close)

  def mousePressEvent(self, event):
    self.close()


class ClickableHistogramEmitter(QObject):
  histogram_click_value = Signal(float)

class ClickableHistogramSeaborn(FigureCanvas):

  @staticmethod
  def fisher_transform(r):
    z = 0.5 * np.log((1 + r) / (1 - r))
    return z

  @staticmethod
  def get_bins(data,fisher_first=False):
    if fisher_first:
      data = ClickableHistogramSeaborn.fisher_transform(data)
    IQR = np.subtract(*np.percentile(data, [75, 25]))
    bin_width = 2 * IQR * len(data) ** (-1/3)
    bins = int((data.max() - data.min()) / bin_width)
    return bins


  def __init__(self, data,bins=None, parent=None, width=5, height=4, dpi=100):
    # Signals
    self.emitter = ClickableHistogramEmitter()

    # Matplotlib
    sns.set_theme(style="white", context="paper")
    self.fig = Figure(figsize=(width, height), dpi=dpi)
    self.axes = self.fig.add_subplot(111)
    FigureCanvas.__init__(self, self.fig)
    self.setParent(parent)
    self.data = data
    self.bins = bins
    if self.bins is None:
      self.bins = self.get_bins(data)
    self.bins = max(self.bins,40) # min number of bins
    self.vline = None
    self.text = None
    self.plot_histogram()

  def plot_histogram(self):
    self.axes.clear()
    sns.histplot(self.data, bins=self.bins, kde=False, ax=self.axes)
    self.draw()

  def mousePressEvent(self, event: QMouseEvent):
    if event.button() == Qt.LeftButton:
      x, y = event.x(), event.y()
      self.calculate_coordinates(x, y)

  def calculate_coordinates(self, x, y):
    canvas_width = self.width()
    canvas_height = self.height()

    # Coordinates relative to the canvas size
    rel_x = x / canvas_width
    rel_y = y / canvas_height

    # Get the axis limits
    x_min, x_max = self.axes.get_xlim()
    y_min, y_max = self.axes.get_ylim()

    # Convert to data coordinates
    data_x = x_min + rel_x * (x_max - x_min)
    y_min + rel_y * (y_max - y_min)

    print(f"Clicked on x = {data_x:.2f}")
    self.emitter.histogram_click_value.emit(data_x)

    if self.vline:  # Remove the previous line if it exists
      self.vline.remove()

    if self.text:  # Remove the previous text if it exists
      self.text.remove()


    self.vline = self.axes.axvline(x=data_x, color='black', linestyle='--')

    # Add text annotation directly above the vline in axis coordinates
    self.text = self.axes.text(data_x, y_max, f"{data_x:.2f}", ha='center', va='bottom',fontsize=14)

    self.draw()



class ClickableHistogramMatplotlib(FigureCanvas):
  def __init__(self, data, parent=None, width=5, height=4, dpi=100):

    self.fig = Figure(figsize=(width, height), dpi=dpi)
    self.axes = self.fig.add_subplot(111)
    FigureCanvas.__init__(self, self.fig)
    self.setParent(parent)
    self.data = data
    self.vline = None
    self.text = None
    self.plot_histogram()

  def plot_histogram(self):
    self.axes.hist(self.data, bins=20, edgecolor='black')
    self.draw()

  def mousePressEvent(self, event: QMouseEvent):
    if event.button() == Qt.LeftButton:
      x, y = event.x(), event.y()
      self.calculate_coordinates(x, y)

  def calculate_coordinates(self, x, y):
    canvas_width = self.width()
    canvas_height = self.height()

    # Coordinates relative to the canvas size
    rel_x = x / canvas_width
    rel_y = y / canvas_height

    # Get the axis limits
    x_min, x_max = self.axes.get_xlim()
    y_min, y_max = self.axes.get_ylim()

    # Convert to data coordinates
    data_x = x_min + rel_x * (x_max - x_min)
    y_min + rel_y * (y_max - y_min)

    print(f"Clicked on x = {data_x:.2f}")

    if self.vline:  # Remove the previous line if it exists
      self.vline.remove()

    if self.text:  # Remove the previous text if it exists
      self.text.remove()


    self.vline = self.axes.axvline(x=data_x, color='black', linestyle='--')

    # Add text annotation directly above the vline in axis coordinates
    self.text = self.axes.text(data_x, y_max, f"{data_x:.2f}", ha='center', va='bottom',fontsize=14)

    self.draw()

class NoCheckComboBox(QComboBox):
  
  def paintEvent(self, event):
    painter = QPainter(self)
    option = QStyleOptionComboBox()

    self.initStyleOption(option)
    option.currentIcon = self.itemIcon(self.currentIndex())
    option.iconSize = self.iconSize()

    # Removes focus rectangle around the icons
    option.state &= ~QStyle.State_HasFocus

    self.style().drawComplexControl(QStyle.CC_ComboBox, option, painter, self)
    self.style().drawControl(QStyle.CE_ComboBoxLabel, option, painter, self)



class NoCheckListView(QListView):
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)

  def drawBranches(self, painter, rect, index):
    pass

  def drawRow(self, painter, option, index):
    # Customize the drawing options here
    option.state &= ~QStyle.State_HasFocus & ~QStyle.State_Selected
    super().drawRow(painter, option, index)



class FloatingVSlider(QSlider):
  def __init__(self):
    super().__init__(Qt.Vertical)
    self.setRange(0, 100)
    self.setValue(50)
    self.setWindowFlags(Qt.Popup)
    self.setMinimumSize(20, 200)  # Min Width, Min Height
    self.setMaximumSize(20, 400)  # Max Width, Max Height


class ISOSlider(QSlider):
  """
  Slider to control the ISO level for a map.
  TODO: Should this be just for active map or individual on each map entry?
  TODO: Separate to a more MVC style
  """
  def __init__(self,map_ref=None):
    super().__init__(Qt.Horizontal)
    self._map_ref = None
    if map_ref is not None:
      self.map_ref = map_ref



    # self.setRange(0, 100)
    # self.setValue(50)
    #self.setWindowFlags(Qt.Popup)
    self.setMinimumSize(300, 20)  # Min Width, Min Height
    self.setMaximumSize(600, 20)  # Max Width, Max Height

  @property
  def map_ref(self):
    return self._map_ref

  @map_ref.setter
  def map_ref(self,map_ref):

    # ISO
    mm = map_ref.map_manager
    a = mm.map_data().as_numpy_array()
    real_min = 0.0
    real_max = np.max(a)
    real_start = 1.0#np.std(a)*3
    self.real_start = real_start
    self.iso_start = real_start
    granularity = 0.01  # Set granularity (the smallest change)


    # Calculate scaled values
    slider_iso_min = int(real_min / granularity)
    slider_iso_max = int(real_max / granularity)
    slider_iso_start = self.real_to_slider_iso(real_start,granularity)
    self.slider_iso_granularity = granularity
    slider_iso_step = int(0.05 / granularity)  # Set small step granularity here
    slider_iso_page_step = int(0.1 / granularity)  # Set large step granularity here

    slider_iso = self
    slider_iso.setMinimum(slider_iso_min)
    slider_iso.setMaximum(slider_iso_max)
    slider_iso.setValue(slider_iso_start)
    slider_iso.setSingleStep(slider_iso_step)
    slider_iso.setPageStep(slider_iso_page_step)
    #slider_iso.setMinimumWidth(250)


  def real_to_slider_iso(self,real_value, granularity):
    return int(real_value / granularity)

  def slider_iso_to_real(self,slider_value, granularity):
    return slider_value * granularity



class ISOWidget(QWidget):
  """
  Widget to hold a label too.
  TODO: Refactor out map_ref from view
  """
  def __init__(self,parent=None,map_ref=None):
    super().__init__(parent)
    layout = QHBoxLayout()
    self._map_ref = None
    if map_ref is not None:
      self.map_ref = map_ref
    self.slider = ISOSlider(map_ref=map_ref)

    layout.addWidget(self.slider)

    self.setLayout(layout)

    # # Set popup flag
    # self.setWindowFlags(Qt.Popup)
    self.setMinimumSize(30, 20)  # Min Width, Min Height



  @property
  def map_ref(self):
    return self._map_ref

  @map_ref.setter
  def map_ref(self,map_ref):
    self.slider.map_ref = map_ref
    self.parent().iso_label.setText(f"ISO: {round(self.slider.real_start,2)}")  # Initialize with the scaled value
    self._map_ref = map_ref






class OpacityWidget(QDialog):
  def __init__(self,parent=None):
    super().__init__(parent)
    layout = QHBoxLayout()
    self.default_opacity = 0.8


    self.slider = QSlider(Qt.Horizontal)
    self.slider.setInvertedAppearance(True)  # Invert visual appearance

    self.slider.setMinimumSize(300, 20)  # Min Width, Min Height
    self.slider.setMaximumSize(600, 20)  # Max Width, Max Height
    self.slider.setMinimum(0)
    self.slider.setMaximum(100)
    self.slider.setValue(int((self.default_opacity)*100))
    self.slider.setSingleStep(1)
    self.slider.setPageStep(5)

    self.label = QLabel(f"Opacity: {round(self.default_opacity,2)}")  # Initialize with the scaled value

    # Add slider and label to layout
    layout.addWidget(self.label)
    layout.addWidget(self.slider)

    self.setLayout(layout)

    # Set popup flag
    #self.setWindowFlags(Qt.Popup)