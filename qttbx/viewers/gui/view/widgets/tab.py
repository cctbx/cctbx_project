"""
These classes are the base classes for all tabs in the gui. 
They implements functionality to show, hide, and focus tabs, as well as functionality
  to drag a tab out into its own window. Which is useful for workspace organization.

Usage: 
  # Create the root tab widget.
  tabs = GUITabWidget()

  # This automaticallly sets the tab bar, tabs.tabBar(), to an instance of DraggableTabBar

  # Add tabs using GUITab or a subclass
  tabs.addTab(GUITab(),"ExampleTab")


"""
from PySide2.QtCore import Qt, QMimeData, QByteArray, QDataStream, QIODevice, Signal, QEvent, QPoint, QCoreApplication
from PySide2.QtGui import QGuiApplication, QDrag, QCursor, QMouseEvent
from PySide2.QtWidgets import QApplication, QMainWindow, QTabWidget, QTabBar, QWidget, QVBoxLayout, QLabel


class GUITab(QTabWidget):
  """
  Base class for all tabs, and therefore the base 'view' for many 
    utilities if they are organized by tab.
  """
  def __init__(self, parent=None, order_index=None):
    self.parent_explicit = parent  # often the implicit parent, self.parent(), is not what you want
    super().__init__(parent)
    self.was_visited = False
    self.order_index = order_index
    self.debug = True

  def log(self,*args):
    if self.debug:
      print(*args)

  def on_first_visit(self):
    pass

  def set_focus_on(self, *args):
    self.set_visible_on()
    name = self.parent_explicit.tabs.findNameByTab(self)
    i = self.parent_explicit.tabs.findIndexByName(name)
    self.parent_explicit.tabs.setCurrentIndex(i)

  def set_visible_on(self, *args):
    name = self.parent_explicit.tabs.findNameByTab(self)
    self.parent_explicit.tabs.toggle_tab_visible(name, show=True)

  def set_visible_off(self, *args):
    name = self.parent_explicit.tabs.findNameByTab(self)
    self.parent_explicit.tabs.toggle_tab_visible(name, show=False)

class ChildWindow(QMainWindow):
  """
  An independent window to hold tabs which have been dragged out of a tab bar 
  """
  def __init__(self, original_tab_widget, tab_widget, tab_index, tab_name, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.setParent(original_tab_widget)  # set parent to close with original window
    self.setWindowFlags(Qt.Window)
    self.original_tab_widget = original_tab_widget
    self.tab_widget = tab_widget
    self.tab_index = tab_index
    self.tab_name = tab_name

    # Get screen geometry
    desktop = QGuiApplication.primaryScreen()
    screen_rect = desktop.availableGeometry()

    # Set main window geometry to take up the left half of the screen
    self.setGeometry(screen_rect.x(), screen_rect.y(), screen_rect.width() // 2, screen_rect.height())

  def closeEvent(self, event):
    # Remove the tab from the child window's QTabWidget
    widget = self.tab_widget.widget(self.tab_index)
    self.tab_widget.removeTab(self.tab_index)

    # Add the tab back to the original QTabWidget
    self.original_tab_widget.insertTab(self.tab_index, widget, self.tab_name)

    event.accept()

class DraggableTabBar(QTabBar):
  """
  A tab bar to hold tabs that can be dragged and dropped into their own window. 
  """
  childSignal = Signal(str)  # Declare the Signal at the class level

  def __init__(self, parent=None):
    super().__init__(parent)
    self.drag_start_index = 0
    self.debug = True

  def log(self,*args):
    if self.debug:
      print(*args)

  def mousePressEvent(self, event):
    drag_start_index = self.tabAt(event.pos())
    if not drag_start_index or drag_start_index<0:
      drag_start_index = 0
    self.drag_start_index = drag_start_index
    self.drag_start_pos = event.pos()
    super().mousePressEvent(event)

  def mouseMoveEvent(self, event):
    if not event.buttons() & Qt.LeftButton:
      return

    if (event.pos() - self.drag_start_pos).manhattanLength() < QApplication.startDragDistance():
      return

    drag = QDrag(self)
    mime_data = QMimeData()

    # Use a custom MIME type to prevent the default file drop behavior
    custom_data = QByteArray()
    stream = QDataStream(custom_data, QIODevice.WriteOnly)
    stream.writeInt32(self.drag_start_index)
    mime_data.setData('application/x-qtabbar-index', custom_data)

    drag.setMimeData(mime_data)
    result = drag.exec_(Qt.CopyAction | Qt.MoveAction)

    if result == Qt.IgnoreAction or result == Qt.MoveAction:
      # The drop was outside the window
      global_pos = QCursor.pos()
      window = self.window()
      window_rect = window.geometry()
      self.log("global_cursor_pos: ",global_pos)

      if not window_rect.contains(global_pos):
        tab_widget = self.parentWidget()
        if tab_widget:
          old_content = tab_widget.widget(self.drag_start_index)
          tab_text = self.tabText(self.drag_start_index)

          # Remove the tab without deleting the content
          tab_widget.removeTab(self.drag_start_index)

          # Create a new QTabWidget
          new_tab_widget = QTabWidget(self)

          # Create a new QMainWindow
          new_window = ChildWindow(self.parentWidget(), new_tab_widget, self.drag_start_index, tab_text)

          # Add the old content to the new QTabWidget
          new_tab_widget.addTab(old_content, tab_text)

          # Make sure it's visible
          old_content.show()

          # Set central widget
          new_window.setCentralWidget(new_tab_widget)

          # Show the new window
          new_window.show()
          new_window.raise_()
          new_window.activateWindow()

          # Persistence
          self.new_window = new_window

          # Notify main window
          self.childSignal.emit("created")

class DraggableTabWidget(QTabWidget):
  """
  A generic base class to implement draggable tabs
  """
  def __init__(self, parent=None):
    super(DraggableTabWidget, self).__init__(parent)
    self.setTabBar(DraggableTabBar())
    self.setMouseTracking(True)
    self.setAcceptDrops(True)

  def dragEnterEvent(self, event):
    self.log("dragEnterEvent on QTabWidget")
    event.acceptProposedAction()

  def dropEvent(self, event):
    self.log("dropEvent on QTabWidget")
    event.ignore()  # Ignore drop events to prevent unintended behavior
    
  def simulate_drag_out(tab_widget, index):
    tab_bar = tab_widget.tabBar()
    drag_start_pos = tab_bar.tabRect(index).center()

    # Calculate the global position of the drag start
    drag_start_global_pos = tab_bar.mapToGlobal(drag_start_pos)

    # Print debug information
    self = tab_widget
    self.log(f"Starting drag simulation for tab at index {index}")
    self.log(f"Drag start position: {drag_start_pos}")
    self.log(f"Drag start global position: {drag_start_global_pos}")

    # Manually set the drag start index and position in the tab bar
    tab_bar.drag_start_index = index
    tab_bar.drag_start_pos = drag_start_pos

    # Directly call the pop-out logic from DraggableTabBar
    global_pos = QPoint(-1000, drag_start_pos.y())  # Simulate a position far outside the window
    window = tab_bar.window()
    window_rect = window.geometry()
    
    if not window_rect.contains(global_pos):
      tab_widget = tab_bar.parentWidget()
      if tab_widget:
        old_content = tab_widget.widget(tab_bar.drag_start_index)
        tab_text = tab_bar.tabText(tab_bar.drag_start_index)

        # Remove the tab without deleting the content
        tab_widget.removeTab(tab_bar.drag_start_index)

        # Create a new QTabWidget
        new_tab_widget = QTabWidget(tab_bar)

        # Create a new QMainWindow
        new_window = ChildWindow(tab_widget, new_tab_widget, tab_bar.drag_start_index, tab_text)

        # Add the old content to the new QTabWidget
        new_tab_widget.addTab(old_content, tab_text)

        # Make sure it's visible
        old_content.show()

        # Set central widget
        new_window.setCentralWidget(new_tab_widget)

        # Show the new window
        new_window.show()
        new_window.raise_()
        new_window.activateWindow()

        # Persistence
        tab_bar.new_window = new_window

        # Notify main window
        tab_bar.childSignal.emit("created")

class GUITabWidget(DraggableTabWidget):
  """
  The top level tab widget, implementing draggable tabs in the GUI
  """
  def __init__(self,parent=None,order_index=None):
    super().__init__(parent=parent)
    self.parent_explicit = parent
    self.order_index = order_index
    self.currentChanged.connect(self.on_tab_changed)
    self.hiddenTabs = {}  # Track hidden tabs as {tabName: widget}
    self.debug = True

  def log(self,*args):
    if self.debug:
      print(*args)


  def on_tab_changed(self,index):
    current_tab_widget = self.widget(index)
    if hasattr(current_tab_widget,"was_visited"):
      if not current_tab_widget.was_visited:
        current_tab_widget.on_first_visit()
        current_tab_widget.was_visited = True

  def toggle_tab_visible(self, tab_name, show=True):
      self.log("toggle_tab_visible: ",tab_name)
      if show:
          if tab_name in self.hiddenTabs:
              # Re-add the tab
              widget = self.hiddenTabs.pop(tab_name)
              i = widget.order_index
              if i is None:
                self.addTab(widget,tab_name)
              else:
                self.insertTab(i,widget, tab_name)
      else:
          index = self.findIndexByName(tab_name)
          if index != -1:
              widget = self.widget(index)
              self.removeTab(index)
              self.hiddenTabs[tab_name] = widget  # Keep track of the widget

  def tabDict(self):
    d = {self.tabText(i):self.widget(i)  for i in range(self.count())}
    d.update(self.hiddenTabs)
    return d

  def findTabByName(self,tab_name):
    return self.tabDict()[tab_name]

  def findIndexByName(self, tab_name):
      for i in range(self.count()):
          if self.tabText(i) ==tab_name:
              return i
      return -1

  def findNameByTab(self,tab_widget):
    d = {v:k for k,v in self.tabDict().items()}
    return d[tab_widget]

  def findNameByIndex(self,tab_index):
    return self.findNameByTab(self.widget(tab_index))

  @property
  def widgets(self):
    tab_widgets = []
    for index in range(self.count()):
      tab_widgets.append(self.widget(index))
    return tab_widgets
