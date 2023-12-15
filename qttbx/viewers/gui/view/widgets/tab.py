from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide2.QtWidgets import QTabWidget, QTabBar, QApplication, QMainWindow
from PySide2.QtCore import Qt, QMimeData, Signal
from PySide2.QtGui import QGuiApplication, QDrag, QCursor

# Keep track of active toasts at module level
_active_toasts = []

class GUITab(QTabWidget):
  def __init__(self,parent=None):
    self.parent_explicit = parent # often the implicit parent, self.parent(), is not what you want
    super().__init__(parent)
    self.was_visited = False

  def on_first_visit(self):
    print("First time visiting: "+str(self))


class ChildWindow(QMainWindow):
  
  def __init__(self, original_tab_widget, tab_widget, tab_index, tab_name,*args, **kwargs):
    super().__init__(*args, **kwargs)
    self.setParent(original_tab_widget) # set parent to close with original window
    self.setWindowFlags(Qt.Window)
    self.original_tab_widget = original_tab_widget
    self.tab_widget = tab_widget
    self.tab_index = tab_index
    self.tab_name = tab_name
    # Get screen geometry
    desktop = QGuiApplication.primaryScreen()
    screen_rect = desktop.availableGeometry()

    # Set main window geometry to take up the left half of the screen
    self.setGeometry(screen_rect.x(), screen_rect.y(),
                            screen_rect.width() // 2, screen_rect.height())

    

  def closeEvent(self, event):
    # Remove the tab from the child window's QTabWidget
    widget = self.tab_widget.widget(self.tab_index)
    self.tab_widget.removeTab(self.tab_index)
    
    # Add the tab back to the original QTabWidget
    self.original_tab_widget.insertTab(self.tab_index, widget, self.tab_name)

    
    event.accept()




class DraggableTabBar(QTabBar):
  childSignal = Signal(str)  # Declare the Signal at the class level

  def __init__(self,parent=None):
    super().__init__(parent)
    
  def mousePressEvent(self, event):
    self.drag_start_index = self.tabAt(event.pos())
    self.drag_start_pos = event.pos()
    super().mousePressEvent(event)

  def mouseMoveEvent(self, event):
    print("Mouse move event triggered.")
    if not event.buttons() & Qt.LeftButton:
        print("Left button not pressed.")
        return

    drag = QDrag(self)
    mime_data = QMimeData()
    tab_name = self.tabText(self.drag_start_index)
    mime_data.setText(tab_name)
    drag.setMimeData(mime_data)

    result = drag.exec_(Qt.CopyAction | Qt.MoveAction)

    if result == Qt.IgnoreAction or result == Qt.MoveAction:
        # The drop was outside the window
        global_pos = QCursor.pos()
        window = self.window()
        window_rect = window.geometry()
        print(f"Global position: {global_pos}")
        print(f"Window rectangle: {window_rect}")

        if not window_rect.contains(global_pos):
            print("Dropped outside window!")
    
        # Your code to handle external drop
        tab_widget = self.parentWidget()
        if tab_widget:
          old_content = tab_widget.widget(self.drag_start_index)
          tab_text = self.tabText(self.drag_start_index)
          
          # Remove the tab without deleting the content
          tab_widget.removeTab(self.drag_start_index)




          # Create a new QTabWidget
          new_tab_widget = QTabWidget(self)

          # Create a new QMainWindow
          new_window = ChildWindow(self.parentWidget(),new_tab_widget,self.drag_start_index,tab_name)
          
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

          # persistence
          self.new_window = new_window

          # notify main window
          self.childSignal.emit("created")
      
class DraggableTabWidget(QTabWidget):
  def __init__(self, parent=None):
    super(DraggableTabWidget, self).__init__(parent)
    self.setTabBar(DraggableTabBar())
    self.setMouseTracking(True)
    self.setAcceptDrops(True)
    

  def dragEnterEvent(self, event):
    print("Enter drag Event")
    event.acceptProposedAction()

  def dropEvent(self, event):
    print("Drop event")
    tab_name = event.mimeData().text()
    # Handle the drop: create a new window, add the tab, etc.
    new_window = QMainWindow()
    new_tab_widget = DraggableTabWidget()
    new_window.setCentralWidget(new_tab_widget)
    new_window.show()

    event.acceptProposedAction()

class GUITabWidget(DraggableTabWidget):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    self.currentChanged.connect(self.on_tab_changed)

  def on_tab_changed(self,index):
    current_tab_widget = self.widget(index)
    if hasattr(current_tab_widget,"was_visited"):
      if not current_tab_widget.was_visited:
        current_tab_widget.on_first_visit()
        current_tab_widget.was_visited = True
  @property
  def widgets(self):
    tab_widgets = []
    for index in range(self.count()):
      tab_widgets.append(self.widget(index))
    return tab_widgets


