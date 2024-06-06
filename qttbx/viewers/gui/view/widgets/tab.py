from PySide2.QtWidgets import QMainWindow, QTabBar, QTabWidget
from PySide2.QtCore import Qt, QMimeData, Signal
from PySide2.QtGui import QGuiApplication, QDrag, QCursor

# Keep track of active toasts at module level
_active_toasts = []

class GUITab(QTabWidget):
  def __init__(self,parent=None,order_index=None):
    self.parent_explicit = parent # often the implicit parent, self.parent(), is not what you want
    super().__init__(parent)
    self.was_visited = False
    self.order_index=order_index

  def on_first_visit(self):
    #print("First time visiting: "+str(self))
     pass

  def set_focus_on(self,*args):
    self.set_visible_on()
    name = self.parent_explicit.tabs.findNameByTab(self)
    i = self.parent_explicit.tabs.findIndexByName(name)
    self.parent_explicit.tabs.setCurrentIndex(i)

  def set_visible_on(self,*args):
    name = self.parent_explicit.tabs.findNameByTab(self)
    self.parent_explicit.tabs.toggle_tab_visible(name,show=True)

  def set_visible_off(self,*args):
    name = self.parent_explicit.tabs.findNameByTab(self)
    self.parent_explicit.tabs.toggle_tab_visible(name,show=False)

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
    self.drag_start_index = 0

  def mousePressEvent(self, event):
    self.drag_start_index = self.tabAt(event.pos())
    self.drag_start_pos = event.pos()
    super().mousePressEvent(event)

  def mouseMoveEvent(self, event):
    #print("Mouse move event triggered.")
    if not event.buttons() & Qt.LeftButton:
        #print("Left button not pressed.")
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
        #print(f"Global position: {global_pos}")
        #print(f"Window rectangle: {window_rect}")

        if not window_rect.contains(global_pos):
            #print("Dropped outside window!")
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
    #print("Enter drag Event")
    event.acceptProposedAction()

  def dropEvent(self, event):
    #print("Drop event")
    event.mimeData().text()
    # Handle the drop: create a new window, add the tab, etc.
    new_window = QMainWindow()
    new_tab_widget = DraggableTabWidget()
    new_window.setCentralWidget(new_tab_widget)
    new_window.show()

    event.acceptProposedAction()

class GUITabWidget(DraggableTabWidget):
  def __init__(self,parent=None,order_index=None):
    super().__init__(parent=parent)
    self.parent_explicit = parent
    self.order_index = order_index
    self.currentChanged.connect(self.on_tab_changed)
    self.hiddenTabs = {}  # Track hidden tabs as {tabName: widget}

  def set_focus_on(self,*args):
    self.set_visible_on()
    name = self.parent_explicit.tabs.findNameByTab(self)
    i = self.parent_explicit.tabs.findIndexByName(name)
    self.parent_explicit.tabs.setCurrentIndex(i)

  def set_visible_on(self,*args):
    # This differs from toggle_tab_visible in that this makes self
    # visible, toggle_tab_visible operates on child tabs
    name = self.parent_explicit.tabs.findNameByTab(self)
    self.parent_explicit.tabs.toggle_tab_visible(name,show=True)

  def on_tab_changed(self,index):
    current_tab_widget = self.widget(index)
    if hasattr(current_tab_widget,"was_visited"):
      if not current_tab_widget.was_visited:
        current_tab_widget.on_first_visit()
        current_tab_widget.was_visited = True

  def toggle_tab_visible(self, tab_name, show=True):
      print("toggle_tab_visible: ",tab_name)
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

  @property
  def widgets(self):
    tab_widgets = []
    for index in range(self.count()):
      tab_widgets.append(self.widget(index))
    return tab_widgets
