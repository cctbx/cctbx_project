from PySide2.QtWidgets import QVBoxLayout
try:
  from qtconsole.rich_jupyter_widget import RichJupyterWidget
  from qtconsole.inprocess import QtInProcessKernelManager
except:
  pass
from PySide2.QtCore import Slot
from PySide2.QtWidgets import QVBoxLayout,QLineEdit, QPlainTextEdit
from PySide2.QtCore import Slot

from ..widgets.tab import GUITab
#from ..widgets import Toast, _active_toasts

class JSConsoleTab(GUITab):
  def __init__(self, parent=None,web_view=None):
    super().__init__(parent=parent)
    self.web_view = web_view
    self.web_view.page().console_message.connect(self.append_js_message)

    self.console_output = QPlainTextEdit()
    self.console_output.setReadOnly(True)

    self.console_input = QLineEdit()
    self.console_input.setPlaceholderText("Execute JavaScript. Access molstar app with 'this.viewer'. Attach new variables to 'window'. Use double quotes. ")
    self.console_input.returnPressed.connect(self.execute_js)

    layout = QVBoxLayout()
    layout.addWidget(self.console_output)
    layout.addWidget(self.console_input)

    self.setLayout(layout)

  @Slot(str, int, str)
  def append_js_message(self, msg, line, sourceID):
    message = msg #f"JS console message: {msg}, line: {line}, sourceID: {sourceID}"
    self.console_output.appendPlainText(message)

  @Slot()
  def execute_js(self):
    js_code = self.console_input.text()
    
    js_code_exec = f"""
    (function() {{
      try {{
        const result = eval('{js_code}');
        if (result === undefined) {{
          return "undefined";
        }}
        return JSON.stringify(result);
      }} catch (error) {{
        if (error instanceof TypeError && error.message.includes('Converting circular structure to JSON')) {{
          return "object";
        }}
        return String(error);
      }}
    }})()
    """
    self.console_output.appendPlainText(f">> {js_code}")
    self.web_view.page().runJavaScript(js_code_exec, 0, self.js_callback)


  def js_callback(self, result):
    self.console_output.appendPlainText(f"[Out]:\n{result}")


# class JSConsoleTab(GUITab):
#   def __init__(self,parent, web_view=None, *args, **kwargs):
#     assert web_view is not None, "Initialize with a WebView object"
#     super().__init__(parent,*args, **kwargs)
#     self.web_view = web_view

#     self.console_input = QTextEdit()
#     self.console_output = QTextEdit()
#     self.console_output.setReadOnly(True)
    
#     execute_button = QPushButton("Execute JS")
#     execute_button.clicked.connect(self.execute_js)

#     layout = QVBoxLayout()
#     layout.addWidget(self.console_input)
#     layout.addWidget(execute_button)
#     layout.addWidget(self.console_output)

#     self.setLayout(layout)

#   @Slot()
#   def execute_js(self):
#     js_code = self.console_input.toPlainText()
#     world_id = 0  # default world
#     self.web_view.page().runJavaScript(js_code, world_id, self.js_callback)

#   def js_callback(self, result):
#     print(result)
#     self.console_output.append(f"Result: {result}")


class JupyterTabWidget(GUITab):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    

    layout = QVBoxLayout()
    
    kernel_manager = QtInProcessKernelManager()
    

    kernel_manager.start_kernel()
    kernel_manager.kernel.shell.banner1 = f"""
    Execute python code here.
    Access the program with 'app'. 

    For example:
      {'app.view':<30} : {'the Qt main window':<30}
      {'app.controller.molstar.viewer':<30} : {'the python interface for the Mol* viewer':<30}
      {'app.state.data_manager':<30} : {'the cctbx data manager':<30}
    """ 
    kernel_manager.kernel.gui = 'qt'

    kernel_client = kernel_manager.client()
    kernel_client.start_channels()

    self.jupyter_widget = RichJupyterWidget()
    self.jupyter_widget.kernel_manager = kernel_manager
    self.jupyter_widget.kernel_client = kernel_client



    layout.addWidget(self.jupyter_widget)
    self.setLayout(layout)
    


  def on_first_visit(self):
    pass
    # msg = "Interact with the program via 'app'\n\nExample:\n \t' app.model.get_atoms() '"#\n\t' import phenix '\n\t' from cctbx.array_family import flex '"
    # self.toast = Toast(msg,duration=4000,parent_widget=self)
    # self.toast.show()
    #QMessageBox.information(None, "Title", "Your message here.")


    
  




