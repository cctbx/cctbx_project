import subprocess

from PySide2.QtCore import QObject, QEvent, QSize
from PySide2.QtWidgets import QWidget,QDialog, QVBoxLayout, QPushButton

class GlobalEventFilter(QObject):
  def eventFilter(self, watched, event):
    if event.type() == QEvent.Enter:
      if isinstance(watched, QWidget):
        print(f"Mouse entered: {watched.objectName()} ({type(watched).__name__})")
    elif event.type() == QEvent.Leave:
      if isinstance(watched, QWidget):
        print(f"Mouse left: {watched.objectName()} ({type(watched).__name__})")
    return super().eventFilter(watched, event)



def check_program_access(programs):
    inaccessible_programs = []

    for program in programs:
        try:
            subprocess.run([program, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except FileNotFoundError:
            inaccessible_programs.append(program)
        except subprocess.CalledProcessError:
            pass

    return inaccessible_programs


class ViewerChoiceDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Choose Viewer")
        self.choice = None
        self.initUI()

    def initUI(self):
        self.setMinimumSize(QSize(200, 100)) 


        layout = QVBoxLayout()
        

        # Buttons for choices
        btnOption1 = QPushButton("Molstar", self)
        btnOption1.clicked.connect(self.select_molstar)

        btnOption2 = QPushButton("ChimeraX", self)
        btnOption2.clicked.connect(self.select_chimerax)

        layout.addWidget(btnOption1)
        layout.addWidget(btnOption2)

        self.setLayout(layout)

    def select_molstar(self):
        self.choice = 'molstar'
        self.accept()

    def select_chimerax(self):
        self.choice = 'chimerax'
        self.accept()