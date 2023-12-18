from PySide2.QtCore import Signal

from ..widgets.tab import GUITabWidget
from ..restraints.bonds import BondsTabView
from ..restraints.angles import AnglesTabView


class RestraintsTopTabView(GUITabWidget):
  trigger_processing = Signal(object)
  def __init__(self,parent=None):
    super().__init__(parent=parent)


    # Bonds
    self.bonds = BondsTabView(parent=self)
    self.addTab(self.bonds, "Bonds")

    # Angles
    self.angles = AnglesTabView(parent=self)
    self.addTab(self.angles, "Angles")
