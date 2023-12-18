from pathlib import Path

from PySide2.QtCore import Signal

from ..widgets.tab import GUITabWidget
from ..restraints_table.bonds import BondsTableTabView
from ..restraints_table.angles import AnglesTableTabView


class RestraintsTableTopTabView(GUITabWidget):
  trigger_processing = Signal(object)
  def __init__(self,parent=None):
    super().__init__(parent=parent)


    # Bonds
    self.bonds = BondsTableTabView(parent=self)
    self.addTab(self.bonds, "Bonds")

    # Angles
    self.angles = AnglesTableTabView(parent=self)
    self.addTab(self.angles, "Angles")
