from pathlib import Path

from PySide2.QtCore import Signal

from ..widgets.tab import GUITabWidget
from ..restraints_table import (
  BondsTableTabView,
  AnglesTableTabView,
  DihedralsTableTabView,
  ChiralsTableTabView,
  PlanarityTableTabView,
  NonbondedTableTabView
  )


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

    # Dihedrals
    self.dihedrals = DihedralsTableTabView(parent=self)
    self.addTab(self.dihedrals, "Dihedrals")

    # Chirals
    self.chirals = ChiralsTableTabView(parent=self)
    self.addTab(self.chirals, "Chirals")

    # Planes
    self.planes = PlanarityTableTabView(parent=self)
    self.addTab(self.planes, "Planes")

    # Nonbonded
    self.nonbonded = NonbondedTableTabView(parent=self)
    self.addTab(self.nonbonded, "NonBonded")