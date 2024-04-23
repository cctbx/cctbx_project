
from ..widgets.tab import GUITabWidget
from .tables import (
  BondsTableTabView,
  AnglesTableTabView,
  DihedralsTableTabView,
  ChiralsTableTabView,
  PlanarityTableTabView,
  NonbondedTableTabView
)

class RestraintsTabView(GUITabWidget): # Top view
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
    self.chirals= ChiralsTableTabView(parent=self)
    self.addTab(self.chirals, "Chirals")

    # Planarity
    self.planes = PlanarityTableTabView(parent=self)
    self.addTab(self.planes, "Planes")

    # Nonbonded
    self.nonbonded = NonbondedTableTabView(parent=self)
    self.addTab(self.nonbonded, "Non-bonded")