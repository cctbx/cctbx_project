
from PySide2.QtCore import Signal
from .restraints_staging import RestraintStagingTabView
from .restraints_table import (
  BondsTableTabView,
  AnglesTableTabView,
  DihedralsTableTabView,
  ChiralsTableTabView,
  PlanarityTableTabView,
  NonbondedTableTabView
)
from ..widgets.tab import GUITabWidget



class RestraintsTabView(GUITabWidget): # Top view
  def __init__(self,parent=None):
    super().__init__(parent=parent)


    # staging
    #self.staging = RestraintStagingTabView(parent=self)
    #self.addTab(self.staging, "Staging")


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