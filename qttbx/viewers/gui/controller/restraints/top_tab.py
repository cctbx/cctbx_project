
from PySide2.QtWidgets import QApplication, QMessageBox
from PySide2 import QtCore
from ..controller import Controller
from ...state.geometry import Geometry
from ...state.ref import GeometryRef
from ...state.ref import SelectionRef
from .tables import (
  BondTableController,
  AngleTableController,
  DihedralTableController,
  ChiralTableController,
  PlanarityTableController,
  NonbondedTableController
)
from pathlib import Path
import pandas as pd



class RestraintTableTopTabController(Controller):
  """
  Simple controller to manage sub-tabs
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.bonds = BondTableController(parent=self,view=view.bonds)
    self.angles = AngleTableController(parent=self,view=view.angles)
    self.dihedrals = DihedralTableController(parent=self,view=view.dihedrals)
    self.chirals = ChiralTableController(parent=self,view=view.chirals)
    self.planes = PlanarityTableController(parent=self,view=view.planes)
    self.nonbonded = NonbondedTableController(parent=self,view=view.nonbonded)

    # Start hidden
    self.view.toggle_tab_visible("Bonds",show=False)
    self.view.toggle_tab_visible("Angles",show=False)
    self.view.toggle_tab_visible("Dihedrals",show=False)
    self.view.toggle_tab_visible("Chirals",show=False)
    self.view.toggle_tab_visible("Planes",show=False)
    self.view.toggle_tab_visible("Non-bonded",show=False)