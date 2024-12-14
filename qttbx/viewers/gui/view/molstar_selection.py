from pathlib import Path
import qtawesome as qta
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QPushButton,
    QVBoxLayout
)

from .widgets.tab import GUITab
from .molstar_base import MolstarTabView
from .viewer_controls_selection import ViewerControlsSelectionView
from .scroll_list import ScrollableListView
from .models import ModelLikeEntryView

class MolstarSelectionTabView(MolstarTabView):
  def __init__(self,parent=None):
    super().__init__(parent)

    self.viewer_controls = ViewerControlsSelectionView()
    self.layout.addWidget(self.viewer_controls)
    
