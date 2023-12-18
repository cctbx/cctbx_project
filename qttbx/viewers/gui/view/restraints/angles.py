
from .restraints import RestraintsTab
from ..widgets.scroll_list import ScrollableListView
from ..widgets.scroll_entry import ScrollEntryView


from ..widgets.tab import GUITab
from pathlib import Path

class AngleEntryView(ScrollEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)


class AngleListView(ScrollableListView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)



class AnglesTabView(RestraintsTab):
  def __init__(self,parent=None,title="Angles"):
    super().__init__(parent=parent,title=title)

    self.angle_list_view = AngleListView(self)
    self.layout.addWidget(self.angle_list_view)
    self.setLayout(self.layout)
