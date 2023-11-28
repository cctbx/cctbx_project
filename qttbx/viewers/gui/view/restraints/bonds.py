from .restraints import RestraintsTab
from ..widgets.scroll_list import ScrollableListView
from ..widgets.scroll_entry import ScrollEntryView


class BondEntryView(ScrollEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)


class BondListView(ScrollableListView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
  


class BondsTabView(RestraintsTab):
  def __init__(self,parent=None,title="Bonds"):
    super().__init__(parent=parent,title=title)
    
    self.entry_list_view = BondListView(self)
    self.layout.addWidget(self.entry_list_view)
    


