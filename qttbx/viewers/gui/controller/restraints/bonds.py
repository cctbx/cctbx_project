

from ..controller import Controller
from ..scroll_entry import ScrollEntryController
from ...view.widgets.scroll_entry import ScrollEntryView
from ..scroll_list import ScrollableListController
from ...state.ref import RestraintRef

class BondEntryController(ScrollEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    # Ref is a restraints ref. Bonds,Angles, etc do not have their own ref
    super().__init__(parent=parent,view=view,ref=ref)


  def toggle_active_func(self,is_checked):
    # TODO: Move this to data tab?
    if is_checked:
      self.state.active_selection_ref = self.ref.selection_ref
    else:
      pass

class BondListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.state.signals.restraints_change.connect(self.update)

  def update(self,event):
    entry_list = self
    refs= [ref for ref in self.state.references.values() if isinstance(ref,RestraintRef)]
    for ref in refs:
      if ref not in entry_list.refs:
        entry_view = ScrollEntryView()
        entry_controller = BondEntryController(parent=self,view=entry_view,ref=ref)
        entry_controller.view.active_toggle.is_checked = True
        entry_list.add_entry(entry_controller)



class BondTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.entry_list = BondListController(parent=self,view=self.view.entry_list_view)
