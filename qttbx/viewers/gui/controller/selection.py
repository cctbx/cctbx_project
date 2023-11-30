from .models import ModelLikeEntryController
from ..view.tabs.selection import SelectionEntryView
from ..view.widgets import InfoDialog
from .scroll_list import ScrollableListController
from .controller import Controller

class SelectionEntryController(ModelLikeEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)
    self.state.signals.selection_change.connect(self.handle_selection_change)
    self.view.button_info.clicked.connect(self.display_info)


  def handle_selection_change(self):
    # Just disable the toggle if self is not active
    if (self.state.active_selection_ref is  None  or
     self.state.active_selection_ref.id != self.ref.id):
      self.view.active_toggle.is_checked = False

  def toggle_active_func(self,is_checked):
    # TODO: Move this to data tab?
    if is_checked:
      self.state.active_selection_ref = self.ref

    else:
      #print("The entry is unchecked.")
      if self.state.active_selection_ref == self.ref:
        self.state.active_selection_ref = None 

  def display_info(self):
    text = f"""
    Reference id: {self.ref.id}
    Model Reference id: {self.ref.model_ref.id}
    External ids:
    {self.model_ref.external_ids}
    
    Phenix string: {self.ref.query.phenix_string}
    """
    dialog = InfoDialog(text, title="Selection Info", parent=self.view)
    dialog.exec_()

class SelectionListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.next_selection_number = 1 # for labeling new selections

    self.state.signals.selection_change.connect(self.update)


  def update(self):
    selection_list = self
    for ref in self.state.references_selection:
      if ref not in selection_list.refs:
        if ref.show_in_list:
          entry_view = SelectionEntryView()
          entry_controller = SelectionEntryController(parent=self,view=entry_view,ref=ref)
          entry_controller.view.active_toggle.is_checked = True
          selection_list.add_entry(entry_controller)
          print(ref.label)
          ref.label = f"Selection {selection_list.next_selection_number}"
          print(ref.label)
          entry_controller.view.label_name.setText(ref.label)
          selection_list.next_selection_number+=1
          # make new selection active
          


class SelectionTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.selection_list_controller = SelectionListController(parent=self,view=self.view.selection_list_view)

