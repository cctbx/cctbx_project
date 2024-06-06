

#from ..view.tabs.cif import add_tabs
from ..cif.cif_browser import CifBrowserController


class RestraintBrowserController(CifBrowserController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    
    #self.state.signals.ciffile_change.connect(self.update_file) # parent already implemented this
    self.state.signals.restraint_activated.connect(self.update_file) # should implement this

    # For restraints, the 'files' dropdown menu will only have one
    self.single_file_ref = None

  def update_file(self,ref):
    # Special callback to set the single file
    self.single_file_ref = ref

    super().update_file(ref)


  @property
  def cif_refs(self):
    return [self.single_file_ref]