


from ..controller import Controller
from .cif_files import CifFileListController
from .cif_browser import CifBrowserController

class CifTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    #self.files = CifFileListController(parent=self,view=self.view.files.list_view)
    self.browser = CifBrowserController(parent=self,view=self.view.browser)
