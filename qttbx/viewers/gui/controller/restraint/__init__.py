from ..controller import Controller
from .restraint_browser import RestraintBrowserController
from .restraint_files import RestraintFileListController

class RestraintTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.files = RestraintFileListController(parent=self,view=self.view.files.list_view)
    self.browser = RestraintBrowserController(parent=self,view=self.view.browser)
