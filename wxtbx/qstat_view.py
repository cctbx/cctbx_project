
from libtbx import sge_utils
import libtbx.load_env
from libtbx.utils import Sorry
import wx
import sys, os, time

if sys.platform in ["linux", "linux2", "darwin"] :
  user = os.environ['USER']
else :
  user = None # is it even worth bothering with Windows here?

job_attrs = ["job_id", "state", "name", "user", "submit", "queue"]
job_labels = ["Job ID", "Status", "Name", "User", "Start time", "Queue"]
status_codes = ["d", "E", "h", "r", "R", "s", "S", "t", "T", "w"]
status_imgs = [3, 4, 0, 1, 0, 0, 0, 0, 0, 2]
col_sizes = [wx.LIST_AUTOSIZE] * 4 + [200,200]
icon_lib = os.path.join("phenix","wxGUI2","images","icons","crystal")

class qsub_list_data (object) :
  def __init__ (self) :
    self._data = []
    self._sortby = None
    self._sort_descending = False

  def Refresh (self) :
    self._data = sge_utils.qstat_parse()
    if self._sortby is not None :
      self.SortItems(self._sortby, swap_order=False)

  def GetItemCount (self) :
    return len(self._data)

  def GetItemText (self, item, col) :
    return getattr(self._data[item], job_attrs[col])

  def GetItemImage (self, item) :
    status = self._data[item].state[-1]
    img_id = status_imgs[status_codes.index(status)]
    return img_id

  def SortItems (self, col, swap_order=True) :
    if swap_order :
      if self._sortby == col :
        self._sort_descending = (not self._sort_descending)
      else :
        self._sort_descending = False
    if col == 0 :
      self._data.sort(lambda x, y: cmp(int(x.job_id), int(y.job_id)))
    elif col == 4 :
      fmt = "%m/%d/%Y %H:%M:%S"
      self._data.sort(lambda x, y: cmp(time.strptime(x.submit, fmt),
        time.strptime(y.submit, fmt)))
    else :
      attr = job_attrs[col]
      self._data.sort(lambda x, y: cmp(getattr(x, attr), getattr(y, attr)))
    if self._sort_descending :
      self._data.reverse()
    self._sortby = col

  def GetOwners (self, job_ids) :
    names = []
    for job in self._data :
      if job.job_id in job_ids :
        names.append(job.user)
    return set(names)

class qsub_list_view (wx.ListCtrl) :
  def __init__ (self, *args, **kwds) :
    wx.ListCtrl.__init__(self, *args, **kwds)
    self.SetupImages()
    self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect, self)
    self.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnDeSelect, self)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick, self)
    self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightClick, self)
    self.Bind(wx.EVT_LIST_COL_CLICK, self.OnSort, self)
    for i, label in enumerate(job_labels) :
      self.InsertColumn(i, label)
      self.SetColumnWidth(i, col_sizes[i]) #wx.LIST_AUTOSIZE)
    self.dataSource = qsub_list_data()
    self.RefreshAllItems()

  def SetupImages (self) :
    icon_dir = libtbx.env.find_in_repositories(relative_path=icon_lib,
      test=os.path.isdir)
    if icon_dir is None :
      return
    il = wx.ImageList(16, 16, True)
    for icon in ["blank", "run", "recur", "stop", "status_unknown"] :
      icon_path = os.path.join(icon_dir, "16x16", "actions", icon + ".png")
      img = wx.Image(icon_path, type=wx.BITMAP_TYPE_ANY, index=-1)
      il.Add(img.ConvertToBitmap())
    self.AssignImageList(il, wx.IMAGE_LIST_SMALL)

  def OnGetItemImage (self, item) :
    return self.dataSource.GetItemImage(item)

  def OnGetItemAttr (self, item) :
    pass

  def OnGetItemText (self, item, col) :
    return self.dataSource.GetItemText(item, col)

  def RefreshAllItems (self) :
    n_items = self.dataSource.GetItemCount()
    self.SetItemCount(n_items)
    self.RefreshItems(0, n_items)

  def GetSelectedJobIDs (self) :
    jobs = []
    item = self.GetFirstSelected()
    while item != -1 :
      jobs.append(self.dataSource.GetItemText(item, 0))
      item = self.GetNextSelected(item)
    return jobs

  def GetOwners (self, job_ids) :
    return self.dataSource.GetOwners(job_ids)

  def OnSelect (self, event) :
    pass

  def OnDeSelect (self, event) :
    pass

  def OnDoubleClick (self, event) :
    pass

  def OnRightClick (self, event) :
    pass

  def OnSort (self, event) :
    col = event.GetColumn()
    self.dataSource.SortItems(col)
    self.RefreshAllItems()

  def Update (self) :
    self.dataSource.Refresh()
    self.RefreshAllItems()

class queue_list_frame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.statusbar = self.CreateStatusBar()
    self.SetupToolbar()
    self.list_ctrl = qsub_list_view(parent=self,
      id=-1,
      size=(800,600),
      style=wx.LC_REPORT|wx.LC_VIRTUAL)
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)
    self.Update()
    self._timer = wx.Timer(owner=self)
    self.Bind(wx.EVT_TIMER, self.OnUpdate)
    self.statusbar.Bind(wx.EVT_LEFT_DCLICK, self.OnUpdate)
    self._timer.Start(10000)

  def SetupToolbar (self) :
    icon_dir = libtbx.env.find_in_repositories(relative_path=icon_lib,
      test=os.path.isdir)
    if icon_dir is None :
      return
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    commands = [
      (["32x32","actions","reload.png"], "OnUpdate", "Refresh list"),
      (["32x32","actions","editdelete.png"], "OnDelete", "Delete selected"),
    ]
    for (icon_paths, fname, label) in commands :
      full_path = os.path.join(icon_dir, *icon_paths)
      img = wx.Image(full_path, type=wx.BITMAP_TYPE_PNG, index=-1)
      tool_button = self.toolbar.AddLabelTool(-1, label, img.ConvertToBitmap(),
        shortHelp=label, kind=wx.ITEM_NORMAL)
      self.Bind(wx.EVT_MENU, getattr(self, fname), tool_button)
    self.toolbar.Realize()

  def OnClose (self, event) :
    self.Destroy()

  def OnDestroy (self, event) :
    pass

  def OnUpdate (self, event) :
    self.Update()

  def OnDelete (self, event) :
    job_ids = self.list_ctrl.GetSelectedJobIDs()
    if len(job_ids) == 0 :
      return
    print "qdel %s" % " ".join(job_ids)
    print self.list_ctrl.GetOwners(job_ids)

  def SetUpdateInterval (self, interval) : # in seconds, not ms
    self._timer.Stop()
    self._timer.Start(interval * 1000)

  def Update (self) :
    self.list_ctrl.Update()
    self.statusbar.SetStatusText("Last updated at %s" % get_time())

def get_time () :
  return time.strftime("%m-%d-%y %H:%M:%S", time.localtime())

#-----------------------------------------------------------------------
def run (args) :
  app = wx.App(0)
  frame = queue_list_frame(None, -1, "SGE Queue Status")
  frame.Show()
  frame.Fit()
  app.MainLoop()

if __name__ == "__main__" :
  run(sys.argv[1:])
