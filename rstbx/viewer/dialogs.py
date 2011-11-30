
from wxtbx.phil_controls import path, ints
from wxtbx import phil_controls
from wxtbx import icons
import wx
from libtbx import object_oriented_patterns as oop

RSTBX_SELECT_IMAGE_IDS = 1

class SelectDatasetDialog (wx.Dialog) :
  def __init__ (self, *args, **kwds) :
    self._datasets = []
    style = wx.CAPTION
    dlg_style = kwds.get("style", 0)
    kwds['style'] = style
    wx.Dialog.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    szr.Add(szr2, 1, wx.ALL, 5)
    szr3 = wx.BoxSizer(wx.HORIZONTAL)
    szr2.Add(szr3)
    bmp = wx.StaticBitmap(self, -1, icons.img_file.GetBitmap())
    szr3.Add(bmp, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    caption = "Please select a dataset to index.  Most common detector " +\
      "file formats are supported (ADSC, R-AXIS, MAR, Pilatus, CBF, etc.)."
    if (dlg_style & RSTBX_SELECT_IMAGE_IDS) :
      caption += "  If you wish you may specify which frames you want to "+ \
        "use; otherwise the program will attempt to pick sensible defaults."
    caption_txt = wx.StaticText(self, -1, caption)
    caption_txt.Wrap(500)
    szr3.Add(caption_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid = wx.FlexGridSizer(cols=2)
    szr.Add(grid, 1, wx.EXPAND|wx.ALL)
    txt1 = wx.StaticText(self, -1, "Directory:")
    grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.dir_ctrl = path.PathCtrl(
      parent=self,
      style=path.WXTBX_PHIL_PATH_DIRECTORY)
    grid.Add(self.dir_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(phil_controls.EVT_PHIL_CONTROL, self.OnChooseDirectory,
      self.dir_ctrl)
    txt2 = wx.StaticText(self, -1, "Image set:")
    grid.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.stack_ctrl = wx.Choice(
      parent=self,
      size=(400,-1))
    grid.Add(self.stack_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnChooseDataset, self.stack_ctrl)
    if (dlg_style & RSTBX_SELECT_IMAGE_IDS) :
      txt3 = wx.StaticText(self, -1, "Use frames:")
      grid.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self.frame_ctrl = ints.IntsCtrl(
        parent=self,
        size=(400,-1))
      grid.Add(self.frame_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    else :
      self.frame_ctrl = oop.null()
    btn_sizer = wx.StdDialogButtonSizer()
    szr.Add(btn_sizer, 0, wx.ALL|wx.ALIGN_RIGHT, 10)
    cancel_btn = wx.Button(self, wx.ID_CANCEL)
    ok_btn = wx.Button(self, wx.ID_OK)
    btn_sizer.Add(cancel_btn, 0, wx.RIGHT, 5)
    btn_sizer.Add(ok_btn)
    szr.Fit(self)

  def OnOkay (self, event) :
    pass

  def OnChooseDirectory (self, event) :
    dir_name = self.dir_ctrl.GetPhilValue()
    if (dir_name is not None) :
      from iotbx.detectors import identify_dataset
      self._datasets = identify_dataset(dir_name)
      choices = [ d.format() for d in self._datasets ]
      self.stack_ctrl.SetItems(choices)

  def OnChooseDataset (self, event) :
    print self.stack_ctrl.GetSelection()

  def GetDataset (self) :
    if (len(self._datasets) == 0) :
      return None, None
    else :
      i = self.stack_ctrl.GetSelection()
      return self._datasets[i], None

def select_dataset (parent=None,
                    title="Select a dataset",
                    pick_frames=False) :
  style = 0
  if (pick_frames) :
    style |= RSTBX_SELECT_IMAGE_IDS
  dlg = SelectDatasetDialog(
    parent=parent,
    title=title,
    style=style)
  dataset = frames = None
  if (dlg.ShowModal() == wx.ID_OK) :
    dataset, frames = dlg.GetDataset()
  wx.CallAfter(dlg.Destroy)
  if (pick_frames) :
    return dataset, frames
  else :
    return dataset

# regression testing
if (__name__ == "__main__") :
  app = wx.App(0)
  dataset = select_dataset()
  print dataset
