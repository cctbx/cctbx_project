
from rstbx.viewer import results_base, indexing, integration
from wxtbx.phil_controls import path, ints
import wx.lib.agw.flatnotebook
import wx
import os
import sys

class ProcessingFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    self.statusbar = self.CreateStatusBar()
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.nb = wx.lib.agw.flatnotebook.FlatNotebook(self)
    self.sizer.Add(self.nb, 1, wx.EXPAND)
    self.nb.SetMinSize((800,40))
    self.start_panel = StartPanel(self.nb)
    self.nb.AddPage(self.start_panel, "Setup")
    self.indexing_panel = indexing.IndexingPanel(self.nb)
    self.nb.AddPage(self.indexing_panel, "Indexing")
    self.integration_panel = integration.IntegrationPanel(self.nb)
    self.nb.AddPage(self.integration_panel, "Integration")
    self.SetSize((800,600))

  def LoadResults (self, dir_name) :
    self.result = results_base.result(dir_name)
    self.indexing_panel.SetIndexingResults(self.result.get_indexing())
    int_results, summaries = self.result.get_integration()
    self.integration_panel.SetResults(int_results, summaries)

  def OnRun (self, evt) :
    pass

class StartPanel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    box = wx.StaticBox(self, -1, "Indexing setup")
    bszr = wx.StaticBoxSizer(box, wx.VERTICAL)
    szr.Add(bszr, 0, wx.ALL|wx.EXPAND, 5)
    grid = wx.FlexGridSizer(cols=2)
    bszr.Add(grid, 0, wx.ALL|wx.EXPAND)
    txt1 = wx.StaticText(self, -1, "First image:")
    grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.file_ctrl = path.PathCtrl(parent=self,
      style=0) #path.WXTBX_PHIL_PATH_VIEW_BUTTON
    grid.Add(self.file_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid.Add((1,1))
    self.info_txt = wx.StaticText(self, -1, "No images loaded.")
    grid.Add(self.info_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt2 = wx.StaticText(self, -1, "Index on images:")
    grid.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.img_ctrl = ints.IntsCtrl(parent=self)
    grid.Add(self.img_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    btn = wx.Button(self, -1, "Run LABELIT...")
    szr.Add(btn, 0, wx.ALL, 10)
    frame = self.GetTopLevelParent()
    frame.Bind(wx.EVT_BUTTON, frame.OnRun, btn)

if (__name__ == "__main__") :
  dir_name = sys.argv[1]
  assert os.path.isdir(dir_name)
  app = wx.App(0)
  frame = ProcessingFrame(None, -1, "LABELIT")
  frame.LoadResults(dir_name)
  frame.Show()
  app.MainLoop()
