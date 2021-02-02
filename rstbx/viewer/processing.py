from __future__ import absolute_import, division, print_function

from rstbx.viewer import dataset, results_base, indexing, integration
from rstbx.viewer.frame import XrayFrame
from wxtbx import process_control, icons
import wxtbx.app
from wxtbx.phil_controls import path
import wx.lib.agw.flatnotebook
import wx
import os
import sys

class ProcessingFrame(wx.Frame):
  def __init__(self, *args, **kwds):
    wx.Frame.__init__(self, *args, **kwds)
    self.viewer = None
    self.toolbar = self.CreateToolBar(style=wx.TB_TEXT)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Image viewer",
      bitmap=icons.hkl_file.GetBitmap(),
      shortHelp="Image viewer",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnLaunchViewer, btn)
    self.toolbar.Realize()
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

  def LoadResults(self, dir_name):
    self.result = results_base.result(dir_name)
    self.indexing_panel.SetIndexingResults(self.result.get_indexing())
    self.integration_panel.SetResults(self.result)
    self.nb.SetSelection(1)

  def OnRunIndexing(self, evt):
    dataset, frames = self.start_panel.GetDataset()
    output_dir = self.start_panel.GetOutputDir()
    result = self.run_indexing(
      dataset=dataset,
      frames=frames,
      output_dir=output_dir)
    self.LoadResults(output_dir)

  def run_indexing(self, **kwds):
    from rstbx.viewer import drivers
    run = drivers.run_indexing(**kwds)
    indexing_result = process_control.run_function_as_process_in_dialog(
      parent=self,
      thread_function=run,
      title="Running LABELIT",
      message="Indexing images and performing simple test integration")
    return indexing_result

  def launch_viewer_frame(self):
    if (self.viewer is None):
      self.viewer = XrayFrame(
        parent=self,
        title="Image viewer")
      self.viewer.Show()
      self.Bind(wx.EVT_CLOSE, self.OnCloseViewer, self.viewer)

  def get_viewer_frame(self):
    self.launch_viewer_frame()
    return self.viewer

  def set_viewer_frame(self, frame):
    assert (self.viewer is None)
    self.viewer = frame

  def OnCloseViewer(self, evt):
    self.viewer.Destroy()
    self.viewer = None

  def OnLaunchViewer(self, evt):
    self.launch_viewer_frame()

class StartPanel(wx.Panel, dataset.SelectDatasetPanelMixin):
  def __init__(self, *args, **kwds):
    wx.Panel.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    box = wx.StaticBox(self, -1, "Indexing setup")
    bszr = wx.StaticBoxSizer(box, wx.VERTICAL)
    szr.Add(bszr, 1, wx.ALL|wx.EXPAND, 5)
    self.draw_dataset_controls(bszr)
    btn = wx.Button(self, -1, "Run LABELIT...")
    szr.Add(btn, 0, wx.ALL, 10)
    frame = self.GetTopLevelParent()
    frame.Bind(wx.EVT_BUTTON, frame.OnRunIndexing, btn)

  def add_controls_to_grid(self, sizer):
    txt = wx.StaticText(self, -1, "Output directory:")
    self.output_ctrl = path.PathCtrl(
      parent=self,
      style=path.WXTBX_PHIL_PATH_DIRECTORY)
    self.output_ctrl.SetValue(os.getcwd())
    sizer.Add(txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    sizer.Add(self.output_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def GetOutputDir(self):
    return self.output_ctrl.GetPhilValue()

if (__name__ == "__main__"):
  app = wxtbx.app.CCTBXApp(0)
  frame = ProcessingFrame(None, -1, "LABELIT")
  frame.Show()
  if (len(sys.argv) > 1) and (os.path.isdir(sys.argv[1])):
    frame.LoadResults(sys.argv[1])
  app.MainLoop()
