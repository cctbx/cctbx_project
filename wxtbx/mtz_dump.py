
from __future__ import absolute_import, division, print_function
import wxtbx.app
from libtbx.utils import Sorry
import wx
import sys

columns = [ "Order", "Label", "#valid", "%valid", "min", "max", "type", ]
column_widths = [ 60, 100, 80, 80, 100, 100, 100 ]

aln_flags = wx.ALL|wx.ALIGN_CENTER_VERTICAL

class MtzInspectionFrame(wx.Frame):
  def __init__(self, *args, **kwds):
    wx.Frame.__init__(self, *args, **kwds)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.panel = MtzContentsPanel(self)
    self.sizer.Add(self.panel, 1, wx.ALL|wx.EXPAND, 0)
    self.sizer.Fit(self.panel)
    self.Fit()

  def SetMtzFile(self, *args, **kwds):
    self.panel.SetMtzFile(*args, **kwds)
    self.sizer.Fit(self.panel)
    self.Fit()

  def OnOpen(self, event):
    from wxtbx import path_dialogs
    file_name = path_dialogs.manager().select_file(
      parent=self,
      message="Choose an MTZ file to view",
      wildcard="MTZ files (*.mtz)|*.mtz")
    if (file_name is not None):
      self.SetMtzFile(file_name)

class MtzContentsPanel(wx.Panel):
  def __init__(self, *args, **kwds):
    wx.Panel.__init__(self, *args, **kwds)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    grid1 = wx.FlexGridSizer(cols=2)
    lbl1 = wx.StaticText(self, label="Title:")
    self.title_txt = wx.StaticText(self, label="", size=(300,-1))
    lbl2 = wx.StaticText(self, label="Space group:")
    self.sg_txt = wx.StaticText(self, label="", size=(300,-1))
    lbl3 = wx.StaticText(self, label="Resolution:")
    self.d_max_min_txt = wx.StaticText(self, label="", size=(300,-1))
    lbl4 = wx.StaticText(self, label="Select dataset:")
    self.dataset_chooser = wx.Choice(self, size=(400,-1))
    self.Bind(wx.EVT_CHOICE, self.OnChooseDataset, self.dataset_chooser)
    grid1.Add(lbl1, 0, aln_flags, 5)
    grid1.Add(self.title_txt, 0, aln_flags, 5)
    grid1.Add(lbl2, 0, aln_flags, 5)
    grid1.Add(self.sg_txt, 0, aln_flags, 5)
    grid1.Add(lbl3, 0, aln_flags, 5)
    grid1.Add(self.d_max_min_txt, 0, aln_flags, 5)
    grid1.Add(lbl4, 0, aln_flags, 5)
    grid1.Add(self.dataset_chooser, 0, aln_flags, 5)
    self.sizer.Add(grid1, 0, wx.ALL|wx.EXPAND)
    self._dataset_panels = []
    self._crystals_and_datasets = []
    self._dataset_labels = []
    self._mtz_obj = None

  def SetMtzFile(self, file_name):
    from iotbx import mtz
    try :
      self._mtz_obj = mtz.object(file_name=file_name)
    except RuntimeError as e :
      raise Sorry(("The file '%s' could not be read as an MTZ file "+
        "(original error: %s)") % (file_name, str(e)))
    self.title_txt.SetLabel(self._mtz_obj.title())
    self.sg_txt.SetLabel(str(self._mtz_obj.space_group_name()))
    self.d_max_min_txt.SetLabel("%g - %g Angstrom" %
      self._mtz_obj.max_min_resolution())
    self._dataset_labels = []
    self._crystals_and_datasets = []
    for i_crystal, crystal in enumerate(self._mtz_obj.crystals()):
      if (crystal.name() == "HKL_base"):
        continue
      for i_dataset, dataset in enumerate(crystal.datasets()):
        label = "/%s/%s" % (crystal.name(), dataset.name())
        self._crystals_and_datasets.append((crystal, dataset))
        self._dataset_labels.append(label)
    self.dataset_chooser.SetItems(self._dataset_labels)
    p = MtzDatasetPanel(self, style=wx.RAISED_BORDER)
    self._dataset_panels.append(p)
    self.sizer.Add(p, 1, wx.ALL|wx.EXPAND, 0)
    if (len(self._dataset_labels) > 0):
      self.dataset_chooser.SetSelection(0)
      self.OnChooseDataset(None)

  def OnChooseDataset(self, event):
    if (len(self._dataset_panels) == 0) : return
    sel = self.dataset_chooser.GetSelection()
    crystal, dataset = self._crystals_and_datasets[sel]
    p = self._dataset_panels[0]
    p.SetMtzDataset(
      crystal=crystal,
      dataset=dataset,
      n_refl=self._mtz_obj.n_reflections())
    self.Layout()
    self.sizer.Fit(self)
    self.Fit()

class MtzDatasetPanel(wx.Panel):
  def __init__(self, *args, **kwds):
    wx.Panel.__init__(self, *args, **kwds)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    grid = wx.FlexGridSizer(cols=2)
    lbl1 = wx.StaticText(self, label="Unit cell:")
    self.uc_txt = wx.StaticText(self, label="", size=(300,-1))
    lbl2 = wx.StaticText(self, label="Wavelength:")
    self.wl_txt = wx.StaticText(self, label="", size=(300,-1))
    grid.Add(lbl1, 0, aln_flags, 5)
    grid.Add(self.uc_txt, 0, aln_flags, 5)
    grid.Add(lbl2, 0, aln_flags, 5)
    grid.Add(self.wl_txt, 0, aln_flags, 5)
    self.sizer.Add(grid, 0, wx.ALL|wx.EXPAND)
    self.lc = MtzColumnList(self, -1, size=(720,240), style=wx.LC_REPORT)
    self.sizer.Add(self.lc, 1, wx.ALL|wx.EXPAND, 5)
    #self.sizer.Fit(self)
    self.Fit()

  def SetMtzDataset(self, crystal, dataset, n_refl):
    self.lc.DeleteAllItems()
    self.lc.SetNReflections(n_refl)
    self.uc_txt.SetLabel("%g %g %g %g %g %g" % crystal.unit_cell().parameters())
    self.wl_txt.SetLabel("%g" % dataset.wavelength())
    self.lc.AddMtzDataset(dataset)
    self.Layout()
    self.Refresh()

class MtzColumnList(wx.ListCtrl):
  def __init__(self, *args, **kwds):
    style = kwds.get('style', 0)
    if (not style & wx.LC_REPORT):
      style &= wx.LC_REPORT
    kwds['style'] = style
    wx.ListCtrl.__init__(self, *args, **kwds)
    self.n_refl = None
    for i, label in enumerate(columns):
      self.InsertColumn(i, label)
      self.SetColumnWidth(i, column_widths[i])

  def SetNReflections(self, n_refl):
    self.n_refl = n_refl

  def AddMtzColumn(self, fields):
    assert (len(fields) == len(columns))
    n = self.GetItemCount() + 1
    item = self.InsertStringItem(sys.maxunicode, str(n))
    for i, field in enumerate(fields[:-2]):
      self.SetStringItem(item, i+1, field)
    self.SetStringItem(item, len(fields)-1, "%s %s" % (fields[-2], fields[-1]))

  def AddMtzDataset(self, dataset):
    assert (self.n_refl is not None)
    for i_col, column in enumerate(dataset.columns()):
      fields = column.format_fields_for_mtz_dump(self.n_refl)
      self.AddMtzColumn(fields)

if (__name__ == "__main__"):
  app = wxtbx.app.CCTBXApp(0)
  frame = MtzInspectionFrame(None, title="Inspect MTZ file contents")
  if (len(sys.argv) == 0):
    frame.OnOpen(None)
  else :
    frame.SetMtzFile(sys.argv[1])
  frame.Show()
  app.MainLoop()
