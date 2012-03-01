
from wxtbx import reports
import wx
import sys

class SymmetryListFrame (wx.Frame, reports.PDBLinkMixin) :
  def __init__ (self, *args, **kwds) :
    super(SymmetryListFrame, self).__init__(*args, **kwds)
    self.statusbar = self.CreateStatusBar()
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    panel = wx.Panel(self)
    szr.Add(panel, 1, wx.EXPAND)
    pszr = wx.BoxSizer(wx.VERTICAL)
    panel.SetSizer(pszr)
    box = wx.StaticBox(panel, label="Symmetry search results")
    boxszr = wx.StaticBoxSizer(box, wx.VERTICAL)
    pszr.Add(boxszr, 1, wx.ALL|wx.EXPAND, 5)
    box2 = wx.BoxSizer(wx.HORIZONTAL)
    boxszr.Add(box2)
    box2.Add(wx.StaticText(panel, -1, "Input symmetry:"), 0, wx.ALL, 5)
    self.symm_txt = wx.StaticText(panel, -1, "", size=(540,-1))
    box2.Add(self.symm_txt, 0, wx.ALL, 5)
    boxszr.Add(wx.StaticText(panel, -1, "Closest matching PDB entries:"), 0,
      wx.ALL, 5)
    self.symm_list = PDBSymmList(
      parent=self,
      style=wx.LC_REPORT|wx.LC_SINGLE_SEL,
      size=(740,200))
    boxszr.Add(self.symm_list, 1, wx.ALL|wx.EXPAND, 5)
    self.create_pdb_buttons(panel, boxszr)
    szr.Fit(panel)
    self.Layout()
    self.Fit()
    self.SetMinSize(self.GetSize())
    self.web_frame = None

  def SetResults (self, symmetry, hits) :
    sg_txt = str(symmetry.space_group_info())
    uc_txt = "%g %g %g %g %g %g" % symmetry.unit_cell().parameters()
    self.symm_txt.SetLabel("%s (%s)" % (sg_txt, uc_txt))
    self.statusbar.SetStatusText("%d results shown" % len(hits))
    self.symm_list.SetResults(hits)

  def get_pdb_id_for_viewing (self) :
    return self.symm_list.GetSelectedID()

class PDBSymmList (wx.ListCtrl) :
  def __init__ (self, *args, **kwds) :
    super(PDBSymmList, self).__init__(*args, **kwds)
    cols = ["PDB ID", "RMSD", "Volume ratio", "Space group", "Unit cell"]
    widths = [100, 80, 120, 120, 300]
    r, l = wx.LIST_FORMAT_RIGHT, wx.LIST_FORMAT_LEFT
    alignments = [l, r, r, r, r]
    for i, col in enumerate(cols) :
      self.InsertColumn(i, col, alignments[i])
      self.SetColumnWidth(i, widths[i])

  def SetResults (self, results) :
    self.DeleteAllItems()
    for result in results :
      i = self.InsertStringItem(sys.maxsize, result.pdb_id)
      self.SetStringItem(i, 1, "%6.3f" % result.rmsd)
      self.SetStringItem(i, 2, "%g" % result.volume_ratio)
      self.SetStringItem(i, 3, str(result.pdb_symmetry.space_group_info()))
      self.SetStringItem(i, 4, "%g %g %g %g %g %g" %
        result.pdb_symmetry.unit_cell().parameters())

  def GetSelectedID (self) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      return self.GetItem(item, 0).GetText()
    return None

def run (args) :
  from mmtbx.command_line import search_pdb_symmetry
  results = search_pdb_symmetry.run(args=args)
  app = wx.App(0)
  frame = SymmetryListFrame(None, -1, "PDB symmetry search")
  frame.SetResults(results.crystal_symmetry, results.hits)
  frame.Show()
  app.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])
