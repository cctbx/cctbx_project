from __future__ import absolute_import, division, print_function

from wxtbx import reports
import wx
import sys

class BlastFrame(wx.Frame, reports.PDBLinkMixin):
  def __init__(self, *args, **kwds):
    super(BlastFrame, self).__init__(*args, **kwds)
    self.statusbar = self.CreateStatusBar()
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    panel = wx.Panel(self)
    szr.Add(panel, 1, wx.EXPAND)
    pszr = wx.BoxSizer(wx.VERTICAL)
    panel.SetSizer(pszr)
    box = wx.StaticBox(panel, label="BLAST results")
    boxszr = wx.StaticBoxSizer(box, wx.VERTICAL)
    pszr.Add(boxszr, 1, wx.ALL|wx.EXPAND, 5)
    boxszr.Add(wx.StaticText(panel, -1, "Query sequence:"))
    self.seq_ctrl = wx.TextCtrl(panel, style=wx.TE_READONLY|wx.TE_MULTILINE,
      size=(600,100))
    boxszr.Add(self.seq_ctrl, 0, wx.ALL, 5)
    boxszr.Add(wx.StaticText(panel, -1, "BLAST hits in PDB:"))
    self.hit_list = BlastList(panel, style=wx.LC_REPORT|wx.LC_SINGLE_SEL,
      size=(600,240))
    #self.hit_list.SetMinSize((600,240))
    boxszr.Add(self.hit_list, 1, wx.ALL|wx.EXPAND, 5)
    self.create_pdb_buttons(panel, boxszr)
    szr.Fit(panel)
    self.Layout()
    self.Fit()
    self.SetMinSize(self.GetSize())
    self.web_frame = None

  def SetResults(self, query, results):
    assert (len(results) > 0)
    self.statusbar.SetStatusText("%d results shown (top match: %s_%s)" %
      (len(results), results[0].pdb_id, results[0].chain_id))
    self.seq_ctrl.SetValue(query)
    self.hit_list.SetResults(results)

  def get_pdb_id_for_viewing(self):
    return self.hit_list.GetSelectedID()

class BlastList(wx.ListCtrl):
  def __init__(self, *args, **kwds):
    super(BlastList, self).__init__(*args, **kwds)
    cols = ["PDB ID", "Chain ID", "E-value", "Length", "% Identity",
            "% Positives"]
    widths = [100, 80, 100, 100, 100, 100]
    r, l = wx.LIST_FORMAT_RIGHT, wx.LIST_FORMAT_LEFT
    alignments = [l, l, r, r, r, r]
    for i, col in enumerate(cols):
      self.InsertColumn(i, col, alignments[i])
      self.SetColumnWidth(i, widths[i])

  def SetResults(self, results):
    self.DeleteAllItems()
    for result in results :
      i = self.InsertStringItem(sys.maxunicode, result.pdb_id)
      self.SetStringItem(i, 1, result.chain_id)
      self.SetStringItem(i, 2, "%g" % result.evalue)
      self.SetStringItem(i, 3, "%d" % result.length)
      self.SetStringItem(i, 4, "%.2f" % result.identity)
      self.SetStringItem(i, 5, "%.2f" % result.positives)

  def GetSelectedID(self):
    item = self.GetFirstSelected()
    if (item >= 0):
      return self.GetItem(item, 0).GetText()
    return None

if (__name__ == "__main__"):
  from iotbx.file_reader import any_file
  from iotbx.bioinformatics.structure import summarize_blast_output
  seq_file = any_file(sys.argv[1], force_type="seq")
  seq_file.check_file_type("seq")
  seq_objects = seq_file.file_object
  assert (len(seq_objects) == 1)
  sequence = seq_objects[0].sequence
  results = summarize_blast_output(blast_file=sys.argv[2])
  app = wx.App(0)
  frame = BlastFrame(None, -1, "BLAST results")
  frame.SetResults(sequence, results)
  frame.Show()
  app.MainLoop()
