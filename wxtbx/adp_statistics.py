
from wxtbx import plots
from libtbx.utils import Sorry, null_out
import wx
import os
import sys

master_phil = """
adp_statistics
  .short_caption = ADP statistics
  .caption = This tool displays tables and plots of statistics for isotropic \
    and anisotropic ADPs in a PDB file.  You can get the same information on \
    the command line by running phenix.pdbtools model_statistics=True \
    model.pdb.
  .style = box auto_align caption_img:icons/custom/phenix.pdbtools.png
{
  pdb_file = None
    .type = path
    .optional = False
    .short_caption = PDB file
    .style = bold file_type:pdb input_file
  cif_file = None
    .type = path
    .multiple = True
    .optional = True
    .short_caption = Restraints (CIF) file
    .style = file_type:cif input_file
}
"""

def run (args=(), params=None, out=None) :
  if (out is None) :
    out = sys.stdout
  if (len(args) > 0) :
    assert (params is None)
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil_string=master_phil,
      pdb_file_def="adp_statistics.pdb_file",
      cif_file_def="adp_statistics.cif_file")
    params = cmdline.work.extract()
    validate_params(params)
  import mmtbx.model
  import mmtbx.restraints
  from mmtbx.monomer_library import pdb_interpretation
  processed_pdb_file = pdb_interpretation.run(
    args=[params.adp_statistics.pdb_file] + params.adp_statistics.cif_file,
    log=null_out())
  geometry = processed_pdb_file.geometry_restraints_manager(show_energies=True)
  restraints_manager = mmtbx.restraints.manager(
    geometry = geometry,
    normalization = True)
  model = mmtbx.model.manager(
    xray_structure     = processed_pdb_file.xray_structure(),
    pdb_hierarchy      = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    restraints_manager = restraints_manager,
    log                = out)
  stats = model.adp_statistics()
  stats.file_name = params.adp_statistics.pdb_file
  stats.show_1(out=out)
  return stats

def validate_params (params) :
  if (params.adp_statistics.pdb_file is None) :
    raise Sorry("A PDB file is required to run this program.")
  elif (not os.path.isfile(params.adp_statistics.pdb_file)) :
    raise Sorry("The path %s is not a file or does not exist.")

#-----------------------------------------------------------------------
# GUI objects
class ADPStatisticsFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    super(ADPStatisticsFrame, self).__init__(*args, **kwds)
    s = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(s)
    p = wx.Panel(self)
    s.Add(p, 1, wx.EXPAND)
    s2 = wx.BoxSizer(wx.VERTICAL)
    p.SetSizer(s2)
    s0 = wx.BoxSizer(wx.HORIZONTAL)
    s2.Add(s0)
    txt0 = wx.StaticText(p, -1, "PDB file:")
    s0.Add(txt0, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.file_text = wx.StaticText(p, -1, "--unknown--")
    s0.Add(self.file_text, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5)
    txt1 = wx.StaticText(p, -1, "Overall statistics:")
    s2.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    if (sys.platform == "darwin") :
      size1 = (720,120)
    else :
      size1 = (720,160) # FIXME check this on Linux
    self.list1 = wx.ListCtrl(p, style=wx.LC_REPORT, size=size1)
    self.list1.InsertColumn(0, "Atom type", width=100)
    self.list1.InsertColumn(1, "# iso", width=60, format=wx.LIST_FORMAT_RIGHT)
    self.list1.InsertColumn(2, "# aniso", width=60, format=wx.LIST_FORMAT_RIGHT)
    self.list1.InsertColumn(3, "Iso min.", width=80,format=wx.LIST_FORMAT_RIGHT)
    self.list1.InsertColumn(4, "Iso max.", width=80,format=wx.LIST_FORMAT_RIGHT)
    self.list1.InsertColumn(5, "Iso mean", width=80,format=wx.LIST_FORMAT_RIGHT)
    self.list1.InsertColumn(6, "Aniso min.", width=80,
      format=wx.LIST_FORMAT_RIGHT)
    self.list1.InsertColumn(7, "Aniso max.", width=80,
      format=wx.LIST_FORMAT_RIGHT)
    self.list1.InsertColumn(8, "Aniso mean", width=80,
      format=wx.LIST_FORMAT_RIGHT)
    for i in range(9) : # FIXME this doesn't work!!!
      continue
      col = self.list1.GetColumn(i)
      font = col.GetFont()
      font.SetPointSize(9)
      col.SetFont(font)
    s2.Add(self.list1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt2 = wx.StaticText(p, -1,
      "Distribution of isotropic (or equivalent) ADP for non-H atoms:")
    s2.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    s3 = wx.BoxSizer(wx.HORIZONTAL)
    s3.Add((5,1))
    s2.Add(s3, 0, wx.EXPAND)
    self.list2 = wx.ListCtrl(p, style=wx.LC_REPORT, size=(300,200))
    self.list2.InsertColumn(0, "Bin", width=50)
    self.list2.InsertColumn(1, "Value range", width=180)
    self.list2.InsertColumn(2, "#atoms", width=60, format=wx.LIST_FORMAT_RIGHT)
    for i in range(3) :
      continue
      col = self.list2.GetColumn(i)
      font = col.GetFont()
      font.SetPointSize(9)
      col.SetFont(font)
    s3.Add(self.list2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5)
    self.plot1 = adp_histogram(p, figure_size=(4,3))
    s3.Add(self.plot1, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5)
    txt3 = wx.StaticText(p, -1, "Distribution of anisotropy:")
    s2.Add(txt3)
    s4 = wx.BoxSizer(wx.HORIZONTAL)
    s4.Add((5,1))
    s2.Add(s4, 0, wx.EXPAND)
    self.list3 = wx.ListCtrl(p, style=wx.LC_REPORT, size=(300,200))
    self.list3.InsertColumn(0, "Bin", width=50)
    self.list3.InsertColumn(1, "Value range", width=180)
    self.list3.InsertColumn(2, "#atoms", width=60, format=wx.LIST_FORMAT_RIGHT)
    for i in range(3) :
      continue
      col = self.list3.GetColumn(i)
      font = col.GetFont()
      font.SetPointSize(9)
      col.SetFont(font)
    s4.Add(self.list3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5)
    self.plot2 = adp_histogram(p, figure_size=(4,3))
    s4.Add(self.plot2, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND, 5)
    s.Fit(p)
    self.Fit()

  def show_statistics (self, stats) :
    file_name = getattr(stats, "file_name", None)
    if (file_name is not None) :
      self.file_text.SetLabel(file_name)
    t1, t2, t3 = stats.format_tables()
    for row in t1 :
      i = self.list1.InsertStringItem(sys.maxint, str(row[0]))
      for j, value in enumerate(row[1:]) :
        self.list1.SetStringItem(i, j+1, str(value))
    for row in t2 :
      i = self.list2.InsertStringItem(sys.maxint, str(row[0]))
      for j, value in enumerate(row[1:]) :
        self.list2.SetStringItem(i, j+1, str(value))
    for row in t3 :
      i = self.list3.InsertStringItem(sys.maxint, str(row[0]))
      for j, value in enumerate(row[1:]) :
        self.list3.SetStringItem(i, j+1, str(value))
    p1, p2 = stats.format_plots()
    y1, yrange1 = p1
    y2, yrange2 = p2
    self.plot1.convert_histogram(y1)
    self.plot2.convert_histogram(y2)

class adp_histogram (plots.histogram) :
  def convert_histogram (self, bins) :
    # XXX flex.histogram has already binned values, so I make new fake values
    # that matplotlib then re-bins
    values = []
    for i, x in enumerate(bins) :
      values.extend([i] * x)
    self.show_histogram(values, len(bins), y_label="# of atoms")

if (__name__ == "__main__") :
  stats = run(sys.argv[1:])
  app = wx.App(0)
  frame = ADPStatisticsFrame(None, -1, "ADP statistics")
  frame.show_statistics(stats)
  frame.Show()
  app.MainLoop()
