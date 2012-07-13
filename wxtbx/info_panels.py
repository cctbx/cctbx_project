
# Assorted info boxes for displaying essential file information

import wxtbx.plots
import wxtbx.icons
import wx
import re
import cStringIO
import os
import sys

class InfoPanelBase (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    kwds['style'] = wx.CAPTION|wx.CLOSE_BOX
    if (wx.Platform == '__WXMSW__') :
      kwds['style'] |= wx.SYSTEM_MENU
    wx.MiniFrame.__init__(self, *args, **kwds)
    self.panel = wx.Panel(self, -1, style=wx.RAISED_BORDER)
    self.panel_sizer = wx.BoxSizer(wx.VERTICAL)
    self.panel.SetSizer(self.panel_sizer)
    self.info_panel = None
    self.file_name = None
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def OnClose (self, event) :
    self.Destroy()

  def OnDestroy (self, event) :
    pass

class ReflectionFileInfo (InfoPanelBase) :
  def __init__ (self, *args, **kwds) :
    super(ReflectionFileInfo, self).__init__(*args, **kwds)
    self.miller_arrays = None
    self._hkl_in = None
    self._current_array = None
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box, 0)
    bmp = wx.StaticBitmap(self.panel, -1, wxtbx.icons.hkl_file.GetBitmap())
    box.Add(bmp, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box.Add(grid, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt1 = wx.StaticText(self.panel, -1, "File name:")
    font = txt1.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    txt1.SetFont(font)
    grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.file_txt = wx.StaticText(self.panel, -1, "(None)")
    grid.Add(self.file_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt3 = wx.StaticText(self.panel, -1, "Data array:")
    txt3.SetFont(font)
    grid.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.array_choice = wx.Choice(self.panel, -1, size=(200,-1))
    grid.Add(self.array_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnChangeArray, self.array_choice)
    self.Centre(wx.BOTH)

  def OnChangeArray (self, event) :
    array_label = self.array_choice.GetStringSelection()
    self.set_miller_array(array_label)

  def set_file (self, file_name) :
    self.file_name = os.path.abspath(file_name)
    from iotbx import file_reader
    self._hkl_in = file_reader.any_file(file_name, force_type="hkl")
    self._hkl_in.assert_file_type("hkl")
    self.SetTitle("Info for %s" % self.file_name)
    self.file_txt.SetLabel(self.file_name)
    self.miller_arrays = self._hkl_in.file_server.miller_arrays
    labels = [ array.info().label_string() for array in self.miller_arrays ]
    self.array_choice.SetItems(labels)
    self.array_choice.SetSelection(0)
    self.set_miller_array(labels[0])

  def set_miller_array (self, array_label) :
    assert (self.file_name is not None) and (self.miller_arrays is not None)
    array = None
    for array_ in self.miller_arrays :
      if (array_.info().label_string() == array_label) :
        array = array_
    assert (array is not None)
    self._current_array = array
    info_out = cStringIO.StringIO()
    array.show_comprehensive_summary(f=info_out)
    array.show_comprehensive_summary()
    info_list = info_out.getvalue().splitlines()
    if (self.info_panel is not None) :
      self.panel_sizer.Detach(self.info_panel)
      self.info_panel.Destroy()
    self.info_panel = wx.Panel(self.panel, -1)
    self.panel_sizer.Add(self.info_panel, 1, wx.ALL|wx.EXPAND)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.info_panel.SetSizer(szr)
    box = wx.StaticBox(self.info_panel, -1, "Array info:")
    box_szr = wx.StaticBoxSizer(box, wx.VERTICAL)
    szr.Add(box_szr, 1, wx.EXPAND|wx.ALL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box_szr.Add(grid, 1, wx.EXPAND)
    for line in info_list[1:] :
      fields = line.split(":")
      label = fields[0]
      value = ":".join(fields[1:])
      txt1 = wx.StaticText(self.info_panel, -1, label + ":")
      font = txt1.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      txt1.SetFont(font)
      grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      txt2 = wx.StaticText(self.info_panel, -1, value)
      font2 = txt2.GetFont()
      font2.SetFamily(wx.FONTFAMILY_MODERN)
      txt2.SetFont(font2)
      grid.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    if (array.is_complex_array()) :
      btn = wx.Button(self.info_panel, -1, "Show map statistics")
      szr.Add(btn, 0, wx.ALL, 5)
      self.Bind(wx.EVT_BUTTON, self.OnShowMapStats, btn)
    szr.Fit(self.info_panel)
    self.panel.Layout()
    self.panel_sizer.Fit(self.panel)
    self.Fit()

  def OnShowMapStats (self, event) :
    map_info = MapCoeffsInfo(self.GetTopLevelParent(), -1,
      "Map statistics for %s" % self._current_array.info().label_string())
    map_info.set_file(self.file_name)
    map_info.set_map_coeffs(self._current_array)
    map_info.Show()

class PDBFileInfo (InfoPanelBase) :
  def __init__ (self, *args, **kwds) :
    super(PDBFileInfo, self).__init__(*args, **kwds)
    self._pdb_in = None
    self._hierarchy = None
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box, 0)
    bmp = wx.StaticBitmap(self.panel, -1, wxtbx.icons.pdb_file.GetBitmap())
    box.Add(bmp, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box.Add(grid, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt1 = wx.StaticText(self.panel, -1, "File name:")
    font = txt1.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    txt1.SetFont(font)
    grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.file_txt = wx.StaticText(self.panel, -1, "(None)")
    grid.Add(self.file_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Centre(wx.BOTH)

  def set_file (self, file_name) :
    self.file_name = os.path.abspath(file_name)
    from iotbx import file_reader
    import iotbx.pdb
    self._pdb_in = file_reader.any_file(file_name, force_type="pdb")
    self._pdb_in.assert_file_type("pdb")
    self._hierarchy = self._pdb_in.file_object.construct_hierarchy()
    info_list = iotbx.pdb.show_file_summary(
      pdb_in=self._pdb_in.file_object,
      hierarchy=self._hierarchy)
    self.SetTitle("Info for %s" % self.file_name)
    self.file_txt.SetLabel(self.file_name)
    if (self.info_panel is not None) :
      self.panel_sizer.Detach(self.info_panel)
      self.info_panel.Destroy()
    self.info_panel = wx.Panel(self.panel, -1)
    self.panel_sizer.Add(self.info_panel, 1, wx.ALL|wx.EXPAND)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.info_panel.SetSizer(szr)
    box = wx.StaticBox(self.info_panel, -1, "File contents:")
    box_szr = wx.StaticBoxSizer(box, wx.VERTICAL)
    szr.Add(box_szr, 1, wx.EXPAND|wx.ALL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box_szr.Add(grid, 1, wx.EXPAND)
    for label, value in info_list :
      txt1 = wx.StaticText(self.info_panel, -1, label + ":")
      font = txt1.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      txt1.SetFont(font)
      grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      str_value = str(value)
      alert = False
      if (str_value.endswith("***")) :
        str_value = re.sub("\s*\*\*\*", "", str_value)
        alert = True
      txt2 = wx.StaticText(self.info_panel, -1, str_value)
      font2 = txt2.GetFont()
      font2.SetFamily(wx.FONTFAMILY_MODERN)
      if (alert) :
        font2.SetWeight(wx.FONTWEIGHT_BOLD)
        txt2.SetForegroundColour((200,0,0))
        txt1.SetForegroundColour((200,0,0))
      txt2.SetFont(font2)
      grid.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    if (len(self._hierarchy.models()) == 1) :
      if (len(self._hierarchy.models()[0].chains()) > 1) :
        btn = wx.Button(self.info_panel, -1, "B-factors by chain...")
        szr.Add(btn, 0, wx.ALL, 5)
        self.Bind(wx.EVT_BUTTON, self.OnShowChainBstats, btn)
    szr.Fit(self.info_panel)
    self.panel.Layout()
    self.panel_sizer.Fit(self.panel)
    self.Fit()

  def OnShowChainBstats (self, event) :
    assert (self._hierarchy is not None)
    b_panel = PDBChainBisoPanel(self)
    b_panel.set_file(
      file_name=self.file_name,
      hierarchy=self._hierarchy)
    b_panel.Show()

class PDBChainBisoPanel (InfoPanelBase) :
  def __init__ (self, *args, **kwds) :
    super(PDBChainBisoPanel, self).__init__(*args, **kwds)
    self._hierarchy = None
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box, 0)
    bmp = wx.StaticBitmap(self.panel, -1, wxtbx.icons.pdb_file.GetBitmap())
    box.Add(bmp, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box.Add(grid, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt1 = wx.StaticText(self.panel, -1, "File name:")
    font = txt1.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    txt1.SetFont(font)
    grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.file_txt = wx.StaticText(self.panel, -1, "(None)")
    grid.Add(self.file_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Centre(wx.BOTH)

  def set_file (self, file_name, hierarchy=None) :
    self.file_name = os.path.abspath(file_name)
    from scitbx.array_family import flex
    if (hierarchy is None) :
      from iotbx import file_reader
      import iotbx.pdb
      pdb_in = file_reader.any_file(file_name, force_type="pdb")
      pdb_in.assert_file_type("pdb")
      hierarchy = pdb_in.file_object.construct_hierarchy()
    if (len(hierarchy.models()) > 1) :
      raise Sorry("Multi-MODEL PDB files not supported.")
    self._hierarchy = hierarchy
    self.SetTitle("B-factors by chain for %s" % self.file_name)
    self.file_txt.SetLabel(self.file_name)
    chain_list = wx.ListCtrl(self.panel, -1, style=wx.LC_REPORT, size=(480,160))
    chain_list.InsertColumn(0, "Chain info")
    chain_list.InsertColumn(1, "Mean B-iso (range)")
    chain_list.SetColumnWidth(0, 260)
    chain_list.SetColumnWidth(1, 200)
    for chain in hierarchy.models()[0].chains() :
      n_res = len(chain.residue_groups())
      chain_atoms = chain.atoms()
      n_atoms = len(chain_atoms)
      main_conf = chain.conformers()[0]
      chain_type = "other"
      if (main_conf.is_protein()) :
        chain_type = "protein"
      elif (main_conf.is_na()) :
        chain_type = "nucleic acid"
      chain_info = "'%s' (%s, %d res., %d atoms)" % (chain.id, chain_type,
        n_res, n_atoms)
      b_iso = chain_atoms.extract_b()
      b_max = flex.max(b_iso)
      b_min = flex.min(b_iso)
      b_mean = flex.mean(b_iso)
      b_info = "%.2f (%.2f - %.2f)" % (b_mean, b_min, b_max)
      item = chain_list.InsertStringItem(sys.maxint, chain_info)
      chain_list.SetStringItem(item, 1, b_info)
    self.panel_sizer.Add(chain_list, 1, wx.EXPAND|wx.ALL, 5)
    self.panel.Layout()
    self.panel_sizer.Fit(self.panel)
    self.Fit()

class ImageFileInfo (InfoPanelBase) :
  def __init__ (self, *args, **kwds) :
    super(ImageFileInfo, self).__init__(*args, **kwds)
    self._img_in = None
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box, 0)
    bmp = wx.StaticBitmap(self.panel, -1, wxtbx.icons.img_file.GetBitmap())
    box.Add(bmp, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box.Add(grid, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt1 = wx.StaticText(self.panel, -1, "File name:")
    font = txt1.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    txt1.SetFont(font)
    grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.file_txt = wx.StaticText(self.panel, -1, "(None)")
    grid.Add(self.file_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Centre(wx.BOTH)

  def set_file (self, file_name) :
    self.file_name = os.path.abspath(file_name)
    from iotbx import file_reader
    img_in = file_reader.any_file(
      file_name=file_name,
      valid_types=["img"],
      force_type="img")
    img_in.assert_file_type("img")
    self._img_in = img_in
    out = cStringIO.StringIO()
    img_in.file_object.show_header()
    img_in.file_object.show_header(out=out)
    self.SetTitle("Info for %s" % self.file_name)
    self.file_txt.SetLabel(self.file_name)
    if (self.info_panel is not None) :
      self.panel_sizer.Detach(self.info_panel)
      self.info_panel.Destroy()
    self.info_panel = wx.Panel(self.panel, -1)
    self.panel_sizer.Add(self.info_panel, 1, wx.ALL|wx.EXPAND)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.info_panel.SetSizer(szr)
    box = wx.StaticBox(self.info_panel, -1, "File contents:")
    box_szr = wx.StaticBoxSizer(box, wx.VERTICAL)
    szr.Add(box_szr, 1, wx.EXPAND|wx.ALL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box_szr.Add(grid, 1, wx.EXPAND)
    for line in out.getvalue().splitlines()[1:] :
      label, value = line.strip().split(":")
      txt1 = wx.StaticText(self.info_panel, -1, label + ":")
      font = txt1.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      txt1.SetFont(font)
      grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      txt2 = wx.StaticText(self.info_panel, -1, value)
      font2 = txt2.GetFont()
      font2.SetFamily(wx.FONTFAMILY_MODERN)
      txt2.SetFont(font2)
      grid.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr.Fit(self.info_panel)
    self.panel.Layout()
    self.panel_sizer.Fit(self.panel)
    self.Fit()

class MapCoeffsInfo (InfoPanelBase) :
  def __init__ (self, *args, **kwds) :
    super(MapCoeffsInfo, self).__init__(*args, **kwds)
    self._current_array = None
    self._current_map = None
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box, 0)
    bmp = wx.StaticBitmap(self.panel, -1, wxtbx.icons.hkl_file.GetBitmap())
    box.Add(bmp, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box.Add(grid, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt1 = wx.StaticText(self.panel, -1, "File name:")
    font = txt1.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    txt1.SetFont(font)
    grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.file_txt = wx.StaticText(self.panel, -1, "(None)")
    grid.Add(self.file_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt3 = wx.StaticText(self.panel, -1, "Map coefficients:")
    txt3.SetFont(font)
    grid.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.map_txt = wx.StaticText(self.panel, -1, "(None)")
    grid.Add(self.map_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def set_file (self, file_name) :
    self.file_txt.SetLabel(file_name)

  def set_map_coeffs (self, array) :
    assert (array.is_complex_array())
    self.map_txt.SetLabel(str(array.info().label_string()))
    self._current_array = array
    map = array.fft_map(resolution_factor=0.25)
    self._current_map = map
    self.info_panel = wx.Panel(self.panel, -1)
    self.panel_sizer.Add(self.info_panel, 1, wx.ALL|wx.EXPAND)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.info_panel.SetSizer(szr)
    box = wx.StaticBox(self.info_panel, -1, "Map info:")
    box_szr = wx.StaticBoxSizer(box, wx.VERTICAL)
    szr.Add(box_szr, 1, wx.EXPAND|wx.ALL, 5)
    grid = wx.FlexGridSizer(cols=2)
    box_szr.Add(grid, 0, wx.EXPAND)
    import iotbx.map_tools
    info_list = iotbx.map_tools.get_map_summary(map)
    for label, value in info_list :
      txt1 = wx.StaticText(self.info_panel, -1, label + ":")
      font = txt1.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      txt1.SetFont(font)
      grid.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      str_value = str(value)
      alert = False
      if (str_value.endswith("***")) :
        str_value = re.sub("\s*\*\*\*", "", str_value)
        alert = True
      txt2 = wx.StaticText(self.info_panel, -1, str_value)
      font2 = txt2.GetFont()
      font2.SetFamily(wx.FONTFAMILY_MODERN)
      if (alert) :
        font2.SetWeight(wx.FONTWEIGHT_BOLD)
        txt2.SetForegroundColour((200,0,0))
        txt1.SetForegroundColour((200,0,0))
      txt2.SetFont(font2)
      grid.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    histogram = wxtbx.plots.histogram(self.info_panel, figure_size=(6,4))
    histogram.show_histogram(
      data=map.real_map(False).as_1d().as_numpy_array(),
      n_bins=20,#n_bins,
      reference_value=0.0,
      x_label="Sigma level",
      y_label="Number of points",
      title="Distribution of grid values (spacing=d_min/4)")
    box_szr.Add(histogram, 1, wx.EXPAND, 5)
    btn = wx.Button(self.info_panel, -1, "Save plot to file...")
    self.Bind(wx.EVT_BUTTON, lambda evt: histogram.save_image(), btn)
    box_szr.Add(btn, 0, wx.TOP, 5)
    szr.Fit(self.info_panel)
    self.panel.Layout()
    self.panel_sizer.Fit(self.panel)
    self.Fit()

if (__name__ == "__main__") :
  import libtbx.load_env
  pdb1 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2C30.pdb",
    test=os.path.isfile)
  hkl1 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/partial_refine_001_map_coeffs.mtz",
    test=os.path.isfile)
  img1 = "/net/sunbird/raid1/sauter/jcsg/1vph/data/jcsg/ssrl1/9_2/20040718/TB0723W/11318/11318_2_017.img"
  assert (not None in [pdb1, hkl1])
  app = wx.App(0)
  frame1 = ReflectionFileInfo(None)
  frame1.set_file(hkl1)
  frame1.Show()
  frame2 = PDBFileInfo(None)
  frame2.set_file(pdb1)
  frame2.Show()
  frame2b = PDBChainBisoPanel(None)
  frame2b.set_file(pdb1)
  frame2b.Show()
  if (os.path.exists(img1)) :
    frame3 = ImageFileInfo(None)
    frame3.set_file(img1)
    frame3.Show()
  app.MainLoop()
