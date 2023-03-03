
"""
Classes for displaying Xtriage results using wxPython (and matplotlib).
"""

from __future__ import absolute_import, division, print_function
from wxtbx import metallicbutton
import wxtbx.misc_dialogs
import wxtbx.path_dialogs
import wxtbx.windows
import wxtbx.tables
import wxtbx.plots
import mmtbx.scaling
import wx.lib.scrolledpanel
from libtbx.utils import Sorry, to_unicode
from libtbx import easy_pickle, easy_run
import wx
import os.path
import sys
from six.moves import range

TEXT_WIDTH = 800

class wx_panel(wx.lib.scrolledpanel.ScrolledPanel):
  def OnChildFocus(self, evt):
    pass

class wx_output_base(mmtbx.scaling.xtriage_output):
  """
  Implements wxPython controls for Xtriage output handler, independent of
  window type.
  """
  gui_output = True
  def __init__(self, *args, **kwds):
    super(wx_output_base, self).__init__(*args, **kwds)
    self._current_panel = None
    self._current_sizer = None
    self._sections = {}
    self._graphs = []
    self._tables = []
    self._in_box = False

  def create_panel(self):
    raise NotImplementedError()

  def add_panel(self, panel, title):
    raise NotImplementedError()

  def show_big_header(self, text) : pass

  def show_header(self, title):
    """
    Creates a new notebook page with the specified title.  This will be the
    output panel for all subsequent method calls until show_header() is called
    again.
    """
    panel = self.create_panel()
    szr = wx.BoxSizer(wx.VERTICAL)
    panel.SetSizer(szr)
    self._current_panel = panel
    self._current_sizer = szr
    self.add_panel(panel, title)

  def show_sub_header(self, title):
    """
    Create a wx.StaticBox
    """
    assert (self._current_panel is not None)
    if (self._in_box):
      self._current_sizer = self._current_panel.GetSizer()
    box = wx.StaticBox(parent=self._current_panel,
      label=title,
      style=wx.NO_BORDER)
    box_szr = wx.StaticBoxSizer(box, wx.VERTICAL)
    self._current_sizer.Add(box_szr, 0, wx.ALL|wx.EXPAND, 5)
    self._current_sizer = box_szr
    self._sections[title] = box
    self._in_box = True

  def show_text(self, text):
    """
    Create wx.StaticText object(s) with automatic wrapping.  Double-newlines
    will be treated as paragraph breaks, otherwise newlines are replaced by
    spaces.
    """
    assert (self._current_panel is not None)
    blocks = text.split("\n\n")
    for block in blocks :
      block = " ".join([ l.strip() for l in block.splitlines() ]).strip()
      wxtxt = wx.StaticText(parent=self._current_panel,
                            label=to_unicode(block))
      wxtxt.Wrap(TEXT_WIDTH)
      self._current_sizer.Add(wxtxt, 0, wx.ALL, 5)

  def show_preformatted_text(self, text):
    """
    Draw monospaced text, preserving the original formatting.
    """
    assert (self._current_panel is not None)
    if text.startswith("\n"):
      text = text[1:]
    wx_txt = wx.StaticText(parent=self._current_panel,
      label=text)
    font = wx_txt.GetFont()
    font2 = wx.Font(font.GetPointSize(), wx.FONTFAMILY_MODERN, wx.NORMAL,
      wx.FONTWEIGHT_NORMAL) #  Python 3 Fix  , face="Courier")
    # FIXME this seems not to work on wxPython 3/Mac OS 10.9
    wx_txt.SetFont(font2)
    self._current_sizer.Add(wx_txt, 0, wx.ALL, 5)

  def warn(self, text):
    """
    Create wx.StaticText object(s) with automatic wrapping.  Font will be
    boldface with red foreground.
    """
    assert (self._current_panel is not None)
    blocks = text.split("\n\n")
    for block in blocks :
      block = " ".join([ l.strip() for l in block.splitlines() ]).strip()
      wxtxt = wx.StaticText(parent=self._current_panel,
        label=block)
      font = wxtxt.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      wxtxt.SetFont(font)
      wxtxt.SetForegroundColour((240,0,0))
      wxtxt.Wrap(TEXT_WIDTH)
      self._current_sizer.Add(wxtxt, 0, wx.ALL, 5)

  def show_lines(self, text):
    assert (self._current_panel is not None)
    wxtxt = wx.StaticText(parent=self._current_panel,
      label=text)
    wxtxt.Wrap(TEXT_WIDTH)
    self._current_sizer.Add(wxtxt, 0, wx.ALL, 5)

  def show_paragraph_header(self, text):
    """
    Draws left-aligned boldface text in a slightly larger font.
    """
    assert (self._current_panel is not None)
    wx_text = wx.StaticText(
      parent=self._current_panel,
      label=text)
    font = wx_text.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    font.SetPointSize(14)
    wx_text.SetFont(font)
    self._current_sizer.Add(wx_text, 0, wx.ALL, 5)

  def show_table(self, table, indent=0, plot_button=False,
      equal_widths=None):
    """
    Draw a wx.ListCtrl and populate from the table.  Can optionally include
    a button to launch a graph viewer window; this is used when the table
    contains more than one graph.
    """
    assert (self._current_panel is not None)
    height = min(400, max(80, 20 + (20*table.n_rows)))
    width = TEXT_WIDTH - 40
    wxtable = wxtbx.tables.TableView(
      parent=self._current_panel,
      style=wx.LC_REPORT|wx.SIMPLE_BORDER,
      size=(width, height))
    wxtable.SetTable(table)
    self._tables.append(wxtable)
    self._current_sizer.Add(wxtable, 0, wx.ALL|wx.EXPAND, 5)
    if (plot_button):
      assert hasattr(table, "get_graph")
      btn = wx.Button(parent=self._current_panel, label="Show graphs")
      self._current_sizer.Add(btn, 0, wx.ALL, 5)
      self.Bind(wx.EVT_BUTTON, wxtable.OnViewGraphs, btn)

  def show_plot(self, table):
    """
    Create an inline matplotlib plot, using the wxtbx.plots.loggraph wrapper.
    """
    assert (self._current_panel is not None)
    graph = wxtbx.plots.iotbx_data_plot_base(
      parent=self._current_panel,
      tables=[table],
      size=(TEXT_WIDTH - 40, min(500, TEXT_WIDTH * 0.75)))
    graph.set_plot(table.only_plot())
    self._graphs.append((graph, table.title))
    self._current_sizer.Add(graph, 0, wx.ALL|wx.EXPAND, 5)

  def show_plots_row(self, tables):
    szr = wx.BoxSizer(wx.HORIZONTAL)
    self._current_sizer.Add(szr, 0, wx.EXPAND)
    plot_w = (TEXT_WIDTH / 3) - 10
    plot_dimensions = (plot_w, min(360, plot_w * 0.9))
    for table in tables :
      graph = wxtbx.plots.small_plot(
        parent=self._current_panel,
        table=table,
        size=plot_dimensions)
      graph.set_plot(table.only_plot())
      self._graphs.append((graph, table.title))
      szr.Add(graph, 0, wx.ALL|wx.EXPAND, 5)

  def show_text_columns(self, rows, indent=0):
    prefix = " "*indent
    n_cols = len(rows[0])
    sizer = wx.FlexGridSizer(cols=n_cols) # Remove rows Python 3 fix
    self._current_sizer.Add(sizer)
    for row in rows :
      assert (len(row) == n_cols)
      txt1 = wx.StaticText(parent=self._current_panel, label=prefix+row[0])
      font = txt1.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      txt1.SetFont(font)
      sizer.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      for cell in row[1:] :
        txt_cell = wx.StaticText(parent=self._current_panel, label=cell)
        sizer.Add(txt_cell, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def newline(self):
    pass

  def OnSaveImage(self, event):
    graph_names = [ title for (graph, title) in self._graphs ]
    selected = wxtbx.misc_dialogs.ChoiceSelector(
      title="Save plot image",
      message="Please choose a plot to save:",
      choices=graph_names,
      defaultChoice=0)
    graph_index = graph_names.index(selected)
    (graph, title) = self._graphs[graph_index]
    graph.save_image()

  def add_change_symmetry_button(self):
    data_file = os.path.join(os.getcwd(), "xtriage_data.pkl")
    if (not os.path.exists(data_file)):
      return
    btn = wx.Button(parent=self._current_panel,
      label="Save data in selected setting...")
    btn.data_file = data_file
    btn.symm_table = self._tables[-1]
    self._current_sizer.Add(btn, 0, wx.ALL, 5)
    self.Bind(wx.EVT_BUTTON, OnChangeSymmetry, btn)

  def show_issues(self, issues):
    """
    Display a list of possible issues with traffic light-like symbols
    indicating severity.  Each item is a button that links to the appropriate
    results section (if defined).
    """
    grid = wx.FlexGridSizer(cols=1) # Remove rows Python 3 fix
    self._current_sizer.Add(grid, 0, wx.ALL, 10)
    for severity, message, linkto in issues :
      ctrl = DrawStatusLightControl(parent=self._current_panel,
        #size=(32,32),
        message=message,
        name=linkto,
        level=severity)#style=wx.SIMPLE_BORDER).SetLevel(severity)
      self.Bind(wx.EVT_BUTTON, self.OnShowSection, ctrl)
      grid.Add(ctrl, 0, wx.ALL|wx.EXPAND, 0)

  def OnShowSection(self, evt):
    """
    Change the page selection in response to clicking on one of the issue
    summaries, and scroll the corresponding results section into view.
    """
    section_name = evt.GetEventObject().GetName()
    if (section_name is None) or (section_name == "button") : return
    section = self._sections.get(section_name)
    assert (section is not None), section_name
    page = section.GetParent()
    self.SetPage(self.GetPageIndex(page))
    page.ScrollChildIntoView(section)

class wx_output(wx_output_base, wxtbx.windows.ChoiceBook):
  """
  Xtriage output implemented as a wxPython notebook.
  """
  def SetupScrolling(self):
    """
    Initialize scrollbars on child panels.
    """
    for i_page in range(self.GetPageCount()):
      panel = self.GetPage(i_page)
      panel.Layout()
      panel.SetupScrolling(scrollToTop=True)

  def create_panel(self):
    return wx_panel(parent=self, style=wx.SUNKEN_BORDER)

  def add_panel(self, panel, title):
    if (title == "Twinning and intensity statistics summary"):
      self.InsertPage(0, panel, title)
    elif (title == "Xtriage summary"):
      self.InsertPage(0, panel, title)
    else :
      self.AddPage(panel, title)

class wx_output_panel(wx_output_base, wx_panel):
  """
  Single panel for displaying an individual Xtriage result.  If more than
  one section is displayed (corresponding to multiple notebook pages in the
  standard window), an error will be raised.
  """
  def create_panel(self):
    assert (self._current_panel is None)
    return self

  def add_panel(self, panel, title):
    pass

def DrawStatusLightControl(parent, message, name, level):
  """
  Create a control (based on wxtbx.metallicbutton) with a traffic light-like
  symbol indicating the severity of a message, with the actual message next to
  it.
  """
  bmp = wx.EmptyBitmap(32, 32)
  dc = wx.MemoryDC()
  dc.SelectObject(bmp)
  gc = wx.GraphicsContext.Create(dc)
  dc.Clear()
  if (level == 0):
    gc.SetBrush(wx.Brush((0,255,0)))
  elif (level == 1):
    gc.SetBrush(wx.Brush((255,255,0)))
  elif (1 < level < 2):
    gc.SetBrush(wx.Brush((255,140,0)))
  else :
    gc.SetBrush(wx.Brush((255,0,0)))
  gc.DrawEllipse(4, 4, 28, 28)

  # delete contexts after drawing is complete and before setting mask
  # (fixes Windows error, not an issue with OS X or Linux)
  del gc
  del dc

  bmp.SetMask(wx.Mask(bmp,wx.WHITE))
  btn = metallicbutton.MetallicButton(
    parent=parent,
    label2=message,
    size=(800,-1),
    caption_size=12,
    bmp=bmp)
  btn.SetFont(parent.GetFont())
  if (name is not None):
    btn.SetName(name)
  return btn

class XtriageFrame(wx.Frame):
  """
  Frame for displaying a wx_output notebook independently (used in AutoSol
  and similar apps).
  """
  def __init__(self, *args, **kwds):
    super(XtriageFrame, self).__init__(*args, **kwds)
    self.output_panel = self.create_panel()
    self.SetupToolbar()
    self._result = None

  def create_panel(self):
    return wx_output(parent=self)

  def SetupToolbar(self):
    self.toolbar = self.CreateToolBar(style=wx.TB_TEXT)
    self.AddAppSpecificButtons()
    bmp = wxtbx.bitmaps.fetch_icon_bitmap("mimetypes", "spreadsheet")
    btn = self.toolbar.AddLabelTool(-1, "Save graph", bmp,
      shortHelp="Save graph", kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.output_panel.OnSaveImage, btn)
    bmp = wxtbx.bitmaps.fetch_icon_bitmap("mimetypes", "txt")
    btn = self.toolbar.AddLabelTool(-1, "View log file", bmp,
      shortHelp="View log file", kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnDisplayLog, btn)
    self.toolbar.Realize()

  def AddAppSpecificButtons(self):
    pass

  def SetResult(self, result):
    self._result = result
    result.show(out=self.output_panel)
    self.output_panel.SetupScrolling()
    if hasattr(self.output_panel, "SetPage"):
      self.output_panel.SetPage(0)

  def OnDisplayLog(self, event):
    if (self._result.log_file_name is not None):
      if (sys.platform == "darwin"):
        easy_run.call("open -a TextEdit %s" % self._result.log_file_name)
      else :
        pass

class XtriageFrameSingleResult(XtriageFrame):
  def create_panel(self):
    return wx_output_panel(parent=self)

def OnChangeSymmetry(event):
  from mmtbx.scaling.xtriage import change_symmetry
  button = event.GetEventObject()
  data_file = button.data_file
  assert (data_file is not None) and os.path.isfile(data_file)
  item = button.symm_table.GetFirstSelected()
  if (item == -1):
    raise Sorry("Please select a symmetry setting first!")
  default_file = os.path.join(os.path.dirname(data_file), "reindexed.mtz")
  file_name = wxtbx.path_dialogs.manager().select_file(
    parent=button,
    message="Output data file in new symmetry",
    style=wx.FD_SAVE,
    wildcard="MTZ files (*.mtz)|*.mtz",
    current_file=default_file)
  data = easy_pickle.load(data_file)
  space_group_symbol = str(button.symm_table.GetItemText(item))
  change_symmetry(
    miller_array=data,
    space_group_symbol=space_group_symbol,
    file_name=str(file_name),
    log=sys.stdout)
  assert os.path.isfile(file_name)
  wx.MessageBox(("The data (as amplitudes) have been saved to %s "+
      "with the selected symmetry setting.") % file_name)
