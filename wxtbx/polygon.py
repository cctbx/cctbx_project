
import wxtbx.polygon_db_viewer
import wxtbx.bitmaps
import wx.lib.colourselect
import wx
import mmtbx.polygon.output
from mmtbx import polygon
from libtbx import adopt_init_args
from math import radians
import sys

class wx_renderer (mmtbx.polygon.output.renderer) :
  def draw_bin (self, out, start, end, angle, color) :
    gc = out
    path = gc.CreatePath()
    path.MoveToPoint(start[0], start[1])
    path.AddLineToPoint(end[0], end[1])
    path.CloseSubpath()
    gc.PushState()
    gc.SetPen(wx.Pen(color, 10))
    gc.StrokePath(path)
    gc.PopState()

  def draw_box (self, out, points, color) :
    gc = out
    path = gc.CreatePath()
    path.MoveToPoint(points[0][0], points[0][1])
    path.AddLineToPoint(points[1][0], points[1][1])
    path.AddLineToPoint(points[2][0], points[2][1])
    path.AddLineToPoint(points[3][0], points[3][1])
    path.AddLineToPoint(points[0][0], points[0][1])
    path.CloseSubpath()
    gc.PushState()
    gc.SetPen(wx.Pen(color, 1)) #TRANSPARENT_PEN)
    gc.SetBrush(wx.Brush(color))
    gc.FillPath(path)
    gc.PopState()

  def draw_solid_line (self, out, start, end, color) :
    gc = out
    line = gc.CreatePath()
    line.MoveToPoint(start[0], start[1])
    line.AddLineToPoint(end[0], end[1])
    line.CloseSubpath()
    gc.PushState()
    if self.color_model == "gray" :
      gc.SetPen(wx.Pen("red", 2))
    else :
      gc.SetPen(wx.Pen("black", 2))
    gc.StrokePath(line)
    gc.PopState()

  def draw_dashed_line (self, out, start, end, color) :
    pass

  def draw_labels (self, out, label, min, max, value, pos, angle) :
    gc = out
    label_font = wx.Font(14, wx.NORMAL, wx.NORMAL, wx.BOLD)
    stat_font = wx.Font(12, wx.MODERN, wx.NORMAL, wx.NORMAL)
    gc.PushState()
    gc.SetPen(wx.Pen("black", 1))
    gc.SetFont(label_font)
    (text_w, text_h) = gc.GetTextExtent(label)
    (anchor_x, anchor_y) = pos
    (text_x, text_y) = self.get_text_position(text_w, text_h, anchor_x,
      anchor_y, angle)
    gc.DrawText(label, text_x, text_y)
    gc.PopState()
    gc.PushState()
    gc.SetFont(gc.CreateFont(stat_font, wx.RED))
    (text_w, text_h) = gc.GetTextExtent(min)
    text_y += text_h + 2
    gc.DrawText(min, text_x, text_y)
    gc.PopState()
    gc.PushState()
    gc.SetFont(gc.CreateFont(stat_font, wx.BLACK))
    (text_w, text_h) = gc.GetTextExtent(value)
    text_y += text_h + 2
    gc.DrawText(value, text_x, text_y)
    gc.PopState()
    gc.PushState()
    gc.SetFont(gc.CreateFont(stat_font, wx.RED))
    (text_w, text_h) = gc.GetTextExtent(max)
    text_y += text_h + 2
    gc.DrawText(max, text_x, text_y)
    gc.PopState()

  def get_text_position (self, w, h, x, y, angle) :
    if angle >= radians(60) and angle < radians(120) :
      text_x = x - (w/2) - 5
      text_y = y - h - 15
    elif angle >= radians(120) and angle < radians(240) :
      text_x = x - w - 15
      text_y = y - (h/2)
    elif angle >= radians(240) and angle < radians(300) :
      text_x = x - (w/2)
      text_y = y
    else : # 300 =< angle < 420
      text_x = x + 5
      text_y = y - (h/2)
    return (text_x, text_y)

class PolygonPanel (wx.Panel) :
  def __init__ (self, parent, renderer) :
    wx.Panel.__init__(self, parent, -1)
    self.renderer = renderer
    self.renderer.resize((640, 640))
    self.SetMinSize((480,480))
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_SIZE, self.OnSize)

  def OnSize (self, event) :
    self.renderer.resize(self.GetSize())
    self.Refresh()

  def OnPaint (self, event) :
    self.renderer.resize(self.GetSize())
    dc = wx.PaintDC(self)
    gc = wx.GraphicsContext.Create(dc)
    self.renderer.draw(gc)

  def draw_color_key (self, dc) :
    gc = wx.GraphicsContext.Create(dc)
    stat_font = wx.Font(12, wx.MODERN, wx.NORMAL, wx.NORMAL)
    x = 40
    y = self.h - 10
    i = 0
    for shade in bin_colors :
      gc.PushState()
      gc.SetPen(wx.Pen(shade, 10))
      path = gc.CreatePath()
      path.MoveToPoint(x, y)
      path.AddLineToPoint(x + 40, y)
      path.CloseSubpath()
      gc.StrokePath(path)
      gc.PopState()
      if i < len(self.cutoffs) :
        gc.PushState()
        gc.SetFont(gc.CreateFont(stat_font, wx.BLACK))
        gc.DrawText(str(self.cutoffs[i]), x + 50, y - 6)
        gc.PopState()
      x += 80
      i += 1

  def OnChar (self, event) :
    keycode = event.GetKeyCode()
    if keycode == 32 :
      self.OnSave()

  def OnSave (self, event=None) :
    rect = self.GetRect()
    bitmap = wx.EmptyBitmap(rect.width, rect.height)
    memory_dc = wx.MemoryDC()
    memory_dc.SelectObject(bitmap)
    memory_dc.SetBackgroundMode(wx.TRANSPARENT)
    gc = wx.GraphicsContext.Create(memory_dc)
    self.renderer.draw(gc)
    output_file = wx.FileSelector("Save image as:",
      default_filename="polygon.png",
      wildcard="PNG image (*.png)|*.png", flags=wx.SAVE)
    if output_file != "" :
      bitmap.SaveFile(output_file, wx.BITMAP_TYPE_PNG)
    if event is not None :
      event.Skip()

  def reset_layout (self) :
    pass

class PolygonFrame (wx.Frame) :
  def __init__ (self, parent, histogram_data, structure_stats) :
    wx.Frame.__init__(self, parent, -1, "POLYGON", size=(1024,720))
    self.SetMinSize((800,500))
    adopt_init_args(self, locals())
    self.renderer = wx_renderer(histogram_data, structure_stats,
      center=(0.5, 0.475))
    self.setup_toolbar()
    main_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.SetSizer(main_sizer)
    self.main_sizer = main_sizer
    self.info_sizer = wx.BoxSizer(wx.VERTICAL)
    self.main_sizer.Add(self.info_sizer, 0, wx.ALL|wx.EXPAND)
    self.draw_top_panel()
    self.label_panel = None
    self.draw_color_key()
    self.polygon_panel = PolygonPanel(self, self.renderer)
    main_sizer.Add(self.polygon_panel, 1, wx.ALL|wx.EXPAND)
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def setup_toolbar (self) :
    self.toolbar = None
    save_icon = wxtbx.bitmaps.fetch_icon_bitmap("actions", "save_all")
    plot_icon = wxtbx.bitmaps.fetch_icon_bitmap("mimetypes", "spreadsheet")
    if (save_icon is not None) and (plot_icon is not None) :
      self.toolbar = wx.ToolBar(self, style=wx.TB_3DBUTTONS|wx.TB_TEXT)
      if sys.platform == "darwin" :
        save_btn = self.toolbar.AddLabelTool(-1, "Save", save_icon,
          kind=wx.ITEM_NORMAL)
      else :
        save_btn = self.toolbar.AddSimpleTool(-1, save_icon, "Save")
      self.Bind(wx.EVT_MENU, self.OnSave, save_btn)
      hist_btn = self.toolbar.AddLabelTool(-1, "Show histograms", plot_icon,
        kind=wx.ITEM_NORMAL)
      self.Bind(wx.EVT_MENU, self.OnDisplayHistogram, hist_btn)
      self.SetToolBar(self.toolbar)
      self.toolbar.Realize()

  def draw_top_panel (self) :
    top_panel = wx.Panel(self, -1, style=wx.SIMPLE_BORDER)
    top_sizer = wx.BoxSizer(wx.VERTICAL)
    top_panel.SetSizer(top_sizer)
    caption = wx.StaticText(top_panel, -1,
"""This graph shows histograms of the distribution of selected statistics \
across %d PDB entries of similar resolution, with the range specified by \
numbers printed in red.  Statistics for the current structure are printed in \
black; the connecting polygon (in black) shows where these values fall in the \
distribution. A typical well-refined structure will have a small and roughly \
equilateral polygon.""" %
      self.renderer.n_pdb)
    caption.Wrap(320)
    top_sizer.Add(caption, 0, wx.ALL, 5)
    mode_sizer = wx.BoxSizer(wx.HORIZONTAL)
    top_sizer.Add(mode_sizer)
    caption1 = wx.StaticText(top_panel, -1, "Color scheme:")
    mode_sizer.Add(caption1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    mode_choice = wx.Choice(top_panel, -1,
      choices=["Original (relative scaling)",
               "Rainbow (by bin size)",
               "Rainbow (relative scaling)",
               "Red to blue (by bin size)",
               "Red to blue (relative scaling)",
               "Red (by bin size)",
               "Blue (by bin size)",
               "Grayscale (by bin size)"])
    mode_choice.SetSelection(1)
    mode_sizer.Add(mode_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnRecolor, mode_choice)
    caption2 = wx.StaticText(top_panel, -1,
      "Citation: Urzhumtseva et al. Acta Cryst. 2009, D65:297-300.")
    caption2.Wrap(320)
    top_sizer.Add(caption2, 0, wx.ALL, 5)
    top_sizer.Fit(top_panel)
    self.info_sizer.Add(top_panel, 1, wx.ALL|wx.EXPAND)

  def draw_color_key (self) :
    if self.label_panel is not None :
      self.main_sizer.Detach(self.label_panel)
      self.label_panel.Destroy()
    lower_panel = wx.Panel(self, -1, style=wx.SIMPLE_BORDER)
    lower_sizer = wx.BoxSizer(wx.VERTICAL)
    lower_panel.SetSizer(lower_sizer)
    if self.renderer.relative_scale_colors :
      caption = wx.StaticText(lower_panel, -1,
"""Histogram bins are colored based on the ratio of the number of structures \
in each bin to the average number per bin:""")
    else :
      caption = wx.StaticText(lower_panel, -1,
"""Histogram bins are colored by the number of structures in each bin.""")
    caption.Wrap(320)
    lower_sizer.Add(caption, 0, wx.ALL, 5)
    key_sizer = wx.BoxSizer(wx.HORIZONTAL)
    key_sizer = wx.FlexGridSizer(rows=0, cols=6)
    lower_sizer.Add(key_sizer)
    colors, cutoffs = self.renderer.get_color_key()
    for i, color in enumerate(colors) :
      color_widget = ColorBox(lower_panel, -1, "", color, size=(24,24))
      key_sizer.Add(color_widget, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      if i < len(cutoffs) :
        if self.renderer.relative_scale_colors :
          label = wx.StaticText(lower_panel, -1, "=< %s" % str(cutoffs[i]))
        else :
          label = wx.StaticText(lower_panel, -1, "= %s" % str(cutoffs[i]))
        key_sizer.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    lower_sizer.Fit(lower_panel)
    self.info_sizer.Add(lower_panel, 1, wx.ALL|wx.EXPAND)
    self.info_sizer.Layout()
    self.Layout()
    self.label_panel = lower_panel

  def OnRecolor (self, event) :
    mode = event.GetEventObject().GetSelection()
    mode_info = [ ("original", True),
                  ("rainbow", False),
                  ("rainbow", True),
                  ("rmb", False),
                  ("rmb", True),
                  ("red", False),
                  ("blue", False),
                  ("gray", False) ]
    (model_name, relative_scaling) = mode_info[mode]
    self.renderer.set_color_model(model_name, relative_scaling)
    self.draw_color_key()
    self.Refresh()

  def OnClose (self, event) :
    wx.CallAfter(self.Destroy)

  def OnDestroy (self, event) :
    if self.parent is not None and hasattr(self.parent, "polygon_frame") :
      self.parent.polygon_frame = None

  def OnSave (self, event) :
    self.polygon_panel.OnSave()

  def OnResize (self, event) :
    #self.panel.OnResize(event)
    self.polygon_panel.Layout()

  def OnDisplayHistogram (self, event) :
    keys = [ key for key, data in self.histogram_data ]
    dlg = wx.SingleChoiceDialog(
      parent=self,
      message="Which statistic do you want to view as a histogram?",
      caption="Select a histogram to display",
      choices=keys)
    if (dlg.ShowModal() == wx.ID_OK) :
      choice = dlg.GetSelection()
    wx.CallAfter(dlg.Destroy)
    if (choice is not None) :
      frame = wxtbx.polygon_db_viewer.HistogramFrame(
        parent=self)
      frame.show_histogram(
        data=self.histogram_data[choice][1],
        n_bins=10,
        reference_value=self.structure_stats[keys[choice]],
        xlabel=mmtbx.polygon.output.stat_names[keys[choice]])
      frame.Show()

class ColorBox (wx.lib.colourselect.ColourSelect) :
  def OnClick (self, event) :
    pass

if (__name__ == "__main__") :
  app = wx.App(0)
  stats = {
    "r_work" : 0.25,
    "r_free" : 0.28,
    "adp_mean_all" : 20.0,
    "bond_rmsd" : 0.02,
    "angle_rmsd" : 1.8,
    "clashscore" : 20.0
  }
  data = mmtbx.polygon.output.get_basic_histogram_data(d_min=2.5)
  frame = PolygonFrame(
    parent=None,
    histogram_data=data,
    structure_stats=stats)
  frame.Show()
  app.MainLoop()
