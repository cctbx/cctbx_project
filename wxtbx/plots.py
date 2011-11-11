
from __future__ import division
from wxtbx import bitmaps
import wx
from libtbx import object_oriented_patterns as oop
from libtbx import adopt_init_args
import math
import sys
if (sys.version_info[2] >= 6) :
  import warnings
  warnings.simplefilter('ignore', DeprecationWarning)

class plot_container (wx.BoxSizer) :
  def __init__ (self,
                parent,
                figure_size=(8,6),
                font_size=12,
                title_font_size=10,
                facecolor='white',
                transparent=False,
                handle_left_click=False,
                show_data_points=True,
                point_types=('o', '^', '+', 's', 'D'),
                title_alignment="right") :
    wx.BoxSizer.__init__(self, wx.VERTICAL)
    adopt_init_args(self, locals())
    self._fonts = {}
    self.disabled = False
    try :
      import matplotlib
      from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
      from matplotlib.backends.backend_wxagg import FigureManager
      import matplotlib.ticker
      import matplotlib.cm
      import matplotlib.figure
      import matplotlib.font_manager
    except ImportError, e :
      print ""
      print "Error loading matplotlib module:"
      print e
      print ""
      self.disabled = True
      self.canvas = oop.null()
      self.figure = oop.null()
      self.text_font = oop.null()
      self.p = oop.null()
      w = int(figure_size[0] * 72)
      h = int(figure_size[1] * 72)
      panel = wx.Panel(parent=parent,
        size=(w,h))
      panel.SetBackgroundColour((150,150,150))
      szr = wx.BoxSizer(wx.VERTICAL)
      panel.SetSizer(szr)
      txt = wx.StaticText(parent=panel,
        label="Plotting disabled due to missing libraries.")
      szr.Add(txt, 1, wx.ALL|wx.ALIGN_CENTER, 10)
      txt.SetForegroundColour((255,0,0))
      font = txt.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      txt.SetFont(font)
      self.Add(panel, 1, wx.EXPAND|wx.ALL)
      self.null_fmt = oop.null()
    else :
      self.figure = matplotlib.figure.Figure(figure_size, 72, linewidth=0,
        facecolor=facecolor)
      if transparent :
        self.figure.figurePatch.set_alpha(0.0)
      self.canvas = FigureCanvasWxAgg(self.parent, -1, self.figure)
      self.canvas.toolbar = oop.null()
      self.figmgr = FigureManager(self.canvas, 1, self)
      self.Add(self.canvas, 1, wx.EXPAND|wx.ALL)
      self.setup_fonts()
      self.null_fmt = matplotlib.ticker.NullFormatter()
      if self.handle_left_click :
        self.canvas.mpl_connect("button_release_event", self.OnClick)
      else :
        self.canvas.Bind(wx.EVT_CONTEXT_MENU, self.OnRightClick, self.canvas)

  def setup_fonts (self) :
    import matplotlib.font_manager
    self._fonts["basic"] = matplotlib.font_manager.FontProperties(
      family = ["Courier", "Monaco", "monospace"],
      weight = "normal",
      size   = 10)
    self._fonts["value_label"] = matplotlib.font_manager.FontProperties(
      family = ["Courier", "Monaco", "monospace"],
      weight = "normal",
      size   = self.font_size)
    self._fonts["axis_label"] = matplotlib.font_manager.FontProperties(
      family = ["Helvetica", "sans-serif"],
      weight = "bold",
      size   = self.font_size)
    self._fonts["text"] = matplotlib.font_manager.FontProperties(
      family = ["Helvetica", "sans-serif"],
      weight = "normal",
      size   = self.font_size)
    self._fonts["title"] = matplotlib.font_manager.FontProperties(
      family = ["Helvetica", "sans-serif"],
      weight = "normal",
      size   = self.title_font_size)

  def get_font (self, font_type) :
    font = self._fonts.get(font_type, None)
    if (font is None) :
      font = self._fonts.get("basic", None)
    return font

  def GetToolBar (self) :
    return None

  def OnRightClick (self, event) :
    pass

  def OnClick (self, event) :
    if not hasattr(event, "button") :
      raise RuntimeError("The OnClick method of plot_container handles "+
        "matplotlib events only.  Please fix your code!")
    if event.button == 1 :
      self.process_mouse_click(event)

  def process_mouse_click (self, mpl_event) :
    pass

  def save_image (self, default_path="") :
    output_file = wx.FileSelector("Saved image name",
      default_path=default_path,
      default_filename="plot.pdf",
      wildcard="Adobe PDF figure (*.pdf)|*.pdf|" + \
               "PNG image (*.png)|*.png|" + \
               "Postscript figure (*.ps)|*.ps", flags=wx.SAVE)
    if output_file != "" :
      if output_file[-3:] == "pdf" :
        self.figure.savefig(output_file, orientation="landscape", format="pdf")
      elif output_file[-2:] == "ps" :
        self.figure.savefig(output_file, orientation="landscape", format="ps")
      else :
        self.figure.savefig(output_file, format="png")

class histogram (plot_container) :
  def show_histogram (self, data, n_bins, reference_value=None, pos=(1,1,1),
      draw_now=True, x_label=None, y_label=None) :
    assert len(pos) == 3
    self.figure.clear()
    p = self.figure.add_subplot(*pos)
    p.hist(data, n_bins, facecolor='blue')
    if reference_value is not None :
      p.axvline(reference_value, color='r')
    if x_label is not None :
      p.set_xlabel(x_label)
    if y_label is not None :
      p.set_ylabel(y_label)
    if draw_now :
      self.canvas.draw()
    return p

def convert_xyz_value_list (values, null_value=0.0) :
  import numpy
  values = sorted(list(values), cmp=lambda x,y: cmp(x[0], y[0]))
  x_rows = [[values[0]]]
  for i, xyz in enumerate(values[1:]) :
    if (xyz[0] != x_rows[-1][-1][0]) :
      x_rows.append([])
    x_rows[-1].append(xyz)
  x_values = [ x_row[0][0] for x_row in x_rows ]
  y_values = sorted([ x_rows[0][n][1] for n in range(len(x_rows[0])) ])
  assert (len(values) == (len(x_values) * len(y_values)))
  z_values = []
  for j in range(len(y_values)) :
    z_values.append([])
  for i, x_row in enumerate(x_rows) :
    x_row = sorted(x_row, cmp=lambda x,y: cmp(x[1], y[1]))
    for j, (x,y,z) in enumerate(x_row) :
      assert (y in y_values)
      if (z is not None) :
        assert (isinstance(z, int) or isinstance(z, float))
        z_values[j].append(z)
      else :
        z_values[j].append(numpy.NaN)
      #  z_values[j].append(null_value)
  return (numpy.array(x_values), numpy.array(y_values), numpy.array(z_values))

class image_plot (plot_container) :
  def show_plot (self,
                 x_data=(),
                 y_data=(),
                 z_data=(),
                 values=(),
                 x_label=None,
                 y_label=None,
                 title=None,
                 cmap=None,
                 interpolation="nearest") :
    if (len(values) > 0) :
      assert (len(x_data) == len(y_data) == len(z_data) == 0)
      (x_data, y_data, z_data) = convert_xyz_value_list(values)
    from matplotlib.image import NonUniformImage
    self.figure.clear()
    ax = self.figure.add_subplot(111)
    im = NonUniformImage(ax, interpolation=interpolation)
    if (cmap is not None) :
      im.set_cmap(get_colormap(cmap))
    im.set_data(x_data, y_data, z_data)
    ax.images.append(im)
    ax.set_xlim(x_data[0], x_data[-1])
    ax.set_ylim(y_data[0], y_data[-1])
    if (x_label is not None) :
      ax.set_xlabel(x_label)
    if (y_label is not None) :
      ax.set_ylabel(y_label)
    if (title is not None) :
      ax.set_title(title)
    self.canvas.draw()

class iotbx_data_plot_base (plot_container) :
  def __init__ (self,
                parent,
                tables,
                size=(640,480),
                **kwds) :
    adopt_init_args(self, locals())
    (x, y, w, h) = tuple(wx.GetClientDisplayRect())
    (width, height) = size
    fig_w = float(width) / 72.0
    fig_h = float(height) / 72.0
    if (fig_w * 72) > (w - 20) :
      fig_w = int(math.floor((w-40) / 72))
    if (fig_h * 72) > (h - 120) :
      fig_h = int(math.floor((h-160) / 72))
    plot_container.__init__(self, parent, (fig_w, fig_h), **kwds)
    self.p = self.figure.add_subplot(111)
    self.plot_type = None

  def set_tables (self, tables) :
    self.tables = tables

  def set_plot (self, graph_name=None, table_name=None, table_index=0) :
    table = None
    if table_name is None :
      table = self.tables[table_index]
    else :
      for t in self.tables :
        if t.title == table_name :
          table = t
    if table is not None :
      self.plot_type = getattr(table, "plot_type", "GRAPH")
      graph = table.get_graph(graph_name)
      self.show_plot(graph)

  def show_plot (self, graph, line_width=1, show_points=True, show_grid=True) :
    if self.disabled :
      return
    self.figure.clear()
    self.graph = graph
    self.p = self.figure.add_subplot(111)
    if graph is None :
      return
    if self.plot_type == "SCATTER" :
      show_points = True
      show_lines = False
    else :
      show_lines = True
    if self.show_data_points :
      point_types = self.point_types #['o', '^', '+', 's', 'D']
    else :
      point_types = [""]
    point_index = 0
    for (x_values, y_values) in self.graph.get_plots() :
      plot_type = ""
      if show_points :
        plot_type = "%s" % point_types[point_index]
        point_index += 1
        if point_index >= len(point_types) :
          point_index = 0
      if show_lines :
        plot_type += "-"
      self.p.plot(x_values, y_values, plot_type, linewidth=line_width)
    self.format_labels()
    if show_grid :
      self.p.get_axes().grid(True, color="0.75")
    self.p.get_axes().set_autoscale_on(True)
    self.p.set_title(graph.name, fontproperties=self.get_font("title"),
      horizontalalignment=self.title_alignment)
    self.canvas.draw()
    self.parent.Refresh()

  def format_x_axis (self) :
    if self.tables[0].x_is_inverse_d_min :
      xdata = self.tables[0].get_x_as_resolution()
      self.p.get_axes().set_xlabel("Resolution",
        fontproperties=self.get_font("axis_label"))
      xticks = self.p.get_axes().get_xticks()
      xticklabels = []
      for x in xticks :
        if (x != 0) :
          x = math.sqrt(1 / x)
        xticklabels.append("%.2f" % x)
      #self.p.get_axes().set_xticks(xticks)
      self.p.get_axes().set_xticklabels(xticklabels)
    else :
      if self.graph.x_axis_label is not None :
        self.p.get_axes().set_xlabel(self.graph.x_axis_label,
          fontproperties=self.get_font("axis_label"))
      else :
        self.p.get_axes().set_xlabel(self.graph.x_label,
          fontproperties=self.get_font("axis_label"))
    for ticklabel in self.p.get_axes().get_xticklabels() :
      ticklabel.set_fontproperties(self.get_font("value_label"))

  def axvline (self, x, **kwargs) :
    axes = self.p#.get_axes()
    if self.tables[0].x_is_inverse_d_min :
      axes.axvline(x=(1.0 / (x**2)), **kwargs)
    else :
      axes.axvline(x=x, **kwargs)
    self.canvas.draw()
    self.parent.Refresh()

  def format_y_axis (self) :
    for ticklabel in self.p.get_axes().get_yticklabels() :
      ticklabel.set_fontproperties(self.get_font("value_label"))
    if self.graph.y_axis_label is not None :
      self.p.get_axes().set_ylabel(self.graph.y_axis_label,
          fontproperties=self.get_font("axis_label"))

  def format_labels (self) :
    self.figure.legend(self.p.lines, self.graph.y_labels,
      prop=self.get_font("text"))
    self.format_x_axis()
    self.format_y_axis()

  def show_grid (self, show=True) :
    if show :
      self.p.get_axes().grid(True, color="0.75")
    else :
      self.p.get_axes().grid(False)
    self.canvas.draw()
    self.parent.Refresh()

  def OnRightClick (self, event) :
    menu = wx.Menu()
    menu_item = menu.Append(-1, "Save image")
    self.parent.Bind(wx.EVT_MENU, self.OnSave, menu_item)
    self.parent.PopupMenu(menu)

  def OnSave (self, event=None) :
    self.save_image()

class small_plot (iotbx_data_plot_base) :
  def __init__ (self, parent, table, size=(320,320)) :
    iotbx_data_plot_base.__init__(self,
      parent=parent,
      tables=[table],
      size=size,
      font_size=9,
      title_font_size=10,
      title_alignment="center",
      point_types=('+'))

class plot_frame (wx.Frame) :
  controls_on_top = True
  show_controls_default = True
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.setup_toolbar()
    self.toolbar.Realize()
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.draw_top_panel()
    self.plot_panel = self.create_plot_panel()
    if self.controls_on_top :
      self.sizer.Add(self.top_panel, 0, wx.EXPAND|wx.ALL)
      self.sizer.Add(self.plot_panel, 1, wx.EXPAND|wx.ALL)
    else :
      self.sizer.Add(self.plot_panel, 1, wx.EXPAND|wx.ALL)
      self.sizer.Add(self.top_panel, 0, wx.EXPAND|wx.ALL)
    self._show_controls = self.show_controls_default
    if not self.show_controls_default :
      self.top_panel.Hide()
    self.SetSizer(self.sizer)
    self.sizer.Layout()
    self.Fit()
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)

  def setup_toolbar (self) :
    tb_buttons = [
      ("Show/hide controls",
       bitmaps.fetch_icon_bitmap("apps", "advancedsettings"),
       self.OnToggleControls),
      ("Save",
       bitmaps.fetch_icon_bitmap("actions", "save_all"),
       self.OnSave),
      #("Print",
      # bitmaps.fetch_icon_bitmap("devices", "printer1"),
      # self.OnPrint),
    ]
    tb = wx.ToolBar(self, style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    tb.SetToolBitmapSize((32,32))
    self.SetToolBar(tb)
    for (name, bitmap, function) in tb_buttons :
      tool_button = tb.AddLabelTool(-1, name, bitmap, kind=wx.ITEM_NORMAL)
      self.Bind(wx.EVT_MENU, function, tool_button)
    self.toolbar = tb

  def draw_top_panel (self) :
    self.top_panel = wx.Panel(parent=self, style=wx.SUNKEN_BORDER)

  def create_plot_panel (self) :
    return (0,0)

  def OnExit (self, event) :
    self.Close()

  def OnClose (self, event) :
    wx.CallAfter(self.Destroy)

  def OnDestroy (self, event) :
    pass

  def OnSave (self, event) :
    self.plot_panel.save_image()

  def OnToggleControls (self, event) :
    if self._show_controls :
      self._show_controls = False
      self.top_panel.Hide()
    else :
      self._show_controls = True
      self.top_panel.Show()
    self.sizer.Layout()
    self.Fit()

class loggraph (plot_frame) :
  plot_type = "loggraph"
  controls_on_top = True
  show_controls_default = True
  table_selection_label = "Table:"
  plot_selection_label = "Plot:"
  def __init__ (self, parent, title, tables=None, file_name=None,
      processed_lines=None) :
    adopt_init_args(self, locals())
    self.tables = []
    self.graph = None
    self.table_frame = None
    self.table_chooser = None
    self.plot_chooser = None
    self.current_table = None
    self.current_plot = None
    if tables is not None :
      self.tables = tables
    elif file_name is not None or processed_lines is not None :
      self.load_log(file_name, processed_lines, update=False)
      if len(self.tables) > 0 :
        self.current_table = self.tables[0].title
        self.current_plot = self.tables[0].graph_names[0]
    plot_frame.__init__(self, parent=parent, title=title)
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Centre()
    self.switch_plot()

  def draw_top_panel (self) :
    self.top_panel = wx.Panel(parent=self, style=wx.SUNKEN_BORDER)
    cp = self.top_panel
    cp_sizer = wx.BoxSizer(wx.VERTICAL)
    cp.SetSizer(cp_sizer)
    grid = wx.FlexGridSizer(cols=2)
    txt = wx.StaticText(self.top_panel, -1, "Table:")
    # TODO: set bold-face
    grid.Add(txt, 0, wx.ALL|wx.EXPAND, 5)
    self.table_chooser = wx.Choice(parent=cp,
      choices=[ t.title for t in self.tables])
    self.Bind(wx.EVT_CHOICE, self.OnSelectTable, self.table_chooser)
    grid.Add(self.table_chooser, 0, wx.ALL|wx.EXPAND, 5)
    txt2 = wx.StaticText(self.top_panel, -1, "Plot:")
    grid.Add(txt2, 0, wx.ALL|wx.EXPAND, 5)
    plot_choices = self.tables[0].graph_names
    self.plot_chooser = wx.Choice(parent=cp,
      choices=plot_choices)
    self.Bind(wx.EVT_CHOICE, self.OnSelectPlot, self.plot_chooser)
    grid.Add(self.plot_chooser, 0, wx.ALL|wx.EXPAND, 5)
    cp_sizer.Add(grid, 0, wx.EXPAND)
    szr = wx.BoxSizer(wx.HORIZONTAL)
    cp_sizer.Add(szr, 0, wx.EXPAND)
    grid_box = wx.CheckBox(parent=cp,
      label="Show grid")
    grid_box.SetValue(True)
    self.Bind(wx.EVT_CHECKBOX, self.OnToggleGrid, grid_box)
    szr.Add(grid_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    point_box = wx.CheckBox(parent=cp,
      label="Show data points")
    point_box.SetValue(True)
    self.Bind(wx.EVT_CHECKBOX, self.OnTogglePoints, point_box)
    szr.Add(point_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    return cp

  def create_plot_panel (self) :
    self.plot = iotbx_data_plot_base(
      parent=self,
      tables=self.tables,
      transparent=False)
    return self.plot

  def load_log (self, file_name=None, processed_lines=None, update=True) :
    from iotbx import data_plots
    if processed_lines :
      self.tables = data_plots.import_ccp4i_logfile(log_lines=processed_lines)
    else :
      self.tables = data_plots.import_ccp4i_logfile(file_name=file_name)
    if update :
      self.update_interface()

  def update_interface (self, update_tables=True) :
    if len(self.tables) > 0 :
      if update_tables :
        table_choices = [t.title for t in self.tables]
        self.table_chooser.SetItems(table_choices)
        self.table_chooser.SetSelection(0)
      self.plot_chooser.SetItems(self.get_table().graph_names)
      self.plot_chooser.SetSelection(0)
      self.top_panel.Layout()
      self.switch_plot()

  def set_current_plot (self, plot_name) :
    for t in self.tables :
      for t_plot_name in t.graph_names :
        if plot_name == t_plot_name :
          self.table_chooser.SetStringSelection(t.title)
          self.plot_chooser.SetStringSelection(plot_name)
          break
    self.switch_plot()

  def switch_plot (self) :
    table = self.table_chooser.GetStringSelection()
    plot = self.plot_chooser.GetStringSelection()
    if table != "" and plot != "" :
      self.plot.set_plot(graph_name=plot, table_name=table)
      #self.Refresh()

  #--- EVENTS
  def OnSelectTable (self, event) :
    table_name = self.table_chooser.GetStringSelection()
    current_plot = self.plot_chooser.GetStringSelection()
    plot_choices = []
    for t in self.tables :
      if t.title == table_name :
        plot_choices = t.graph_names
    self.plot_chooser.SetItems(plot_choices)
    if current_plot in plot_choices :
      self.plot_chooser.SetStringSelection(current_plot)
    self.switch_plot()

  def OnSelectPlot (self, event) :
    self.switch_plot()

  def OnPrint (self, event) :
    pass

  def OnToggleGrid (self, event) :
    show = event.GetEventObject().GetValue()
    self.plot.show_grid(show)

  def OnTogglePoints (self, event) :
    self.plot.show_data_points = event.GetEventObject().GetValue()
    self.switch_plot()

standard_colormaps = [ ("jet", "Rainbow"),
                       ("Greys", "Greyscale"),
                       ("Reds", "Red"),
                       ("YlGn", "Yellow/Green"),
                       ("Greens", "Green"),
                       ("Blues", "Blue"),
                     ]
colormap_names = [ cm_name for cm_id, cm_name in standard_colormaps ]
colormap_id_dict = dict([ (name, id) for id, name in standard_colormaps ])

def get_colormap (cm_name) :
  import matplotlib.cm
  cm = None
  for cm_id, cm_name2 in standard_colormaps :
    if cm_name2 == cm_name :
      cm = getattr(matplotlib.cm, cm_id, None)
      break
  else :
    cm = getattr(matplotlib.cm, cm_name, None)
  if cm is None :
    cm = matplotlib.cm.jet
  return cm

def exercise () :
  values = [
    (1.5, 3.0, 40.2),
    (2.5, 3.5, 35.1),
    (1.5, 3.5, 45.6),
    (2.0, 3.5, 54.2),
    (1.5, 4.5, 48.2),
    (2.5, 4.5, 28.4),
    (2.5, 4.0, 37.6),
    (2.0, 2.5, None),
    (2.0, 3.0, 29.8),
    (2.0, 4.0, 36.7),
    (2.0, 4.5, 33.2),
    (2.5, 2.5, None),
    (1.5, 2.5, 39.5),
    (2.5, 3.0, None),
    (1.5, 4.0, 50.9),
  ]
  sys.path.pop(0)
  from iotbx import data_plots
  loggraph1 = """\
$TABLE: Resolution shell statistics
$GRAPHS
:R-free vs. resolution
:A:1,3:
:FOM vs. resolution
:A:1,4:
$$
1/resol^2  Nrefl      R-free     FOM       $$
$$
0.02       2004       0.25       0.89
0.04       2084       0.23       0.88
0.06       2037       0.27       0.83
0.08       1949       0.28       0.75
0.1        1783       0.38       *
$$
"""
  app = wx.App(0)
  frame = loggraph(parent=None,
    title="Loggraph test",
    tables=None,
    processed_lines=loggraph1.splitlines())
  frame.Show()
  frame2 = wx.Frame(None, -1, "Heat map")
  panel = wx.Panel(frame2, -1)
  sizer2 = wx.BoxSizer(wx.VERTICAL)
  frame2.SetSizer(sizer2)
  sizer2.Add(panel, 1, wx.EXPAND)
  implot = image_plot(panel)
  implot.show_plot(values=values, x_label="RMSD", y_label="resolution",
    cmap="autumn")
  sizer2.Layout()
  frame2.Fit()
  frame2.Show()
  app.MainLoop()

if __name__ == "__main__" :
  exercise()
