
from wxtbx import bitmaps
import wx
from libtbx import object_oriented_patterns as oop
from libtbx import adopt_init_args
from libtbx.utils import Sorry
import sys
import os
if (sys.version_info[2] >= 6) :
  import warnings
  warnings.simplefilter('ignore', DeprecationWarning)

class plot_container (wx.BoxSizer) :
  def __init__ (self,
                parent,
                figure_size=(8,6),
                facecolor='white',
                transparent=False,
                handle_left_click=False) :
    wx.BoxSizer.__init__(self, wx.VERTICAL)
    from matplotlib.backends.backend_wxagg import Toolbar
    from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, FigureManager
    import matplotlib.ticker
    import matplotlib.cm
    import matplotlib.figure
    import matplotlib.font_manager
    self.parent = parent
    self.figure = matplotlib.figure.Figure(figure_size, 72, linewidth=0,
      facecolor=facecolor)
    if transparent :
      self.figure.figurePatch.set_alpha(0.0)
    self.canvas = FigureCanvasWxAgg(self.parent, -1, self.figure)
    self.canvas.toolbar = oop.null()
    self.figmgr = FigureManager(self.canvas, 1, self)
    self.Add(self.canvas, 1, wx.EXPAND|wx.ALL)
    self.font = matplotlib.font_manager.FontProperties(
      family = ["Courier", "Monaco", "monospace"],
      weight = "normal",
      size   = 10
    )
    self.label_font = matplotlib.font_manager.FontProperties(
      family = ["Helvetica", "sans-serif"],
      weight = "bold",
      size = 11
    )
    self.legend_font = matplotlib.font_manager.FontProperties(
      family = ["Helvetica", "sans-serif"],
      weight = "normal",
      size = 10
    )
    self.value_label_font = matplotlib.font_manager.FontProperties(
      family = ["Courier", "Monaco", "monospace"],
      weight = "normal",
      size   = self.font_size
    )
    self.axis_label_font = matplotlib.font_manager.FontProperties(
      family = ["Helvetica", "sans-serif"],
      weight = "bold",
      size   = self.font_size
    )
    self.text_font = matplotlib.font_manager.FontProperties(
      family = ["Helvetica", "sans-serif"],
      weight = "normal",
      size   = self.font_size
    )
    self.title_font = matplotlib.font_manager.FontProperties(
      family = ["Helvetica", "sans-serif"],
      weight = "normal",
      size   = self.title_font_size
    )
    self.null_fmt = matplotlib.ticker.NullFormatter()
    if self.handle_left_click :
      self.canvas.mpl_connect("button_release_event", self.OnClick)
    else :
      self.canvas.Bind(wx.EVT_CONTEXT_MENU, self.OnRightClick, self.canvas)

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

class iotbx_data_plot_base (plot_container) :
  def __init__ (self,
                parent,
                tables,
                size=(640,480),
                randomize_colors=False,
                use_points=True,
                transparent=True,
                handle_left_click=False) :
    adopt_init_args(self, locals())
    (x, y, w, h) = tuple(wx.GetClientDisplayRect())
    (width, height) = size
    fig_w = float(width) / 72.0
    fig_h = float(height) / 72.0
    if (fig_w * 72) > (w - 20) :
      fig_w = int(math.floor((w-40) / 72))
    if (fig_h * 72) > (h - 120) :
      fig_h = int(math.floor((h-160) / 72))
    plot_container.__init__(self,
      parent=parent,
      figure_size=(fig_w, fig_h),
      transparent=transparent)
    self.p = self.figure.add_subplot(111)
    self.plot_type = None

  def set_plot (self, graph_name=None, table_name=None) :
    table = None
    if table_name is None :
      table = self.tables[0]
    else :
      for t in self.tables :
        if t.title == table_name :
          table = t
    if table is not None :
      self.plot_type = getattr(table, "plot_type", "GRAPH")
      graph = table.get_graph(graph_name)
      self.show_plot(graph)

  def show_plot (self, graph, line_width=1, show_points=True, show_grid=True) :
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
    if self.use_points :
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
    self.p.set_title(graph.name, fontproperties=self.title_font,
      horizontalalignment=self.title_alignment)
    self.canvas.draw()
    self.parent.Refresh()

  def format_x_axis (self) :
    if self.tables[0].x_is_inverse_d_min :
      xdata = self.tables[0].get_x_as_resolution()
      self.p.get_axes().set_xlabel("Resolution",
        fontproperties=self.axis_label_font)
      marks = [5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.25, 1.1,
                1, 0.9, 0.8, 0.7, 0.6, 0.5]
      xticks = []
      xticklabels = []
      start = 0
      xticklabels.append("%.2f" % xdata[0])
      xticks.append(1.0 / (xdata[0]**2))
      for x in marks :
        if x >= xdata[-1] :
          xticklabels.append(str(x))
          xticks.append(1.0 / (x**2))
      self.p.get_axes().set_xticks(xticks)
      self.p.get_axes().set_xticklabels(xticklabels)
    else :
      if self.graph.x_axis_label is not None :
        self.p.get_axes().set_xlabel(self.graph.x_axis_label,
          fontproperties=self.axis_label_font)
      else :
        self.p.get_axes().set_xlabel(self.graph.x_label,
          fontproperties=self.axis_label_font)
    for ticklabel in self.p.get_axes().get_xticklabels() :
      ticklabel.set_fontproperties(self.value_label_font)

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
      ticklabel.set_fontproperties(self.value_label_font)
    if self.graph.y_axis_label is not None :
      self.p.get_axes().set_ylabel(self.graph.y_axis_label,
          fontproperties=self.axis_label_font)

  def format_labels (self) :
    self.figure.legend(self.p.lines, self.graph.y_labels,
      prop=self.text_font)
    self.format_x_axis()
    self.format_y_axis()

  def show_grid (self, show=True) :
    if show :
      self.p.get_axes().grid(True, color="0.75")
    else :
      self.p.get_axes().grid(False)
    self.canvas.draw()
    self.Refresh()

  def OnRightClick (self, event) :
    menu = wx.Menu()
    menu.add_command("Save image", self.OnSave)
    self.parent.PopupMenu(menu)

  def OnSave (self, event=None) :
    self.save_image()

class plot_frame (wx.Frame) :
  controls_on_top = True
  show_controls_default = True
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.setup_toolbar()
    self.toolbar.Realize()
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.top_panel = wx.Panel(parent=self, style=wx.SUNKEN_BORDER)
    self.draw_top_panel()
    self.plot_panel = self.create_plot_panel()
    if self.controls_on_top :
      self.sizer.Add(self.top_panel, 0, wx.EXPAND|wx.ALL)
      self.sizer.Add(self.plot_panel, 1, wx.EXPAND|wx.ALL)
    else :
      self.sizer.Add(self.plot_panel, 1, wx.EXPAND|wx.ALL)
      self.sizer.Add(self.top_panel, 0, wx.EXPAND|wx.ALL)
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
       bitmaps.fetch_icon_bitmap("apps/advancedsettings.png"),
       self.OnToggleControls),
      ("Save",
       bitmaps.fetch_icon_bitmap("actions/save_all.png"),
       self.OnSave),
      #("Print",
      # bitmaps.fetch_icon_bitmap("devices/printer1.png"),
      # self.OnPrint),
    ]
    tb = wx.ToolBar(self, style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    tb.SetToolBitmapSize((32,32))
    self.SetToolBar(tb)
    for (button_name, button_bitmap, button_function) in tb_buttons :
      tool_button = tb.AddLabelTool(-1, button_name, bmp,
          shortHelp=help_string, kind=wx.ITEM_NORMAL)
      self.Bind(wx.EVT_MENU, button_function, tool_button)
    tb.Realize()

  def draw_top_panel (self) :
    pass

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
    cp.Add(grid, 0, wx.EXPAND)
    szr = wx.BoxSizer(wx.HORIZONTAL)
    cp.Add(szr, 0, wx.EXPAND)
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
      self.Refresh()

  #--- EVENTS
  def OnShowTables (self, event) :
    if self.table_frame is None :
      self.table_frame = TableFrame(self, "Tables", self.tables)
    self.table_frame.Show()
    self.table_frame.Raise()

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
    self.plot.use_points = event.GetEventObject().GetValue()
    self.switch_plot()
