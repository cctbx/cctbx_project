from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 07/21/2017
Last Changed: 02/12/2018
Description : IOTA image-tracking GUI module
'''

import os
import wx
from wxtbx import bitmaps
import wx.lib.agw.ultimatelistctrl as ulc
import wx.lib.mixins.listctrl as listmix
import numpy as np

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

from iotbx import phil as ip

from iota.components.iota_dialogs import DIALSSpfDialog
from iota.components.iota_utils import InputFinder
from iota.components.iota_dials import phil_scope, IOTADialsProcessor
import iota.components.iota_threads as thr
import iota.components.iota_controls as ct

import time
assert time # going to keep time around for testing

# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  norm_font_size = 10
  button_font_size = 12
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = 'python'
elif wx.Platform == '__WXMAC__':
  norm_font_size = 12
  button_font_size = 14
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = "Python"
elif (wx.Platform == '__WXMSW__'):
  norm_font_size = 9
  button_font_size = 11
  LABEL_SIZE = 11
  CAPTION_SIZE = 9

ginp = InputFinder()
user = os.getlogin()

default_target = '\n'.join(['verbosity=10',
                            'spotfinder {',
                            '  threshold {',
                            '    dispersion {',
                            '      gain = 1',
                            '      sigma_strong = 3',
                            '      global_threshold = 0',
                            '      }',
                            '   }',
                            '}',
                            'geometry {',
                            '  detector {',
                            '    distance = None',
                            '    slow_fast_beam_centre = None',
                            '  }',
                            '}',
                            'indexing {',
                            '  known_symmetry {',
                            '    space_group = None',
                            '    unit_cell = None',
                            '  }',
                            '  refinement_protocol {',
                            '    n_macro_cycles = 1',
                            '    d_min_start = 2.0',
                            '  }',
                            '  basis_vector_combinations.max_combinations = 10',
                            '  stills { ',
                            '    indexer = stills',
                            '    method_list = fft1d real_space_grid_search',
                            '  }',
                            '}',
                            'refinement {',
                            '  parameterisation {',
                            '    beam.fix=all',
                            '    detector  {',
                            '      fix=all',
                            '      hierarchy_level=0',
                            '    }',
                            '    auto_reduction {',
                            '      action=fix',
                            '      min_nref_per_parameter=1',
                            '    }',
                            '  }',
                            '  reflections {',
                            '    outlier.algorithm=null',
                            '    weighting_strategy  {',
                            '      override=stills',
                            '      delpsi_constant=1000000',
                            '    }',
                            '  }',
                            '}',
                            'integration {',
                            '  integrator=stills',
                            '  profile.fitting=False',
                            '  background {',
                            '    simple {',
                            '      outlier {',
                            '        algorithm = null',
                            '      }',
                            '    }',
                            '  }',
                            '}',
                            'profile {',
                            '  gaussian_rs {',
                            '    min_spots.overall = 0',
                            '  }',
                            '}'
                            ])

class TrackChart(wx.Panel):
  def __init__(self, parent, main_window):
    wx.Panel.__init__(self, parent, size=(100, 100))
    self.main_window = main_window

    self.main_box = wx.StaticBox(self, label='Spotfinding Chart')
    self.main_fig_sizer = wx.StaticBoxSizer(self.main_box, wx.VERTICAL)
    self.SetSizer(self.main_fig_sizer)

    self.track_figure = Figure()
    self.track_axes = self.track_figure.add_subplot(111)
    self.track_axes.set_ylabel('Found Spots')
    self.track_axes.set_xlabel('Frame')

    self.track_figure.set_tight_layout(True)
    self.track_canvas = FigureCanvas(self, -1, self.track_figure)
    self.track_axes.patch.set_visible(False)

    self.plot_sb = wx.Slider(self, minValue=0, maxValue=1)
    self.plot_sb.Hide()

    self.main_fig_sizer.Add(self.track_canvas, 1, wx.EXPAND)
    self.main_fig_sizer.Add(self.plot_sb, flag=wx.EXPAND)

    # Scroll bar binding
    self.Bind(wx.EVT_SCROLL, self.onScroll, self.plot_sb)

    # Plot bindings
    self.track_figure.canvas.mpl_connect('button_press_event', self.onPress)
    # self.track_figure.canvas.mpl_connect('scroll_event', self.onScroll)

    self.reset_chart()



  def onSelect(self, xmin, xmax):
    ''' Called when SpanSelector is used (i.e. click-drag-release); passes on
        the boundaries of the span to tracker window for selection and
        display of the selected images '''
    if self.selector == 'select':
      self.select_span.set_visible(True)
      self.patch_x = int(xmin)
      self.patch_x_last = int(xmax) + 1
      self.patch_width = self.patch_x_last - self.patch_x
      self.bracket_set = True


      self.main_window.update_image_list()
      self.main_window.tracker_panel.image_list_panel.Show()
      self.main_window.tracker_panel.Layout()

    elif self.selector == 'zoom':
      if (int(xmax) - int(xmin) >= 5):
        self.x_min = int(xmin)
        self.x_max = int(xmax)
        self.plot_zoom = True
        self.max_lock = False
        self.chart_range = int(self.x_max - self.x_min)
        self.main_window.tracker_panel.chart_window.toggle.SetValue(True)
        self.main_window.tracker_panel.chart_window.toggle_boxes(flag_on=True)
        self.main_window.tracker_panel.chart_window.ctr.SetValue(self.chart_range)
        sb_center = self.x_min + self.chart_range / 2

        self.plot_sb.SetValue(sb_center)
        self.plot_sb.Show()
        self.draw_plot()

        if self.bracket_set:
          self.bracket_set = False
          self.main_window.tracker_panel.image_list_panel.Hide()
          self.main_window.tracker_panel.Layout()

  def onScroll(self, e):
    sb_center = self.plot_sb.GetValue()
    half_span = (self.x_max - self.x_min) / 2
    if sb_center - half_span == 0:
      self.x_min = 0
      self.x_max = half_span * 2
    else:
      self.x_min = sb_center - half_span
      self.x_max = sb_center + half_span

    if self.plot_sb.GetValue() == self.plot_sb.GetMax():
      self.max_lock = True
    else:
      self.max_lock = False

    self.draw_plot()


  def onPress(self, e):
    ''' If left mouse button is pressed, activates the SpanSelector;
    otherwise, makes the span invisible and sets the toggle that clears the
    image list; if shift key is held, does this for the Selection Span,
    otherwise does this for the Zoom Span '''
    if e.button != 1:
      self.zoom_span.set_visible(False)
      self.select_span.set_visible(False)
      self.bracket_set = False
      self.plot_zoom = False
      self.plot_sb.Hide()
      self.draw_plot()

      # Hide list of images
      self.main_window.tracker_panel.image_list_panel.Hide()
      self.main_window.tracker_panel.chart_window.toggle.SetValue(False)
      self.main_window.tracker_panel.chart_window.toggle_boxes(flag_on=False)
      self.main_window.tracker_panel.Layout()
    else:
      if e.key == 'shift':
        self.selector = 'select'
        self.zoom_span.set_visible(False)
        self.select_span.set_visible(True)
      else:
        self.selector = 'zoom'
        self.zoom_span.set_visible(True)
        self.select_span.set_visible(False)

  def reset_chart(self):
    self.track_axes.clear()
    self.track_figure.patch.set_visible(False)
    self.track_axes.patch.set_visible(False)

    self.xdata = []
    self.ydata = []
    self.x_min = 0
    self.x_max = 1
    self.y_max = 1
    self.bracket_set = False
    self.button_hold = False
    self.plot_zoom = False
    self.chart_range = None
    self.selector = None
    self.max_lock = True
    self.patch_x = 0
    self.patch_x_last = 1
    self.patch_width = 1
    self.start_edge = 0
    self.end_edge = 1

    self.acc_plot = self.track_axes.plot([], [], 'o', color='#4575b4')[0]
    self.rej_plot = self.track_axes.plot([], [], 'o', color='#d73027')[0]
    self.bragg_line = self.track_axes.axhline(0, c='#4575b4', ls=':', alpha=0)
    self.highlight = self.track_axes.axvspan(0.5, 0.5, ls='--', alpha=0,
                                             fc='#deebf7', ec='#2171b5')
    self.track_axes.set_autoscaley_on(True)

    self.select_span = SpanSelector(ax=self.track_axes, onselect=self.onSelect,
                                    direction='horizontal',
                                    rectprops=dict(alpha=0.5, ls = ':',
                                            fc='#deebf7', ec='#2171b5'))
    self.select_span.set_active(False)

    self.zoom_span = SpanSelector(ax=self.track_axes, onselect=self.onSelect,
                                  direction='horizontal',
                                  rectprops=dict(alpha=0.5, ls = ':',
                                                 fc='#ffffd4', ec='#8c2d04'))
    self.zoom_span.set_active(False)

  def draw_bragg_line(self):
    min_bragg = self.main_window.tracker_panel.min_bragg.ctr.GetValue()
    if min_bragg > 0:
     self.bragg_line.set_alpha(1)
    else:
      self.bragg_line.set_alpha(0)
    self.bragg_line.set_ydata(min_bragg)

  def draw_plot(self, new_x=None, new_y=None):
    min_bragg = self.main_window.tracker_panel.min_bragg.ctr.GetValue()

    if new_x is None:
      new_x = []
    if new_y is None:
      new_y = []

    nref_x = np.append(self.xdata, np.array(new_x).astype(np.double))
    nref_y = np.append(self.ydata, np.array(new_y).astype(np.double))
    self.xdata = nref_x
    self.ydata = nref_y
    nref_xy = zip(nref_x, nref_y)
    all_acc = [i[0] for i in nref_xy if i[1] >= min_bragg]
    all_rej = [i[0] for i in nref_xy if i[1] < min_bragg]

    if nref_x != [] and nref_y != []:
      if self.plot_zoom:
        if self.max_lock:
          self.x_max = np.max(nref_x)
          self.x_min = self.x_max - self.chart_range
      else:
        self.x_min = -1
        self.x_max = np.max(nref_x) + 1

      if min_bragg > np.max(nref_y):
        self.y_max = min_bragg + int(0.1 * min_bragg)
      else:
        self.y_max = np.max(nref_y) + int(0.1 * np.max(nref_y))

      self.track_axes.set_xlim(self.x_min, self.x_max)
      self.track_axes.set_ylim(0, self.y_max)

    else:
      self.x_min = -1
      self.x_max = 1

    acc = [i for i in all_acc if i > self.x_min and i < self.x_max]
    rej = [i for i in all_rej if i > self.x_min and i < self.x_max]

    self.acc_plot.set_xdata(nref_x)
    self.rej_plot.set_xdata(nref_x)
    self.acc_plot.set_ydata(nref_y)
    self.rej_plot.set_ydata(nref_y)
    self.acc_plot.set_markevery(acc)
    self.rej_plot.set_markevery(rej)

    self.Layout()

    count = '{}'.format(len([i for i in nref_xy if i[1] >= min_bragg]))
    self.main_window.tracker_panel.count_txt.SetLabel(count)
    self.main_window.tracker_panel.status_sizer.Layout()

    # Set up scroll bar
    if len(self.xdata) > 0:
      self.plot_sb.SetMax(np.max(nref_x))
      if self.max_lock:
        self.plot_sb.SetValue(self.plot_sb.GetMax())

    # Draw extended plots
    self.track_axes.draw_artist(self.acc_plot)
    self.track_axes.draw_artist(self.rej_plot)

class ImageList(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent, size=(450, -1))

    self.main_box = wx.StaticBox(self, label='Hits')
    self.main_sizer = wx.StaticBoxSizer(self.main_box, wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    #self.image_list = FileListCtrl(self)
    self.image_list = VirtualListCtrl(self)
    self.main_sizer.Add(self.image_list, 1, flag=wx.EXPAND)


# class FileListCtrl(ct.CustomImageListCtrl, listmix.ColumnSorterMixin):
#   def __init__(self, parent):
#     ct.CustomImageListCtrl.__init__(self, parent=parent)

class FileListCtrl(ct.CustomImageListCtrl, listmix.ColumnSorterMixin):
  def __init__(self, parent):
    ct.CustomImageListCtrl.__init__(self, parent=parent)
    self.parent = parent
    self.tracker_panel = parent.GetParent()
    self.main_window = parent.GetParent().GetParent().GetParent()


    # Generate columns
    self.ctr.InsertColumn(0, "#", width=50)
    self.ctr.InsertColumn(1, "File", width=150)
    self.ctr.InsertColumn(2, "No. Spots", width=3)
    #self.ctr.setResizeColumn(1)
    #self.ctr.setResizeColumn(2)

    # Bindings
    self.Bind(ulc.EVT_LIST_ITEM_SELECTED, self.OnItemSelected)
    self.Bind(ulc.EVT_LIST_ITEM_DESELECTED, self.OnItemDeselected)


  def OnItemSelected(self, e):
    self.main_window.toolbar.EnableTool(
      self.main_window.tb_btn_view.GetId(), True)

  def OnItemDeselected(self, e):
    ctr = e.GetEventObject()
    if ctr.GetSelectedItemCount() <= 0:
      self.main_window.toolbar.EnableTool(
        self.main_window.tb_btn_view.GetId(), False)

  def GetListCtrl(self):
    return self.ctr

  def OnColClick(self, e):
    print "column clicked"
    e.Skip()

  def instantiate_sorting(self, data):
    self.itemDataMap = data
    listmix.ColumnSorterMixin.__init__(self, 3)
    self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.ctr)

  def add_item(self, img_idx, img, n_obs):
    item = ct.FileListItem(path=img,
                           items={
                             'n_obs':n_obs,
                             'img_idx':img_idx
                           })

    # Insert list item
    idx = self.ctr.InsertStringItem(self.ctr.GetItemCount() + 1,
                                    str(item.img_idx))
    self.ctr.SetStringItem(idx, 1, os.path.basename(item.path))
    self.ctr.SetStringItem(idx, 2, str(item.n_obs))

    # Record index in item data
    item.id = idx

    # Attach data object to item
    self.ctr.SetItemData(idx, item)

    # Resize columns to fit content
    #self.ctr.SetColumnWidth(0, width=50)
    #self.ctr.SetColumnWidth(2, width=150)
    #self.ctr.SetColumnWidth(1, width=-3)

    self.tracker_panel.Layout()

  def delete_item(self, index):
    item = self.ctr.GetItemData(index)
    self.ctr.DeleteItem(index)

    # Refresh widget and list item indices
    for i in range(self.ctr.GetItemCount()):
      item_data = self.ctr.GetItemData(i)
      item_data.id = i
      self.ctr.SetItemData(i, item_data)

  def delete_selected(self, idxs):
    counter = 0
    for idx in idxs:
      self.delete_item(idx - counter)
      counter += 1

  def delete_all(self):
    for idx in range(self.ctr.GetItemCount()):
      self.delete_item(index=0)

class VirtualListCtrl(ct.VirtualImageListCtrl):
  def __init__(self, parent):
    ct.VirtualImageListCtrl.__init__(self, parent=parent)
    self.parent = parent
    self.tracker_panel = parent.GetParent()
    self.main_window = parent.GetParent().GetParent().GetParent()

    # Generate columns
    self.ctr.InsertColumn(0, "#")
    self.ctr.InsertColumn(1, "File")
    self.ctr.InsertColumn(2, "No. Spots")
    self.ctr.setResizeColumn(1)
    self.ctr.setResizeColumn(2)

    # Bindings
    self.Bind(ulc.EVT_LIST_CACHE_HINT, self.OnSelection, self.ctr)
    self.Bind(wx.EVT_LIST_COL_CLICK, self.ctr.OnColClick, self.ctr)

  def initialize_data_map(self, data):
    self.ctr.InitializeDataMap(data)
    self.ctr.SetColumnWidth(0, width=50)
    self.ctr.SetColumnWidth(2, width=150)
    self.ctr.SetColumnWidth(1, width=-3)

  def OnSelection(self, e):
    ctr = e.GetEventObject()
    if ctr.GetSelectedItemCount() <= 0:
      self.main_window.tracker_panel.btn_view_sel.Disable()
    else:
      self.main_window.tracker_panel.btn_view_sel.Enable()

class TrackerPanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent)
    self.parent = parent

    self.main_sizer = wx.GridBagSizer(10, 10)
    self.SetSizer(self.main_sizer)


    # Status box
    self.status_panel = wx.Panel(self)
    self.status_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    self.status_sizer.AddGrowableCol(0)
    self.status_box = wx.StaticBox(self.status_panel, label='Status')
    self.status_box_sizer = wx.StaticBoxSizer(self.status_box, wx.HORIZONTAL)
    self.status_txt = wx.StaticText(self.status_panel, label='')
    self.status_box_sizer.Add(self.status_txt, flag=wx.ALL | wx.ALIGN_CENTER,
                              border=10)
    self.status_sizer.Add(self.status_box_sizer, flag=wx.EXPAND)
    self.status_panel.SetSizer(self.status_sizer)
    self.count_box = wx.StaticBox(self.status_panel, label='Hit Count')
    self.count_box_sizer = wx.StaticBoxSizer(self.count_box, wx.HORIZONTAL)
    self.count_txt = wx.StaticText(self.status_panel, label='')
    self.count_box_sizer.Add(self.count_txt, flag=wx.ALL | wx.ALIGN_CENTER,
                             border=10)
    font = wx.Font(20, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    self.count_txt.SetFont(font)

    self.status_sizer.Add(self.count_box_sizer, flag=wx.EXPAND)
    self.main_sizer.Add(self.status_panel, pos=(0, 0), span=(1, 2),
                        flag=wx.EXPAND | wx.ALL, border=10)

    # Put in chart
    self.graph_panel = wx.Panel(self)
    self.graph_sizer = wx.GridBagSizer(5, 5)

    self.chart = TrackChart(self.graph_panel, main_window=parent)
    self.min_bragg = ct.SpinCtrl(self.graph_panel,
                                 label='Min. Bragg spots',
                                 label_size=wx.DefaultSize,
                                 ctrl_size=(100, -1),
                                 ctrl_value=10,
                                 ctrl_min=0)
    self.chart_window = ct.SpinCtrl(self.graph_panel,
                                    label_size=wx.DefaultSize,
                                    checkbox=True,
                                    checkbox_state=False,
                                    checkbox_label='Finite chart window',
                                    ctrl_size=(100, -1),
                                    ctrl_value=100,
                                    ctrl_min=10,
                                    ctrl_step=10)
    self.spf_options = wx.Button(self.graph_panel,
                                 label='Spotfinding Options...')
    self.graph_sizer.Add(self.chart, flag=wx.EXPAND, pos=(0, 0), span=(1, 3))
    self.graph_sizer.Add(self.min_bragg, flag=wx.ALIGN_LEFT, pos=(1, 0))
    self.graph_sizer.Add(self.chart_window, flag=wx.ALIGN_CENTER, pos=(1, 1))
    self.graph_sizer.Add(self.spf_options, flag=wx.ALIGN_RIGHT, pos=(1, 2))

    self.graph_sizer.AddGrowableRow(0)
    self.graph_sizer.AddGrowableCol(1)
    self.graph_panel.SetSizer(self.graph_sizer)

    # List of images
    self.image_list_panel = wx.Panel(self)
    self.image_list_sizer = wx.BoxSizer(wx.VERTICAL)
    self.image_list = ImageList(self.image_list_panel)
    self.btn_view_sel = wx.Button(self.image_list_panel, label='View Selected')
    self.btn_view_all = wx.Button(self.image_list_panel, label='View All')
    self.btn_wrt_file = wx.Button(self.image_list_panel, label='Write to File')
    self.btn_view_sel.Disable()
    self.btn_wrt_file.Disable()
    self.btn_view_all.Disable()

    btn_sizer = wx.FlexGridSizer(1, 3, 0, 5)
    btn_sizer.Add(self.btn_wrt_file, flag=wx.ALIGN_LEFT)
    btn_sizer.Add(self.btn_view_all, flag=wx.ALIGN_RIGHT)
    btn_sizer.Add(self.btn_view_sel, flag=wx.ALIGN_RIGHT)
    btn_sizer.AddGrowableCol(0)

    self.image_list_sizer.Add(self.image_list, 1, flag=wx.EXPAND)
    self.image_list_sizer.Add(btn_sizer, flag=wx.ALIGN_RIGHT | wx.EXPAND)
    self.image_list_panel.SetSizer(self.image_list_sizer)

    self.main_sizer.Add(self.graph_panel, pos=(1, 0),
                        flag=wx.EXPAND | wx.LEFT | wx.BOTTOM, border=10)
    self.main_sizer.Add(self.image_list_panel, pos=(1, 1),
                        flag=wx.EXPAND | wx.LEFT | wx.BOTTOM | wx.RIGHT,
                        border=10)
    self.main_sizer.AddGrowableCol(0)
    self.main_sizer.AddGrowableRow(1)
    self.image_list_panel.Hide()

class TrackerWindow(wx.Frame):
  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(1500, 600))
    self.parent = parent
    self.term_file = os.path.join(os.curdir, '.terminate_image_tracker')
    self.info_file = os.path.join(os.curdir, '.spotfinding_info')
    self.folder_file = os.path.join(os.curdir, '.data_info')

    self.reset_spotfinder()

    # Setup main sizer
    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    # Setup toolbar
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS | wx.TB_TEXT)
    quit_bmp = bitmaps.fetch_icon_bitmap('actions', 'exit')
    self.tb_btn_quit = self.toolbar.AddLabelTool(wx.ID_EXIT, label='Quit',
                                                 bitmap=quit_bmp,
                                                 shortHelp='Quit',
                                                 longHelp='Quit image tracker')
    self.toolbar.AddSeparator()
    open_bmp = bitmaps.fetch_icon_bitmap('actions', 'open')
    self.tb_btn_open = self.toolbar.AddLabelTool(wx.ID_ANY, label='Open',
                                                bitmap=open_bmp,
                                                shortHelp='Open',
                                                longHelp='Open folder')
    rec_bmp = bitmaps.fetch_icon_bitmap('actions', 'quick_restart')
    self.tb_btn_restore = self.toolbar.AddLabelTool(wx.ID_ANY, label='Restore',
                                                bitmap=rec_bmp,
                                                shortHelp='Restore',
                                                longHelp='Restore aborted run')
    run_bmp = bitmaps.fetch_icon_bitmap('actions', 'run')
    self.tb_btn_run = self.toolbar.AddLabelTool(wx.ID_ANY, label='Run',
                                                bitmap=run_bmp,
                                                shortHelp='Run',
                                                longHelp='Run Spotfinding')
    stop_bmp = bitmaps.fetch_icon_bitmap('actions', 'stop')
    self.tb_btn_stop = self.toolbar.AddLabelTool(wx.ID_ANY, label='Stop',
                                                bitmap=stop_bmp,
                                                shortHelp='Stop',
                                                longHelp='Stop Spotfinding')
    # self.toolbar.AddSeparator()
    # span_view = bitmaps.fetch_custom_icon_bitmap('span_view')
    # self.tb_btn_view = self.toolbar.AddLabelTool(wx.ID_ANY, label='View',
    #                                              bitmap=span_view,
    #                                              kind=wx.ITEM_CHECK,
    #                                              shortHelp='Select to View',
    #                                              longHelp='Select images to view')
    # span_zoom = bitmaps.fetch_custom_icon_bitmap('span_zoom')
    # self.tb_btn_zoom = self.toolbar.AddLabelTool(wx.ID_ANY, label='Zoom In',
    #                                              bitmap=span_zoom,
    #                                              kind=wx.ITEM_CHECK,
    #                                              shortHelp='Zoom In',
    #                                              longHelp='Zoom in on chart')
    # self.toolbar.ToggleTool(self.tb_btn_view.GetId(), True)

    # view_bmp = bitmaps.fetch_custom_icon_bitmap('image_viewer32')
    # self.tb_btn_view = self.toolbar.AddLabelTool(wx.ID_ANY, label='View',
    #                                             bitmap=view_bmp,
    #                                             shortHelp='View images',
    #                                             longHelp='View selected images')
    # self.toolbar.EnableTool(self.tb_btn_view.GetId(), False)

    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), False)

    if os.path.isfile(self.folder_file) and os.path.isfile(self.info_file):
      self.toolbar.EnableTool(self.tb_btn_restore.GetId(), True)
    else:
      self.toolbar.EnableTool(self.tb_btn_restore.GetId(), False)

    self.toolbar.Realize()

    # Setup timer
    self.timer = wx.Timer(self)

    self.tracker_panel = TrackerPanel(self)
    self.data_dict = self.tracker_panel.image_list.image_list.ctr.data.copy()
    self.img_list_initialized = False

    self.main_sizer.Add(self.tracker_panel, 1, wx.EXPAND)

    # Generate default PHIL file and instantiate DIALS stills processor
    default_phil = ip.parse(default_target)
    self.phil = phil_scope.fetch(source=default_phil)
    self.params = self.phil.extract()
    self.params.output.strong_filename = None
    self.phil = self.phil.format(python_object=self.params)

    # Bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onGetImages, self.tb_btn_open)
    self.Bind(wx.EVT_TOOL, self.onRunSpotfinding, self.tb_btn_run)
    self.Bind(wx.EVT_TOOL, self.onRestoreRun, self.tb_btn_restore)
    self.Bind(wx.EVT_TOOL, self.onStop, self.tb_btn_stop)
    self.Bind(wx.EVT_BUTTON, self.onSelView, self.tracker_panel.btn_view_sel)
    self.Bind(wx.EVT_BUTTON, self.onWrtFile, self.tracker_panel.btn_wrt_file)
    self.Bind(wx.EVT_BUTTON, self.onAllView, self.tracker_panel.btn_view_all)

    # Spotfinder / timer bindings
    self.Bind(thr.EVT_SPFDONE, self.onSpfOneDone)
    self.Bind(thr.EVT_SPFALLDONE, self.onSpfAllDone)
    self.Bind(thr.EVT_SPFTERM, self.onSpfTerminated)
    self.Bind(wx.EVT_TIMER, self.onTimer, id=self.timer.GetId())

    # Settings bindings
    self.Bind(wx.EVT_BUTTON, self.onSpfOptions, self.tracker_panel.spf_options)
    self.Bind(wx.EVT_SPINCTRL, self.onMinBragg,
              self.tracker_panel.min_bragg.ctr)
    self.Bind(wx.EVT_SPINCTRL, self.onChartRange,
              self.tracker_panel.chart_window.ctr)

  def reset_spotfinder(self):
    self.frame_count = []
    self.obs_counts = []
    self.done_list = []
    self.data_list = []
    self.new_frames = []
    self.new_counts = []
    self.spotfinding_info = []
    self.all_info = []
    self.current_min_bragg = 0
    self.waiting = False
    self.terminated = False

  def onWrtFile(self, e):
    idxs = []
    listctrl = self.tracker_panel.image_list.image_list.ctr
    if listctrl.GetSelectedItemCount() == 0:
      for index in range(listctrl.GetItemCount()):
        idxs.append(index)
    else:
      index = listctrl.GetFirstSelected()
      idxs.append(index)
      while len(idxs) != listctrl.GetSelectedItemCount():
        index = listctrl.GetNextSelected(index)
        idxs.append(index)
    self.write_images_to_file(idxs=idxs)

  def onSelView(self, e):
    idxs = []
    listctrl = self.tracker_panel.image_list.image_list.ctr
    if listctrl.GetSelectedItemCount() == 0:
      return

    index = listctrl.GetFirstSelected()
    idxs.append(index)
    while len(idxs) != listctrl.GetSelectedItemCount():
      index = listctrl.GetNextSelected(index)
      idxs.append(index)
    self.view_images(idxs=idxs)

  def onAllView(self, e):
    listctrl = self.tracker_panel.image_list.image_list.ctr
    idxs = range(listctrl.GetItemCount())
    self.view_images(idxs=idxs)

  def write_images_to_file(self, idxs):
    # Determine param filepath
    save_dlg = wx.FileDialog(self,
                             message="Save Image Paths to File",
                             defaultDir=os.curdir,
                             defaultFile="*.lst",
                             wildcard="*",
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                             )
    if save_dlg.ShowModal() == wx.ID_OK:
      script_filepath = save_dlg.GetPath()
      file_list = [self.data_dict[idx][1] for idx in idxs]
      with open(script_filepath, 'w') as img_file:
        file_list_string = '\n'.join(file_list)
        img_file.write(file_list_string)

  def view_images(self, idxs):
    file_list = [self.data_dict[idx][1] for idx in idxs]
    file_string = ' '.join(file_list)

    viewer = thr.ImageViewerThread(self,
                                   file_string=file_string)
    viewer.start()

  def onStop(self, e):
    self.terminated = True
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), False)
    with open(self.term_file, 'w') as tf:
      tf.write('')
    self.msg = 'Stopping...'

  def remove_term_file(self):
    try:
      os.remove(self.term_file)
    except Exception:
      pass

  def onGetImages(self, e):
    ''' Select folder to watch for incoming images '''
    open_dlg = wx.DirDialog(self, "Choose the data folder:",
                            style=wx.DD_DEFAULT_STYLE)
    if open_dlg.ShowModal() == wx.ID_OK:
      self.data_folder = open_dlg.GetPath()
      open_dlg.Destroy()
      self.remove_term_file()
      self.reset_spotfinder()

      self.tracker_panel.chart.reset_chart()
      self.tracker_panel.spf_options.Enable()
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), True)
      # self.toolbar.EnableTool(self.tb_btn_calc.GetId(), True)
      timer_txt = '[ ------ ]'
      self.msg = 'Ready to track images in {}'.format(self.data_folder)
      self.tracker_panel.status_txt.SetLabel('{} {}'.format(timer_txt, self.msg))

      self.toolbar.EnableTool(self.tb_btn_restore.GetId(), False)
      try:
        os.remove(self.folder_file)
        os.remove(self.info_file)
      except Exception, e:
        pass
    else:
      open_dlg.Destroy()
      return

  def onMinBragg(self, e):
    self.tracker_panel.chart.draw_bragg_line()

  def onChartRange(self, e):
    chart_range = self.tracker_panel.chart_window.ctr.GetValue()
    self.tracker_panel.chart.chart_range = chart_range
    self.tracker_panel.chart.max_lock = True

  def onSpfOptions(self, e):
    spf_dlg = DIALSSpfDialog(self,
                             phil=self.phil,
                             title='DIALS Spotfinding Options')
    if (spf_dlg.ShowModal() == wx.ID_OK):
      self.phil = self.phil.fetch(source=spf_dlg.spf_phil)
    spf_dlg.Destroy()


  def onRunSpotfinding(self, e):
    self.terminated = False
    self.start_spotfinding()
    self.tracker_panel.chart.draw_bragg_line()
    self.tracker_panel.chart.select_span.set_active(True)
    self.tracker_panel.chart.zoom_span.set_active(True)

  def onRestoreRun(self, e):
    self.remove_term_file()
    self.spotfinding_info = []
    with open(self.info_file, 'r') as f:
      contents = f.readlines()

    with open(self.folder_file, 'r') as ff:
      self.data_folder = ff.read()

    for item in contents:
      items = item.replace('\n', '').split(',')
      self.spotfinding_info.append(items)

    os.remove(self.info_file)
    os.remove(self.folder_file)

    self.new_counts = [int(i[1]) for i in self.spotfinding_info]
    self.new_frames = [int(i[0]) for i in self.spotfinding_info]
    self.done_list = [i[2] for i in self.spotfinding_info]
    self.all_info = self.spotfinding_info

    self.plot_results()
    self.start_spotfinding()

  def start_spotfinding(self):
    ''' Start timer and perform spotfinding on found images '''
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), True)
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_open.GetId(), False)
    self.params = self.phil.extract()
    self.processor = IOTADialsProcessor(params=self.params)
    self.tracker_panel.spf_options.Disable()


    with open(self.folder_file, 'w') as f:
      f.write(self.data_folder)

    self.timer.Start(1000)
    self.spin_update = 0
    self.find_new_images()

  def stop_run(self):
    timer_txt = '[ xxxxxx ]'
    self.msg = 'STOPPED SPOTFINDING!'
    self.tracker_panel.status_txt.SetLabel('{} {}'
                                           ''.format(timer_txt, self.msg))
    self.timer.Stop()

  def run_spotfinding(self):
    ''' Generate the spot-finder thread and run it '''
    self.spf_thread = thr.SpotFinderThread(self,
                                           data_list=self.data_list,
                                           term_file=self.term_file,
                                           processor=self.processor)
    self.spf_thread.start()

  def onSpfOneDone(self, e):
    ''' Occurs on every wx.PostEvent instance; updates lists of images with
    spotfinding results '''
    if not self.terminated:
      if e.GetValue() is not None:
        info = e.GetValue()
        idx = info[0] + len(self.done_list)
        obs_count = info[1]
        img_path = info[2]
        self.obs_counts[idx] = obs_count
        self.new_frames.append(idx)
        self.new_counts.append(obs_count)
        self.spotfinding_info.append([idx, obs_count, img_path])
        self.all_info.append([idx, obs_count, img_path])
      #self.plot_results()

  def onSpfAllDone(self, e):
    self.done_list.extend(e.GetValue())
    self.find_new_images()
    # if e.GetValue() == []:
    #   self.stop_run()
    # else:
    #   self.done_list.extend(e.GetValue())
    #   self.find_new_images()

  def onSpfTerminated(self, e):
    self.stop_run()
    self.toolbar.EnableTool(self.tb_btn_open.GetId(), True)

  def find_new_images(self):
    if self.done_list != []:
      last_file = self.done_list[-1]
    else:
      last_file = None
    found_files = ginp.make_input_list([self.data_folder],
                                       filter=True,
                                       filter_type='image',
                                       last=last_file)

    if found_files != []:
      print 'DEBUG: FIRST FILE - {}'.format(found_files[0])

    self.data_list = list(set(found_files) - set(self.done_list))
    self.data_list = [i for i in self.data_list if not 'tmp' in i]

    print 'DEBUG: FOUND {} FILES\n'.format(len(self.data_list))

    if len(self.data_list) == 0:
      self.msg = 'Waiting for new images in {} ...'.format(self.data_folder)
      self.waiting = True
    else:
      self.waiting = False
      self.msg = 'Tracking new images in {} ...'.format(self.data_folder)
      if self.obs_counts == [] and self.frame_count == []:
        self.obs_counts = [-1] * len(self.data_list)
        self.frame_count = range(len(self.data_list))
      else:
        self.obs_counts = self.obs_counts + ([-1] * len(self.data_list))
        self.frame_count = self.frame_count + range(len(self.data_list))
      self.run_spotfinding()

  def onTimer(self, e):
    ''' Every second, update spotfinding chart (for some odd reason, chart is
    not updated when the wx.PostEvent happens) '''
    if self.spin_update == 5:
      self.spin_update = 0
    else:
      self.spin_update += 1
    tick = ['-o----', '--o---', '---o--', '----o-', '---o--', '--o---']
    timer_txt = '[ {0} ]'.format(tick[self.spin_update])
    self.tracker_panel.status_txt.SetLabel('{} {}'.format(timer_txt, self.msg))

    if not self.terminated:
      self.plot_results()
      if len(self.data_list) == 0:
        self.find_new_images()
    else:
      self.stop_run()

  def update_image_list(self):
    listctrl = self.tracker_panel.image_list.image_list.ctr
    if self.tracker_panel.chart.bracket_set:
      try:
        first_img = self.tracker_panel.chart.patch_x
        last_img = self.tracker_panel.chart.patch_x_last

        new_data_dict = {}
        sel_img_list = self.all_info[first_img:last_img]
        for img in sel_img_list:
          idx = sel_img_list.index(img)
          new_data_dict[idx] = (img[0], img[2], img[1])
        self.data_dict = new_data_dict
        listctrl.InitializeDataMap(self.data_dict)
        self.tracker_panel.btn_view_all.Enable()
        self.tracker_panel.btn_wrt_file.Enable()

      except TypeError:
        pass

    else:
      self.data_dict = {}
      self.tracker_panel.btn_view_all.Disable()
      self.tracker_panel.btn_view_sel.Disable()
      listctrl.InitializeDataMap(self.data_dict)

  def plot_results(self):
    self.tracker_panel.chart.draw_plot(new_x=self.new_frames,
                                       new_y=self.new_counts)
    # # Save results in a text file
    # with open(self.info_file, 'a') as f:
    #   for item in self.spotfinding_info:
    #     f.write('{},{},{}\n'.format(item[0], item[1], item[2]))

    self.spotfinding_info = []
    self.new_frames = []
    self.new_counts = []

  def onQuit(self, e):
    self.timer.Stop()
    with open(self.term_file, 'w') as tf:
      tf.write('')
    self.Close()
