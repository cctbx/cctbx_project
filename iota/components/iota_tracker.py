from __future__ import absolute_import, division, print_function
from six.moves import range, zip

'''
Author      : Lyubimov, A.Y.
Created     : 07/21/2017
Last Changed: 01/30/2019
Description : IOTA image-tracking GUI module
'''

import os
import wx
import argparse

from wxtbx import bitmaps
import wx.lib.agw.ultimatelistctrl as ulc
import wx.lib.mixins.listctrl as listmix
import numpy as np

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

from iotbx import phil as ip

from iota import iota_version
from iota.components.iota_ui_dialogs import DIALSSpfDialog
from iota.components.iota_utils import InputFinder
from iota.components.iota_processing import phil_scope
import iota.components.iota_threads as thr
import iota.components.iota_ui_controls as ct
import iota.components.iota_utils as util

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

help_message = ''' This is a tracking module for collected images '''

def parse_command_args(help_message):
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog = 'iota',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=(help_message),
            epilog=('\n{:-^70}\n'.format('')))
  parser.add_argument('path', type=str, nargs = '?', default = None,
            help = 'Path to data or file with IOTA parameters')
  parser.add_argument('--version', action = 'version',
            version = 'IOTA {}'.format(iota_version),
            help = 'Prints version info of IOTA')
  parser.add_argument('-s', '--start', action = 'store_true',
            help='Automatically start from first image')
  parser.add_argument('-p', '--proceed', action = 'store_true',
            help='Automatically start from latest image')
  parser.add_argument('-t', '--time', type=int, nargs=1, default=0,
            help='Automatically start from image collected n seconds ago')
  parser.add_argument('-b', '--backend', type=str, default='mosflm',
            help='Specify backend for spotfinding / indexing')
  parser.add_argument('-a', '--action', type=str, default='spotfind',
            help='Specify backend for spotfinding / indexing')
  parser.add_argument('-n', type=int, nargs='?', default=0, dest='nproc',
            help = 'Specify a number of cores for a multiprocessor run"')
  parser.add_argument('-f', '--file', type=str, default=None,
            help='Read spotfinding info from file instead of calculating')
  parser.add_argument('--reorder', action = 'store_true', default=True,
            help='Assign order number as images are read in (from file only)')
  parser.add_argument('--bragg', type=int, default=10,
            help = 'Specify the minimum number of Bragg spots for a "hit"')
  parser.add_argument('--paramfile', type=str, default=None,
                      help='Parameter file for processing')

  return parser

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

    self.reset_chart()

  def onSelect(self, xmin, xmax):
    """ Called when SpanSelector is used (i.e. click-drag-release); passes on
        the boundaries of the span to tracker window for selection and
        display of the selected images """
    if self.selector == 'select':
      self.select_span.set_visible(True)
      self.patch_x = int(xmin)
      self.patch_x_last = int(xmax) + 1
      self.patch_width = self.patch_x_last - self.patch_x
      self.bracket_set = True
      self.main_window.update_image_list()
      gp = self.main_window.tracker_panel.graph_panel
      ip = self.main_window.tracker_panel.image_list_panel
      sp = self.main_window.tracker_panel.chart_sash_position
      if sp == 0:
        sp = int(self.main_window.GetSize()[0] * 0.70)
      self.main_window.tracker_panel.chart_splitter.SplitVertically(gp, ip, sp)
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
          self.main_window.tracker_panel.chart_sash_position = \
            self.main_window.tracker_panel.chart_splitter.GetSashPosition()
          self.main_window.tracker_panel.chart_splitter.Unsplit()
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
    """ If left mouse button is pressed, activates the SpanSelector;
    otherwise, makes the span invisible and sets the toggle that clears the
    image list; if shift key is held, does this for the Selection Span,
    otherwise does this for the Zoom Span """
    if e.button != 1:
      self.zoom_span.set_visible(False)
      self.select_span.set_visible(False)
      self.bracket_set = False
      self.plot_zoom = False
      self.plot_sb.Hide()
      self.draw_plot()

      # Hide list of images
      self.main_window.tracker_panel.chart_sash_position = \
        self.main_window.tracker_panel.chart_splitter.GetSashPosition()
      self.main_window.tracker_panel.chart_splitter.Unsplit()
      self.main_window.tracker_panel.chart_window.toggle.SetValue(False)
      self.main_window.tracker_panel.chart_window.toggle_boxes(flag_on=False)
      self.main_window.tracker_panel.Layout()
    else:
      if self.main_window.tb_btn_view.IsToggled():
        self.selector = 'select'
        self.zoom_span.set_visible(False)
        self.select_span.set_visible(True)
      elif self.main_window.tb_btn_zoom.IsToggled():
        self.selector = 'zoom'
        self.zoom_span.set_visible(True)
        self.select_span.set_visible(False)

  def reset_chart(self):
    self.track_axes.clear()
    self.track_figure.patch.set_visible(False)
    self.track_axes.patch.set_visible(False)

    self.xdata = []
    self.ydata = []
    self.idata = []
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
    self.idx_plot = self.track_axes.plot([], [], 'wo', ms=2)[0]
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
    try:
      self.draw_plot()
    except AttributeError:
      pass

  def draw_plot(self, new_x=None, new_y=None, new_i=None, new_p=None):
    min_bragg = self.main_window.tracker_panel.min_bragg.ctr.GetValue()

    if new_x is None:
      new_x = []
    if new_y is None:
      new_y = []
    if new_i is None:
      new_i = []

    nref_x = np.append(self.xdata, np.array(new_x).astype(np.double))
    nref_y = np.append(self.ydata, np.array(new_y).astype(np.double))
    nref_i = np.append(self.idata, np.array(new_i).astype(np.double))
    self.xdata = nref_x
    self.ydata = nref_y
    self.idata = nref_i

    nref_xy = list(zip(nref_x, nref_y))
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

    acc = [int(i) for i in all_acc if self.x_min < i < self.x_max]
    rej = [int(i) for i in all_rej if self.x_min < i < self.x_max]

    self.acc_plot.set_xdata(nref_x)
    self.rej_plot.set_xdata(nref_x)
    self.idx_plot.set_xdata(nref_x)
    self.acc_plot.set_ydata(nref_y)
    self.rej_plot.set_ydata(nref_y)
    self.idx_plot.set_ydata(nref_i)
    self.acc_plot.set_markevery(acc)
    self.rej_plot.set_markevery(rej)

    self.Layout()

    count = '{}'.format(len([i for i in nref_xy if i[1] >= min_bragg]))
    idx_count = '{}'.format(len(nref_i[~np.isnan(nref_i)]))

    self.main_window.tracker_panel.count_txt.SetLabel(count)
    self.main_window.tracker_panel.idx_count_txt.SetLabel(idx_count)
    self.main_window.tracker_panel.info_sizer.Layout()

    # Set up scroll bar
    if len(self.xdata) > 0:
      self.plot_sb.SetMax(np.max(nref_x))
      if self.max_lock:
        self.plot_sb.SetValue(self.plot_sb.GetMax())

    # Draw extended plots
    self.track_axes.draw_artist(self.acc_plot)
    self.track_axes.draw_artist(self.rej_plot)

    # If any new folders are found, place marker at switch
    if new_p is not None:
      for p in new_p:
        self.track_axes.axvline(p[0], ymin=-15, c='red', ls='--')
        label = os.path.basename(p[1])
        self.track_axes.annotate(label, xy=(p[0], 1), ha='left', va='top',
                                 xycoords=('data', 'figure fraction'))


class ImageList(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent, size=(450, -1))

    self.main_box = wx.StaticBox(self, label='Hits')
    self.main_sizer = wx.StaticBoxSizer(self.main_box, wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    self.image_list = VirtualListCtrl(self)
    self.main_sizer.Add(self.image_list, 1, flag=wx.EXPAND)


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
    print ("column clicked")
    e.Skip()

  def instantiate_sorting(self, data):
    self.itemDataMap = data
    listmix.ColumnSorterMixin.__init__(self, 3)
    self.Bind(wx.EVT_LIST_COL_CLICK, self.OnColClick, self.ctr)

  def add_item(self, img_idx, img, n_obs, idxd=False, intd=False):
    item = ct.FileListItem(path=img,
                           items={
                             'n_obs':n_obs,
                             'img_idx':img_idx,
                             'idxd':idxd,
                             'intd':intd
                           })

    # Insert list item
    idx = self.ctr.InsertStringItem(self.ctr.GetItemCount() + 1,
                                    str(item.img_idx))
    self.ctr.SetStringItem(idx, 1, os.path.basename(item.path))
    self.ctr.SetStringItem(idx, 2, str(item.n_obs))
    if idxd:
      self.ctr.SetStringItem(idx, 3, 'X')
    if intd:
      self.ctr.SetStringItem(idx, 4, 'X')

    # Record index in item data
    item.id = idx

    # Attach data object to item
    self.ctr.SetItemData(idx, item)

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
    self.tracker_panel = parent.GetParent().GetParent().GetParent()

    # Generate columns
    self.ctr.InsertColumn(0, "#")
    self.ctr.InsertColumn(1, "File")
    self.ctr.InsertColumn(2, "No. Spots")
    self.ctr.InsertColumn(3, "IDX")
    self.ctr.setResizeColumn(1)
    self.ctr.setResizeColumn(2)

    # Bindings
    self.Bind(ulc.EVT_LIST_CACHE_HINT, self.OnSelection, self.ctr)
    self.Bind(wx.EVT_LIST_COL_CLICK, self.ctr.OnColClick, self.ctr)

  def initialize_data_map(self, data):
    self.ctr.InitializeDataMap(data)
    self.ctr.SetColumnWidth(0, width=30)
    self.ctr.SetColumnWidth(2, width=100)
    self.ctr.SetColumnWidth(3, width=15)
    self.ctr.SetColumnWidth(1, width=-3)

  def OnSelection(self, e):
    ctr = e.GetEventObject()
    if ctr.GetSelectedItemCount() <= 0:
      self.tracker_panel.btn_view_sel.Disable()
    else:
      self.tracker_panel.btn_view_sel.Enable()

class TrackerPanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent)
    self.parent = parent
    self.chart_sash_position = 0

    self.main_sizer = wx.GridBagSizer(10, 10)

    # Status box
    self.info_panel = wx.Panel(self)
    self.info_sizer = wx.FlexGridSizer(1, 4, 0, 10)
    self.info_sizer.AddGrowableCol(3)
    self.info_panel.SetSizer(self.info_sizer)

    self.count_box = wx.StaticBox(self.info_panel, label='Hits')
    self.count_box_sizer = wx.StaticBoxSizer(self.count_box, wx.HORIZONTAL)
    self.count_txt = wx.StaticText(self.info_panel, label='')
    self.count_box_sizer.Add(self.count_txt, flag=wx.ALL | wx.ALIGN_CENTER,
                             border=10)

    self.idx_count_box = wx.StaticBox(self.info_panel, label='Indexed')
    self.idx_count_box_sizer = wx.StaticBoxSizer(self.idx_count_box, wx.HORIZONTAL)
    self.idx_count_txt = wx.StaticText(self.info_panel, label='')
    self.idx_count_box_sizer.Add(self.idx_count_txt, flag=wx.ALL | wx.ALIGN_CENTER,
                          border=10)

    self.pg_box = wx.StaticBox(self.info_panel, label='Best Lattice')
    self.pg_box_sizer = wx.StaticBoxSizer(self.pg_box, wx.HORIZONTAL)
    self.pg_txt = wx.StaticText(self.info_panel, label='')
    self.pg_box_sizer.Add(self.pg_txt, flag=wx.ALL | wx.ALIGN_CENTER,
                          border=10)

    self.uc_box = wx.StaticBox(self.info_panel, label='Best Unit Cell')
    self.uc_box_sizer = wx.StaticBoxSizer(self.uc_box, wx.HORIZONTAL)
    self.uc_txt = wx.StaticText(self.info_panel, label='')
    self.uc_box_sizer.Add(self.uc_txt, flag=wx.ALL | wx.ALIGN_CENTER,
                          border=10)


    font = wx.Font(20, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    self.count_txt.SetFont(font)
    self.idx_count_txt.SetFont(font)
    font = wx.Font(18, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    self.pg_txt.SetFont(font)
    self.uc_txt.SetFont(font)

    self.info_sizer.Add(self.count_box_sizer, flag=wx.EXPAND)
    self.info_sizer.Add(self.idx_count_box_sizer, flag=wx.EXPAND)
    self.info_sizer.Add(self.pg_box_sizer, flag=wx.EXPAND)
    self.info_sizer.Add(self.uc_box_sizer, flag=wx.EXPAND)

    # Put chart and image list into splitter window
    self.chart_splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE |
                                                        wx.SP_3DSASH |
                                                        wx.SP_NOBORDER)

    # Put in chart
    self.graph_panel = wx.Panel(self.chart_splitter)
    self.graph_sizer = wx.GridBagSizer(5, 5)

    self.chart = TrackChart(self.graph_panel, main_window=parent)
    self.min_bragg = ct.SpinCtrl(self.graph_panel, label='Min. Bragg spots',
                                 ctrl_size=(100, -1), ctrl_value=10)
    self.chart_window = ct.SpinCtrl(self.graph_panel, checkbox=True,
                                    checkbox_label='Finite chart window',
                                    ctrl_size=(100, -1), ctrl_value=100,
                                    ctrl_min=10, ctrl_step=10)

    self.graph_sizer.Add(self.chart, flag=wx.EXPAND, pos=(0, 0), span=(1, 3))
    self.graph_sizer.Add(self.min_bragg, flag=wx.ALIGN_LEFT, pos=(1, 0))
    self.graph_sizer.Add(self.chart_window, flag=wx.ALIGN_CENTER, pos=(1, 1))

    self.graph_sizer.AddGrowableRow(0)
    self.graph_sizer.AddGrowableCol(1)
    self.graph_panel.SetSizer(self.graph_sizer)

    # List of images
    self.image_list_panel = wx.Panel(self.chart_splitter)
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

    # Add all to main sizer
    self.main_sizer.Add(self.info_panel, pos=(0, 0),
                        flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.chart_splitter, pos=(1, 0),
                        flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.AddGrowableCol(0)
    self.main_sizer.AddGrowableRow(1)
    self.SetSizer(self.main_sizer)
    self.chart_splitter.SplitVertically(self.graph_panel, self.image_list_panel)
    self.chart_splitter.Unsplit()

class TrackerWindow(wx.Frame):
  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(1500, 600))
    self.parent = parent
    self.term_file = os.path.join(os.curdir, '.terminate_image_tracker')
    self.spf_backend = 'mosflm'
    self.run_indexing = False
    self.run_integration = False
    self.running_clustering = False

    # Status bar
    self.sb = self.CreateStatusBar()
    self.sb.SetFieldsCount(2)
    self.sb.SetStatusWidths([100, -1])

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
    # pref_bmp = bitmaps.fetch_icon_bitmap('apps', 'advancedsettings')
    # self.tb_btn_prefs = self.toolbar.AddLabelTool(wx.ID_ANY,
    #                                               label='Preferences',
    #                                               bitmap=pref_bmp,
    #                                               shortHelp='Preferences',
    #                                               longHelp='IOTA image tracker preferences')
    # self.toolbar.AddSeparator()
    open_bmp = bitmaps.fetch_icon_bitmap('actions', 'open')
    self.tb_btn_open = self.toolbar.AddLabelTool(wx.ID_ANY, label='Open',
                                                bitmap=open_bmp,
                                                shortHelp='Open',
                                                longHelp='Open folder')
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
    self.toolbar.AddSeparator()
    span_view = bitmaps.fetch_custom_icon_bitmap('zoom_list')
    self.tb_btn_view = self.toolbar.AddLabelTool(wx.ID_ANY, label='View',
                                                 bitmap=span_view,
                                                 kind=wx.ITEM_RADIO,
                                                 shortHelp='Select to View',
                                                 longHelp='Select images to view')
    span_zoom = bitmaps.fetch_custom_icon_bitmap('zoom_view')
    self.tb_btn_zoom = self.toolbar.AddLabelTool(wx.ID_ANY, label='Zoom In',
                                                 bitmap=span_zoom,
                                                 kind=wx.ITEM_RADIO,
                                                 shortHelp='Zoom In',
                                                 longHelp='Zoom in on chart')
    self.toolbar.ToggleTool(self.tb_btn_zoom.GetId(), True)
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), False)
    self.toolbar.Realize()

    # Setup timers
    self.spf_timer = wx.Timer(self)
    self.uc_timer = wx.Timer(self)
    self.ff_timer = wx.Timer(self)

    self.tracker_panel = TrackerPanel(self)
    self.data_dict = self.tracker_panel.image_list.image_list.ctr.data.copy()
    self.img_list_initialized = False

    self.main_sizer.Add(self.tracker_panel, 1, wx.EXPAND)

    # Bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onGetImages, self.tb_btn_open)
    self.Bind(wx.EVT_TOOL, self.onRunSpotfinding, self.tb_btn_run)
    self.Bind(wx.EVT_TOOL, self.onStop, self.tb_btn_stop)
    self.Bind(wx.EVT_BUTTON, self.onSelView, self.tracker_panel.btn_view_sel)
    self.Bind(wx.EVT_BUTTON, self.onWrtFile, self.tracker_panel.btn_wrt_file)
    self.Bind(wx.EVT_BUTTON, self.onAllView, self.tracker_panel.btn_view_all)
    self.Bind(wx.EVT_TOOL, self.onZoom, self.tb_btn_zoom)
    self.Bind(wx.EVT_TOOL, self.onList, self.tb_btn_view)

    # Spotfinder / timer bindings
    self.Bind(thr.EVT_SPFDONE, self.onSpfOneDone)

    # Settings bindings
    self.Bind(wx.EVT_SPINCTRL, self.onMinBragg,
              self.tracker_panel.min_bragg.ctr)
    self.Bind(wx.EVT_SPINCTRL, self.onChartRange,
              self.tracker_panel.chart_window.ctr)
    self.Bind(wx.EVT_CHECKBOX, self.onChartRange,
              self.tracker_panel.chart_window.toggle)

    # Initialize tracker
    self.initialize_spotfinder()

  def onZoom(self, e):
    if self.tb_btn_zoom.IsToggled():
      self.toolbar.ToggleTool(self.tb_btn_view.GetId(), False)

  def onList(self, e):
    if self.tb_btn_view.IsToggled():
      self.toolbar.ToggleTool(self.tb_btn_zoom.GetId(), False)

  def initialize_spotfinder(self):
    self.data_folder = None
    self.done_list = []
    self.data_list = []
    self.spotfinding_info = []
    self.plot_idx = 0
    self.bookmark = 0
    self.all_info = []
    self.current_min_bragg = 0
    self.waiting = False
    self.submit_new_images = False
    self.terminated = False

    # Read arguments if any
    self.args, self.phil_args = parse_command_args('').parse_known_args()

    # Generate DIALS PHIL file
    if self.args.paramfile is None:
      default_phil = ip.parse(default_target)
      self.phil = phil_scope.fetch(source=default_phil)
    else:
      with open(self.args.paramfile, 'r') as phil_file:
        phil_string = phil_file.read()
      user_phil = ip.parse(phil_string)
      self.phil = phil_scope.fetch(source=user_phil)
    self.params = self.phil.extract()

    # Set backend
    self.spf_backend = self.args.backend

    # Determine how far the DIALS processing will go
    if 'index' in self.args.action:
      self.run_indexing = True
    elif 'int' in self.args.action:
      self.run_indexing = True
      self.run_integration = True
    self.tracker_panel.min_bragg.ctr.SetValue(self.args.bragg)

    # Determine how the tracker will track images: from file output by
    # iota.single_image, or by turning over actual files. If the latter,
    # determine at what point the tracker will start the tracking
    auto_start = True
    min_back = None
    if self.args.file is not None:
      self.results_file = self.args.file
    elif self.args.path is not None:
      path = os.path.abspath(self.args.path)
      self.open_images_and_get_ready(path=path)
      if self.args.start:
        print ('IMAGE_TRACKER: STARTING FROM FIRST RECORDED IMAGE')
      elif self.args.proceed:
        print ('IMAGE_TRACKER: STARTING FROM IMAGE RECORDED 1 MIN AGO')
        min_back = -1
      elif self.args.time > 0:
        min_back = -self.args.time[0]
        print ('IMAGE_TRACKER: STARTING FROM IMAGE RECORDED {} MIN AGO' \
              ''.format(min_back))
      else:
        auto_start = False
    else:
      auto_start = False

    # Initialize processing thread
    if self.args.file is None:
      self.proc_thread = thr.InterceptorThread(self,
                                               data_folder=self.data_folder,
                                               term_file=self.term_file,
                                               proc_params=self.params,
                                               backend=self.args.backend,
                                               n_proc=self.args.nproc,
                                               min_back=min_back,
                                               run_indexing=self.run_indexing,
                                               run_integration=self.run_integration)
    else:
      self.proc_thread = thr.InterceptorFileThread(self,
                                                   results_file=self.args.file,
                                                   reorder=self.args.reorder)

    if auto_start:
      self.start_spotfinding()


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
    self.stop_run()

  def remove_term_file(self):
    try:
      os.remove(self.term_file)
    except Exception:
      pass

  def onGetImages(self, e):
    """ Select folder to watch for incoming images """
    open_dlg = wx.DirDialog(self, "Choose the data folder:",
                            style=wx.DD_DEFAULT_STYLE)
    if open_dlg.ShowModal() == wx.ID_OK:
      self.data_folder = open_dlg.GetPath()
      open_dlg.Destroy()
      self.open_images_and_get_ready()
    else:
      open_dlg.Destroy()
      return

  def open_images_and_get_ready(self, path=None):
    if path is not None:
      self.data_folder = path

    self.remove_term_file()
    # self.initialize_spotfinder()

    self.tracker_panel.chart.reset_chart()
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), True)
    timer_txt = '[ ------ ]'
    self.msg = 'Ready to track images in {}'.format(self.data_folder)
    self.sb.SetStatusText('{} {}'.format(timer_txt, self.msg), 1)

  def onMinBragg(self, e):
    self.tracker_panel.chart.draw_bragg_line()

  def onChartRange(self, e):
    if self.tracker_panel.chart_window.toggle.GetValue():
      chart_range = self.tracker_panel.chart_window.ctr.GetValue()
      self.tracker_panel.chart.plot_zoom = True
      self.tracker_panel.chart.chart_range = chart_range
      self.tracker_panel.chart.max_lock = True
    else:
      self.tracker_panel.chart.plot_zoom = False
    self.tracker_panel.chart.draw_plot()

  def onSpfOptions(self, e):
    spf_dlg = DIALSSpfDialog(self,
                             phil=self.phil,
                             title='DIALS Spotfinding Options')
    if (spf_dlg.ShowModal() == wx.ID_OK):
      self.phil = self.phil.fetch(source=spf_dlg.spf_phil)
    spf_dlg.Destroy()


  def onRunSpotfinding(self, e):
    self.start_spotfinding()

  def start_spotfinding(self):
    """ Start timer and perform spotfinding on found images """
    self.terminated = False
    self.tracker_panel.chart.draw_bragg_line()
    self.tracker_panel.chart.select_span.set_active(True)
    self.tracker_panel.chart.zoom_span.set_active(True)
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), True)
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_open.GetId(), False)
    self.params = self.phil.extract()

    self.spin_update = 0
    if self.args.file is not None:
      self.sb.SetStatusText('{}'.format('FILE'), 0)
    else:
      self.sb.SetStatusText('{}'.format(self.spf_backend.upper()), 0)

    self.proc_thread.start()
    self.proc_thread.prc_timer.Start(1000)
    self.proc_thread.cls_timer.Start(15000)

  def stop_run(self):
    timer_txt = '[ xxxxxx ]'
    self.msg = 'STOPPED SPOTFINDING!'
    self.sb.SetStatusText('{} {}'.format(timer_txt, self.msg), 1)

    if self.proc_thread is not None:
      if self.proc_thread.prc_timer.IsRunning():
        self.proc_thread.prc_timer.Stop()
      if self.proc_thread.cls_timer.IsRunning():
        self.proc_thread.cls_timer.Stop()
      self.proc_thread.terminate_thread()


  def onSpfOneDone(self, e):
    """ Occurs on every wx.PostEvent instance; updates lists of images with
    spotfinding results """
    if not self.terminated:
      self.msg, self.spotfinding_info, self.cluster_info = e.GetValue()
      self.msg = self.proc_thread.msg
      self.update_spinner()
      self.plot_results()

  def update_spinner(self):
    """ Update spotfinding chart """
    if self.spin_update == 4:
      self.spin_update = 0
    else:
      self.spin_update += 1
    tick = ['-o-----', '--o----', '---o---', '----o--', '-----o-']
    timer_txt = '[ {0} ]'.format(tick[self.spin_update])
    self.sb.SetStatusText('{} {}'.format(timer_txt, self.msg), 1)

  def update_image_list(self):
    listctrl = self.tracker_panel.image_list.image_list.ctr
    if self.tracker_panel.chart.bracket_set:
      try:
        first_img = self.tracker_panel.chart.patch_x
        last_img = self.tracker_panel.chart.patch_x_last

        new_data_dict = {}
        sel_img_list = self.spotfinding_info[first_img:last_img]
        for img in sel_img_list:
          idx = sel_img_list.index(img)
          if img[3] is not None:
            idxd = 'X'
          else:
            idxd = ''
          new_data_dict[idx] = (img[0], img[2], img[1], idxd)
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
    """ Plot processing results; if indexing was turned on, will also provide
    unit cell clustering results """

    # Plot spotfinding / indexing results
    if self.spotfinding_info is not None and len(self.spotfinding_info) > 0:
      if self.plot_idx <= len(self.spotfinding_info):
        new_frames = [i[0] for i in self.spotfinding_info[self.plot_idx:]]
        new_counts = [i[1] for i in self.spotfinding_info[self.plot_idx:]]
        idx_counts = [i[1] if i[3] is not None else np.nan
                           for i in self.spotfinding_info[self.plot_idx:]]

        # Try finding where data comes from a new folder
        new_paths = []
        for info in self.spotfinding_info[self.plot_idx:]:
          cur_dir = os.path.dirname(info[2])
          cur_idx = self.spotfinding_info.index(info)
          if cur_idx > 0:
            prev_idx = cur_idx - 1
            prev_dir = os.path.dirname(self.spotfinding_info[prev_idx][2])
            if cur_dir != prev_dir:
              new_paths.append((info[0], cur_dir))
          else:
            new_paths.append((info[0], cur_dir))

        self.tracker_panel.chart.draw_plot(new_x=new_frames,
                                           new_y=new_counts,
                                           new_i=idx_counts,
                                           new_p=new_paths)
        self.plot_idx = self.spotfinding_info.index(self.spotfinding_info[-1]) + 1

    # Plot clustering results
    if self.cluster_info:
      clusters = sorted(self.cluster_info, key=lambda i: i['number'],
                        reverse=True)
      uc_dims = clusters[0]['uc'].rsplit()
      u = util.UnicodeCharacters()
      uc_line = "{} = {}, {} = {}, {} = {}, {} = {}, {} = {}, {} =  {} " \
                "".format('a', uc_dims[0], 'b', uc_dims[1], 'c', uc_dims[2],
                          u.alpha, uc_dims[3], u.beta, uc_dims[4], u.gamma,
                          uc_dims[5])

      self.tracker_panel.pg_txt.SetLabel(clusters[0]['pg'])
      self.tracker_panel.uc_txt.SetLabel(uc_line)

  def onQuit(self, e):
    self.Close()
    self.stop_run()
