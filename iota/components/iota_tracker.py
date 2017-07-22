from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 07/21/2017
Last Changed: 07/21/2017
Description : IOTA image-tracking GUI module
'''

import os
import wx
from wxtbx import bitmaps
import numpy as np

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

from iotbx import phil as ip

from iota.components.iota_dialogs import DIALSSpfDialog
from iota.components.iota_utils import InputFinder
from iota.components.iota_dials import phil_scope, IOTADialsProcessor
import iota.components.iota_threads as thr
import iota.components.iota_controls as ct

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
                            '    xds {',
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

    self.track_figure.patch.set_visible(False)
    self.track_axes.patch.set_visible(False)

    self.xdata = []
    self.ydata = []
    self.acc_plot = self.track_axes.plot([], [], 'o', color='#4575b4')[0]
    self.rej_plot = self.track_axes.plot([], [], 'o', color='#d73027')[0]
    self.bragg_line = self.track_axes.axhline(0, c='#4575b4', ls=':', alpha=0)

    self.track_axes.set_autoscaley_on(True)

    self.track_figure.set_tight_layout(True)
    self.track_canvas = FigureCanvas(self, -1, self.track_figure)

    self.main_fig_sizer.Add(self.track_canvas, 1, wx.EXPAND)
    self.track_figure.canvas.mpl_connect('scroll_event', self.onScroll)

  def onScroll(self, event):
    step = int(event.step)
    cur_bragg = self.main_window.tracker_panel.min_bragg.ctr.GetValue()
    new_bragg = cur_bragg + step
    self.main_window.tracker_panel.min_bragg.ctr.SetValue(new_bragg)
    self.draw_plot(min_bragg=new_bragg)

  def clear_all(self):
    self.track_axes.clear()
    self.acc_plot = self.track_axes.plot([], [], 'o', color='#4575b4')[0]
    self.rej_plot = self.track_axes.plot([], [], 'o', color='#d73027')[0]
    self.bragg_line = self.track_axes.axhline(0, c='#4575b4', ls=':', alpha=0)
    self.track_axes.set_autoscaley_on(True)
    self.track_canvas.flush_events()

  def refresh(self):
    self.clear_all()
    self.track_canvas.flush_events()

  def draw_plot(self, min_bragg=0, new_x=None, new_y=None):
    self.track_axes.patch.set_visible(False)

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

    if min_bragg > 0:
     self.bragg_line.set_alpha(1)
    else:
      self.bragg_line.set_alpha(0)

    if nref_x != [] and nref_y != []:
      if self.main_window.tracker_panel.chart_window.toggle.GetValue():
        window_size = self.main_window.tracker_panel.chart_window.ctr.GetValue()
        if window_size == '':
          x_min = -1
          x_max = np.max(nref_x) + 1
        elif np.max(nref_x) > int(window_size):
          x_min = np.max(nref_x) - int(window_size) - 1
          x_max = np.max(nref_x)
        else:
          x_min = -1
          x_max = int(window_size) + 1

      else:
        x_min = -1
        x_max = np.max(nref_x) + 1

      if min_bragg > np.max(nref_y):
        y_max = min_bragg + int(0.1 * min_bragg)
      else:
        y_max = np.max(nref_y) + int(0.1 * np.max(nref_y))

      self.track_axes.set_xlim(x_min, x_max)
      self.track_axes.set_ylim(0, y_max)

    else:
      x_min = -1
      x_max = 1

    acc = [i for i in all_acc if i > x_min and i < x_max]
    rej = [i for i in all_rej if i > x_min and i < x_max]


    self.acc_plot.set_xdata(nref_x)
    self.rej_plot.set_xdata(nref_x)
    self.acc_plot.set_ydata(nref_y)
    self.rej_plot.set_ydata(nref_y)
    self.bragg_line.set_ydata(min_bragg)
    self.acc_plot.set_markevery(acc)
    self.rej_plot.set_markevery(rej)

    self.Layout()

    count = '{}'.format(len([i for i in nref_xy if i[1] >= min_bragg]))
    self.main_window.tracker_panel.count_txt.SetLabel(count)
    self.main_window.tracker_panel.status_sizer.Layout()

    try:
      self.track_axes.draw_artist(self.acc_plot)
      self.track_axes.draw_artist(self.rej_plot)
    except ValueError, e:
      print 'DEBUG: VALUE ERROR!!'
      self.refresh()

class ImageList(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, size=(450, -1))

    self.main_box = wx.StaticBox(self, label='Hits')
    self.main_sizer = wx.StaticBoxSizer(self.main_box, wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    self.image_list = FileListCtrl(self)
    self.main_sizer.Add(self.image_list, 1, flag=wx.EXPAND)

class FileListItem(object):
  ''' Class that will contain all the elements of a file list entry '''
  def __init__(self, path, buttons, n_obs=0):
    self.id = None
    self.n_obs = n_obs
    self.path = path
    self.buttons = buttons

class FileListCtrl(ct.CustomListCtrl):
  def __init__(self, parent):
    ct.CustomListCtrl.__init__(self, parent=parent)

    self.parent = parent
    self.main_window = parent.GetParent()

    # Generate columns
    self.ctr.InsertColumn(0, "File")
    self.ctr.InsertColumn(1, "")
    self.ctr.setResizeColumn(0)

  def add_item(self, img, n_obs):
    item = FileListItem(path=img,
                        n_obs=n_obs,
                        buttons=ct.MiniButtonBoxInput(self.ctr))
    item.buttons.btn_mag.Bind(wx.EVT_BUTTON, self.onMagButton)
    item.buttons.btn_delete.Bind(wx.EVT_BUTTON, self.onDelButton)
    item.buttons.btn_info.Bind(wx.EVT_BUTTON, self.onInfoButton)

    view_bmp = bitmaps.fetch_custom_icon_bitmap('image_viewer16')
    item.buttons.btn_mag.SetBitmapLabel(view_bmp)

    # Insert list item
    idx = self.ctr.InsertStringItem(self.ctr.GetItemCount() + 1,
                                    os.path.basename(item.path))
    self.ctr.SetItemWindow(idx, 1, item.buttons, expand=True)

    # Record index in all relevant places
    item.id = idx
    item.buttons.index = idx

    # Resize columns to fit content
    self.ctr.SetColumnWidth(1, width=-1)
    self.ctr.SetColumnWidth(0, width=-3)

    # Attach data object to item
    self.ctr.SetItemData(item.id, item)

    self.main_window.Layout()

  def delete_item(self, index):
    self.ctr.DeleteItem(index)

    # Refresh widget and list item indices
    for i in range(self.ctr.GetItemCount()):
      item_data = self.ctr.GetItemData(i)
      item_data.id = i
      item_data.buttons.index = i
      self.ctr.SetItemData(i, item_data)

  def delete_all(self):
    for idx in range(self.ctr.GetItemCount()):
      self.delete_item(index=0)

  def onMagButton(self, e):
    idx = e.GetEventObject().GetParent().index
    item_obj = self.ctr.GetItemData(idx)
    path = item_obj.path
    viewer = thr.ImageViewerThread(self,
                                   backend='dials',
                                   file_string=path)
    viewer.start()

  def onDelButton(self):
    pass

  def onInfoButton(self):
    pass


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
    self.image_list_sizer.Add(self.image_list, 1, flag=wx.EXPAND)
    self.image_list_panel.SetSizer(self.image_list_sizer)

    self.main_sizer.Add(self.graph_panel, pos=(1, 0),
                        flag=wx.EXPAND | wx.LEFT | wx.BOTTOM, border=10)
    self.main_sizer.Add(self.image_list_panel, pos=(1, 1),
                        flag=wx.EXPAND | wx.LEFT | wx.BOTTOM | wx.RIGHT,
                        border=10)
    self.main_sizer.AddGrowableCol(0)
    self.main_sizer.AddGrowableRow(1)


class TrackerWindow(wx.Frame):
  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(1500, 600))
    self.parent = parent
    self.term_file = os.path.join(os.curdir, '.terminate_image_tracker')
    self.info_file = os.path.join(os.curdir, '.spotfinding_info')
    self.folder_file = os.path.join(os.curdir, '.data_info')

    self.frame_count = []
    self.obs_counts = []
    self.done_list = []
    self.new_frames = []
    self.new_counts = []
    self.spotfinding_info = []
    self.all_info = []
    self.hits = []
    self.refresh_chart = False
    self.current_min_bragg = 0
    self.waiting = False

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
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), False)
    # self.toolbar.EnableTool(self.tb_btn_calc.GetId(), False)
    if os.path.isfile(self.folder_file) and os.path.isfile(self.info_file):
      self.toolbar.EnableTool(self.tb_btn_restore.GetId(), True)
    else:
      self.toolbar.EnableTool(self.tb_btn_restore.GetId(), False)

    self.toolbar.Realize()

    # Setup timer
    self.timer = wx.Timer(self)

    self.tracker_panel = TrackerPanel(self)
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

    # Spotfinder / timer bindings
    self.Bind(thr.EVT_SPFDONE, self.onSpfOneDone)
    self.Bind(thr.EVT_SPFALLDONE, self.onSpfAllDone)
    self.Bind(wx.EVT_TIMER, self.onTimer, id=self.timer.GetId())

    # Settings bindings
    self.Bind(wx.EVT_BUTTON, self.onSpfOptions, self.tracker_panel.spf_options)
    self.Bind(wx.EVT_SPINCTRL, self.onMinBragg,
              self.tracker_panel.min_bragg.ctr)

  def onStop(self, e):
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
    min_bragg = self.tracker_panel.min_bragg.ctr.GetValue()
    self.tracker_panel.chart.draw_plot(min_bragg=min_bragg)

  def onSpfOptions(self, e):
    spf_dlg = DIALSSpfDialog(self,
                             phil=self.phil,
                             title='DIALS Spotfinding Options')
    if (spf_dlg.ShowModal() == wx.ID_OK):
      self.phil = self.phil.fetch(source=spf_dlg.spf_phil)
    spf_dlg.Destroy()


  def onRunSpotfinding(self, e):
    self.start_spotfinding()

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
    self.params = self.phil.extract()
    self.processor = IOTADialsProcessor(params=self.params)

    self.tracker_panel.spf_options.Disable()

    #self.data_list = ginp.make_input_list([self.data_folder])
    #self.obs_counts = [-1] * len(self.data_list)
    #self.frame_count = range(len(self.data_list))

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
    self.new_frames = []
    self.new_counts = []
    self.spotfinding_info = []
    self.results = []
    self.data_list = []
    self.done_list = []
    self.tracker_panel.chart.clear_all()
    self.tracker_panel.spf_options.Enable()

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

  def onSpfAllDone(self, e):
    if e.GetValue() == []:
      self.stop_run()
    else:
      self.done_list.extend(e.GetValue())
      self.find_new_images()

  def find_new_images(self):
    if self.done_list != []:
      last_file = self.done_list[-1]
    else:
      last_file = None
    found_files = ginp.make_input_list([self.data_folder], last=last_file)

    self.data_list = [i for i in found_files if i not in self.done_list]
    self.data_list = [i for i in self.data_list if not 'tmp' in i]
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

    if not os.path.isfile(self.term_file):
      self.update_image_list()
      self.plot_results()
      if len(self.data_list) == 0:
        self.find_new_images()
    else:
      if self.waiting:
        self.stop_run()

  def update_image_list(self):
    min_bragg = self.tracker_panel.min_bragg.ctr.GetValue()

    if self.current_min_bragg != min_bragg:
      self.current_min_bragg = min_bragg
      self.tracker_panel.image_list.image_list.delete_all()
      updated_hits = [i for i in self.all_info if i[1] >= min_bragg]
      if len(updated_hits) > 0:
       for hit in updated_hits:
         self.tracker_panel.image_list.image_list.add_item(img=hit[2],
                                                           n_obs=hit[1])

    new_hits = [i for i in self.spotfinding_info if i[1] >= min_bragg]
    if len(new_hits) > 0:
      for hit in new_hits:
        if hit[1] >= min_bragg:
          self.tracker_panel.image_list.image_list.add_item(img=hit[2],
                                                          n_obs=hit[1])
    self.hits.extend(new_hits)


  def plot_results(self):
    min_bragg = self.tracker_panel.min_bragg.ctr.GetValue()
    self.tracker_panel.chart.draw_plot(min_bragg=min_bragg,
                                       new_x=self.new_frames,
                                       new_y=self.new_counts)
    # Save results in a text file
    with open(self.info_file, 'a') as f:
      for item in self.spotfinding_info:
        f.write('{},{},{}\n'.format(item[0], item[1], item[2]))

    self.spotfinding_info = []
    self.new_frames = []
    self.new_counts = []

  def onQuit(self, e):
    self.timer.Stop()
    with open(self.term_file, 'w') as tf:
      tf.write('')
    self.Close()