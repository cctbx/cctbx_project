from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.track_images
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 04/19/2017
Last Changed: 04/26/2017
Description : IOTA image-tracking GUI module
'''

import os
import wx
from wxtbx import bitmaps
import numpy as np

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

# from dxtbx.datablock import DataBlockFactory
# from dials.algorithms.background import RadialAverage
from iotbx import phil as ip
from libtbx import easy_run

from iota import iota_version
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
    wx.Panel.__init__(self, parent)
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
    cur_bragg = self.main_window.tracker_panel.options.min_bragg.ctr.GetValue()
    new_bragg = cur_bragg + step
    self.main_window.tracker_panel.options.min_bragg.ctr.SetValue(new_bragg)
    self.draw_plot(min_bragg=new_bragg)

  def draw_plot(self, min_bragg=0, new_x=None, new_y=None):
    self.track_axes.patch.set_visible(False)

    if new_x is None:
      new_x = []
    if new_y is None:
      new_y = []

    nref_x = np.append(self.acc_plot.get_xdata(),
                       np.array(new_x).astype(np.double))
    nref_y = np.append(self.acc_plot.get_ydata(),
                       np.array(new_y).astype(np.double))
    nref_xy = zip(nref_x, nref_y)
    acc = [i[0] for i in nref_xy if i[1] >= min_bragg]
    rej = [i[0] for i in nref_xy if i[1] < min_bragg]

    self.acc_plot.set_xdata(nref_x)
    self.rej_plot.set_xdata(nref_x)
    self.acc_plot.set_ydata(nref_y)
    self.rej_plot.set_ydata(nref_y)
    self.bragg_line.set_ydata(min_bragg)
    self.acc_plot.set_markevery(acc)
    self.rej_plot.set_markevery(rej)

    if min_bragg > 0:
     self.bragg_line.set_alpha(1)
    else:
      self.bragg_line.set_alpha(0)

    # self.track_axes.relim()
    # self.track_axes.autoscale_view()

    if nref_x != [] and nref_y != []:
      if self.main_window.tracker_panel.options.chart_window.toggle:
        window_size = self.main_window.tracker_panel.options.chart_window.window.GetValue()
        if window_size == '':
          x_min = -0.5
          x_max = np.max(nref_x) + 0.5
        elif np.max(nref_x) > int(window_size):
          x_min = np.max(nref_x) - int(window_size) - 0.5
          x_max = np.max(nref_x)
        else:
          x_min = -0.5
          x_max = int(window_size) + 0.5
      else:
        x_min = -0.5
        x_max = np.max(nref_x) + 0.5
      self.track_axes.set_xlim(x_min, x_max)
      self.track_axes.set_ylim(0, np.max(nref_y) + int(0.1 * np.max(nref_y)))

    self.Layout()

    self.track_axes.draw_artist(self.acc_plot)
    self.track_axes.draw_artist(self.rej_plot)

    # self.track_canvas.flush_events()


class SpotfinderSettings(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    # Minimum Bragg spots cutoff (maybe make a slider?)
    common_size = (150, -1)
    self.display = wx.Panel(self)
    self.display_box = wx.StaticBox(self.display, label='Display Options')
    self.display_sizer = wx.StaticBoxSizer(self.display_box, wx.VERTICAL)
    self.display.SetSizer(self.display_sizer)
    self.min_bragg = ct.SpinCtrl(self.display,
                                 label='Min. Bragg spots',
                                 label_size=common_size,
                                 ctrl_size=(100, -1),
                                 ctrl_value=10,
                                 ctrl_min = 0)
    self.display_sizer.Add(self.min_bragg, flag=wx.BOTTOM, border=10)

    self.chart_window = ct.OptionCtrl(self.display,
                                      items=[('window', '100')],
                                      checkbox=True,
                                      checkbox_label='Finite chart window',
                                      checkbox_state=True,
                                      label_size=common_size,
                                      ctrl_size=(100, -1))
    self.display_sizer.Add(self.chart_window, flag=wx.BOTTOM, border=10)

    #Spotfinding settings
    self.options = wx.Panel(self)
    self.options_box = wx.StaticBox(self.options, label='Spotfinding Options')
    self.options_sizer = wx.StaticBoxSizer(self.options_box, wx.VERTICAL)
    self.options.SetSizer(self.options_sizer)
    self.sigma_background = ct.OptionCtrl(self.options,
                                          items=[('s_bkg', 6)],
                                          label='Sigma background',
                                          label_size=common_size,
                                          ctrl_size=(100, -1))
    self.options_sizer.Add(self.sigma_background, flag=wx.BOTTOM, border=10)

    self.sigma_strong = ct.OptionCtrl(self.options,
                                      items=[('s_strong', 3)],
                                      label='Sigma strong',
                                      label_size=common_size,
                                      ctrl_size=(100, -1))
    self.options_sizer.Add(self.sigma_strong, flag=wx.BOTTOM, border=10)

    self.global_threshold = ct.OptionCtrl(self.options,
                                          items=[('threshold', 0)],
                                          label='Global threshold',
                                          label_size=common_size,
                                          ctrl_size=(100, -1))
    self.options_sizer.Add(self.global_threshold, flag=wx.BOTTOM, border=10)

    self.min_local = ct.OptionCtrl(self.options,
                                   items=[('min_local', 2)],
                                   label='Min. local',
                                   label_size=common_size,
                                   ctrl_size=(100, -1))
    self.options_sizer.Add(self.min_local, flag=wx.BOTTOM, border=10)

    self.gain = ct.OptionCtrl(self.options,
                              items=[('gain', 1.0)],
                              label='Detector gain',
                              label_size=common_size,
                              ctrl_size=(100, -1))
    self.options_sizer.Add(self.gain, flag=wx.BOTTOM, border=10)

    self.kernel_size = ct.OptionCtrl(self.options,
                                     items=[('kernel', '3 3')],
                                     label='Kernel size',
                                     label_size=common_size,
                                     ctrl_size=(100, -1))
    self.options_sizer.Add(self.kernel_size, flag=wx.BOTTOM, border=10)

    self.main_sizer.Add(self.display, flag=wx.EXPAND | wx.BOTTOM, border=10)
    self.main_sizer.Add(self.options, 1, wx.EXPAND)

    self.Fit()

class TrackerPanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent)
    self.parent = parent

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)


    # Status box
    self.status_panel = wx.Panel(self)
    self.status_sizer = wx.BoxSizer(wx.VERTICAL)
    self.status_box = wx.StaticBox(self.status_panel, label='Status')
    self.status_box_sizer = wx.StaticBoxSizer(self.status_box, wx.HORIZONTAL)
    self.status_txt = wx.StaticText(self.status_panel, label='')
    self.status_box_sizer.Add(self.status_txt, flag=wx.ALL | wx.ALIGN_CENTER,
                              border=10)
    self.status_sizer.Add(self.status_box_sizer,
                          flag=wx.EXPAND | wx.ALL, border=3)
    self.status_panel.SetSizer(self.status_sizer)
    self.main_sizer.Add(self.status_panel, flag=wx.EXPAND|wx.ALL, border=10)

    self.graph_panel = wx.Panel(self)
    self.graph_sizer = wx.FlexGridSizer(1, 2, 0, 0)
    self.graph_sizer.AddGrowableRow(0)
    self.graph_sizer.AddGrowableCol(1)

    # Setting controls
    self.options = SpotfinderSettings(self.graph_panel)
    self.graph_sizer.Add(self.options, 1, flag=wx.EXPAND)

    # Put in chart
    self.chart = TrackChart(self.graph_panel, main_window=parent)
    self.graph_sizer.Add(self.chart, 1, flag=wx.EXPAND|wx.LEFT, border=10)

    self.graph_panel.SetSizer(self.graph_sizer)
    self.main_sizer.Add(self.graph_panel, 1, flag=wx.EXPAND|wx.ALL, border=10)


class TrackerWindow(wx.Frame):
  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(1200, 500))
    self.parent = parent
    self.term_file = os.path.join(os.curdir, '.terminate_image_tracker')

    self.obj_counts = None
    self.frame_count = None
    self.done_list = []
    self.new_frames = []
    self.new_counts = []

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
                                                longHelp='Open Images')
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
    self.Bind(wx.EVT_TOOL, self.onStop, self.tb_btn_stop)

    # Spotfinder / timer bindings
    self.Bind(thr.EVT_SPFDONE, self.onSpfOneDone)
    self.Bind(thr.EVT_SPFALLDONE, self.onSpfAllDone)
    self.Bind(wx.EVT_TIMER, self.onTimer, id=self.timer.GetId())

    # Settings bindings
    self.Bind(wx.EVT_SPINCTRL, self.onMinBragg,
              self.tracker_panel.options.min_bragg.ctr)

  def onStop(self, e):
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), False)
    with open(self.term_file, 'w') as tf:
      tf.write('')
    self.msg = 'Stopping...'

  def onResume(self, e):
    os.remove(self.term_file)

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
      timer_txt = '[ ------ ]'
      self.msg = 'Ready to track images in {}'.format(self.data_folder)
      self.tracker_panel.status_txt.SetLabel('{} {}'.format(timer_txt, self.msg))
    else:
      open_dlg.Destroy()
      return

  def onMinBragg(self, e):
    min_bragg = self.tracker_panel.options.min_bragg.ctr.GetValue()
    self.tracker_panel.chart.draw_plot(min_bragg=min_bragg)

  def onRunSpotfinding(self, e):
    ''' Start timer and perform spotfinding on found images '''
    self.toolbar.EnableTool(self.tb_btn_stop.GetId(), True)
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)

    s_bkg = self.tracker_panel.options.sigma_background.s_bkg.GetValue()
    s_str = self.tracker_panel.options.sigma_strong.s_strong.GetValue()
    thresh = self.tracker_panel.options.global_threshold.threshold.GetValue()
    min_local = self.tracker_panel.options.min_local.min_local.GetValue()
    gain = self.tracker_panel.options.gain.gain.GetValue()
    kernel = self.tracker_panel.options.kernel_size.kernel.GetValue()

    phil_string = '\n'.join(['spotfinder {',
                             '  threshold {',
                             '    xds {',
                             '      gain = {}'.format(gain),
                             '      kernel_size = {}'.format(kernel),
                             '      sigma_background = {}'.format(s_bkg),
                             '      sigma_strong = {}'.format(s_str),
                             '      min_local = {}'.format(min_local),
                             '      global_threshold = {}'.format(thresh),
                             '    }',
                             '  }',
                             '}'
                             ])
    spf_phil = ip.parse(phil_string)
    self.phil = self.phil.fetch(source=spf_phil)
    self.params = self.phil.extract()
    self.processor = IOTADialsProcessor(params=self.params)

    self.tracker_panel.options.options.Disable()

    self.data_list = ginp.make_input_list([self.data_folder])
    self.obs_counts = [-1] * len(self.data_list)
    self.frame_count = range(len(self.data_list))

    self.timer.Start(1000)
    self.spin_update = 0

    if len(self.data_list) == 0:
      self.msg = 'Waiting for new images in {} ...'.format(self.data_folder)
    else:
      self.msg = 'Tracking images in {} ...'.format(self.data_folder)
      self.run_spotfinding()

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
      self.obs_counts[idx] = obs_count
      self.new_frames.append(idx)
      self.new_counts.append(obs_count)

  def onSpfAllDone(self, e):
    if e.GetValue() == []:
      timer_txt = '[ xxxxxx ]'
      self.msg = 'STOPPED SPOTFINDING!'
      self.tracker_panel.status_txt.SetLabel('{} {}'
                                             ''.format(timer_txt, self.msg))
      self.timer.Stop()
    else:
      self.done_list.extend(e.GetValue())
      self.find_new_images()

  def find_new_images(self):
    if self.done_list != []:
      last_file = self.done_list[-1]
    else:
      last_file = None
    found_files = ginp.make_input_list([self.data_folder], last=last_file)

    print 'DEBUG: ACQUIRED {} FILES!'.format(len(found_files))

    self.data_list = [i for i in found_files if i not in self.done_list]
    if len(self.data_list) == 0:
      self.msg = 'Waiting for new images in {} ...'.format(self.data_folder)
    else:
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
    not updated when the wx.PostEvent happens'''
    if self.spin_update == 5:
      self.spin_update = 0
    else:
      self.spin_update += 1
    tick = ['-o----', '--o---', '---o--', '----o-', '---o--', '--o---']
    timer_txt = '[ {0} ]'.format(tick[self.spin_update])
    self.tracker_panel.status_txt.SetLabel('{} {}'.format(timer_txt, self.msg))

    if not os.path.isfile(self.term_file)  :
      self.plot_results()
      if len(self.data_list) == 0:
        self.find_new_images()

  def plot_results(self):
    min_bragg = self.tracker_panel.options.min_bragg.ctr.GetValue()
    self.tracker_panel.chart.draw_plot(min_bragg=min_bragg,
                                       new_x=self.new_frames,
                                       new_y=self.new_counts)
    self.new_frames = []
    self.new_counts = []

  def onQuit(self, e):
    self.timer.Stop()
    self.Close()

class MainApp(wx.App):
  ''' App for the main GUI window  '''

  def OnInit(self):
    self.frame = TrackerWindow(None, -1, title='IOTA IMAGE TRACKER v.{}'
                                               ''.format(iota_version))
    self.frame.SetMinSize(self.frame.GetEffectiveMinSize())
    self.frame.SetPosition((150, 150))
    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

# ---------------------------------------------------------------------------- #

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()


# ---------------------------------------------------------------------------- #
# Below are un-used functions for radial average calculations. I probably
# won't use them now, but maybe later...

# def do_calculation(self, img):
#   # Create datablock
#   datablock = DataBlockFactory.from_filenames([img])[0]
#   self.make_radial_average(datablock=datablock)
#
# def make_radial_average(self, datablock):
#   imageset = datablock.extract_imagesets()[0]
#   beam = imageset.get_beam()
#   detector = imageset.get_detector()
#   scan_range = (0, len(imageset))
#
#   summed_data = None
#   summed_mask = None
#   for i in range(*scan_range):
#     data = imageset.get_raw_data(i)
#     mask = imageset.get_mask(i)
#     assert isinstance(data, tuple)
#     assert isinstance(mask, tuple)
#     if summed_data is None:
#       summed_mask = mask
#       summed_data = data
#     else:
#       summed_data = [ sd + d for sd, d in zip(summed_data, data) ]
#       summed_mask = [ sm & m for sm, m in zip(summed_mask, mask) ]
#   num_bins = sum(sum(p.get_image_size()) for p in detector)
#   vmin = (1.0 / 20) ** 2  #0
#   vmax = (1.0 / 1.5) ** 2 #detector.get_max_resolution(beam.get_s0())) ** 2
#
#   # Compute the radial average
#   radial_average = RadialAverage(beam, detector, vmin, vmax, num_bins)
#   for d, m in zip(summed_data, summed_mask):
#     radial_average.add(d.as_double() / (scan_range[1] - scan_range[0]), m)
#   mean = radial_average.mean()
#   reso = radial_average.inv_d2()
