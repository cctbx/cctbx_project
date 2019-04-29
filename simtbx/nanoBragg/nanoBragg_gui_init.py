from __future__ import division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 12/12/2017
Last Changed: 10/30/2018
Description : SIMTBX (nanoBragg) GUI Initialization module
'''

import os
import wx
from wxtbx import bitmaps
import numpy as np

from simtbx.nanoBragg import nanoBragg_gui_frames as frm
from simtbx.nanoBragg import nanoBragg_input as inp


# ------------------------------- Main Window -------------------------------- #

class MainWindow(wx.Frame):
  ''' Frame housing the entire app; all windows open from this one '''

  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title)
    self.parent = parent
    self.simtbx_phil = inp.generate_simtbx_phil()
    self.sparams = self.simtbx_phil.extract()

    #self.simtbx_phil.show()

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
    run_bmp = bitmaps.fetch_icon_bitmap('actions', 'run')
    self.tb_btn_run = self.toolbar.AddLabelTool(wx.ID_ANY, label='Run',
                                                bitmap=run_bmp,
                                                shortHelp='Run',
                                                longHelp='Run Spotfinding')
    self.toolbar.Realize()

    #Processing status bar
    self.sb = self.CreateStatusBar()
    self.sb.SetFieldsCount(2)
    self.sb.SetStatusWidths([-1, -2])

    # Output gauge in status bar
    self.gauge_process = wx.Gauge(self.sb, -1, style=wx.GA_HORIZONTAL | wx.GA_SMOOTH)
    rect = self.sb.GetFieldRect(0)
    self.gauge_process.SetPosition((rect.x + 2, rect.y + 2))
    self.gauge_process.SetSize((rect.width - 4, rect.height - 4))
    self.gauge_process.Hide()

    # Instantiate windows
    self.top_window = frm.TopPanel(self)

    # Single input window
    self.main_sizer.Add(self.top_window, 1,
                        flag=wx.ALL | wx.EXPAND,
                        border=10)
    self.main_sizer.Add((-1, 20))

    # Bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_run)
    self.Bind(wx.EVT_BUTTON, self.onRunPreview,
              self.top_window.opt_btn_prev)

  def onQuit(self, e):
    self.Close()

  def onRun(self, e):
    self.init_settings()
    self.top_window.sparams = self.sparams

    init = InitAll(self, self.sparams)
    good_to_go = init.run()

    if good_to_go:
      self.top_window.run_simulator(init=init)

  def onRunPreview(self, e):
    self.init_settings()
    self.top_window.sparams = self.sparams
    self.top_window.run_preview()


  def init_settings(self):

    # Pull settings from input window
    # TODO: once tabs are made, pull user-phil files from each tab
    self.top_window.generate_phil()
    self.simtbx_phil = self.simtbx_phil.fetch(self.top_window.input_phil)

    # Generate params python object
    self.sparams = self.simtbx_phil.extract()


# ------------------------------- Initialize --------------------------------- #

class InitAll(object):
  """ Class to initialize current nanoBragg run

  """

  def __init__(self, parent, params):
    self.parent = parent
    self.params = params
    self.img_base = None

  def make_input_list(self):
    if self.img_base is not None:
      start = self.params.dataset.start_phi
      finish = self.params.dataset.finish_phi
      step = self.params.dataset.oscillation
      angles = list(np.arange(start, finish, step))

      # generate list of image filenames
      img_prefix = self.params.image_prefix
      img_format = self.params.image_format
      inputs = []
      for phi in angles:
        idx = angles.index(phi) + 1
        input = {'index':idx,
                 'path':'{}_{:05d}.{}'.format(img_prefix, idx, img_format),
                 'phi':phi, 'osc':step}
        inputs.append(input)
      return inputs
    else:
      return None



  def run(self):
    ''' Initialize simtbx UI'''

    # Create folder for synthetic images
    self.img_base = os.path.join(os.path.abspath(self.params.output), 'images')

    # Generate list of inputs
    self.input_list = self.make_input_list()

    if self.input_list is not None:
      for i in self.input_list: print(i)
      return True
    else:
      return False

