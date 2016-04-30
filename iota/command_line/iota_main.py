from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 04/29/2016
Description : IOTA main module. Version 1.0.003
'''

import os
import wx
from threading import Thread

import difflib
import math
import numpy as np

import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
from matplotlib.figure import Figure

from libtbx.easy_mp import parallel_map
from libtbx import easy_pickle as ep

from iota.components.iota_analysis import Analyzer
import iota.components.iota_gui_init as ginit
import iota.components.iota_input as inp
import iota.components.iota_misc as misc
import iota.components.iota_image as img

iota_version = '1.0.003'
description = '''The integration optimization, triage and analysis (IOTA) toolkit for the processing of serial diffraction data.

Reference: Lyubimov, et al., J Appl Cryst, submitted
'''
license = ''' IOTA is distributed under open source license '''
icons = os.path.join(os.path.dirname(os.path.abspath(ginit.__file__)), 'icons/')

# Threading
# Set up events for finishing one cycle and for finishing all cycles
tp_EVT_ONEDONE = wx.NewEventType()
EVT_ONEDONE = wx.PyEventBinder(tp_EVT_ONEDONE, 1)
tp_EVT_ALLDONE = wx.NewEventType()
EVT_ALLDONE = wx.PyEventBinder(tp_EVT_ALLDONE, 1)

# Class for one cycle
class OneDone(wx.PyCommandEvent):
  def __init__(self, etype, eid, result=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.result = result
  def GetValue(self):
    return self.result

# Class for all cycles
class AllDone(wx.PyCommandEvent):
  def __init__(self, etype, eid, img_objects=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.image_objects = img_objects
  def GetValue(self):
    return self.image_objects

# Wrapper class (instantiate inside the importer class)
class ImportImage():
  def __init__(self, init, input_entry):
    self.init = init
    self.input_entry = input_entry
  def run(self):
    img_object = img.SingleImage(self.input_entry, self.init)
    imp_image = img_object.import_image()
    return imp_image

# Wrapper class (instantiate inside the importer class)
class ImportImageObject():
  def __init__(self, init, input_entry):
    self.init = init
    self.input_entry = input_entry
  def run(self):
    img_object = self.input_entry[2]
    imp_object = img_object.import_int_file(self.init)
    return imp_object

# Wrapper class (instantiate inside the importer class)
class IntegrateImage():
  def __init__(self, init, input_entry):
    self.init = init
    self.input_entry = input_entry
  def run(self):
    img_object = self.input_entry[2]
    proc_image = img_object.process()
    return proc_image

class ProcessImage():
  ''' Wrapper class to do one-fell-swoop processing of an image '''
  def __init__(self, init, input_entry, input_type = 'image'):
    self.init = init
    self.input_entry = input_entry
    self.input_type = input_type
  def run(self):
    if self.input_type == 'image':
      img_object = img.SingleImage(self.input_entry, self.init)
      img_object.import_image()
    elif self.input_type == 'object':
      img_object = self.input_entry[2]
      img_object.import_int_file(self.init)

    if self.init.params.image_conversion.convert_only:
      return img_object
    else:
      img_object.process()
      return img_object

# Work class (thread)
class ProcThread(Thread):
  ''' Worker thread; generated so that the GUI is not locked up when processing is running '''
  def __init__(self, parent, init, iterable, stage='import', input_type='image'):
    Thread.__init__(self)
    self.parent = parent
    self.init = init
    self.iterable = iterable
    self.stage = stage
    self.type = input_type

  def run(self):
    img_objects = parallel_map(iterable=self.iterable,
                               #func=self.proc_wrapper,
                               func = self.full_proc_wrapper,
                               callback = self.callback,
                               processes=self.init.params.n_processors)
    evt = AllDone(tp_EVT_ALLDONE, -1, img_objects)
    wx.PostEvent(self.parent, evt)

  def callback(self, result):
    evt = OneDone(tp_EVT_ONEDONE, -1, result)
    wx.PostEvent(self.parent, evt)

  def full_proc_wrapper(self, input_entry):
    proc_image_instance = ProcessImage(self.init, input_entry, self.type)
    proc_image = proc_image_instance.run()
    return proc_image

  def proc_wrapper(self, input_entry):
    if self.stage == 'import':
      imp_image_instance = ImportImage(self.init, input_entry)
      imp_image = imp_image_instance.run()
      return imp_image
    elif self.stage == 'gs_import':
      imp_object_instance = ImportImageObject(self.init, input_entry)
      imp_object = imp_object_instance.run()
      return imp_object
    elif self.stage == 'process':
      int_image_instance = IntegrateImage(self.init, input_entry)
      int_image = int_image_instance.run()
      return int_image


# ------------------------------------- Main Window ------------------------------------- #

class MainWindow(wx.Frame):
  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title)

    # Menu bar
    menubar = wx.MenuBar()

    # Status bar
    self.sb = self.CreateStatusBar()
    self.sb.SetFieldsCount(3)
    self.sb.SetStatusWidths([320, 200, -2])

    # Help menu item with the about dialog
    m_help = wx.Menu()
    m_file = wx.Menu()
    self.mb_load_script = m_file.Append(wx.ID_OPEN, '&Load Script...')
    self.mb_save_script = m_file.Append(wx.ID_SAVE, '&Save Script...')
    self.mb_about = m_help.Append(wx.ID_ANY, '&About')
    menubar.Append(m_file, '&File')
    menubar.Append(m_help, '&Help')

    self.SetMenuBar(menubar)

    main_box = wx.BoxSizer(wx.VERTICAL)

    # Toolbar
    self.toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.tb_btn_quit = self.toolbar.AddLabelTool(wx.ID_EXIT, label='Quit',
                                                 bitmap=wx.Bitmap('{}/32x32/exit.png'.format(icons)),
                                                 shortHelp='Quit',
                                                 longHelp='Quit IOTA')

    self.toolbar.AddSeparator()
    self.tb_btn_run = self.toolbar.AddLabelTool(wx.ID_ANY, label='Run',
                                                bitmap=wx.Bitmap('{}/32x32/run.png'.format(icons)),
                                                shortHelp='Run', longHelp='Run all stages of refinement')
    self.tb_btn_imglist = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                    label='Write Image List',
                                                    bitmap=wx.Bitmap('{}/32x32/list.png'.format(icons)),
                                                    shortHelp='Write Image List',
                                                    longHelp='Collect list of raw images and output to file')
    self.tb_btn_convert = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                    label='Convert Images',
                                                    bitmap=wx.Bitmap('{}/32x32/convert.png'.format(icons)),
                                                    shortHelp='Convert Image Files',
                                                    longHelp='Convert raw images to image pickle format')

    # These buttons will be disabled until input path is provided
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_imglist.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_convert.GetId(), False)
    self.toolbar.Realize()

    # Instantiate windows
    self.input_window = ginit.InputWindow(self)
    self.gparams = self.input_window.gparams

    # Single input window
    main_box.Add(self.input_window, flag=wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, border=10)
    main_box.Add((-1, 20))

    # Draw the main window sizer
    self.SetSizer(main_box)

    # Toolbar button bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_BUTTON, self.onInput, self.input_window.inp_btn_browse)
    self.Bind(wx.EVT_TEXT, self.onInput, self.input_window.inp_ctr)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_run)
    self.Bind(wx.EVT_TOOL, self.onWriteImageList, self.tb_btn_imglist)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_convert)

    # Menubar button bindings
    self.Bind(wx.EVT_MENU, self.OnAboutBox, self.mb_about)
    self.Bind(wx.EVT_MENU, self.onOutputScript, self.mb_save_script)
    self.Bind(wx.EVT_MENU, self.onLoadScript, self.mb_load_script)

  def init_settings(self):
    # Grab params from main window class
    self.gparams = self.input_window.gparams
    int_index = self.input_window.int_cbx_choice.GetCurrentSelection()
    self.gparams.advanced.integrate_with = str(self.input_window.int_cbx_choice.GetString(int_index)[:5]).lower()
    self.gparams.description = self.input_window.inp_des_ctr.GetValue()
    self.gparams.input = [self.input_window.inp_ctr.GetValue()]
    self.gparams.output = self.input_window.out_ctr.GetValue()
    self.gparams.advanced.random_sample.flag_on = self.input_window.opt_chk_random.GetValue()
    self.gparams.advanced.random_sample.number = self.input_window.opt_spc_random.GetValue()
    self.gparams.n_processors = self.input_window.opt_spc_nprocs.GetValue()

  def OnAboutBox(self, e):
    ''' About dialog '''
    info = wx.AboutDialogInfo()
    info.SetName('IOTA')
    info.SetVersion(iota_version)
    info.SetDescription(description)
    info.SetWebSite('http://cci.lbl.gov/xfel')
    info.SetLicense(license)
    info.AddDeveloper('Art Lyubimov')
    info.AddDeveloper('Monarin Uervirojnangkoorn')
    info.AddDeveloper('Aaron Brewster')
    info.AddDeveloper('Nick Sauter')
    info.AddDeveloper('Axel Brunger')
    info.AddDocWriter('Art Lyubimov')
    info.AddTranslator('Art Lyubimov')
    wx.AboutBox(info)

  def onInput(self, e):
    if self.input_window.inp_ctr.GetValue() != '':
      self.toolbar.EnableTool(self.tb_btn_imglist.GetId(), True)
      self.toolbar.EnableTool(self.tb_btn_convert.GetId(), True)
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), True)
    else:
      self.toolbar.EnableTool(self.tb_btn_imglist.GetId(), False)
      self.toolbar.EnableTool(self.tb_btn_convert.GetId(), False)
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)

  def onRun(self, e):
    # Run full processing
    if e.GetId() == self.tb_btn_convert.GetId():
      self.gparams.image_conversion.convert_only = True
      title = 'Image Conversion'
    else:
      title = 'Image Processing'
    self.init_settings()
    self.proc_window = ProcWindow(self, -1, title=title, params=self.gparams)
    if self.proc_window.good_to_go:
      self.proc_window.Show(True)

  def onOutputScript(self, e):
    self.init_settings()

    # Generate text of params
    final_phil = inp.master_phil.format(python_object=self.gparams)
    with misc.Capturing() as txt_output:
      final_phil.show()
    txt_out = ''
    for one_output in txt_output:
      txt_out += one_output + '\n'

    # Save param file
    save_dlg = wx.FileDialog(self,
                             message="Save IOTA Script",
                             defaultDir=os.curdir,
                             defaultFile="*.param",
                             wildcard="*",
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                             )
    if save_dlg.ShowModal() == wx.ID_OK:
      with open(save_dlg.GetPath(), 'w') as savefile:
        savefile.write(txt_out)

  def onLoadScript(self, e):
    load_dlg = wx.FileDialog(self,
                             message="Load script file",
                             defaultDir=os.curdir,
                             defaultFile="*.param",
                             wildcard="*",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if load_dlg.ShowModal() == wx.ID_OK:
      filepath = load_dlg.GetPaths()[0]

      # Extract params from file
      import iotbx.phil as ip
      user_phil = ip.parse(open(filepath).read())
      working_phil = inp.master_phil.fetch(sources=[user_phil])
      self.gparams = working_phil.extract()

      # Update all params within sub-windows
      self.input_window.gparams = self.gparams

      # Set input window params
      self.input_window.inp_ctr.SetValue(str(self.gparams.input[1]))
      self.input_window.out_ctr.SetValue(str(self.gparams.output))
      self.input_window.inp_des_ctr.SetValue(str(self.gparams.description))
      self.input_window.opt_chk_random.SetValue(self.gparams.advanced.random_sample.flag_on)
      if self.gparams.advanced.random_sample.flag_on:
        self.input_window.opt_txt_random.Enable()
        self.input_window.opt_spc_random.Enable()
      else:
        self.input_window.opt_txt_random.Disable()
        self.input_window.opt_spc_random.Disable()
      self.input_window.opt_spc_random.SetValue(int(self.gparams.advanced.random_sample.number))
      self.input_window.opt_spc_nprocs.SetValue(int(self.gparams.n_processors))
      if str(self.gparams.advanced.integrate_with).lower() == 'cctbx':
        self.input_window.int_cbx_choice.SetSelection(0)
      elif str(self.gparams.advanced.integrate_with).lower() == 'dials':
        self.input_window.int_cbx_choice.SetSelection(1)

  def onWriteImageList(self, e):
    img_list_dlg = wx.FileDialog(self,
                                 message="Save List of Images",
                                 defaultDir=os.curdir,
                                 defaultFile="*.lst",
                                 wildcard="*",
                                 style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                                 )

    if img_list_dlg.ShowModal() == wx.ID_OK:
      self.init_settings()
      self.init = ginit.InitAll(iota_version)
      good_init = self.init.run(self.gparams, list_file = img_list_dlg.GetPath())

      if good_init:
        wx.MessageBox('List of images saved in {}'.format(img_list_dlg.GetPath()),
                      'Info', wx.OK | wx.ICON_INFORMATION)

  def onQuit(self, e):
    self.Close()

# ============================= Processing Window ============================= #

class LogPage(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

    self.log_sizer = wx.BoxSizer(wx.VERTICAL)
    self.log_window = wx.TextCtrl(self,
                                  style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_DONTWRAP)
    self.log_window.SetFont(wx.Font(9, wx.TELETYPE, wx.NORMAL, wx.NORMAL, False))
    self.log_sizer.Add(self.log_window, proportion=1, flag= wx.EXPAND | wx.ALL, border=10)
    self.SetSizer(self.log_sizer)

class GraphPage(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent)
    self.proc_sizer = wx.BoxSizer(wx.VERTICAL)
    self.proc_figure = Figure()
    gsp = gridspec.GridSpec(4, 4)
    self.int_axes = self.proc_figure.add_subplot(gsp[2:, 2:])
    self.bxy_axes = self.proc_figure.add_subplot(gsp[2:, :2])
    self.nsref_axes = self.proc_figure.add_subplot(gsp[0, :])
    self.nsref_axes.set_ylabel('Strong Spots')
    self.res_axes = self.proc_figure.add_subplot(gsp[1, :])
    self.res_axes.set_xlabel('Frame')
    self.res_axes.set_ylabel('Resolution')
    self.proc_figure.set_tight_layout(True)
    self.canvas = FigureCanvas(self, -1, self.proc_figure)
    self.proc_sizer.Add(self.canvas, proportion=1, flag =wx.EXPAND | wx.ALL, border=10)
    self.SetSizer(self.proc_sizer)


class ProcWindow(wx.Frame):
  def __init__(self, parent, id, title, params, test=False):
    wx.Frame.__init__(self, parent, id, title, size=(800, 800),
                      style= wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER)

    self.logtext = ''
    self.objects_in_progress = []
    self.obj_counter = 0
    self.gparams = params

    # Toolbar
    self.proc_toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.tb_btn_abort = self.proc_toolbar.AddLabelTool(wx.ID_ANY, label='Abort',
                                                bitmap=wx.Bitmap('{}/32x32/stop.png'.format(icons)),
                                                shortHelp='Abort')

    # Tabbed output window(s)
    self.proc_panel = wx.Panel(self)
    self.proc_nb = wx.Notebook(self.proc_panel, style=0)
    self.graph_tab = GraphPage(self.proc_nb)
    self.log_tab = LogPage(self.proc_nb)
    self.proc_nb.AddPage(self.log_tab, 'Log')
    self.proc_nb.AddPage(self.graph_tab, 'Graphs')
    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.main_sizer.Add(self.proc_nb, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.proc_panel.SetSizer(self.main_sizer)

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

    # Output polling timer
    self.timer = wx.Timer(self)

    # Event bindings
    self.Bind(EVT_ONEDONE, self.onFinishedOneTask)
    self.Bind(EVT_ALLDONE, self.onFinishedProcess)
    self.sb.Bind(wx.EVT_SIZE, self.onStatusBarResize)
    self.Bind(wx.EVT_TIMER, self.onTimer, id=self.timer.GetId())

    # Button bindings
    self.Bind(wx.EVT_TOOL, self.onAbort, self.tb_btn_abort)


    if not test:
      self.run()


  def onStatusBarResize(self, e):
    rect = self.sb.GetFieldRect(0)
    self.gauge_process.SetPosition((rect.x + 2, rect.y + 2))
    self.gauge_process.SetSize((rect.width - 4, rect.height - 4))


  def onAbort(self, e):
    with open(self.tmp_abort_file, 'w') as af:
      af.write('')
    self.sb.SetStatusText('Aborting...', 1)
    self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)


  def run(self):
    # Initialize IOTA parameters and log
    self.init = ginit.InitAll(iota_version)
    good_init = self.init.run(self.gparams)

    # Start process
    if good_init:
      self.tmp_abort_file = os.path.join(self.init.int_base, '.abort.tmp')
      #self.stage = 'import'
      #self.import_images()
      self.process_images()
      self.good_to_go = True
      self.timer.Start(5000)
    else:
      self.good_to_go = False


  # def import_images(self):
  #   if self.init.params.cctbx.selection.select_only.flag_on:
  #     img_list = [[i, len(self.init.gs_img_objects) + 1, j] for i, j in enumerate(self.init.gs_img_objects, 1)]
  #     img_import = ProcThread(self, self.init, img_list, stage='gs_import')
  #   else:
  #     img_list = [[i, len(self.init.input_list) + 1, j] for i, j in enumerate(self.init.input_list, 1)]
  #     img_import = ProcThread(self, self.init, img_list, stage='import')
  #   self.status_summary = [0] * len(img_list)
  #   self.objects_in_progress = []
  #   img_import.start()
  #   self.gauge_process.SetRange(len(img_list))
  #
  # def integrate_images(self):
  #   # Remove rejected images from image object list
  #   acc_img_objects = [i.fail for i in self.img_objects if i.fail == None]
  #
  #   # Exit if none of the images have diffraction
  #   if str(self.init.params.image_triage.type).lower() != 'none':
  #     if len(acc_img_objects) == 0:
  #       misc.main_log(self.init.logfile, 'No images have diffraction!', True)
  #       return False
  #     else:
  #       misc.main_log(self.init.logfile, "{} out of {} images have diffraction (at " \
  #                                       "least {} Bragg peaks)" \
  #                                       "".format(len(acc_img_objects),
  #                                                 len(self.img_objects),
  #                                                 self.init.params.image_triage.min_Bragg_peaks))
  #   if self.init.params.image_conversion.convert_only:
  #     return False
  #   else:
  #     img_list = [[i, len(self.img_objects) + 1, j] for i, j in enumerate(self.img_objects, 1)]
  #     self.nref_list = [0] * len(self.img_objects)
  #     self.nref_xaxis = [i.img_index for i in self.img_objects]
  #     self.res_list = [0] * len(self.img_objects)
  #     self.objects_in_progress = []
  #     img_process = ProcThread(self, self.init, img_list, stage='process')
  #     img_process.start()
  #     self.gauge_process.SetRange(len(img_list))
  #     return True


  def process_images(self):
    ''' One-fell-swoop importing / triaging / integration of images '''

    if self.init.params.cctbx.selection.select_only.flag_on:
      self.img_list = [[i, len(self.init.gs_img_objects) + 1, j] for i, j in enumerate(self.init.gs_img_objects, 1)]
      type = 'object'
    else:
      self.img_list = [[i, len(self.init.input_list) + 1, j] for i, j in enumerate(self.init.input_list, 1)]
      type = 'image'

    self.gauge_process.SetRange(len(self.img_list))
    self.status_summary = [0] * len(self.img_list)
    self.nref_list = [0] * len(self.img_list)
    self.nref_xaxis = [i[0] for i in self.img_list]
    self.res_list = [0] * len(self.img_list)

    img_process = ProcThread(self, self.init, self.img_list, input_type=type)
    img_process.start()


  def analyze_results(self):
    # Analysis of integration results
    if os.path.isfile(self.tmp_abort_file):
      font = self.sb.GetFont()
      font.SetWeight(wx.BOLD)
      self.sb.SetFont(font)
      self.sb.SetStatusText('ABORTED BY USER', 1)
      return
    if len(self.final_objects) == 0:
      self.sb.SetStatusText('No images successfully integrated', 1)
    elif not self.gparams.image_conversion.convert_only:
      analysis = Analyzer(self.init, self.img_objects, iota_version)
      analysis.print_results()
      analysis.unit_cell_analysis()
      analysis.print_summary()
      analysis.make_prime_input()

    self.display_log()
    self.plot_integration()
    self.timer.Stop()


  def onFinishedOneTask(self, e):
    ''' Called by event at the end of each process
    @param e: event object
    '''
    obj = e.GetValue()
    self.objects_in_progress.append(obj)
    self.nref_list[obj.img_index - 1] = obj.final['strong']
    self.res_list[obj.img_index - 1] = obj.final['res']

  # def plot_import(self):
  #
  #   names_numbers = [
  #     ['no diffraction', len([i for i in self.objects_in_progress if i.fail == 'failed triage'])],
  #     ['have diffraction', len([i for i in self.objects_in_progress if i.status == 'imported' and i.fail == None])]
  #   ]
  #   names_numbers.append(['not imported', len(self.init.input_list) - sum([i[1] for i in names_numbers[:2]])])
  #   names = [i[0] for i in names_numbers if i[1] > 0]
  #   numbers = [i[1] for i in names_numbers if i[1] > 0]
  #   colors = ['skyblue', 'paleturquoise', 'royalblue']
  #   self.graph_tab.imp_axes.clear()
  #   self.graph_tab.imp_axes.pie(numbers, autopct='%.0f%%', colors=colors, radius=0.75)
  #   self.graph_tab.imp_axes.legend(names, loc='upper left', bbox_to_anchor=(-0.25, 1.25))
  #   self.graph_tab.imp_axes.axis('equal')
  #   self.graph_tab.canvas.draw()


  def display_log(self):
    ''' Display log output '''
    with open(self.init.logfile, 'r') as logfile:
      self.old_logtext = self.logtext
      self.logtext = logfile.readlines()

    ins_pt = self.log_tab.log_window.GetInsertionPoint()
    diff = difflib.Differ()
    for i in diff.compare(self.old_logtext, self.logtext):
      if i.startswith('+ '):
        self.log_tab.log_window.AppendText(i.replace('+', '', 1))
        self.log_tab.log_window.SetInsertionPoint(ins_pt)


  def plot_integration(self):
    try:
      # Summary pie chart
      names_numbers = [
        ['have diffraction', len([i for i in self.objects_in_progress if i.status == 'imported' and i.fail == None])],
        ['failed triage', len([i for i in self.objects_in_progress if i.fail == 'failed triage'])],
        ['failed indexing / integration', len([i for i in self.objects_in_progress if i.fail == 'failed grid search'])],
        ['failed prefilter', len([i for i in self.objects_in_progress if i.fail == 'failed prefilter'])],
        ['failed spotfinding', len([i for i in self.objects_in_progress if i.fail == 'failed spotfinding'])],
        ['failed indexing', len([i for i in self.objects_in_progress if i.fail == 'failed indexing'])],
        ['failed integration', len([i for i in self.objects_in_progress if i.fail == 'failed integration'])],
        ['integrated', len([i for i in self.objects_in_progress if i.status == 'final' and i.final['final'] != None])]
      ]
      names_numbers.append(['not processed', len(self.img_list) - sum([i[1] for i in names_numbers])])
      names = [i[0] for i in names_numbers if i[1] > 0]
      numbers = [i[1] for i in names_numbers if i[1] > 0]
      colors = ['crimson', 'darkorchid', 'plum', 'salmon', 'coral', 'sandybrown', 'gold']
      self.graph_tab.int_axes.clear()
      self.graph_tab.int_axes.pie(numbers, autopct='%.0f%%', colors=colors)
      self.graph_tab.int_axes.legend(names, loc='lower left', fontsize=9, fancybox=True)
      self.graph_tab.int_axes.axis('equal')

      # Strong reflections per frame
      self.graph_tab.nsref_axes.clear()
      self.graph_tab.nsref_axes.grid(True)
      self.graph_tab.nsref_axes.scatter(self.nref_xaxis, self.nref_list, marker='o', color='red', picker=True)
      self.graph_tab.nsref_axes.set_xlim(0, np.max(self.nref_xaxis) + 1)
      self.graph_tab.nsref_axes.set_ylim(ymin=0)
      self.graph_tab.nsref_axes.set_ylabel('Ref. (I / sigI > {})' \
                                           ''.format(self.gparams.cctbx.selection.min_sigma))

      # Resolution per frame
      self.graph_tab.res_axes.clear()
      self.graph_tab.res_axes.grid(True)
      self.graph_tab.res_axes.scatter(self.nref_xaxis, self.res_list, marker='d', color='green', picker=True)
      self.graph_tab.res_axes.set_xlim(0, np.max(self.nref_xaxis) + 1)
      self.graph_tab.res_axes.set_ylim(ymin=0)
      self.graph_tab.res_axes.set_xlabel('Frame')
      self.graph_tab.res_axes.set_ylabel('Resolution')

      # Beam XY (cumulative)
      info = []
      wavelengths = []
      distances = []
      cells = []

      # Import relevant info
      for root, dirs, files in os.walk(self.init.fin_base):
        for filename in files:
          found_file = os.path.join(root, filename)
          if found_file.endswith(('pickle')):
            beam = ep.load(found_file)
            info.append([i, beam['xbeam'], beam['ybeam']])
            wavelengths.append(beam['wavelength'])
            distances.append(beam['distance'])
            cells.append(beam['observations'][0].unit_cell().parameters())

      # Calculate beam center coordinates and distances
      beamX = [i[1] for i in info]
      beamY = [j[2] for j in info]
      beam_dist = [math.hypot(i[1] - np.median(beamX), i[2] - np.median(beamY)) for i in info]

      wavelength = np.median(wavelengths)
      det_distance = np.median(distances)
      a = np.median([i[0] for i in cells])
      b = np.median([i[1] for i in cells])
      c = np.median([i[2] for i in cells])

      # Calculate predicted L +/- 1 misindexing distance for each cell edge
      aD = det_distance * math.tan(2 * math.asin(wavelength / (2 * a)))
      bD = det_distance * math.tan(2 * math.asin(wavelength / (2 * b)))
      cD = det_distance * math.tan(2 * math.asin(wavelength / (2 * c)))

      # Calculate axis limits of beam center scatter plot
      beamxy_delta = np.ceil(np.max(beam_dist))
      xmax = round(np.median(beamX) + beamxy_delta)
      xmin = round(np.median(beamX) - beamxy_delta)
      ymax = round(np.median(beamY) + beamxy_delta)
      ymin = round(np.median(beamY) - beamxy_delta)

      # Plot beam center scatter plot
      self.graph_tab.bxy_axes.axis([xmin, xmax, ymin, ymax])
      self.graph_tab.bxy_axes.scatter(beamX, beamY, alpha=1, s=20, c='grey', lw=1)
      self.graph_tab.bxy_axes.plot(np.median(beamX), np.median(beamY), markersize=8, marker='o', c='yellow', lw=2)

      # Plot projected mis-indexing limits for all three axes
      circle_a = plt.Circle((np.median(beamX), np.median(beamY)), radius=aD, color='r', fill=False, clip_on=True)
      circle_b = plt.Circle((np.median(beamX), np.median(beamY)), radius=bD, color='g', fill=False, clip_on=True)
      circle_c = plt.Circle((np.median(beamX), np.median(beamY)), radius=cD, color='b', fill=False, clip_on=True)
      self.graph_tab.bxy_axes.add_patch(circle_a)
      self.graph_tab.bxy_axes.add_patch(circle_b)
      self.graph_tab.bxy_axes.add_patch(circle_c)
      self.graph_tab.bxy_axes.set_xlabel('BeamX (mm)', fontsize=15)
      self.graph_tab.bxy_axes.set_ylabel('BeamY (mm)', fontsize=15)
      self.graph_tab.bxy_axes.set_title('Beam Center Coordinates')

      self.graph_tab.canvas.draw()

    except ValueError, e:
      print self.nref_list
      print self.res_list


  def onTimer(self, e):
    # During image processing, plot stats on screen

    # if self.stage in ('import', 'gs_import') and len(self.objects_in_progress) > 0:
    #   self.graph_tab.imp_axes.set_visible(True)
    #   self.plot_import()
    #
    # elif self.stage == 'process' and len(self.objects_in_progress) > 0:
    #   self.graph_tab.int_axes.set_visible(True)
    #   self.graph_tab.nsref_axes.set_visible(True)
    #   self.graph_tab.res_axes.set_visible(True)
    #   self.plot_integration()

    if len(self.objects_in_progress) > self.obj_counter:
      self.plot_integration()
      self.obj_counter = len(self.objects_in_progress)


    # Update gauge
    self.gauge_process.Show()
    self.gauge_process.SetValue(len(self.objects_in_progress))

    # Update status bar
    if os.path.isfile(self.tmp_abort_file):
      self.sb.SetStatusText('Aborting......', 1)
    # elif self.stage == 'import':
    #   img_with_diffraction = [i for i in self.objects_in_progress if i.status == 'imported' and i.fail == None]
    #   self.sb.SetStatusText('{} of {} images imported, {} have diffraction'\
    #                         ''.format(len(self.objects_in_progress), len(self.init.input_list),
    #                                   len(img_with_diffraction)), 1)
    # elif self.stage == 'gs_import':
    #   self.sb.SetStatusText(
    #     '{} of {} objects imported'.format(len(self.objects_in_progress), len(self.init.gs_img_objects)), 1)
    #elif self.stage == 'process':
    elif self.gparams.image_conversion.convert_only:
      img_with_diffraction = [i for i in self.objects_in_progress if i.status == 'imported' and i.fail == None]
      self.sb.SetStatusText('{} of {} images imported, {} have diffraction'\
                            ''.format(len(self.objects_in_progress), len(self.init.input_list),
                                      len(img_with_diffraction)), 1)
    else:
      processed_images = [i for i in self.objects_in_progress if i.status == 'final']
      self.sb.SetStatusText('{} of {} images processed, {} successfully integrated' \
                            ''.format(len(self.objects_in_progress), len(self.img_list), len(processed_images)), 1)

    # Update log
    self.display_log()


  def onFinishedProcess(self, e):
    if os.path.isfile(self.tmp_abort_file):
      font = self.sb.GetFont()
      font.SetWeight(wx.BOLD)
      self.sb.SetFont(font)
      self.sb.SetStatusText('ABORTED BY USER', 1)
      return
    # if self.stage == 'import':
    #   self.img_objects = e.GetValue()
    #   self.plot_import()
    #   if self.integrate_images():
    #     self.sb.SetStatusText('Imported {} images'.format(len(self.img_objects)), 1)
    #     self.stage = 'process'
    #   else:
    #     self.gauge_process.Hide()
    #     self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)
    #     self.sb.SetStatusText('Converted {} images'.format(len(self.img_objects)), 1)
    #     return
    # elif self.stage == 'process':
    else:
      self.img_objects = e.GetValue()
      self.final_objects = [i for i in self.img_objects if i.fail == None]
      self.gauge_process.Hide()
      self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)
      self.sb.SetStatusText('{} of {} images successfully integrated'\
                            ''.format(len(self.final_objects), len(self.img_objects)), 1)
      self.plot_integration()
      self.analyze_results()


# ====================================================================================== #

class InputApp(wx.App):
  ''' App for the input window; will only collect user input and generate parameters for the IOTA run;
      The processing itself will be run in a separate ProcessingApp class'''
  def OnInit(self):
    self.frame = MainWindow(None, -1, title='IOTA v.{}'.format(iota_version))
    self.frame.Fit()
    self.frame.SetPosition((150, 150))
    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

if __name__ == '__main__':
  app = InputApp(0)
  app.MainLoop()
