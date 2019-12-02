from __future__ import division, print_function, absolute_import
from past.builtins import range
from six.moves import range

'''
Author      : Lyubimov, A.Y.
Created     : 01/17/2017
Last Changed: 12/02/2019
Description : IOTA GUI Windows / frames
'''

import os
import numpy as np
import time
import warnings
import multiprocessing

import wx
import wx.lib.buttons as btn
from wx import richtext as rt
from wx.aui import AuiNotebook

import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib import colors
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D

from wxtbx import bitmaps
from libtbx import easy_run, easy_pickle as ep
from libtbx.utils import to_unicode, Sorry

from prime.postrefine.mod_mx import mx_handler
import prime.postrefine.mod_plotter as ppl

from iota import iota_version, gui_description, gui_license
from iota.components import iota_init as init
from iota.components.iota_base import ProcInfo
from iota.components.gui.base import IOTABaseFrame, IOTABasePanel, \
  IOTABaseScrolledPanel
from iota.components.iota_analysis import Analyzer
from iota.components.gui.plotter import Plotter, PlotWindow
from iota.components.gui import make_phil_index
import iota.components.iota_input as inp
import iota.components.gui.controls as ct
import iota.components.iota_threads as thr
import iota.components.gui.dialogs as d
import iota.components.iota_utils as ut

import iota.components.gui.phil_controls as pct

assert Axes3D
f = ut.WxFlags()

# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  plot_font_size = 9
  norm_font_size = 9
  button_font_size = 10
  big_button_font_size = 12
  LABEL_SIZE = 12
  CAPTION_SIZE = 10
  python = 'python'
elif wx.Platform == '__WXMAC__':
  plot_font_size = 9
  norm_font_size = 12
  button_font_size = 14
  big_button_font_size = 16
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = "Python"
elif wx.Platform == '__WXMSW__':
  plot_font_size = 9
  norm_font_size = 9
  button_font_size = 11
  big_button_font_size = 15
  LABEL_SIZE = 11
  CAPTION_SIZE = 9
  python = "Python"  # TODO: make sure it's right!


# ------------------------------ Input Window -------------------------------- #

class InputWindow(pct.PHILDialogPanel):
  """ This window will enable all input handling """
  def __init__(self, parent, scope, *args, **kwargs):
    super(InputWindow, self).__init__(parent, scope, *args, **kwargs)


class MainWindow(IOTABaseFrame):
  """ Frame housing the entire app; all windows open from this one """

  def __init__(self, parent, id, title, input_dict=None, phil=None, msg=None):
    IOTABaseFrame.__init__(self, parent, id, title, size=(800, 500))
    self.parent = parent
    self.iota_phil = phil
    self.target_phil = None

    # Toolbar
    self.initialize_toolbar()
    self.tb_btn_quit = self.add_tool(label='Quit',
                                     bitmap=('actions', 'exit'),
                                     shortHelp='Quit')
    self.tb_btn_prefs = self.add_tool(label='Preferences',
                                      bitmap=('apps', 'advancedsettings'),
                                      shortHelp='Preferences')
    self.add_toolbar_separator(stretch=False)
    self.tb_btn_load = self.add_tool(label='Load Script',
                                     bitmap=('actions', 'download'),
                                     shortHelp='Load Script')
    self.tb_btn_save = self.add_tool(label='Save Script',
                                     bitmap=('actions', 'save_all'),
                                     shortHelp='Save Script')
    self.tb_btn_reset = self.add_tool(label='Reset',
                                      bitmap=('actions', 'reload'),
                                      shortHelp='Reset Settings')
    self.add_toolbar_separator(stretch=False)
    self.tb_btn_analysis = self.add_tool(label='Recover',
                                         bitmap=('actions', 'list'),
                                         shortHelp='Recover')
    self.add_toolbar_separator(stretch=True)
    self.tb_dlg_test = self.add_tool(label='Test',
                                     bitmap=('actions', 'utilities'),
                                     shortHelp='Test Dialog')
    self.toolbar.RemoveTool(self.tb_dlg_test.GetId())
    self.Bind(wx.EVT_TOOL, self.onTest, self.tb_dlg_test)

    # These buttons will be disabled until input path is provided
    # self.set_tool_state(self.tb_btn_run, False)
    self.realize_toolbar()

    # Toolbar button bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onPreferences, self.tb_btn_prefs)
    self.Bind(wx.EVT_TOOL, self.onRecovery, self.tb_btn_analysis)
    self.Bind(wx.EVT_TOOL, self.onLoadScript, self.tb_btn_load)
    self.Bind(wx.EVT_TOOL, self.onOutputScript, self.tb_btn_save)
    self.Bind(wx.EVT_TOOL, self.onReset, self.tb_btn_reset)

    # Add panels (make separate FlexGridSizer to make sure the sizing works)
    self.panel_sizer = wx.FlexGridSizer(2, 1, 0, 0)
    self.panel_sizer.AddGrowableRow(0)
    self.panel_sizer.AddGrowableCol(0)

    # Instantiate PHIL indices
    self.initialize_IOTA_index()
    self.update_IOTA_index(phil=self.iota_phil,
                           input_dict=input_dict,
                           messages=msg)
    self.initialize_backend_index()
    self.update_backend_index(phil=self.target_phil)

    # Instantiate input window
    main_window_scopes = ['description', 'output', 'input']
    input_scope = self.iota_index.get_scopes(include=main_window_scopes)
    self.input_window = InputWindow(self,
                                    size=(600, -1),
                                    scope=input_scope,
                                    phil_index=self.iota_index,
                                    path_extras={
                                      "file_types": ['image pickle file',
                                                     'raw image file',
                                                     'hdf5 image file',
                                                     'image file',
                                                     'hdf5 image list',
                                                     'raw image list',
                                                     'image pickle list',
                                                     'image list',
                                                     'mixed input list'],
                                      "folder_types": ['raw image folder',
                                                       'image pickle folder',
                                                       'hdf5 image folder',
                                                       'mixed input folder',
                                                       'image folder'],
                                      'data_types': ['image'],
                                      'input_filter':'image'
                                    })

    # Front options panel
    self.bottom_sizer = wx.GridBagSizer(5, 5)
    line = wx.StaticLine(self, style=wx.LI_HORIZONTAL, size=(-1, -1))
    self.bottom_sizer.Add(line, pos=(0, 0), span=(1, 3),
                          flag=wx.EXPAND| wx.TOP | wx.BOTTOM, border=5)
    opt_scope = self.iota_index.get_scopes(include=['data_selection'])
    self.option_panel = pct.PHILDialogPanel(self, scope=opt_scope,
                                            phil_index=self.iota_index)
    self.bottom_sizer.Add(self.option_panel, pos=(1, 0), span=(2, 1),
                          flag=wx.EXPAND)
    line = wx.StaticLine(self, style=wx.LI_VERTICAL, size=(-1, -1))
    self.bottom_sizer.Add(line, pos=(1, 1), span=(2, 1),
                          flag=wx.EXPAND| wx.RIGHT | wx.LEFT, border=15)
    self.bottom_sizer.AddGrowableCol(0)
    self.bottom_sizer.AddGrowableRow(0)

    # Front options
    btn_box = wx.StaticBox(self, label='')
    btn_sizer = wx.StaticBoxSizer(btn_box, wx.VERTICAL)
    self.opt_btn_iota = wx.Button(self, label='IOTA Settings...')
    self.opt_btn_proc = wx.Button(self, label='Processor Options...')
    btn_sizer.Add(self.opt_btn_iota, flag=wx.EXPAND| wx.ALL, border=5)
    btn_sizer.Add(self.opt_btn_proc, flag=wx.EXPAND| wx.ALL, border=5)
    self.bottom_sizer.Add(btn_sizer, pos=(1, 2), flag=wx.EXPAND)

    run_bmp = bitmaps.fetch_icon_bitmap('actions', 'run')
    self.btn_run = ct.GradButton(self,
                                 size=(140, 40),
                                 start_color=(161,217,155),
                                 highlight_color=(199,233,192),
                                 gradient_percent=-25,
                                 bmp=run_bmp,
                                 label='   RUN IOTA',
                                 label_size=big_button_font_size)
    self.bottom_sizer.Add(self.btn_run, pos=(2, 2),
                          flag=wx.SHAPED | wx.ALIGN_BOTTOM)

    # do not disable if there's already pre-loaded
    if not self.gparams.input:
      self.btn_run.Disable()

    # Option button bindings
    self.Bind(wx.EVT_BUTTON, self.onIOTAOptions, self.opt_btn_iota)
    self.Bind(wx.EVT_BUTTON, self.onProcOptions, self.opt_btn_proc)
    self.Bind(wx.EVT_BUTTON, self.onRun, self.btn_run)

    # Add panels to sizers and lay out the window
    self.panel_sizer.Add(self.input_window, 1, flag=wx.EXPAND)
    self.panel_sizer.Add(self.bottom_sizer, flag=wx.EXPAND )
    self.main_sizer.Add(self.panel_sizer, 1, flag=wx.EXPAND | wx.ALL, border=15)
    self.main_sizer.Add((0, 0), flag=wx.EXPAND | wx.ALL, border=5)
    self.Fit()
    self.Layout()

  def onTest(self, e):
    self.test_index = getattr(self, 'test_index', None)
    if not self.test_index:
      self.test_index = make_phil_index(master_phil=pct.get_test_phil(),
                                        working_phil=pct.get_test_phil(),
                                        fetch_new=True)

    scope_names = ['string_definition',
                   'multi_string',
                   'child_scope_alpha',
                   'child_scope_beta']
    test_scope = self.test_index.get_scopes(include=scope_names)
    test_dlg = pct.PHILDialog(self, scope=test_scope,
                              phil_index=self.test_index, name='test')

    if test_dlg.ShowModal() == wx.ID_OK:
      print ('DEBUG: TEST RESULTS: ')
      self.test_index.working_phil.show()
    test_dlg.Destroy()

  def onReset(self, e):
    self.reset_settings()

  def open_options_dialog(self, phil_index, include=None, name=None,
                          title=None):
    if not name:
      if isinstance(include, list):
        name = include[0]
      else:
        name = include
      if not name or name.isspace():
        name = 'options'

    phil_scope = phil_index.get_scopes(include=include)
    phil_dlg = pct.PHILDialog(self,
                              scope=phil_scope,
                              phil_index=phil_index,
                              name=name,
                              title=title)

    if phil_dlg.ShowModal() == wx.ID_OK:
      OK = True
    else:
      OK = False

    phil_dlg.Destroy()
    return OK

  def onPreferences(self, e):
    opt_OK = self.open_options_dialog(phil_index=self.iota_index,
                                      include=['mp', 'gui'],
                                      title='GUI Preferences')
    if opt_OK:
      self.gparams = self.iota_index.get_python_object(make_copy=True)

  def onIOTAOptions(self, e):
    iota_scopes = ['image_import', 'cctbx_xfel', 'analysis', 'advanced']
    opt_OK = self.open_options_dialog(phil_index=self.iota_index,
                                      include=iota_scopes,
                                      title='IOTA Settings')
    if opt_OK:
      self.gparams = self.iota_index.get_python_object(make_copy=True)
      if self.gparams.cctbx_xfel.target is not None:
        try:
          with open(self.gparams.cctbx_xfel.target, 'r') as pf:
            self.target_phil = pf.read()
        except Exception:
          self.target_phil = None
        else:
          try:
            self.bknd_index.update_phil(phil_string=self.target_phil)
          except Exception as e:
            raise Sorry('IOTA PHIL ERROR: ', e)

  def onProcOptions(self, e):
    command_list = [
      ('Spotfinder...',
       lambda evt: self.open_options_dialog(
         phil_index=self.bknd_index, include=['spotfinder'],
         title='Backend Spotfinder Options')),
      ('Indexing...',
       lambda evt: self.open_options_dialog(
         phil_index=self.bknd_index, include=['indexing'],
         title='Backend Indexing Options')),
      ('Refinement...',
       lambda evt: self.open_options_dialog(
         phil_index=self.bknd_index, include=['refinement'],
         title='Backend Refinement Options')),
      ('Integration...',
       lambda evt: self.open_options_dialog(
         phil_index=self.bknd_index, include=['integration'],
         title='Backend Integration Options')),
      ('Advanced...',
       lambda evt: self.open_options_dialog(
         phil_index=self.bknd_index,
         include=['geometry', 'profile', 'prediction', 'significance_filter'],
         title='Advanced Backend Options'))
    ]
    browse_menu = ct.Menu(self)
    browse_menu.add_commands(command_list)
    self.PopupMenu(browse_menu)
    browse_menu.Destroy()

  def init_settings(self):
    ''' Consolidate settings from the main window and initiate run '''

    input_phil_string = self.input_window.GetPHIL(expand=True)
    sel_phil_string = self.option_panel.GetPHIL(expand=True)
    phil_string = input_phil_string + '\n' + sel_phil_string

    self.iota_index.update_phil(phil_string=phil_string)
    self.gparams = self.iota_index.get_python_object(make_copy=True)
    self.iota_phil = self.iota_index.working_phil

    ok_init, self.info, msg = init.initialize_new_run(
      phil=self.iota_phil,
      target_phil=self.bknd_index.get_diff().as_str()
    )
    if not ok_init and msg:
      msg_string = 'IOTA Reported the following error(s):{}\n\n' \
                   ''.format(msg)
      wx.MessageBox(caption='IOTA Errors!',
                    message=msg_string, style=wx.OK | wx.ICON_ERROR)

    return ok_init

  def OnAboutBox(self, e):
    """ About dialog """
    from wx import adv

    info = adv.AboutDialogInfo()
    info.SetName('IOTA')
    info.SetVersion(iota_version)
    info.SetDescription(gui_description)
    info.SetWebSite('http://cci.lbl.gov/xfel')
    info.SetLicense(gui_license)
    info.AddDeveloper('Art Lyubimov')
    info.AddDeveloper('Monarin Uervirojnangkoorn')
    info.AddDeveloper('Aaron Brewster')
    info.AddDeveloper('Nick Sauter')
    info.AddDeveloper('Axel Brunger')
    info.AddDocWriter('Art Lyubimov')
    adv.AboutBox(info)

  def onRecovery(self, e):
    ''' Recover finished / terminated runs '''

    # Find integration folder
    int_folder = os.path.abspath('{}/integration'.format(os.curdir))
    if not os.path.isdir(int_folder):
      open_dlg = wx.DirDialog(self, "Choose the integration folder:",
                              style=wx.DD_DEFAULT_STYLE)
      if open_dlg.ShowModal() == wx.ID_OK:
        int_folder = open_dlg.GetPath()
        open_dlg.Destroy()
      else:
        open_dlg.Destroy()
        return

    # Find paths of runs
    paths = [os.path.join(int_folder, p) for p in os.listdir(int_folder)]
    paths = [p for p in paths if os.path.isdir(p)]

    # Launch recovery dialog
    path_dlg = d.RecoveryDialog(self)
    path_dlg.insert_paths(paths)

    # Recover run
    if path_dlg.ShowModal() == wx.ID_OK:
      selected = path_dlg.selected
      recovery_mode = path_dlg.recovery_mode
      int_path = selected[1]

      # Recover info object
      info = ProcInfo.from_folder(path=int_path)

      if hasattr(info, 'paramfile'):
        self.load_script(filepath=info.paramfile)
      else:
        wx.MessageBox(caption='Run Recovery Error!',
                      message='Cannot recover the parameters for this run',
                      style=wx.OK | wx.ICON_ERROR)
        return

      # Re-open processing window with results of the run
      if recovery_mode == 0:
        title = 'Image Processing Run {}'.format(selected[2])
        self.proc_window = ProcWindow(self, -1, title=title,
                                      phil=self.iota_phil)
        self.proc_window.recover(int_path=info.int_base,
                                 status=selected[0],
                                 params=self.gparams,
                                 info=info)
        self.proc_window.place_and_size(set_by='parent')
        self.proc_window.Show(True)

  def onRun(self, e):
    # Start full run timer
    self.start_time = time.time()

    # Run full processing
    if self.init_settings():
      title = 'Image Processing Run {}'.format(self.info.run_number)
      self.proc_window = ProcWindow(self, -1, title=title, phil=self.iota_phil)
      self.proc_window.info = self.info
      self.proc_window.start_processing()
      self.proc_window.place_and_size(set_by='parent')
      self.proc_window.Show(True)

  def onOutputScript(self, e):

    # Determine param filepath
    save_dlg = wx.FileDialog(self,
                             message="Save IOTA Script",
                             defaultDir=os.curdir,
                             defaultFile="*.param",
                             wildcard="*",
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                             )
    if save_dlg.ShowModal() == wx.ID_OK:
      script_filepath = save_dlg.GetPath()
      self.save_script(script_filepath)

  def save_script(self, filepath):
    ''' Finalize IOTA settings and save to file '''
    self.init_settings()
    with open(filepath, 'w') as sf:
      sf.write(self.iota_phil.as_str())

  def onLoadScript(self, e):
    """
    Widget event either for Load Script menu item or toolbar button
    :param e: event object
    :return:
    """
    load_dlg = wx.FileDialog(self,
                             message="Load script file",
                             defaultDir=os.curdir,
                             defaultFile="*.param",
                             wildcard="*.param",
                             style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if load_dlg.ShowModal() == wx.ID_OK:
      self.load_script(load_dlg.GetPaths()[0])

  def load_script(self, filepath):
    '''Clears settings and loads new settings from IOTA param file

    :param filepath: path to IOTA parameter file
    :return: None
    '''

    self.iota_phil, _ = inp.get_input_phil(paramfile=filepath, gui=True)
    self.update_IOTA_index(phil=self.iota_phil)

    # Pass on target PHIL (if found) to input window
    target = self.gparams.cctbx_xfel.target
    if target:
      try:
        with open(target, 'r') as pf:
          self.target_phil = pf.read()
      except Exception:
        self.target_phil = None
    else:
      self.target_phil = None
    self.update_backend_index(phil=self.target_phil)

    self.input_window.redraw_panel()
    self.option_panel.redraw_panel()

  def reset_settings(self):
    """ Clear all controls in input window """

    # Reset index PHILs to defaults
    self.reset_IOTA_index()
    self.reset_backend_index()

    # Redraw PHIL panels
    self.input_window.redraw_panel(reset=True, exempt=['description'])
    self.option_panel.redraw_panel(reset=True)
    self.gparams = self.iota_index.working_phil.extract()

  def initialize_IOTA_index(self):
    iota_master_phil, _ = inp.get_input_phil(gui=True)
    self.iota_index = make_phil_index(master_phil=iota_master_phil)
    self.gparams = self.iota_index.get_python_object(make_copy=True)

  def initialize_backend_index(self):
    from dials.command_line.stills_process import phil_scope
    self.bknd_index = make_phil_index(master_phil=phil_scope)
    method = self.gparams.advanced.processing_backend
    target_phil, _ = inp.write_defaults(method=method,
                                        write_target_file=False,
                                        write_param_file=False)
    self.bknd_index.update_phil(phil_string=target_phil.as_str())

  def reset_IOTA_index(self):
    description = self.gparams.description
    n_processors = self.gparams.mp.n_processors
    self.iota_index.reset_phil(reindex=True)
    self.gparams = self.iota_index.working_phil.extract()
    self.gparams.description = description
    self.gparams.mp.n_processors = n_processors
    self.iota_index.update_from_python(python_object=self.gparams)

  def reset_backend_index(self):
    self.bknd_index.reset_phil(reindex=True)

  def update_IOTA_index(self, phil=None, input_dict=None, messages=None):
    """ Update IOTA PHIL with params passed on via command line args
    :param phil: A PHIL object or string
    :param input_dict: Dictionary with input paths
    :param messages: Error message list
    """

    # update from existing PHIL first, THEN append the contents of input_dict
    if phil:
      if not isinstance(phil, str):
        try:
          phil = phil.as_str()
        except Exception as e:
          raise Sorry('IOTA GUI ERROR: Cannot read PHIL object! ', e)
      self.iota_index.update_phil(phil_string=phil)

    # Update gparams here, to include new stuff from PHIL
    self.gparams = self.iota_index.get_python_object(make_copy=True)

    # Pass input from IOTA PHIL through Input Finder to resolve wildcards, etc.
    if self.gparams.input is not None:
      inputs = [
        self.gparams.input.pop(i) for i in range(len(self.gparams.input))]
      phil_input_dict = ut.ginp.process_mixed_input(paths=inputs)
      if input_dict:
        input_dict.update(phil_input_dict)
      else:
        input_dict = phil_input_dict

    # Read input OR update from window params
    if input_dict:
      if messages:
        msg_string = 'IOTA Reported the following error(s):{}\n\n' \
                     ''.format('\n'.join(messages))
        wx.MessageBox(caption='IOTA Errors!',
                      message=msg_string, style=wx.OK | wx.ICON_ERROR)

      if len(input_dict['imagepaths']) > 10:
        n_input = len([i for i in os.listdir(self.gparams.output)
                       if i.startswith('input_')])
        input_fn = "input_{:04d}.lst".format(n_input + 1)
        input_list_file = os.path.join(self.gparams.output, input_fn)
        with open(input_list_file, 'w') as lf:
          for f in input_dict['imagefiles']:
            lf.write('{}\n'.format(f))
        self.gparams.input = input_list_file
      else:
        for path in input_dict['imagepaths']:
          if path not in self.gparams.input:
            self.gparams.input.append(path)

    # update n_processors
    if self.gparams.mp.n_processors <= 1:
      self.gparams.mp.n_processors = int(multiprocessing.cpu_count() * 0.75)

    self.iota_index.update_from_python(python_object=self.gparams)

  def update_backend_index(self, phil=None):
    if phil:
      if not isinstance(phil, str):
        try:
          phil = phil.as_str()
        except Exception as e:
          raise Sorry('IOTA GUI ERROR: Cannot read PHIL object! ', e)
      self.bknd_index.update_phil(phil_string=phil)

  def onQuit(self, e):
    # Need a try block in case the C++ portion of the proc window doesn't exist
    try:
      # Check if processing window has been launched
      if hasattr(self, "proc_window"):

        # Export info object with "aborted" status
        info = self.proc_window.info
        info.status = 'aborted'
        info.export_json()

        # Check if proc_thread exists
        if hasattr(self.proc_window, 'proc_thread'):

          # Check if proc_thread is running
          if self.proc_window.proc_thread.is_alive():
            tmp_aborted_file = os.path.join(self.proc_window.tmp_aborted_file)
            with open(tmp_aborted_file, 'w') as tf:
              tf.write('')
            self.proc_window.proc_thread.abort()

            # Close window only when thread is dead
            while self.proc_window.proc_thread.is_alive():
              continue

        import shutil

        try:
          # Clean up temp folders and files if necessary
          if os.path.isdir(self.proc_window.info.tmp_base):
            shutil.rmtree(self.proc_window.info.tmp_base)

          dials_log_dir = os.path.join(self.proc_window.info.log_base,
                                       'logs/dials_logs')
          if os.path.isdir(dials_log_dir):
            shutil.rmtree(dials_log_dir)
        except Exception:
          pass
        print('JOB TERMINATED!')
    except Exception:
      pass

    print ('Closing IOTA...')

    self.Close()

  def onClose(self, e):

    self.DestroyChildren()
    e.Skip()
# ----------------------------  Processing Window ---------------------------  #

class LogTab(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

    self.log_sizer = wx.BoxSizer(wx.VERTICAL)
    self.log_window = rt.RichTextCtrl(self,
                                      style=rt.RE_MULTILINE |
                                            rt.RE_READONLY |
                                            wx.TE_DONTWRAP)
    self.log_window.SetFont(
      wx.Font(norm_font_size, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL,
              wx.FONTWEIGHT_NORMAL, False))
    self.log_sizer.Add(self.log_window, proportion=1, flag=wx.EXPAND | wx.ALL,
                       border=10)

    buttons = {'forward': (-1, 'Forward'),
               'reverse': (-1, 'Reverse')}
    self.find_string = ct.TextCtrlWithButtons(self, buttons=buttons,
                                              ctrl_label='Find String:')
    self.log_sizer.Add(self.find_string, flag=wx.EXPAND | wx.ALL, border=10)
    self.SetSizer(self.log_sizer)

    self.Bind(wx.EVT_BUTTON, self.onSearchForward, self.find_string.btn_forward)
    self.Bind(wx.EVT_BUTTON, self.onSearchReverse, self.find_string.btn_reverse)
    self.Bind(wx.EVT_TEXT_ENTER, self.onEnter, self.find_string.txt_ctrl)

  def onEnter(self, e):
    self.search_forward()

  def onSearchForward(self, e):
    self.search_forward()

  def search_forward(self):
    if self.log_window.GetCaretPosition() == -1:
      self.log_window.SetCaretPosition(0)
    pos = self.log_window.GetCaretPosition()
    search_string = self.find_string.txt_ctrl.GetValue().lower()
    log_string = self.log_window.GetValue()[pos:-1].lower()
    if search_string.replace(' ', '') not in ('', None):
      found_pos = log_string.find(search_string)
      if found_pos == -1:
        if pos > 0:
          msg_text = 'String Not Found! Search From the Top?'
          msg = wx.MessageDialog(None, msg_text, 'Not Found!',
                                 wx.YES_NO | wx.ICON_QUESTION)
          if msg.ShowModal() == wx.ID_YES:
            log_string = self.log_window.GetValue()[0:pos].lower()
            found_pos = log_string.find(search_string)
            if found_pos == -1:
              wx.MessageBox('String Not Found!', 'Not Found!',
                            wx.OK | wx.ICON_EXCLAMATION)
              return
          else:
            return
        else:
          wx.MessageBox('String Not Found!', 'Not Found!',
                        wx.OK | wx.ICON_EXCLAMATION)
          return
      else:
        found_pos += pos
      sel_range = (found_pos, found_pos + len(search_string))
      self.log_window.SetSelectionRange(sel_range)
    else:
      found_pos = 0
    self.log_window.SetCaretPosition(found_pos + len(search_string))
    if not self.log_window.IsPositionVisible(found_pos):
      self.log_window.ShowPosition(found_pos)

  def onSearchReverse(self, e):
    self.search_reverse()

  def search_reverse(self):
    if self.log_window.GetCaretPosition() == -1:
      self.log_window.SetCaretPosition(0)
    pos = self.log_window.GetCaretPosition()
    search_string = self.find_string.txt_ctrl.GetValue().lower()
    full_log = self.log_window.GetValue()
    log_string = full_log[0:pos].lower()
    log_end = len(full_log)
    if search_string.replace(' ', '') not in ('', None):
      found_pos = log_string.rfind(search_string)
      if found_pos == -1:
        if pos < log_end:
          msg_text = 'String Not Found! Search From the Bottom?'
          msg = wx.MessageDialog(None, msg_text, 'Not Found!',
                                 wx.YES_NO | wx.ICON_QUESTION)
          if msg.ShowModal() == wx.ID_YES:
            log_string = full_log[pos:-1].lower()
            found_pos = log_string.rfind(search_string)
            if found_pos == -1:
              wx.MessageBox('String Not Found!', 'Not Found!',
                            wx.OK | wx.ICON_EXCLAMATION)
              return
            else:
              found_pos += pos
          else:
            return
        else:
          wx.MessageBox('String Not Found!', 'Not Found!',
                        wx.OK | wx.ICON_EXCLAMATION)
          return
      sel_range = (found_pos, found_pos + len(search_string))
      self.log_window.SetSelectionRange(sel_range)
    else:
      found_pos = 0
    self.log_window.SetCaretPosition(found_pos - len(search_string))
    if not self.log_window.IsPositionVisible(found_pos):
      self.log_window.ShowPosition(found_pos)


# class ProcessingTab(d.ScrolledPanel):
class ProcessingTab(IOTABasePanel):
  def __init__(self, parent, gparams, *args, **kwargs):
    super(ProcessingTab, self).__init__(parent, *args, **kwargs)
    self.gparams = gparams
    self.info = None
    self.proc_stats = None
    self.parent = parent

    self.hkl_view_axis = 'l'
    self.user_sg = 'P1'
    self.pick = {'image': None, 'index': 0, 'axis': None, 'picked': False}
    self.proc_fnames = None

    self.processed = []
    self.nsref_x = []
    self.nsref_y = []
    self.res_x = []
    self.res_y = []

    self.dblclick = False

    self.main_fig_sizer = wx.GridBagSizer(0, 0)

    # Set regular font
    mpl.rc('font', family='sans-serif', size=plot_font_size)
    mpl.rc('mathtext', default='regular')

    # Integration figure (resolution & reflections / frame)
    self.int_panel = wx.Panel(self)
    int_sizer = wx.BoxSizer(wx.VERTICAL)
    self.int_panel.SetSizer(int_sizer)

    # Image info sizer
    self.info_sizer = wx.GridBagSizer(0, 5)
    self.info_txt = wx.TextCtrl(self.int_panel, style=wx.TE_READONLY)
    view_bmp = bitmaps.fetch_custom_icon_bitmap('image_viewer16')
    r_bmp = bitmaps.fetch_icon_bitmap('actions', '1rightarrow', size=16)
    l_bmp = bitmaps.fetch_icon_bitmap('actions', '1leftarrow', size=16)
    self.btn_right = btn.GenBitmapButton(self, bitmap=r_bmp)
    self.btn_left = btn.GenBitmapButton(self, bitmap=l_bmp)
    self.btn_viewer = btn.GenBitmapButton(self, bitmap=view_bmp)
    self.info_sizer.Add(self.info_txt, pos=(0, 1), flag=wx.EXPAND)
    self.info_sizer.Add(self.btn_left, pos=(0, 2))
    self.info_sizer.Add(self.btn_right, pos=(0, 3))
    self.info_sizer.Add(self.btn_viewer, pos=(0, 4))
    self.info_sizer.AddGrowableCol(1)
    int_sizer.Add(self.info_sizer, flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                  border=10)
    self.info_txt.Hide()
    self.btn_right.Hide()
    self.btn_left.Hide()
    self.btn_viewer.Hide()

    self.btn_right.Disable()
    self.btn_left.Disable()
    self.btn_viewer.Disable()

    self.Bind(wx.EVT_BUTTON, self.onImageView, self.btn_viewer)
    self.Bind(wx.EVT_BUTTON, self.onArrow, self.btn_right)
    self.Bind(wx.EVT_BUTTON, self.onArrow, self.btn_left)

    # Charts
    # For some reason MatPlotLib 2.2.3 on GTK 6 does not create transparent
    # patches (tried set_visible(False), set_facecolor("none"), set_alpha(0),
    # to no avail). Thus, for GTK only, setting facecolor to background
    # color; otherwide to 'none'. This may change if other tests reveal the
    # same problem in other systems.
    if wx.Platform == '__WXGTK__':
      bg_color = [i/255 for i in self.GetBackgroundColour()]
    else:
      bg_color = 'none'

    # Integration plot
    self.int_figure = Figure(figsize=(1, 2.5), facecolor=bg_color)
    int_gsp = gridspec.GridSpec(2, 1, wspace=0, hspace=0)

    # Resolution / No. strong reflections plot
    self.res_axes = self.int_figure.add_subplot(int_gsp[0])
    self.res_axes.set_ylabel('Resolution ({})'.format(r'$\AA$'))
    self.res_axes.tick_params(labelbottom=False)

    self.nsref_axes = self.int_figure.add_subplot(int_gsp[1])
    self.nsref_axes.set_xlabel('Frame')
    self.nsref_axes.set_ylabel(
      'Spots (I/{0}(I)>{1})'
      ''.format(r'$\sigma$',
                self.gparams.data_selection.image_triage.strong_sigma))

    # Initialize blank charts, medians, and picks
    self.res_chart = self.res_axes.plot([], [], 'o', c='#0571b0', zorder=1,
                                        mec='black', picker=5)[0]
    self.res_med = self.res_axes.axhline(0, zorder=0, c='#0571b0', ls=':')
    self.res_pick, = self.res_axes.plot(0, 0, 'o', zorder=2, ms=12, alpha=0.5,
                                        mec='black', c='yellow', visible=False)
    self.nsref_chart = self.nsref_axes.plot([], [], 'o', c='#ca0020', zorder=1,
                                            mec='black', picker=5)[0]
    self.nsref_med = self.nsref_axes.axhline(0, zorder=0, c='#ca0020', ls=':')
    self.nsref_pick, = self.nsref_axes.plot(0, 0, 'o', ms=12, alpha=0.5,
                                            mec='black', zorder=2, c='yellow',
                                            visible=False)
    self.res_axes.set_autoscaley_on(True)
    self.nsref_axes.set_autoscaley_on(True)

    self.info_txt.Show()
    self.btn_right.Show()
    self.btn_left.Show()
    self.btn_viewer.Show()

    self.int_canvas = FigureCanvas(self.int_panel, -1, self.int_figure)
    int_sizer.Add(self.int_canvas, 1, flag=wx.EXPAND)
    self.int_canvas.draw()
    self.int_figure.set_tight_layout(True)

    # Wilson (<I> vs. res) plot
    self.wp_panel = wx.Panel(self)
    wp_sizer = wx.BoxSizer(wx.VERTICAL)
    self.wp_panel.SetSizer(wp_sizer)

    self.wp_figure = Figure(figsize=(0.3, 0.6), facecolor=bg_color)
    self.wp_axes = self.wp_figure.add_subplot(111)
    self.wp_axes.set_ylabel("<I>")

    self.wp_figure.set_tight_layout(True)
    self.wp_canvas = FigureCanvas(self.wp_panel, -1, self.wp_figure)
    wp_sizer.Add(self.wp_canvas, 1, flag=wx.EXPAND)

    # HKL (or slice) plot
    self.hkl_panel = wx.Panel(self)
    hkl_sizer = wx.BoxSizer(wx.VERTICAL)
    self.hkl_panel.SetSizer(hkl_sizer)

    self.hkl_figure = Figure(figsize=(0.3, 0.3), facecolor=bg_color)
    self.hkl_axes = self.hkl_figure.add_subplot(111)
    self.hkl_axes.set_xticks([])
    self.hkl_axes.set_yticks([])
    self.hkl_axes.set_frame_on(False)

    self.hkl_canvas = FigureCanvas(self.hkl_panel, -1, self.hkl_figure)
    hkl_sizer.Add(self.hkl_canvas, 1, flag=wx.EXPAND)

    self.hkl_sg = ct.OptionCtrl(self.hkl_panel, items=[('sg', 'P1')],
                                checkbox=True, checkbox_label='Space Group: ',
                                label_size=wx.DefaultSize,
                                ctrl_size=wx.DefaultSize)
    hkl_sizer.Add(self.hkl_sg, flag=wx.ALIGN_CENTER)

    # Proc Tab bindings
    self.Bind(wx.EVT_TEXT_ENTER, self.onSGTextEnter, self.hkl_sg.sg)
    self.Bind(wx.EVT_CHECKBOX, self.onSGCheckbox, self.hkl_sg.toggle)

    # Processing summary figure
    self.proc_panel = wx.Panel(self, size=(-1, 120))
    proc_sizer = wx.BoxSizer(wx.VERTICAL)
    self.proc_panel.SetSizer(proc_sizer)

    self.proc_figure = Figure(figsize=(0.5, 2.5), facecolor=bg_color)
    self.sum_axes = self.proc_figure.add_subplot(111)
    self.sum_axes.axis('off')

    self.bracket = mpatches.Rectangle((0, 0), 1, 1, fc='white',
                                      ec='black', lw=2, visible=False)

    self.proc_canvas = FigureCanvas(self.proc_panel, -1, self.proc_figure)
    self.proc_canvas.draw()
    self.proc_figure.set_tight_layout(True)
    proc_sizer.Add(self.proc_canvas, flag=wx.EXPAND | wx.BOTTOM, border=10)

    self.main_fig_sizer.Add(self.int_panel, pos=(0, 0), span=(2, 6),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.wp_panel, pos=(2, 0), span=(2, 4),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.hkl_panel, pos=(2, 4), span=(2, 2),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.proc_panel, pos=(4, 0), span=(1, 6),
                            flag=wx.EXPAND)

    self.main_fig_sizer.AddGrowableCol(0)
    self.main_fig_sizer.AddGrowableCol(1)
    self.main_fig_sizer.AddGrowableCol(2)
    self.main_fig_sizer.AddGrowableCol(3)
    self.main_fig_sizer.AddGrowableCol(4)
    self.main_fig_sizer.AddGrowableCol(5)
    self.main_fig_sizer.AddGrowableRow(1)
    self.main_fig_sizer.AddGrowableRow(2)
    self.main_fig_sizer.AddGrowableRow(3)

    cid = self.int_canvas.mpl_connect('pick_event', self.on_pick)
    sid = self.proc_canvas.mpl_connect('pick_event', self.on_bar_pick)
    xid = self.int_canvas.mpl_connect('button_press_event',
                                      self.on_button_press)
    xid = self.hkl_canvas.mpl_connect('button_press_event', self.on_hkl_press)
    xid = self.proc_canvas.mpl_connect('button_press_event',
                                       self.on_button_press)
    xid = self.proc_canvas.mpl_connect('button_release_event',
                                       self.on_button_release)

    self.SetSizer(self.main_fig_sizer)

  def _update_canvas(self, canvas, draw_idle=True):
    """ Update a canvas (passed as arg)
    :param canvas: A canvas to be updated via draw_idle
    """
    # Draw_idle is useful for regular updating of the chart; straight-up draw
    # without flush_events() will have to be used when buttons are clicked to
    # avoid recursive calling of wxYield
    if draw_idle:
      canvas.draw_idle()
      try:
        canvas.flush_events()
      except (NotImplementedError, AssertionError):
        pass
    else:
      canvas.draw()
    canvas.Refresh()

  def onSGTextEnter(self, e):
    user_sg = str(self.hkl_sg.sg.GetValue())

    # Attempt to create a space group object; if symbol or number are
    # invalid, this will fail
    from cctbx import crystal
    if user_sg:
      try:
        sym = crystal.symmetry(space_group_symbol=str(user_sg))
      except RuntimeError:
        pass
      else:
        self.user_sg = str(sym.space_group_info())
    self.hkl_sg.sg.SetValue(self.user_sg)
    self.draw_measured_indices()

  def onSGCheckbox(self, e):
    if e.IsChecked():
      if self.gparams.cctbx_xfel.target_space_group:
        self.hkl_sg.sg.SetValue(
          str(self.gparams.cctbx_xfel.target_space_group))
      self.user_sg = str(self.hkl_sg.sg.GetValue())
    else:
      self.user_sg = 'P1'
    self.draw_measured_indices()

  def draw_summary(self):

    try:
      n_img = len(self.info.categories['total'][0])
      nonzero = []
      names = []
      patches = []
      cat_keys = ('failed_triage', 'failed_spotfinding', 'failed_indexing',
                  'failed_grid_search', 'failed_integration', 'failed_filter',
                  'integrated', 'not_processed')

      for key in cat_keys:
        name = self.info.categories[key][1]
        color = self.info.categories[key][3]
        n_in_cat = len(self.info.categories[key][0])
        if color is not None and n_in_cat > 0:
          nonzero.append([name, color, n_in_cat])
          names.append(name)
          patches.append(None)

      self.sum_axes.clear()
      self.sum_axes.axis([0, 1.01 * n_img, -0.5, 0.75])
      self.sum_axes.axis('off')

      for i in range(len(nonzero)):
        percent = np.round(nonzero[i][2] / n_img * 100)
        previous = [j[2] for j in nonzero[:i]]
        lf = np.sum(previous)
        barh = self.sum_axes.barh(0, nonzero[i][2],
                                  left=lf,
                                  color=nonzero[i][1],
                                  align='center',
                                  label='{}%'.format(percent),
                                  picker=True)
        patch = barh[0]
        bl = patch.get_xy()
        patches[i] = ((bl[0], bl[0] + patch.get_width()),
                      (bl[1], bl[1] + patch.get_height()))

        x = 0.5 * patch.get_width() + bl[0]
        y = 0.5 * patch.get_height() + bl[1]
        self.sum_axes.text(x, y, "{}%".format(percent),
                           ha='center', va='center')
      self.sum_axes.legend(names,
                           ncol=len(nonzero),
                           bbox_to_anchor=(0, -1, 1, 1),
                           loc='upper center', fontsize=8,
                           frameon=False, handlelength=1,
                           mode='expand')

      if self.bracket.get_visible():
        idx = self.pick['index']
        p_idx = [patches.index(i) for i in patches if
                 i[0][0] <= idx <= i[0][1]][0]
        cat = names[p_idx].replace(' ', '_')
        self.proc_fnames = self.info.categories[cat][0]

        px = patches[p_idx][0][0]
        pw = patches[p_idx][0][1] - patches[p_idx][0][0]
        py = patches[p_idx][1][0]
        ph = patches[p_idx][1][1] - patches[p_idx][1][0]
        self.bracket.set_bounds(px, py, pw, ph)
        self.sum_axes.add_patch(self.bracket)

      self._update_canvas(canvas=self.proc_canvas)

    except ValueError as e:
      print('SUMMARY PLOT ERROR: ', e)
      return

  def draw_plots(self):
    self.draw_integration_plots()
    self.draw_b_factors()
    self.draw_measured_indices()

  def draw_integration_plots(self):
    # Extract data and make arrays:
    try:
      if self.info.stats:
        idx = None
        filenames = None
        img_idx = 0
        spt = None
        res = None

        # Strong reflections
        if self.info.stats['strong']['lst']:
          idx, filenames, img_idx, spt = zip(*self.info.stats['strong']['lst'])
          self.nsref_x = np.append(self.nsref_x,
                                   np.array(idx).astype(np.double))
          self.nsref_y = np.append(self.nsref_y,
                                   np.array(spt).astype(np.double))

        # Resolution
        if self.info.stats['res']['lst']:
          idx, filenames, img_idx, res = zip(*self.info.stats['res']['lst'])
          self.res_x = np.append(self.res_x, np.array(idx).astype(np.double))
          self.res_y = np.append(self.res_y, np.array(res).astype(np.double))

        # Update arrays
        if (idx and filenames and spt and res):
          self.processed.extend(zip(idx, filenames, img_idx, spt, res))
          res_median = np.median(self.res_y)
          nsref_median = np.median(self.nsref_y)
        else:
          return
    except ValueError as e:
      print('IOTA PLOTTING (PROC) ERROR: ', e)
    else:
      # Resolution per frame
      res_m = np.isfinite(self.res_y)
      self.res_chart.set_xdata(self.res_x[res_m])
      self.res_chart.set_ydata(self.res_y[res_m])
      self.res_med.set_ydata(res_median)

      self.res_axes.set_xlim(0, np.nanmax(self.res_x) + 2)
      res_ymax = np.nanmax(self.res_y) * 1.1
      res_ymin = np.nanmin(self.res_y) * 0.9
      if res_ymin == res_ymax:
        res_ymax = res_ymin + 1
      self.res_axes.set_ylim(ymin=res_ymin, ymax=res_ymax)
      self.res_axes.draw_artist(self.res_chart)

      # Plot strong reflections per frame
      self.nsref_chart.set_xdata(self.nsref_x)
      self.nsref_chart.set_ydata(self.nsref_y)
      self.nsref_med.set_ydata(nsref_median)

      self.nsref_axes.set_xlim(0, np.nanmax(self.nsref_x) + 2)
      nsref_ymax = np.nanmax(self.nsref_y) * 1.25 + 10
      if nsref_ymax == 0:
        nsref_ymax = 100
      self.nsref_axes.set_ylim(ymin=0, ymax=nsref_ymax)
      self.nsref_axes.draw_artist(self.nsref_chart)

      self._update_canvas(canvas=self.int_canvas)

  def draw_b_factors(self):
    self.wp_axes.clear()
    self.wp_axes.set_xlabel('B-factor')
    self.wp_axes.set_ylabel('Count')
    self.wp_axes.set_title('Wilson B-factor Histogram')
    if self.info.b_factors:
      self.wp_axes.hist(self.info.b_factors, 50, density=False,
                        facecolor='#4575b4', histtype='stepfilled')

    self._update_canvas(canvas=self.wp_canvas)

  def draw_measured_indices(self):
    # Draw a h0, k0, or l0 slice of merged data so far
    self.hkl_axes.clear()
    self.hkl_axes.axhline(0, lw=0.5, c='black', ls='-')
    self.hkl_axes.axvline(0, lw=0.5, c='black', ls='-')
    self.hkl_axes.set_xticks([])
    self.hkl_axes.set_yticks([])

    try:
      self.hkl_colorbar.remove()
    except AttributeError:
      pass
    except KeyError:
      pass

    try:
      hkl_slice = self.info.get_hkl_slice(sg=self.user_sg,
                                          axis=self.hkl_view_axis)
    except AttributeError as e:
      print ('IOTA HKL PLOT ERROR: ', e)
      return

    try:
      hkl, freq = zip(*hkl_slice)
    except ValueError:
      hkl = [(0, 0, 0)]
      freq = 1
    except TypeError:
      return

    h, k, l = zip(*hkl)
    if self.hkl_view_axis == 'l':
      x = h
      y = k
      self.hkl_axes.set_xlabel('h', weight='bold')
      self.hkl_axes.set_ylabel('k', weight='bold', rotation='horizontal')
    elif self.hkl_view_axis == 'k':
      x = h
      y = l
      self.hkl_axes.set_xlabel('h', weight='bold')
      self.hkl_axes.set_ylabel('l', weight='bold', rotation='horizontal')
    elif self.hkl_view_axis == 'h':
      x = k
      y = l
      self.hkl_axes.set_xlabel('k', weight='bold')
      self.hkl_axes.set_ylabel('l', weight='bold', rotation='horizontal')
    else:
      # if for some reason this goes to plotting without any indices
      x = np.zeros(500)
      y = np.zeros(500)
      freq = np.zeros(500)

    # Format plot
    pt_size = int(max(self.hkl_panel.GetSize()) / 100)
    if pt_size == 0:
      pt_size = 1

    hkl_scatter = self.hkl_axes.scatter(x, y, c=freq, cmap="jet",
                                        s=pt_size, edgecolor='none')
    self.hkl_axes.xaxis.set_label_coords(x=1, y=0.5)
    self.hkl_axes.yaxis.set_label_coords(x=0.5, y=1)

    try:
      xmax = abs(max(x, key=abs))
      ymax = abs(max(y, key=abs))

      # Zero values will result in a matplotlib UserWarning; set to an
      # arbitrarily small number to avoid this
      if xmax == 0:
        xmax += 0.00001
      if ymax == 0:
        ymax += 0.00001

      self.hkl_axes.set_xlim(xmin=-xmax, xmax=xmax)
      self.hkl_axes.set_ylim(ymin=-ymax, ymax=ymax)
    except ValueError:
      pass

    vmax = 2 if np.max(freq) <= 2 else np.max(freq)
    norm = colors.Normalize(vmin=1, vmax=vmax)
    self.hkl_colorbar = self.hkl_figure.colorbar(hkl_scatter,
                                                 ax=self.hkl_axes,
                                                 cmap='jet',
                                                 norm=norm,
                                                 orientation='vertical',
                                                 aspect=40)

    self._update_canvas(canvas=self.hkl_canvas)

  def onImageView(self, e):
    # filepath = self.info_txt.GetValue().split(' ')[0]
    # viewer = self.gparams.gui.image_viewer
    # if os.path.isfile(filepath):
    #   viewer = thr.ImageViewerThread(self,
    #                                  viewer=viewer,
    #                                  file_string=filepath)
    #   viewer.start()

    idx = str(self.pick['index'])
    exp = self.info.pointers[idx]['experiments']
    ref = self.info.pointers[idx]['reflections']
    file_string = 'experiments={} reflections={}'.format(exp, ref)
    if os.path.isfile(exp) and os.path.isfile(ref):
      viewer = thr.ImageViewerThread(self,
                                     viewer=self.gparams.gui.image_viewer,
                                     file_string=file_string)
      viewer.start()

  def view_proc_images(self):
    filenames = None
    if len(self.proc_fnames) > 5:
      view_warning = d.ViewerWarning(self.parent, len(self.proc_fnames))
      if view_warning.ShowModal() == wx.ID_OK:
        # parse 'other' entry
        img_no_string = str(view_warning.no_images).split(',')
        filenames = []
        for n in img_no_string:
          if '-' in n:
            img_limits = n.split('-')
            start = int(min(img_limits))
            end = int(max(img_limits))
            if start <= len(self.proc_fnames) and end <= len(self.proc_fnames):
              filenames.extend(self.proc_fnames[start:end])
          else:
            if int(n) <= len(self.proc_fnames):
              filenames.append(self.proc_fnames[int(n)])
    elif 0 < len(self.proc_fnames) <= 5:
      filenames = self.proc_fnames

    # handle filepath tuples
    if isinstance(filenames, list) or isinstance(filenames, tuple):
      filenames = zip(*filenames)[1]

    if filenames:
      file_string = ' '.join(filenames)
      viewer = self.gparams.gui.image_viewer
      viewer = thr.ImageViewerThread(self,
                                     viewer=viewer,
                                     file_string=file_string)
      viewer.start()

  def onArrow(self, e):
    idx = self.pick['index']
    self.pick['image'] = None
    if e.GetId() == self.btn_right.GetId():
      direction = 1
    elif e.GetId() == self.btn_left.GetId():
      direction = -1

    search = True
    while search:
      idx += direction
      if idx > len(self.processed) or idx < 0:
        break
      else:
        entry = [i for i in self.processed if i[0] == idx]
        if entry:
          search = False
          point_idx, img, img_idx, spt, res = entry[0]
          self.info_txt.SetValue(img)
          self.pick['index'] = idx
          self.pick['image'] = img
          self.pick['image_index'] = img_idx
          self.nsref_pick.set_data(idx, spt)
          self.res_pick.set_data(idx, res)
        else:
          search = True

    self._update_canvas(canvas=self.int_canvas, draw_idle=False)

  def on_pick(self, event):
    self.bracket.set_visible(False)
    self.nsref_pick.set_visible(True)
    self.res_pick.set_visible(True)
    self._update_canvas(canvas=self.proc_canvas, draw_idle=False)

    idx = int(round(event.mouseevent.xdata))
    entry = [i for i in self.processed if i[0] == idx]
    if entry:
      point_idx, img, img_idx, spt, res = entry[0]
      self.pick['picked'] = True
      if event.mouseevent.inaxes == self.nsref_axes:
        self.pick['axis'] = 'nsref'
      elif event.mouseevent.inaxes == self.res_axes:
        self.pick['axis'] = 'res'
      self.pick['image'] = img
      self.pick['index'] = idx
      self.pick['img_idx'] = img_idx
      self.nsref_pick.set_data(point_idx, spt)
      self.res_pick.set_data(point_idx, res)
      txt = '{} ({})'.format(img, img_idx)
      self.toggle_pick(enabled=True, txt=txt)

    self._update_canvas(canvas=self.int_canvas, draw_idle=False)

  def on_bar_pick(self, event):
    self.show_image_group(e=event.mouseevent)

  def show_image_group(self, e):
    self.pick['picked'] = True
    if e.inaxes == self.sum_axes:
      self.pick['axis'] = 'summary'
      self.pick['index'] = e.xdata
      self.bracket.set_visible(True)
      self.draw_summary()
    self.toggle_pick(enabled=False)

  def toggle_pick(self, enabled=False, txt=''):
    self.nsref_pick.set_visible(enabled)
    self.res_pick.set_visible(enabled)
    self.info_txt.SetValue(txt)

    if enabled:
      self.btn_left.Enable()
      self.btn_right.Enable()
      self.btn_viewer.Enable()
    else:
      self.btn_left.Disable()
      self.btn_right.Disable()
      self.btn_viewer.Disable()

    self._update_canvas(canvas=self.int_canvas, draw_idle=False)
    self.Refresh()

  def on_hkl_press(self, event):
    if event.inaxes == self.hkl_axes:
      if self.hkl_view_axis == 'h':
        self.hkl_view_axis = 'k'
      elif self.hkl_view_axis == 'k':
        self.hkl_view_axis = 'l'
      elif self.hkl_view_axis == 'l':
        self.hkl_view_axis = 'h'
      self.draw_measured_indices()
    self.Refresh()

  def on_button_press(self, event):
    if event.button != 1:
      self.pick['picked'] = False
      if event.inaxes == self.sum_axes:
        self.bracket.set_visible(False)
        self._update_canvas(canvas=self.proc_canvas, draw_idle=False)
      elif event.inaxes in (self.nsref_axes, self.res_axes):
        self.nsref_pick.set_visible(False)
        self.res_pick.set_visible(False)
        self.info_txt.SetValue('')
        self.btn_left.Disable()
        self.btn_right.Disable()
        self.btn_viewer.Disable()
        self._update_canvas(canvas=self.int_canvas, draw_idle=False)

    if event.dblclick:
      self.dblclick = True
    else:
      self.dblclick = False

    self.Refresh()

  def on_button_release(self, event):
    if event.button == 1 and self.dblclick:
      self.dblclick = False
      if not self.bracket.get_visible():
        self.show_image_group(e=event)
      self.view_proc_images()


class LiveAnalysisTab(IOTABaseScrolledPanel):
  def __init__(self, parent, gparams=None, *args, **kwargs):
    self.parent = parent
    self.info = None
    self.gparams = gparams

    super(LiveAnalysisTab, self).__init__(parent, *args, **kwargs)
    self.main_fig_sizer = wx.GridBagSizer(0, 0)

    # Set regular font
    mpl.rc('font', family='sans-serif', size=plot_font_size)
    mpl.rc('mathtext', default='regular')

    # For some reason MatPlotLib 2.2.3 on GTK 6 does not create transparent
    # patches (tried set_visible(False), set_facecolor("none"), set_alpha(0),
    # to no avail). Thus, for GTK only, setting facecolor to background
    # color; otherwide to 'none'. This may change if other tests reveal the
    # same problem in other systems.
    if wx.Platform == '__WXGTK__':
      bg_color = [i/255 for i in self.GetBackgroundColour()]
    else:
      bg_color = 'none'

    # UC Histogram / cluster figure
    self.uc_panel = wx.Panel(self)
    uc_box = wx.StaticBox(self.uc_panel, label='Unit Cell Histograms')
    uc_sizer = wx.StaticBoxSizer(uc_box, wx.VERTICAL)
    self.uc_panel.SetSizer(uc_sizer)
    self.uc_figure = Figure(figsize=(1, 2.5), facecolor=bg_color)

    uc_gsub = gridspec.GridSpec(2, 3, wspace=0, hspace=0)
    self.a_axes = self.uc_figure.add_subplot(uc_gsub[0])
    self.a_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.a_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.b_axes = self.uc_figure.add_subplot(uc_gsub[1], sharey=self.a_axes)
    self.b_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.b_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.b_axes.set_yticklabels(list('' * 5), visible=False)
    self.c_axes = self.uc_figure.add_subplot(uc_gsub[2], sharey=self.a_axes)
    self.c_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.c_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.c_axes.set_yticklabels(list('' * 5), visible=False)
    self.alpha_axes = self.uc_figure.add_subplot(uc_gsub[3])
    self.alpha_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.alpha_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.beta_axes = self.uc_figure.add_subplot(uc_gsub[4],
                                                sharey=self.alpha_axes)
    self.beta_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.beta_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.beta_axes.set_yticklabels(list('' * 5), visible=False)
    self.gamma_axes = self.uc_figure.add_subplot(uc_gsub[5],
                                                 sharey=self.alpha_axes)
    self.gamma_axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    self.gamma_axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.gamma_axes.set_yticklabels(list('' * 5), visible=False)

    self.uc_canvas = FigureCanvas(self.uc_panel, -1, self.uc_figure)
    self.uc_figure.set_tight_layout(True)
    uc_sizer.Add(self.uc_canvas, 1, flag=wx.EXPAND)

    # UC Clustering Result
    self.cluster_panel = wx.Panel(self)
    cluster_box = wx.StaticBox(self.cluster_panel, label='Unit Cell Clustering')
    cluster_box_sizer = wx.StaticBoxSizer(cluster_box, wx.VERTICAL)
    self.cluster_panel.SetSizer(cluster_box_sizer)
    self.cluster_list = ct.CustomListCtrl(self.cluster_panel, size=(-1, 100),
                                          content_font_size=norm_font_size)
    font = wx.Font(norm_font_size, wx.FONTFAMILY_TELETYPE,
                              wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
    self.cluster_list.SetFont(font)
    self.cluster_list.ctr.InsertColumn(0, "#")
    self.cluster_list.ctr.InsertColumn(1, "Lattice")
    self.cluster_list.ctr.InsertColumn(2, "Unit Cell", width=200)
    self.cluster_list.ctr.setResizeColumn(3)
    cluster_box_sizer.Add(self.cluster_list, proportion=1, flag=wx.EXPAND)

    # PRIME result
    self.tb1_panel = wx.Panel(self)
    tb1_box = wx.StaticBox(self.tb1_panel, label='LivePRIME Results')
    tb1_sizer = wx.StaticBoxSizer(tb1_box, wx.VERTICAL)
    self.tb1_panel.SetSizer(tb1_sizer)
    self.tb1_plot = ppl.Plotter(self.tb1_panel, params=None, info=None)
    self.tb1_plot.initialize_figure(figsize=(0.1, 0.1))
    tb1_sizer.Add(self.tb1_plot, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Analysis Options
    self.opt_panel = wx.Panel(self)
    self.opt_sizer = wx.GridBagSizer(0, 0)
    self.opt_panel.SetSizer(self.opt_sizer)
    self.pg_uc = ct.OptionCtrl(self.opt_panel,
                               items=[('pg', ''), ('uc', '')],
                               sub_labels=('symmetry', 'cell'),
                               grid=(2, 2),
                               expand_cols=(1),
                               label_size=wx.DefaultSize,
                               ctrl_size=wx.DefaultSize)
    self.btn_run_analysis = wx.Button(self.opt_panel, label='Run Analysis')
    self.opt_sizer.Add(self.pg_uc, pos=(0, 0), span=(2, 1),
                       flag=wx.EXPAND | wx.ALL, border=5)
    self.opt_sizer.Add(self.btn_run_analysis, pos=(0, 1), span=(1, 1),
                       flag=wx.ALL, border=5)
    self.opt_sizer.AddGrowableCol(0)

    self.main_fig_sizer.Add(self.uc_panel, pos=(0, 0), span=(2, 4),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.cluster_panel, pos=(2, 0), span=(3, 2),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.tb1_panel, pos=(2, 2), span=(2, 2),
                            flag=wx.EXPAND)
    self.main_fig_sizer.Add(self.opt_panel, pos=(4, 2), span=(1, 2),
                            flag=wx.EXPAND)

    self.main_fig_sizer.AddGrowableCol(0)
    self.main_fig_sizer.AddGrowableCol(1)
    self.main_fig_sizer.AddGrowableCol(2)
    self.main_fig_sizer.AddGrowableRow(1)
    self.main_fig_sizer.AddGrowableRow(2)
    self.main_fig_sizer.AddGrowableRow(3)

    self.SetSizer(self.main_fig_sizer)

  def draw_plots(self):
    if self.info and self.info.init_proc:
      if self.info.categories['integrated'][0]:
        self.draw_uc_histograms()
      if hasattr(self.info, 'clusters') and self.info.clusters:
        self.report_clustering_results()
      if hasattr(self.info, 'prime_info') and self.info.prime_info:
        self.report_prime_results()

      self.SetupScrolling()
      self.Refresh()

  def report_clustering_results(self):
    self.cluster_list.ctr.DeleteAllItems()
    for c in self.info.clusters:
      i = self.info.clusters.index(c)
      idx = self.cluster_list.ctr.InsertStringItem(i, str(c['number']))
      self.cluster_list.ctr.SetStringItem(idx, 1, str(c['pg']))
      self.cluster_list.ctr.SetStringItem(idx, 2, str(c['uc']))

    self.cluster_list.ctr.SetColumnWidth(0, width=-1)
    self.cluster_list.ctr.SetColumnWidth(1, width=-1)
    self.cluster_list.ctr.SetColumnWidth(2, width=-1)


  def report_prime_results(self):

    try:
      self.tb1_plot.table.remove()
    except Exception:
      pass

    try:
      self.tb1_plot.info = self.info.prime_info
      self.tb1_plot.table_one_figure()
    except Exception as e:
      import traceback
      traceback.print_exc()
      print('PRIME PLOTTER ERROR: ', e)


  def calculate_uc_histogram(self, a, axes, xticks_loc='top', set_ylim=False):
    # n, bins = np.histogram(a, 50)
    # left = np.array(bins[:-1])
    # right = np.array(bins[1:])
    # bottom = np.zeros(len(left))
    # top = bottom + n
    # XY = np.array([[left, left, right, right], [bottom, top, top, bottom]]).T
    # barpath = path.Path.make_compound_path_from_polys(XY)
    # patch = mpatches.PathPatch(barpath, fc='#4575b4', lw=0, alpha=0.75)
    # axes.add_patch(patch)

    # axes.set_xlim(left[0], right[-1])
    # if set_ylim:
    #   axes.set_ylim(bottom.min(), 1.05 * top.max())

    axes.hist(a, 50, density=False, facecolor='#4575b4', histtype='stepfilled')
    axes.xaxis.get_major_ticks()[0].label1.set_visible(False)
    axes.xaxis.get_major_ticks()[-1].label1.set_visible(False)

    if xticks_loc == 'top':
      axes.xaxis.tick_top()
    elif xticks_loc == 'bottom':
      axes.xaxis.tick_bottom()

  def draw_uc_histograms(self):
    try:
      # Unit cell histograms
      self.a_axes.clear()
      self.b_axes.clear()
      self.c_axes.clear()
      self.alpha_axes.clear()
      self.beta_axes.clear()
      self.gamma_axes.clear()

      a, b, c, alpha, beta, gamma, sg = zip(*self.info.cluster_iterable)

      self.calculate_uc_histogram(a, self.a_axes, set_ylim=True)
      edge_ylabel = 'a, b, c ({})'.format(r'$\AA$')
      self.a_axes.set_ylabel(edge_ylabel)

      self.calculate_uc_histogram(b, self.b_axes)
      self.b_axes.set_yticklabels(list('' * 5), visible=False)

      self.calculate_uc_histogram(c, self.c_axes)
      self.c_axes.set_yticklabels(list('' * 5), visible=False)

      self.calculate_uc_histogram(alpha, self.alpha_axes,
                                  xticks_loc='bottom', set_ylim=True)
      ang_ylabel = '{}, {}, {} ({})'.format(r'$\alpha$', r'$\beta$',
                                            r'$\gamma$', r'$^\circ$')
      self.alpha_axes.set_ylabel(ang_ylabel)

      self.calculate_uc_histogram(beta, self.beta_axes, xticks_loc='bottom')
      self.beta_axes.set_yticklabels(list('' * 5), visible=False)

      self.calculate_uc_histogram(gamma, self.gamma_axes, xticks_loc='bottom')
      self.gamma_axes.set_yticklabels(list('' * 5), visible=False)

      self.uc_canvas.draw()
      self.uc_canvas.Refresh()

    except ValueError as e:
      print('UC HISTOGRAM ERROR: ', e)


class SummaryTab(IOTABaseScrolledPanel):
  def __init__(self, parent, info, gparams, *args, **kwargs):
    super(SummaryTab, self).__init__(parent, *args, **kwargs)

    self.parent = parent
    self.info = info
    self.gparams = gparams

    summary_sizer = wx.BoxSizer(wx.VERTICAL)

    sfont = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,
                    wx.FONTWEIGHT_NORMAL)
    bfont = wx.Font(norm_font_size, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL,
                    wx.FONTWEIGHT_BOLD)
    self.SetFont(bfont)

    # Run information
    run_box = wx.StaticBox(self, label='Run Information')
    run_box.SetFont(sfont)
    run_box_sizer = wx.StaticBoxSizer(run_box, wx.VERTICAL)
    run_box_grid = wx.FlexGridSizer(3, 2, 5, 20)
    self.title_txt = wx.StaticText(self, label='')
    self.title_txt.SetFont(sfont)
    self.folder_txt = wx.StaticText(self, label='')
    self.folder_txt.SetFont(sfont)

    run_box_grid.AddMany([(wx.StaticText(self, label='Title')),
                          (self.title_txt, 1, wx.EXPAND),
                          (wx.StaticText(self, label='Directory')),
                          (self.folder_txt, 1, wx.EXPAND)])

    run_box_grid.AddGrowableCol(1, 1)
    run_box_sizer.Add(run_box_grid, flag=wx.EXPAND | wx.ALL, border=10)

    summary_sizer.Add(run_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Integration summary
    int_box = wx.StaticBox(self, label='Analysis of Integration')
    int_box.SetFont(sfont)
    int_box_sizer = wx.StaticBoxSizer(int_box, wx.HORIZONTAL)
    self.int_box_grid = wx.GridBagSizer(5, 20)

    # Button & binding for heatmap display
    heatmap_bmp = bitmaps.fetch_custom_icon_bitmap('heatmap24')
    self.int_heatmap = ct.GradButton(self,
                                     bmp=heatmap_bmp,
                                     label='  Spotfinding Heatmap',
                                     size=(250, -1))
    self.int_box_grid.Add(self.int_heatmap, pos=(0, 1), flag=wx.ALIGN_RIGHT)
    self.Bind(wx.EVT_BUTTON, self.onPlotHeatmap, self.int_heatmap)
    self.int_box_grid.AddGrowableCol(0)

    # Insert into sizers
    int_box_sizer.Add(self.int_box_grid, 1, flag=wx.ALL | wx.EXPAND, border=10)
    summary_sizer.Add(int_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    if self.gparams.advanced.processing_backend != 'ha14':
      summary_sizer.Hide(int_box_sizer, recursive=True)
      self.Layout()

    # Dataset Info
    dat_box = wx.StaticBox(self, label='Dataset Information')
    dat_box.SetFont(sfont)
    dat_box_sizer = wx.StaticBoxSizer(dat_box, wx.HORIZONTAL)
    self.dat_box_grid = wx.GridBagSizer(5, 20)

    # Buttons for res. histogram and beam xy plot
    hist_bmp = bitmaps.fetch_icon_bitmap('mimetypes', 'spreadsheet', size=32,
                                         scale=(24, 24))
    self.dat_reshist = ct.GradButton(self,
                                     bmp=hist_bmp,
                                     label='  Resolution Histogram',
                                     size=(250, -1))
    beamXY_bmp = bitmaps.fetch_custom_icon_bitmap('scatter_plot_24')
    self.dat_beamxy = ct.GradButton(self,
                                    bmp=beamXY_bmp,
                                    label='  Beam XY Plot', size=(250, -1))
    self.dat_beam3D = ct.GradButton(self,
                                    bmp=beamXY_bmp,
                                    label='  Beam XYZ Plot', size=(250, -1))
    self.dat_box_grid.Add(self.dat_reshist, pos=(0, 1), flag=wx.ALIGN_RIGHT)
    self.dat_box_grid.Add(self.dat_beamxy, pos=(1, 1), flag=wx.ALIGN_RIGHT)
    self.dat_box_grid.Add(self.dat_beam3D, pos=(2, 1), flag=wx.ALIGN_RIGHT)
    self.Bind(wx.EVT_BUTTON, self.onPlotBeamXY, self.dat_beamxy)
    self.Bind(wx.EVT_BUTTON, self.onPlotBeam3D, self.dat_beam3D)
    self.Bind(wx.EVT_BUTTON, self.onPlotResHist, self.dat_reshist)

    # Insert into sizers
    dat_box_sizer.Add(self.dat_box_grid, 1, flag=wx.ALL | wx.EXPAND, border=10)
    summary_sizer.Add(dat_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Clustering info
    self.cluster_panel = wx.Panel(self)
    cluster_box = wx.StaticBox(self.cluster_panel, label='Unit Cell Clustering')
    cluster_box_sizer = wx.StaticBoxSizer(cluster_box, wx.VERTICAL)
    self.cluster_panel.SetSizer(cluster_box_sizer)
    self.cluster_info = ct.CustomListCtrl(self.cluster_panel, size=(-1, 150))
    self.cluster_info.ctr.InsertColumn(0, "#")
    self.cluster_info.ctr.InsertColumn(1, "Lattice")
    self.cluster_info.ctr.InsertColumn(2, "Unit Cell")
    self.cluster_info.ctr.InsertColumn(3, "Filename")
    self.cluster_info.ctr.setResizeColumn(3)
    cluster_box_sizer.Add(self.cluster_info, proportion=1, flag=wx.EXPAND)
    summary_sizer.Add(self.cluster_panel, flag=wx.EXPAND | wx.ALL, border=10)

    # Hide if not done
    if not self.gparams.analysis.clustering.flag_on:
      self.cluster_panel.Hide()

    # Summary
    smr_box = wx.StaticBox(self, label='Run Summary')
    smr_box.SetFont(sfont)
    smr_box_sizer = wx.StaticBoxSizer(smr_box, wx.HORIZONTAL)
    self.smr_box_grid = wx.GridBagSizer(5, 20)

    prime_bmp = bitmaps.fetch_custom_icon_bitmap('prime32', scale=(24, 24))
    self.smr_runprime = ct.GradButton(self,
                                      bmp=prime_bmp,
                                      label='  Run PRIME', size=(250, -1))
    cluster_bmp = bitmaps.fetch_custom_icon_bitmap('distance_difference',
                                                   scale=(24, 24))
    self.smr_runcluster = ct.GradButton(self,
                                        bmp=cluster_bmp,
                                        label='  Run CLUSTER', size=(250, -1))

    self.smr_box_grid.Add(self.smr_runcluster, pos=(0, 1), flag=wx.ALIGN_RIGHT)
    self.smr_box_grid.Add(self.smr_runprime, pos=(1, 1), flag=wx.ALIGN_RIGHT)
    self.Bind(wx.EVT_BUTTON, self.onPRIME, self.smr_runprime)
    self.Bind(wx.EVT_BUTTON, self.onCLUSTER, self.smr_runcluster)

    smr_box_sizer.Add(self.smr_box_grid, 1, flag=wx.ALL | wx.EXPAND, border=10)
    summary_sizer.Add(smr_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    self.SetFont(sfont)
    self.SetSizer(summary_sizer)
    self.SetupScrolling()

  def onPRIME(self, e):
    from prime.postrefine.mod_gui_frames import PRIMEWindow

    self.prime_window = PRIMEWindow(None, -1, title='PRIME',
                                    prefix=self.gparams.advanced.prime_prefix)
    self.prime_window.load_script(out_dir=self.info.int_base)
    self.prime_window.place_and_size(set_size='v_default', set_by='mouse',
                                     center=True)

    os.chdir(self.info.int_base)
    self.prime_window.Show(True)

  def onCLUSTER(self, e):
    cluster_dlg = d.ClusterDialog(self)
    cluster_dlg.write_files.SetValue(
      self.gparams.analysis.clustering.write_files)
    cluster_dlg.cluster_threshold.ctr.SetValue(
      self.gparams.analysis.clustering.threshold)
    cluster_dlg.cluster_limit.ctr.SetValue(
      self.gparams.analysis.clustering.limit)
    if self.gparams.analysis.clustering.n_images > 0:
      cluster_dlg.cluster_n_images.ctr.SetValue(
        self.gparams.analysis.clustering.n_images)

    if (cluster_dlg.ShowModal() == wx.ID_OK):
      self.cluster_panel.Show()
      self.Layout()
      self.gparams.analysis.clustering.flag_on = True
      self.gparams.analysis.clustering.write_files = \
        cluster_dlg.write_files.GetValue()
      self.gparams.analysis.clustering.threshold = \
        cluster_dlg.cluster_threshold.ctr.GetValue()
      self.gparams.analysis.clustering.limit = \
        cluster_dlg.cluster_limit.ctr.GetValue()
      if cluster_dlg.cluster_n_images.toggle.GetValue():
        self.gparams.analysis.clustering.n_images = int(
          cluster_dlg.cluster_n_images.ctr.GetValue())
      else:
        self.gparams.analysis.clustering.n_images = 0

      analysis = Analyzer(info=self.info, params=self.gparams, gui_mode=True)
      clusters = analysis.unit_cell_analysis()
      if clusters:
        self.report_clustering_results(clusters=clusters)

  def report_clustering_results(self, clusters):
    self.cluster_info.ctr.DeleteAllItems()
    clusters = sorted(clusters, key=lambda i: i['number'], reverse=True)
    for c in clusters:
      i = clusters.index(c)
      idx = self.cluster_info.ctr.InsertStringItem(i, str(c['number']))
      self.cluster_info.ctr.SetStringItem(idx, 1, str(c['pg']))
      self.cluster_info.ctr.SetStringItem(idx, 2, str(c['uc']))
      if 'filename' in c and c['filename'] not in ('*', None):
        self.cluster_info.ctr.SetStringItem(idx, 3, c['filename'])
    self.Refresh()
    self.SetupScrolling()

  def update(self):

    # Add title and directory information
    description = self.info.description if hasattr(self.info, 'description') \
      else self.gparams.description
    self.title_txt.SetLabel(description)
    int_base = self.info.int_base if hasattr(self.info, 'int_base') else \
      self.gparams.output
    self.folder_txt.SetLabel(int_base)

    # # Add grid search stats if exist
    # if hasattr(self.info, 'gs_stats'):
    #   gs_stats = self.info.gs_stats
    #   rkeys = ['s', 'h', 'a']
    #   clabels = ['min', 'max', 'avg', 'std']
    #   rlabels = [gs_stats[k]['label'] for k in rkeys]
    #   contents = [[gs_stats[k]['min'], gs_stats[k]['max'], gs_stats[k]['mean'],
    #                gs_stats[k]['std']] for k in rkeys]
    #   gs_table = ct.TableCtrl(self,
    #                           clabels=clabels,
    #                           rlabels=rlabels,
    #                           contents=contents)
    #   self.int_box_grid.Add(gs_table, span=(len(rlabels), 1),
    #                         pos=(0, 0), flag=wx.EXPAND)

    # Add dataset information
    if hasattr(self.info, 'stats'):
      lres = self.info.stats['lres']['mean']
      hres = self.info.stats['res']['mean']
      res = to_unicode(u'{:.2f} - {:.2f} {}'.format(lres, hres, u'\u212B'))

      if type(self.info.best_uc) in (list, tuple):
        a, b, c, alpha, beta, gamma = self.info.best_uc
        uc = '{:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}' \
             ''.format(a, b, c, alpha, beta, gamma)
      else:
        uc = str(self.info.best_uc)

      beamxy = 'X = {:.2f}, Y = {:.2f}'.format(self.info.stats['beamX']['mean'],
                                               self.info.stats['beamY']['mean'])
      data = zip(
        ['Bravais lattice: ', 'Unit cell: ','Resolution: ', 'Beam XY (''mm): '],
        [self.info.best_pg, uc, res, beamxy]
      )

      # Make dataset stat plot
      stat_plot = Plotter(self, info=self.info)
      stat_plot.initialize_figure(figsize=(0.1, 0.1))
      self.dat_box_grid.Add(stat_plot, span=(4, 1), pos=(0, 0), flag=wx.EXPAND)
      stat_plot.plot_table(data=data)
      self.dat_box_grid.AddGrowableCol(0)
      self.dat_box_grid.AddGrowableRow(3)

    # Add processing summary
    if hasattr(self.info, 'categories'):
      ckeys = ['total', 'have_diffraction', 'failed_triage',
               'failed_spotfinding', 'failed_indexing', 'failed_grid_search',
               'failed_integration', 'failed_filter', 'integrated']
      rlabels = []
      contents = []
      for k in ckeys:
        if self.info.categories[k][0]:
          contents.append(str(len(self.info.categories[k][0])))
          rlabels.append(self.info.categories[k][1])
      proc_data = zip(rlabels, contents)

      # Make plot
      proc_plot = Plotter(self, info=self.info)
      proc_plot.initialize_figure(figsize=(0.1, 0.1))
      self.smr_box_grid.Add(proc_plot, span=(5, 1), pos=(0, 0), flag=wx.EXPAND)
      proc_plot.plot_table(data=proc_data)
      self.smr_box_grid.AddGrowableCol(0)
      self.smr_box_grid.AddGrowableRow(4)

    # Add clustering info
    if self.info.clusters:
      self.report_clustering_results(clusters=self.info.clusters)

    self.Layout()
    self.SetupScrolling()

  def initialize_standalone_plot(self, figsize=(8, 8)):
    self.plot_window = PlotWindow(self, -1, title='IOTA Plot')
    self.plot = Plotter(self.plot_window, params=self.gparams, info=self.info)
    self.plot_window.plot_panel = self.plot
    self.plot.initialize_figure(figsize=figsize)

  def show_plot(self):
    self.plot_window.add_plot_to_window()
    self.plot_window.place_and_size(set_by='parent', set_size=True,
                                    position=(25, 25))
    self.plot_window.Show()

  def onPlotHeatmap(self, e):
    if self.info.final_objects is not None:
      self.initialize_standalone_plot()
      self.plot.plot_spotfinding_heatmap()
      self.show_plot()

  def onPlotBeamXY(self, e):
    if self.info.final_objects is not None:
      self.initialize_standalone_plot()
      self.plot.plot_beam_xy()
      self.show_plot()

  def onPlotBeam3D(self, e):
    if self.info.final_objects is not None:
      self.initialize_standalone_plot()
      self.plot.plot_beam_xy(threeD=True)
      self.show_plot()

  def onPlotResHist(self, e):
    if self.info.final_objects is not None:
      self.initialize_standalone_plot()
      self.plot.plot_res_histogram()
      self.show_plot()


class ProcWindow(IOTABaseFrame):
  """ New frame that will show processing info """

  def __init__(self, parent, id, title, phil):
    IOTABaseFrame.__init__(self, parent, id, title,
                           size=(800, 900),
                           style=wx.SYSTEM_MENU | wx.CAPTION |
                                 wx.CLOSE_BOX | wx.RESIZE_BORDER)

    self.parent = parent
    self.info = None
    self.bookmark = 0
    self.gparams = phil.extract()

    self.state = 'process'
    self.recovery = False

    self.abort_initiated = False
    self.run_aborted = False
    self.monitor_mode = False
    self.monitor_mode_timeout = None
    self.timeout_start = None
    self.find_new_images = False
    self.finding_images = False
    self.start_object_finder = True
    self.plotter_time = []
    self.obj_reader_time = []
    self.message = 'Running...'

    self.run_cluster = False
    self.running_cluster = False
    self.run_prime = False
    self.running_prime = False
    self.running_manual_analysis = False
    self.job_id = None

    self.finished_objects = []
    self.read_object_files = []
    self.new_images = []
    self.indices = []
    self.b_factors = []
    self.clusters = None
    self.pparams = None
    self.prime_info = None

    # Initialize mx handler
    self.mxh = mx_handler()

    # Toolbar
    self.initialize_toolbar()
    self.tb_btn_abort = self.add_tool(label='Abort',
                                      bitmap=('actions', 'stop'),
                                      shortHelp='Quit')

    resume_bmp = bitmaps.fetch_icon_bitmap('actions', 'quick_restart')
    self.tb_btn_resume = self.add_tool(label='Resume',
                                       bitmap=('actions', 'quick_restart'),
                                       shortHelp='Resume aborted run')
    self.add_toolbar_separator()
    watch_bmp = bitmaps.fetch_icon_bitmap('apps', 'search')
    self.tb_btn_monitor = self.add_tool(label='Monitor',
                                        kind=wx.ITEM_CHECK,
                                        bitmap=('apps', 'search'),
                                        shortHelp='Monitor Mode')
    self.tb_btn_analysis = self.add_tool(label='Analysis',
                                         kind=wx.ITEM_CHECK,
                                         bitmap=(
                                         'mimetypes', 'spreadsheet', 32),
                                         shortHelp='Toggle Runtime Analysis Tab')
    self.set_tool_state(self.tb_btn_resume, False)
    self.realize_toolbar()

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

    # Tabbed output window(s)
    self.proc_panel = wx.Panel(self)
    self.proc_nb = AuiNotebook(self.proc_panel, style=wx.aui.AUI_NB_TOP)
    self.proc_tab = ProcessingTab(self.proc_nb, gparams=self.gparams)
    self.log_tab = LogTab(self.proc_nb)
    self.chart_tab = LiveAnalysisTab(self.proc_nb)
    self.proc_nb.AddPage(self.log_tab, 'Log')
    self.proc_nb.AddPage(self.proc_tab, 'Processing')
    self.proc_nb.AddPage(self.chart_tab, 'Analysis')
    self.proc_nb.RemovePage(2)
    self.proc_nb.SetSelection(1)
    self.proc_sizer = wx.BoxSizer(wx.VERTICAL)
    self.proc_sizer.Add(self.proc_nb, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.proc_panel.SetSizer(self.proc_sizer)

    self.main_sizer.Add(self.status_panel, flag=wx.EXPAND | wx.ALL, border=3)
    self.main_sizer.Add(self.proc_panel, 1, flag=wx.EXPAND | wx.ALL, border=3)

    # Processing status bar
    self.sb.SetFieldsCount(2)
    self.sb.SetStatusWidths([-1, -2])

    # Output polling timer
    self.proc_timer = wx.Timer(self)
    self.chart_timer = wx.Timer(self)

    # PostEvent bindings
    self.Bind(thr.EVT_ALLDONE, self.onFinishedProcess)
    self.Bind(thr.EVT_IMGDONE, self.onFinishedImageFinder)
    self.Bind(thr.EVT_OBJDONE, self.onFinishedObjectReader)
    self.Bind(thr.EVT_CLUSTERDONE, self.onFinishedCluster)
    self.Bind(thr.EVT_PRIMEDONE, self.onFinishedPRIME)

    # Event bindings
    self.Bind(wx.EVT_TIMER, self.onProcTimer, id=self.proc_timer.GetId())
    self.Bind(wx.EVT_TIMER, self.onChartTimer, id=self.chart_timer.GetId())

    # Button bindings
    self.Bind(wx.EVT_TOOL, self.onAbort, self.tb_btn_abort)
    self.Bind(wx.EVT_TOOL, self.onResume, self.tb_btn_resume)
    self.Bind(wx.EVT_TOOL, self.onMonitor, self.tb_btn_monitor)
    self.Bind(wx.EVT_TOOL, self.onAnalysis, self.tb_btn_analysis)
    self.Bind(wx.EVT_BUTTON, self.onAnalysisManualRun,
              self.chart_tab.btn_run_analysis)

    # Determine if monitor mode was previously selected
    if self.gparams.gui.monitor_mode:
      self.toolbar.ToggleTool(self.tb_btn_monitor.GetId(), True)
      self.monitor_mode = True
      self.find_new_images = True
      if self.gparams.gui.monitor_mode_timeout:
        if self.gparams.gui.monitor_mode_timeout_length is None:
          self.monitor_mode_timeout = 30
        else:
          self.monitor_mode_timeout = self.gparams.gui.monitor_mode_timeout_length

  def onAnalysisManualRun(self, e):
    # Run analysis calculations when Run Analysis button is pressed; this
    # will enable analysis run even when IOTA isn't running anymore
    if not (self.running_cluster or self.running_prime):
      self.running_manual_analysis = True
      self.chart_tab.btn_run_analysis.Disable()
      self.run_clustering_thread()

  def onAnalysis(self, e):
    if self.toolbar.GetToolState(self.tb_btn_analysis.GetId()):

      # Start timer only if the proc timer is running
      if self.proc_timer.IsRunning():
        self.chart_timer.Start(15000)

      # Insert Analysis page
      if wx.__version__[0] == '4':
        self.proc_nb.InsertPage(page_idx=2, page=self.chart_tab,
                                caption='Analysis', select=True)
      else:
        self.proc_nb.InsertPage(n=2, page=self.chart_tab,
                                text='Analysis')
      self.proc_nb.SetSelection(2)
      self.plot_live_analysis(force_plot=True)
      self.chart_tab.Show()   # Need to show when re-opening page (see below)
    else:

      # Stop chart timer
      if self.chart_timer.IsRunning():
        self.chart_timer.Stop()

      # Move selection to different page
      if self.proc_nb.GetSelection() == 2:
        self.proc_nb.SetSelection(1)

      # Remove Analysis tab; since removal of the tab doesn't hide the
      # contents, have to do that first. Cannot use .DeletePage() because I
      # then will have to re-create the entire chart_tab.
      self.chart_tab.Hide()
      self.proc_nb.RemovePage(2)

    self.proc_nb.Refresh()


  def onMonitor(self, e):
    if self.toolbar.GetToolState(self.tb_btn_monitor.GetId()):
      self.monitor_mode = True
      if self.gparams.gui.monitor_mode_timeout:
        if self.gparams.gui.monitor_mode_timeout_length is None:
          self.monitor_mode_timeout = 30
        else:
          self.monitor_mode_timeout = self.gparams.gui.monitor_mode_timeout_length
    elif not self.toolbar.GetToolState(self.tb_btn_monitor.GetId()):
      self.monitor_mode = False
      self.monitor_mode_timeout = None
    self.find_new_images = self.monitor_mode

  def onAbort(self, e):
    if self.gparams.mp.method == 'lsf':
      kill_command = 'bkill -J {}'.format(self.job_id)
      easy_run.fully_buffered(kill_command)
    self.abort_initiated = True
    self.set_tool_state(self.tb_btn_abort, False)

  def onResume(self, e):
    """ Restarts an aborted run if the processing window is still open.
    Basically goes through self.finished_objects, extracts the raw image
    names and regenerates the self.img_list to only have those image paths;
    then finds any 'new' images (which includes unprocessed images as well as
    any images that may have been added during the abort pause) and runs
    processing """

    # Remove abort signal file(s)
    if os.path.isfile(self.tmp_abort_file):
      os.remove(self.tmp_abort_file)
    if os.path.isfile(self.tmp_aborted_file):
      os.remove(self.tmp_aborted_file)
    self.abort_initiated = False
    self.run_aborted = False

    # Reset tmp folder (if necessary)
    if not os.path.isdir(self.info.tmp_base):
      os.makedirs(self.info.tmp_base)

    # Reconstitute object list (if necessary)
    with open(self.info.obj_list_file, 'w') as olf:
      if not os.path.isfile(self.info.obj_list_file):
        finished_objects = self.info.get_finished_objects()
        for obj in finished_objects:
          olf.write('{}\n'.format(obj.obj_file))
      olf.seek(0, 2)
      self.info.bookmark = olf.tell()

    # Reset toolbar buttons
    self.set_tool_states([(self.tb_btn_abort, True),
                          (self.tb_btn_resume, False),
                          (self.tb_btn_monitor, True)])

    # Determine state
    if len(self.info.categories['not_processed']) <= 0:
      self.state = 'finished'
    else:
      self.info.reset_input_list()
      self.state = 'resume'

    # Set status font properties
    font = self.sb.GetFont()
    font.SetWeight(wx.FONTWEIGHT_NORMAL)
    self.status_txt.SetFont(font)
    self.status_txt.SetForegroundColour('black')
    self.status_txt.SetLabel(
      'Resuming run #{} ...'.format(self.info.run_number))

    # Run processing, etc.
    self.proc_timer.Start(3000)
    self.process_images()

  def recover(self, int_path, status, info, params):
    self.recovery = True
    self.gparams = params
    self.tmp_abort_file = os.path.join(int_path, '.abort.tmp')
    self.tmp_aborted_file = os.path.join(int_path, '.aborted.tmp')

    # Get info and apply necessary type conversions
    self.info = info

    # unicode-to-str in cluster iterable
    for item in self.info.cluster_iterable:
      item[6] = str(item[6])

    self.status_txt.SetLabel('Searching in {} ...'.format(int_path))
    self.state = status
    self.finish_process()

  def start_processing(self):
    ''' Run processing '''
    self.tmp_abort_file = os.path.join(self.info.int_base, '.abort.tmp')
    self.tmp_aborted_file = os.path.join(self.info.int_base, '.aborted.tmp')
    self.status_txt.SetForegroundColour('black')

    self.state = 'start'
    self.process_images()
    self.proc_timer.Start(3000)

  def process_images(self):
    """ Initiate processing in a new or resumed run """

    # Finish up if resuming a run that was already finished
    if self.state == 'finished':
      self.finish_process()

    # Instantiate and start submit thread
    self.proc_thread = thr.JobSubmitThread(self, params=self.gparams)
    self.proc_thread.name = 'IOTAJobSubmitThread'
    self.proc_thread.start()

  def display_log(self):
    if os.path.isfile(self.info.logfile):
      with open(self.info.logfile, 'r') as out:
        out.seek(self.bookmark)
        output = out.read()
        self.bookmark = out.tell()

      ins_pt = self.log_tab.log_window.GetInsertionPoint()
      self.log_tab.log_window.AppendText(output)
      self.log_tab.log_window.SetInsertionPoint(ins_pt)

  def plot_integration(self, force_plot=False):
    """ This function will plot fast-drawing runtime processing charts on the
        "Processing" tab """
    sw = wx.StopWatch()
    if self.proc_nb.GetSelection() == 1 or force_plot:
      if hasattr(self.info, 'unplotted_stats'):
        s = time.time()
        self.proc_tab.info = self.info
        self.proc_tab.draw_plots()
        self.info.unplotted_stats = {}
        self.proc_tab.draw_summary()
      self.plotter_time.append(sw.Time())

  def plot_live_analysis(self, force_plot=False):
    """ This function will plot in-depth analysis that will (importantly)
        involve expensive and slow-drawing charts on the Live Analysis tab """

    if self.proc_nb.GetSelection() == 2 or force_plot:
      self.chart_tab.info = self.info
      self.chart_tab.draw_plots()

  def onChartTimer(self, e):
    self.plot_live_analysis()
    if not (self.running_cluster or self.running_prime):
      self.run_clustering_thread()

  def run_clustering_thread(self):
    self.chart_tab.btn_run_analysis.Disable()
    self.running_cluster = True
    self.cluster_thread = thr.ClusterThread(self,
                                            iterable=self.info.cluster_iterable)
    self.cluster_thread.start()

  def onFinishedCluster(self, e):
    self.info.clusters, errors = e.GetValue()
    self.running_cluster = False

    if self.info.clusters:
      self.info.best_pg = self.info.clusters[0]['pg']
      self.info.best_uc = self.info.clusters[0]['uc']

    if errors:
      for err in errors:
        print ('IOTA ERROR (CLUSTERING): ', err)
      self.info.errors.extend(errors)

    # Output cluster results
    ep.dump(self.info.cluster_info_file, obj=self.clusters)

    if not (hasattr(self, 'prime_thread') and self.prime_thread.isAlive()):
      self.pparams = None
      self.run_prime_thread()

  def run_prime_thread(self):
    # Run PRIME (basic merge only)
    self.chart_tab.btn_run_analysis.Disable()
    self.running_prime = True
    pg = ut.makenone(self.chart_tab.pg_uc.pg.GetValue())
    uc = ut.makenone(self.chart_tab.pg_uc.uc.GetValue())
    self.prime_thread = thr.PRIMEThread(self, self.info, self.gparams, pg, uc)
    self.prime_thread.start()

  def onFinishedPRIME(self, e):
    self.prime_info = e.GetValue()
    self.running_prime = False
    if self.running_manual_analysis:
      self.info.prime_info = self.prime_info
      self.plot_live_analysis(force_plot=True)
      self.running_manual_analysis = False
    self.chart_tab.btn_run_analysis.Enable()

  def adjust_timer(self, new_interval):
    self.proc_timer.Stop()
    self.proc_timer.Start(new_interval)
    print("IOTA: TIMER SET TO {} SECONDS".format(new_interval / 1000))

  def onProcTimer(self, e):

    # Check if processing has been initialized
    if not self.info.init_proc:
      self.status_txt.SetLabel('Initializing IOTA run ...')
      self.info = ProcInfo.from_json(self.info.info_file)
      return

    # Catch mean of empty slice warning (for cosmetic reasons)
    with warnings.catch_warnings():
      warnings.simplefilter("ignore", category=RuntimeWarning)
      sw_means = np.mean(self.plotter_time) + np.mean(self.obj_reader_time)
    timer_adjustment = np.ceil(sw_means * 5 / 1000) * 1000
    if not (
            np.isnan(timer_adjustment) or
            timer_adjustment == self.proc_timer.GetInterval()
    ):
      self.adjust_timer(new_interval=timer_adjustment)

    # Plot integration and display lob
    self.info.export_json()
    self.plot_integration()
    self.display_log()

    if self.abort_initiated:
      self.status_txt.SetForegroundColour('red')
      self.status_txt.SetLabel('Aborting...')
      self.run_aborted = not (hasattr(self, 'proc_thread') and
                              self.proc_thread.is_alive())
      if not self.run_aborted:
        self.proc_thread.abort()
        if self.running_cluster:
          self.cluster_thread.abort()
        if self.running_prime:
          self.prime_thread.abort()

      if self.run_aborted:
        self.status_txt.SetLabel('ABORTED BY USER!')
        print("IOTA: RUN ABORTED!")
        self.finish_process()
    else:
      self.run_aborted = False

      # If all current images processed, check for new ones and/or finish process
      if len(self.info.categories['not_processed'][0]) <= 0:
        self.monitor_filesystem()
      else:
        prcd = len(self.info.categories['total'][0]) - \
               len(self.info.categories['not_processed'][0])
        intd = len(self.info.categories['integrated'][0])
        update = '- {} processed ({} integrated)'.format(prcd, intd)
        self.status_txt.SetLabel('Processing {} images {}'
                                 ''.format(self.info.n_input_images, update))

        # Instantiate and start processing thread
        if not (hasattr(self,
                        'object_reader') and self.object_reader.is_alive()):
          self.obj_sw = wx.StopWatch()
          self.object_reader = thr.ObjectReaderThread(self, info=self.info)
          self.object_reader.name = 'IOTAObjectReader'
          self.object_reader.start()

  def monitor_filesystem(self):
    # Check if all images have been looked at; if yes, finish process
    if self.monitor_mode:
      if self.new_images:
        self.status_txt.SetLabel(
          'Found {} new images'.format(len(self.new_images)))
        self.timeout_start = None
        self.state = 'new images'
        self.info.update_input_list(new_input=self.new_images)
        self.new_images = []
        self.info.export_json()
        self.process_images()
      else:
        if self.monitor_mode_timeout:
          if self.timeout_start is None:
            self.timeout_start = time.time()
          else:
            interval = time.time() - self.timeout_start
            if interval >= self.monitor_mode_timeout:
              self.status_txt.SetLabel('Timed out. Finishing...')
              self.finish_process()
            else:
              timeout_msg = 'No new images found! Timing out in {} seconds' \
                            ''.format(
                int(self.monitor_mode_timeout - interval))
              self.status_txt.SetLabel(timeout_msg)
        else:
          self.status_txt.SetLabel('No new images found! Waiting ...')

        if self.find_new_images:
          self.get_images_from_filesystem()

    else:
      self.status_txt.SetLabel('Wrapping up ...')
      # interrupt long-running PRIME thread
      try:
        if self.running_cluster:
          self.cluster_thread.abort()
        if self.running_prime:
          self.prime_thread.abort()
      except Exception:
        pass
      self.finish_process()

  def get_images_from_filesystem(self):
    self.find_new_images = False
    img_finder = thr.ImageFinderThread(
      self,
      input=self.gparams.input,
      input_list=self.info.categories['total'][0])
    img_finder.start()

  def onFinishedProcess(self, e):
    self.finding_images = False
    if not self.monitor_mode:
      self.finish_process()

  def onFinishedImageFinder(self, e):
    self.new_images = e.GetValue()
    self.find_new_images = True

  def onFinishedObjectReader(self, info):
    # Read info object
    self.info = info.GetValue()

    # Update w/ latest Cluster and PRIME results
    self.info.prime_info = self.prime_info
    self.info.cluster_info = self.clusters

    if not self.recovery:
      self.obj_reader_time.append(self.obj_sw.Time())

  def finish_process(self):
    # Stop timer
    self.proc_timer.Stop()
    self.chart_timer.Stop()

    # Aborted run
    if self.run_aborted:
      with open(self.tmp_aborted_file, 'w') as tf:
        tf.write('')
      self.info.status = 'aborted'
      end_color = 'red'
      end_msg = 'ABORTED BY USER'
      self.set_tool_state(self.tb_btn_resume, True)

    # Recovered run display
    elif self.recovery:
      end_color = 'black'
      end_msg = 'Run #{} Loaded!'.format(self.info.run_number)
      self.set_tool_states([(self.tb_btn_abort, False),
                            (self.tb_btn_resume, True),
                            (self.tb_btn_monitor, False, False)])

    # Normal finish
    else:
      self.info.status = 'finished'
      self.set_tool_states([(self.tb_btn_abort, False),
                            (self.tb_btn_monitor, False, False)])
      successfully_processed = len(self.info.categories['integrated'][0])
      end_color = 'blue'
      if successfully_processed > 0:
        from iota.components.iota_analysis import Analyzer
        analyzer = Analyzer(info=self.info, params=self.gparams, gui_mode=True)
        self.info = analyzer.run_all(get_results=False)
        end_msg = 'DONE'
      else:
        end_msg = 'NO IMAGES PROCESSED'

    # Final analysis and display summary
    if (
            not self.run_aborted and
            len(self.info.categories['integrated'][0]) > 0 and
            len(self.info.categories['not_processed'][0]) == 0
    ):
      self.set_tool_state(tool=self.tb_btn_resume, enable=False)
      self.summary_tab = SummaryTab(self.proc_nb, self.info, self.gparams)
      self.proc_nb.AddPage(self.summary_tab, 'Summary', select=True)
      self.proc_nb.SetSelection(self.proc_nb.GetPageCount())
      self.summary_tab.update()

      import shutil
      try:
        shutil.rmtree(self.info.tmp_base)
        shutil.rmtree(self.info.dials_log_base)
      except Exception:
        pass

    # Export final info file
    self.info.export_json()

    # Signal end of run
    font = self.sb.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    self.status_txt.SetFont(font)
    self.status_txt.SetForegroundColour(end_color)
    self.status_txt.SetLabel(end_msg)

    # Finish up
    self.plot_integration(force_plot=True)
    self.plot_live_analysis(force_plot=True)

    # Report run time
    start_time = getattr(self.parent, 'start_time', None)
    if start_time:
      import datetime
      runtime = datetime.timedelta(seconds=int(time.time() - start_time))
      hours, minutes, seconds = str(runtime).split(':')
      ut.main_log(self.info.logfile,
                  "Total run time: {} hours, {} minutes, {} seconds"
                  "".format(hours, minutes, seconds),
                  print_tag=True)
    self.display_log()

# -- end
