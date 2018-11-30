from __future__ import division, print_function, absolute_import
from past.builtins import range

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 11/29/2018
Description : IOTA GUI Initialization module
'''


import os
import wx

import wx.lib.agw.ultimatelistctrl as ulc
import multiprocessing
import argparse

from iotbx import phil as ip
from libtbx import easy_pickle as ep
from cctbx import miller
assert miller

from iota import iota_version, gui_description, gui_license
import iota.components.iota_input as inp
import iota.components.iota_utils as util
import iota.components.iota_ui_frames as frm
import iota.components.iota_ui_dialogs as dlg
from iota.components.iota_ui_base import IOTABaseFrame

ginp = util.InputFinder()
pid = os.getpid()

try:
  user = os.getlogin()
except OSError:
  user = 'iota'

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
  python = "Python"    #TODO: make sure it's right!


# --------------------------- Command-line Parser ---------------------------- #

def parse_command_args(help_message):
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog = 'iota',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=(help_message),
            epilog=('\n{:-^70}\n'.format('')))
  parser.add_argument('path', type=str, nargs = '*', default = None,
            help = 'Path to data or file with IOTA parameters')
  parser.add_argument('--version', action = 'version',
            version = 'IOTA {}'.format(iota_version),
            help = 'Prints version info of IOTA')
  parser.add_argument('-w', type=int, nargs=1, default=0, dest='watch',
            help = 'Run IOTA in watch mode - check for new images')
  parser.add_argument('-random', type=int, nargs=1, default=0,
            help = 'Size of randomized subset, e.g. "--random 10"')
  parser.add_argument('--range', type=str, nargs='?', default=None,
            help = 'Range of images, e.g."--range 1-5,25,200-250"')
  parser.add_argument('-n', type=int, nargs='?', default=0, dest='nproc',
            help = 'Specify a number of cores for a multiprocessor run"')
  parser.add_argument('--tmp', type=str, nargs = 1, default = None,
            help = 'Path to temp folder')
  return parser


# ------------------------------- Main Window -------------------------------- #

class MainWindow(IOTABaseFrame):
  """ Frame housing the entire app; all windows open from this one """

  def __init__(self, parent, id, title):
    IOTABaseFrame.__init__(self, parent, id, title, size=(800, 500))
    self.parent = parent

    self.prefs_phil = None
    self.target_phil = None
    self.args = None
    self.phil_args = None

    # Menu bar
    menubar = wx.MenuBar()

    # Status bar
    self.sb = self.CreateStatusBar()

    # Help menu item with the about dialog
    m_help = wx.Menu()
    m_file = wx.Menu()
    self.mb_load_script = m_file.Append(wx.ID_OPEN, '&Load Script...')
    self.mb_save_script = m_file.Append(wx.ID_SAVE, '&Save Script...')
    m_file.AppendSeparator()
    self.mb_reset = m_file.Append(wx.ID_ANY, '&Reset Settings')
    self.mb_about = m_help.Append(wx.ID_ANY, '&About')
    menubar.Append(m_file, '&File')
    menubar.Append(m_help, '&Help')

    self.SetMenuBar(menubar)

    # Toolbar
    self.initialize_toolbar()
    self.tb_btn_quit = self.add_tool(label='Quit',
                                     bitmap=('actions', 'exit'),
                                     shortHelp='Quit')

    self.tb_btn_prefs = self.add_tool(label='Preferences',
                                      bitmap=('apps', 'advancedsettings'),
                                      shortHelp='Preferences')
    self.toolbar.AddSeparator()
    self.tb_btn_load = self.add_tool(label='Load Script',
                                     bitmap=('actions', 'download'),
                                     shortHelp='Load Script')
    self.tb_btn_save = self.add_tool(label='Save Script',
                                     bitmap=('actions', 'save_all'),
                                     shortHelp='Save Script')
    self.tb_btn_reset = self.add_tool(label='Reset',
                                      bitmap=('actions', 'reload'),
                                      shortHelp='Reset Settings')
    self.toolbar.AddSeparator()
    self.tb_btn_analysis = self.add_tool(label='Recover',
                                         bitmap=('actions', 'list'),
                                         shortHelp='Recover')
    self.tb_btn_run = self.add_tool(label='Run',
                                    bitmap=('actions', 'run'),
                                    shortHelp='Run')

    # These buttons will be disabled until input path is provided
    self.set_tool_state(self.tb_btn_run, False)
    self.realize_toolbar()

    # Instantiate settings
    self.on_startup()

    # Instantiate input windows
    self.input_window = frm.InputWindow(self, phil=self.iota_phil)
    self.main_sizer.Add(self.input_window, 1,
                        flag=wx.ALL | wx.EXPAND,
                        border=10)
    self.main_sizer.Add((-1, 20))
    self.update_input_window()

    # button bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onPreferences, self.tb_btn_prefs)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_run)
    self.Bind(wx.EVT_TOOL, self.onRecovery, self.tb_btn_analysis)
    self.Bind(wx.EVT_TOOL, self.onLoadScript, self.tb_btn_load)
    self.Bind(wx.EVT_TOOL, self.onOutputScript, self.tb_btn_save)
    self.Bind(wx.EVT_TOOL, self.onReset, self.tb_btn_reset)

    # Menubar button bindings
    self.Bind(wx.EVT_MENU, self.OnAboutBox, self.mb_about)
    self.Bind(wx.EVT_MENU, self.onOutputScript, self.mb_save_script)
    self.Bind(wx.EVT_MENU, self.onLoadScript, self.mb_load_script)
    self.Bind(wx.EVT_MENU, self.onReset, self.mb_reset)

    # Bindings to Input Window
    self.Bind(wx.EVT_BUTTON, self.onImportOptions,
              self.input_window.opt_btn_import)
    self.Bind(wx.EVT_BUTTON, self.onProcessOptions,
              self.input_window.opt_btn_process)
    self.Bind(wx.EVT_BUTTON, self.onAnalysisOptions,
              self.input_window.opt_btn_analysis)

    # File list control bindings
    self.Bind(ulc.EVT_LIST_INSERT_ITEM, self.onItemInserted,
              self.input_window.input)

    from platform import python_version
    print('Python  : ', python_version())
    print('wxPython: ', wx.__version__)

  def on_startup(self):
    ''' Read command-line args, create PHIL, apply command-line args
    :return: inputs: list of input entries to update input window
    '''

    help_msg = '''This command will run the IOTA GUI '''
    self.args, self.phil_args = parse_command_args(help_msg).parse_known_args()

    # Understand path(s)
    inputs = []
    messages = []
    paramfile = None
    if self.args.path and len(self.args.path) > 0:
      for carg in self.args.path:
        msg = None
        if os.path.isfile(carg):
          ptype = ginp.get_file_type(carg)
          if 'image' in ptype:
            inputs.append(carg)
          elif 'IOTA settings' in ptype:
            paramfile = carg
          else:
            if 'settings' in ptype:
              msg = '- Cannot recognize settings file {}'.format(carg)
            elif 'file' in ptype:
              msg = '- Format for {} not recognized'.format(carg)
        elif os.path.isdir(carg):
          ptype = ginp.get_folder_type(carg)
          if 'image' in ptype:
            inputs.append(carg)
          else:  #TODO: empty folder check
            msg = '- Invalid input in {} ({})'.format(carg, ptype)
        else:
          msg = '- Error accessing {}: no such file or directory'.format(carg)
        if msg:
          messages.append(msg)

    # Load script or generate IOTA PHIL:
    if paramfile:
      self.load_script(filepath=paramfile, reset=False,
                       update_input_window=False)
    else:
      self.iota_phil = inp.process_ui_input(args=self.args,
                                            phil_args=self.phil_args,
                                            input_source=None)
      self.gparams = self.iota_phil.extract()

    # Add any inputs found in command-line args
    # TODO: add better handling of wildcards
    if inputs:
      for path in inputs:
        self.gparams.input.append(os.path.abspath(path))

    # Report errors if any (can get REALLY annoying if a lot of text)
    if messages:
      msg_string = 'IOTA Reported the following error(s):{}\n\n' \
                   ''.format('\n'.join(messages))
      wx.MessageBox(caption='Command-line Argument Error!',
                    message=msg_string, style=wx.OK | wx.ICON_ERROR)

  def onItemInserted(self, e):
    print (self.input_window.input.all_data_images)

  def onReset(self, e):
    self.reset_settings()

  def onPreferences(self, e):
    """ Opens dialog for IOTA preferences
    :param e: event object for self.tb_btn_prefs
    :return: modifies self.iota_phil with updated parameters
    """
    prefs = dlg.IOTAPreferences(self, phil=self.iota_phil)
    prefs.set_choices()

    if prefs.ShowModal() == wx.ID_OK:
      full_phil = self.iota_phil.fetch(source=prefs.prefs_phil)
      fixer = inp.PHILFixer()  # Backend switching requires the PHIL fixer
      self.iota_phil = fixer.run(old_phil=full_phil)
    prefs.Destroy()

    self.input_window.input_phil = self.iota_phil
    self.gparams = self.iota_phil.extract()
    self.input_window.gparams = self.gparams


  def onImportOptions(self, e):
    """ Opens dialog for image import options
    :param e: event object for self.input_window.opt_btn_import
    :return: modifies self.iota_phil with updated parameters
    """
    if self.gparams.advanced.processing_backend == 'cctbx.xfel':
      imp_dialog = dlg.ImportWindow(self,
                                    phil=self.iota_phil,
                                    title='Import Options',
                                    style=wx.DEFAULT_DIALOG_STYLE |
                                          wx.STAY_ON_TOP)
      imp_dialog.Fit()
    elif self.gparams.advanced.processing_backend == 'ha14':
      imp_dialog = dlg.HA14ImportWindow(self,
                                        phil=self.iota_phil,
                                        title='Import Options',
                                        style=wx.DEFAULT_DIALOG_STYLE |
                                              wx.STAY_ON_TOP)
      imp_dialog.Fit()
    else:
      imp_dialog = None

    if (imp_dialog and imp_dialog.ShowModal() == wx.ID_OK):
      self.iota_phil = self.iota_phil.fetch(source=imp_dialog.import_phil)
    imp_dialog.Destroy()

  def onProcessOptions(self, e):
    """ Opens dialog for image processing options, either for cctbx or DIALS
    depending on user selection.
    :param e: event object for self.input_window.opt_btn_process
    :return: modifies self.iota_phil with updated parameters
    """

    # For current cctbx.xfel options
    if self.gparams.advanced.processing_backend == 'cctbx.xfel':
      int_dialog = dlg.BackendOptions(self,
                                      phil=self.iota_phil,
                                      target=self.target_phil,
                                      title='cctbx.xfel Options',
                                      style=wx.DEFAULT_DIALOG_STYLE |
                                          wx.STAY_ON_TOP | wx.RESIZE_BORDER)
      int_dialog.SetMinSize((600, -1))
      int_dialog.Fit()

    # For deprecated cctbx.xfel HA14 options
    elif self.gparams.advanced.processing_backend == 'ha14':
      int_dialog = dlg.HA14BackendOptions(self,
                                          phil=self.iota_phil,
                                          target=self.target_phil,
                                          title='cctbx.xfel HA14 Options',
                                          style=wx.DEFAULT_DIALOG_STYLE |
                                          wx.STAY_ON_TOP | wx.RESIZE_BORDER)
      int_dialog.SetMinSize((600, -1))
      int_dialog.Fit()
    else:
      int_dialog = None

    # Get values and set parameters
    if (int_dialog and int_dialog.ShowModal() == wx.ID_OK):
      self.iota_phil = self.iota_phil.fetch(source=int_dialog.proc_phil)
      self.target_phil = int_dialog.target_phil
      int_dialog.Destroy()

  def onAnalysisOptions(self, e):
    """ Opens dialog for integrated dataset analysis options
    :param e: event object for self.input_window.opt_btn_analysis
    :return: modifies self.iota_phil with updated parameters
    """

    if self.gparams.advanced.processing_backend == 'cctbx.xfel':
      an_dialog = dlg.AnalysisWindow(self,
                                     phil=self.iota_phil,
                                     title='Dataset Analysis Options',
                                     style=wx.DEFAULT_DIALOG_STYLE |
                                           wx.STAY_ON_TOP | wx.RESIZE_BORDER)
    elif self.gparams.advanced.processing_backend == 'ha14':
      an_dialog = dlg.HA14AnalysisWindow(self,
                                         phil=self.iota_phil,
                                         title='Dataset Analysis Options',
                                         style=wx.DEFAULT_DIALOG_STYLE |
                                               wx.STAY_ON_TOP | wx.RESIZE_BORDER)
    else:
      an_dialog = None

    if an_dialog:
      an_dialog.SetMinSize((600, -1))
      an_dialog.Fit()

      # Get values and set parameters
      if (an_dialog.ShowModal() == wx.ID_OK):
        self.iota_phil = self.iota_phil.fetch(source=an_dialog.viz_phil)
      an_dialog.Destroy()

  def init_settings(self):
    # Grab params from main window class

    # Get list of inputs from input window
    idxs = self.input_window.input.ctr.GetItemCount()
    inputs = [self.input_window.input.ctr.GetItemData(i).path for i in range(idxs)]

    # Set all main window params (including inputs)
    self.gparams = self.iota_phil.extract()
    self.gparams.input = inputs
    self.gparams.description = util.noneset(
      self.input_window.project_title.ctr.GetValue())
    self.gparams.output = self.input_window.project_folder.ctr.GetValue()

    # Format main IOTA PHIL
    self.iota_phil = self.iota_phil.format(python_object=self.gparams)

  def OnAboutBox(self, e):
    """ About dialog """
    info = wx.adv.AboutDialogInfo()
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
    info.AddTranslator('Art Lyubimov')
    wx.adv.AboutBox(info)

  def onRecovery(self, e):
    # Find finished runs and display results
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

    paths = [os.path.join(int_folder, p) for p in os.listdir(int_folder)]
    paths = [p for p in paths if os.path.isdir(p)]

    path_dlg = dlg.RecoveryDialog(self)
    path_dlg.insert_paths(paths)

    if path_dlg.ShowModal() == wx.ID_OK:
      self.reset_settings()
      selected = path_dlg.selected
      recovery_mode = path_dlg.recovery_mode
      int_path = selected[1]
      init_file = os.path.join(int_path, 'init.cfg')

      rec_init = None
      if os.path.isfile(init_file):
        try:
          rec_init = ep.load(init_file)
          tmp_phil = rec_init.iota_phil
        except Exception:
          pass
        else:
          fixer = inp.PHILFixer()
          tmp_phil.show()
          self.iota_phil = fixer.run(old_phil=tmp_phil)

      if not rec_init:
        rec_init = UIInitAll()
        rec_init.int_base = int_path
        rec_init.obj_base = os.path.join(int_path, 'image_objects')
        rec_init.fin_base = os.path.join(int_path, 'final')
        rec_init.log_base = os.path.join(int_path, 'logs')
        rec_init.viz_base = os.path.join(int_path, 'visualization')
        rec_init.logfile = os.path.join(int_path, 'iota.log')
        with open(rec_init.logfile, 'r') as lf:
          lines = lf.readlines()[4:86]
          log_phil = ip.parse(''.join(lines))
        self.iota_phil = self.iota_phil.fetch(source=log_phil)
        rec_init.params = self.iota_phil.extract()
        input_entries = [i for i in rec_init.params.input if i is not None]
        rec_init.input_list = ginp.make_input_list(input_entries)

      self.gparams = self.iota_phil.extract()

      # Re-populate input window with settings from read-in run (check that
      # nothing has been moved)
      rec_target_phil_file = os.path.join(rec_init.int_base, 'target.phil')
      with open(rec_target_phil_file, 'r') as pf:
        rec_target_phil = pf.read()
      self.target_phil = rec_target_phil
      self.update_input_window()

      # Re-open processing window with results of the run
      if recovery_mode == 0:
        title = 'Image Processing Run {}'.format(selected[2])
        self.proc_window = frm.ProcWindow(self, -1, title=title,
                                          target_phil=rec_target_phil,
                                          phil=self.iota_phil)
        self.proc_window.recover(int_path=rec_init.int_base,
                                 init=rec_init,
                                 status=selected[0],
                                 params=self.gparams)
        self.proc_window.place_and_size(set_by='parent')
        self.proc_window.Show(True)

  def onRun(self, e):
    # Run full processing

    self.init_settings()
    input_list = []
    input_items = self.input_window.input.all_data_images

    for key, imageset in input_items.items():
      input_list.extend(imageset)

    self.proc_window = frm.ProcWindow(self, -1, title='',
                                      target_phil=self.target_phil,
                                      phil=self.iota_phil)
    init = UIInitAll()

    init.iota_phil = self.iota_phil
    if self.target_phil:
      init.target_phil = self.target_phil
    init.input_list = input_list
    self.proc_window.run(init)

    if self.proc_window.good_to_go:
      self.proc_window.place_and_size(set_by='parent')
      self.proc_window.SetTitle('Image Processing Run {}' \
                                ''.format(int(os.path.basename(init.int_base))))
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

      # Finalize settings
      self.init_settings()
      self.gparams = self.iota_phil.extract()

      # Save target PHIL file, if a PHIL script exists
      if self.target_phil is not None:
        fn = '{}.phil'.format(self.gparams.advanced.processing_backend.replace('.', '_'))
        phil_filepath = os.path.join(os.path.dirname(script_filepath), fn)
        self.gparams.cctbx_xfel.target = phil_filepath

      # Generate text of params
      final_phil = self.iota_phil.format(python_object=self.gparams)
      final_phil_txt = util.convert_phil_to_text(final_phil, script_filepath)

      if self.target_phil is not None:
        with open(phil_filepath, 'w') as phil_file:
          phil_file.write(self.target_phil)

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

  def load_script(self, filepath, reset=True, update_input_window=True):
    '''Clears settings and loads new settings from IOTA param file

    :param filepath: path to IOTA parameter file
    :param update_input_window: set to True to update input window with settings
    :return: None
    '''

    if reset:
      self.reset_settings()

    self.iota_phil = inp.process_ui_input(input_source=filepath, args=self.args,
                                          phil_args=self.phil_args, mode='file')
    self.prefs_phil = self.iota_phil
    self.gparams = self.iota_phil.extract()

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

    # Update input window with loaded parameters
    if update_input_window:
      self.update_input_window()

  def reset_settings(self):
    """ Clear all controls in input window """

    # Reset inputs
    self.input_window.input.delete_all()

    # Reset IOTA PHIL to defaults
    self.iota_phil = inp.master_phil
    self.gparams = self.iota_phil.extract()

    # Reset target phil to blank
    self.target_phil = None

    # Reset input window with default values
    self.gparams.description = ''
    self.gparams.output = os.path.abspath(os.curdir)
    self.gparams.mp.n_processors = multiprocessing.cpu_count()
    self.update_input_window()

  def update_input_window(self):
    """ Update input window with parameters in PHIL """

    # Description
    if self.gparams.description is not None:
      self.input_window.project_title.ctr.SetValue(self.gparams.description)
    else:
      self.input_window.project_title.ctr.SetValue('')

    # Output folder
    if self.gparams.output is not None:
      self.input_window.project_folder.ctr.SetValue(self.gparams.output)
    else:
      self.input_window.project_folder.ctr.SetValue(os.path.abspath(os.curdir))

    # Inputs
    for inp_path in self.gparams.input:
      if inp_path is not None:
        self.input_window.input.add_item(inp_path)

    # PHILs
    self.input_window.input_phil = self.iota_phil
    self.input_window.target_phil = self.target_phil

  def onQuit(self, e):

    # Check if processing window has been launched
    if hasattr(self, "proc_window"):

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
        shutil.rmtree(self.proc_window.init.tmp_base)
      except Exception:
        pass
      print('JOB TERMINATED!')

    self.Close()

# ------------------------------ Initialization  ----------------------------- #

from iota.components.iota_base import InitBase

class UIInitAll(InitBase):
  def __init__(self):
    InitBase.__init__(self)

  def sanity_check(self):
    """ Check for conditions necessary to starting the run
    @return: True if passed, False if failed
    """
    msg = ''
    # Check for data not found
    if len(self.input_list) == 0:
      msg = 'IOTA_UI_INIT_ERROR: Data Not Found!'
      wx.MessageBox('ERROR: Data Not Found!', 'ERROR', wx.OK | wx.ICON_ERROR)
      return False, msg

    return True, msg

  def initialize_interface(self):
    """ Override with UI-specific init procedure

    :param gparams: IOTA parameters from UI and script
    :param input_list: list of input filepaths
    :return: True = initialization OK; False = initialization Failed
    """

    # Generate default backend settings if none are specified
    if self.target_phil is None:
      self.target_phil, _ = inp.write_defaults(current_path=self.params.output,
                                               write_param_file=False,
                                               method='cctbx.xfel')

    # If input list for some reason isn't transmitted from main window, make it
    if self.input_list is None:
      self.input_list = self.make_input_list()

    # In case any sub-setting of the input is anticipated, create a list of
    # all found input for future reference (i.e. run recovery)
    self.all_input = self.input_list

    # Select range of images if turned on
    if self.params.advanced.image_range.flag_on:
      self.input_list = self.select_image_range(self.input_list)

    # Select a random subset of images if turned on
    if self.params.advanced.random_sample.flag_on and \
            self.params.advanced.random_sample.number < len(self.input_list):
      self.input_list = self.select_random_subset(self.input_list)

    # Run the sanity check procedure & return result
    is_good, msg = self.sanity_check()

    return is_good, msg
