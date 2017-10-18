from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 10/18/2017
Description : IOTA GUI Initialization module
'''

import os
import wx
import wx.lib.agw.ultimatelistctrl as ulc
from wxtbx import bitmaps
import multiprocessing

from iotbx import phil as ip
from libtbx import easy_pickle as ep
from cctbx import miller
assert miller

from iota import iota_version, gui_description, gui_license
import iota.components.iota_input as inp
import iota.components.iota_misc as misc
import iota.components.iota_frames as frm
import iota.components.iota_dialogs as dlg
from iota.components.iota_utils import InputFinder

ginp = InputFinder()
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


# ------------------------------- Main Window -------------------------------- #

class MainWindow(wx.Frame):
  ''' Frame housing the entire app; all windows open from this one '''

  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(800, 500))

    # TODO: Allow GUI to be loaded with command line args, parse into PHIL
    self.iota_phil = inp.master_phil
    self.prefs_phil = None
    self.target_phil = None
    self.term_file = None

    # Create some defaults on startup
    # Figure out temp folder
    self.iparams = self.iota_phil.extract()
    tmp_folder = '/tmp/{}_{}'.format(user, pid)
    self.iparams.advanced.temporary_output_folder = tmp_folder

    self.iota_phil = self.iota_phil.format(python_object=self.iparams)

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

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    # Toolbar
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS | wx.TB_TEXT)
    quit_bmp = bitmaps.fetch_icon_bitmap('actions', 'exit')
    self.tb_btn_quit = self.toolbar.AddLabelTool(wx.ID_EXIT, label='Quit',
                                                 bitmap=quit_bmp,
                                                 shortHelp='Quit',
                                                 longHelp='Quit IOTA')
    pref_bmp = bitmaps.fetch_icon_bitmap('apps', 'advancedsettings')
    self.tb_btn_prefs = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                  label='Preferences',
                                                  bitmap=pref_bmp,
                                                  shortHelp='Preferences',
                                                  longHelp='IOTA Preferences')
    self.toolbar.AddSeparator()
    load_bmp = bitmaps.fetch_icon_bitmap('actions', 'open')
    self.tb_btn_load = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                 label='Load Script',
                                                 bitmap=load_bmp,
                                                 shortHelp='Load Script',
                                                 longHelp='Load IOTA Script')
    save_bmp = bitmaps.fetch_icon_bitmap('actions', 'save')
    self.tb_btn_save = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                 label='Save Script',
                                                 bitmap=save_bmp,
                                                 shortHelp='Save Script',
                                                 longHelp='Save IOTA Script')
    reset_bmp = bitmaps.fetch_icon_bitmap('actions', 'reload')
    self.tb_btn_reset = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                  label='Reset',
                                                  bitmap=reset_bmp,
                                                  shortHelp='Reset Settings',
                                                  longHelp='Reset IOTA settings with defaults')
    self.toolbar.AddSeparator()
    analyze_bmp = bitmaps.fetch_icon_bitmap('mimetypes', 'text-x-generic-2')
    self.tb_btn_analysis = self.toolbar.AddLabelTool(wx.ID_ANY, label='Recover',
                                                     bitmap=analyze_bmp,
                                                     shortHelp='Recover',
                                                     longHelp='Recover run, show statistics and restart if aborted ')
    run_bmp = bitmaps.fetch_icon_bitmap('actions', 'run')
    self.tb_btn_run = self.toolbar.AddLabelTool(wx.ID_ANY, label='Run',
                                                bitmap=run_bmp,
                                                shortHelp='Run',
                                                longHelp='Run all stages of refinement')

    # Test buttons for test windows - comment out when not needed
    # self.toolbar.AddSeparator()
    test_bmp = bitmaps.fetch_icon_bitmap('actions', 'utilities')
    self.tb_btn_test = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                 label='Test',
                                                 bitmap=test_bmp)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_test)
    self.toolbar.RemoveTool(self.tb_btn_test.GetId())

    # These buttons will be disabled until input path is provided
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.Realize()

    # Instantiate windows
    self.input_window = frm.InputWindow(self, phil=self.iota_phil)

    # Single input window
    self.main_sizer.Add(self.input_window, 1,
                        flag=wx.ALL | wx.EXPAND,
                        border=10)
    self.main_sizer.Add((-1, 20))

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

  def onItemInserted(self, e):
    print self.input_window.input.all_data_images

  def onReset(self, e):
    self.reset_settings()

  def onPreferences(self, e):
    ''' Opens dialog for IOTA preferences
    :param e: event object for self.tb_btn_prefs
    :return: modifies self.iota_phil with updated parameters
    '''
    prefs = dlg.IOTAPreferences(self, phil=self.iota_phil)
    prefs.set_choices()

    if prefs.ShowModal() == wx.ID_OK:
      self.iota_phil = self.iota_phil.fetch(source=prefs.prefs_phil)
    prefs.Destroy()

    self.input_window.input_phil = self.iota_phil
    self.input_window.gparams = self.iota_phil.extract()

  def onImportOptions(self, e):
    ''' Opens dialog for image import options
    :param e: event object for self.input_window.opt_btn_import
    :return: modifies self.iota_phil with updated parameters
    '''
    imp_dialog = dlg.ImportWindow(self,
                                  phil=self.iota_phil,
                                  title='Import Options',
                                  style=wx.DEFAULT_DIALOG_STYLE | wx.STAY_ON_TOP)
    imp_dialog.Fit()

    if (imp_dialog.ShowModal() == wx.ID_OK):
      self.iota_phil = self.iota_phil.fetch(source=imp_dialog.import_phil)
    imp_dialog.Destroy()

  def onProcessOptions(self, e):
    ''' Opens dialog for image processing options, either for cctbx or DIALS
    depending on user selection.
    :param e: event object for self.input_window.opt_btn_process
    :return: modifies self.iota_phil with updated parameters
    '''

    # For cctbx.xfel options
    if self.input_window.int_box.ctr.GetCurrentSelection() == 0:
      int_dialog = dlg.CCTBXOptions(self,
                                    phil=self.iota_phil,
                                    target=self.target_phil,
                                    title='cctbx.xfel Options',
                                    style=wx.DEFAULT_DIALOG_STYLE |
                                          wx.STAY_ON_TOP | wx.RESIZE_BORDER)
      int_dialog.SetMinSize((600, -1))
      int_dialog.Fit()

      # Get values and set parameters
      if (int_dialog.ShowModal() == wx.ID_OK):
        self.iota_phil = self.iota_phil.fetch(sources=[int_dialog.proc_phil])
        self.target_phil = int_dialog.target_phil
      int_dialog.Destroy()

    # For DIALS options
    elif self.input_window.int_box.ctr.GetCurrentSelection() == 1:
      int_dialog = dlg.DIALSOptions(self,
                                    phil=self.iota_phil,
                                    target=self.target_phil,
                                    title='DIALS Options',
                                    style=wx.DEFAULT_DIALOG_STYLE |
                                          wx.STAY_ON_TOP | wx.RESIZE_BORDER)
      int_dialog.SetMinSize((600, -1))
      int_dialog.Fit()

      # Get values and set parameters
      if (int_dialog.ShowModal() == wx.ID_OK):
        self.iota_phil = self.iota_phil.fetch(source=int_dialog.proc_phil)
        self.target_phil = int_dialog.target_phil
      int_dialog.Destroy()

  def onAnalysisOptions(self, e):
    ''' Opens dialog for integrated dataset analysis options
    :param e: event object for self.input_window.opt_btn_analysis
    :return: modifies self.iota_phil with updated parameters
    '''

    an_dialog = dlg.AnalysisWindow(self,
                                   phil=self.iota_phil,
                                   title='Dataset Analysis Options',
                                   style=wx.DEFAULT_DIALOG_STYLE |
                                         wx.STAY_ON_TOP | wx.RESIZE_BORDER)
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
    self.gparams.description = misc.noneset(
      self.input_window.project_title.ctr.GetValue())
    self.gparams.output = self.input_window.project_folder.ctr.GetValue()
    self.gparams.n_processors = self.input_window.opt_spc_nprocs.ctr.GetValue()
    self.gparams.advanced.integrate_with = \
      str(self.input_window.int_box.ctr.GetString(
        self.input_window.int_box.ctr.GetSelection())).lower()

    # Format main IOTA PHIL
    self.iota_phil = self.iota_phil.format(python_object=self.gparams)

  def OnAboutBox(self, e):
    ''' About dialog '''
    info = wx.AboutDialogInfo()
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
    wx.AboutBox(info)

  def onInput(self, e):
    if self.input_window.inp_box.ctr.GetValue() != '':
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), True)
    else:
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)

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
      int_path = selected[1]

      init_file = os.path.join(int_path, 'init.cfg')

      if os.path.isfile(init_file):
        rec_init = ep.load(init_file)
        tmp_phil = inp.master_phil.format(python_object=rec_init.params)
        self.iota_phil = self.iota_phil.fetch(source=tmp_phil)
      else:
        rec_init = InitAll(iver=iota_version)
        rec_init.int_base = int_path
        rec_init.obj_base = os.path.join(int_path, 'image_objects')
        rec_init.fin_base = os.path.join(int_path, 'final')
        rec_init.log_base = os.path.join(int_path, 'logs')
        rec_init.viz_base = os.path.join(int_path, 'visualization')
        rec_init.logfile = os.path.join(int_path, 'iota.log')
        with open(rec_init.logfile, 'r') as lf:
          log_phil = ip.parse(''.join(lf.readlines()[4:86]))
        self.iota_phil = self.iota_phil.fetch(source=log_phil)
        rec_init.params = self.iota_phil.extract()
        input_entries = [i for i in rec_init.params.input if i != None]
        rec_init.input_list = ginp.make_input_list(input_entries)

      self.gparams = self.iota_phil.extract()

      # # Test / fix input folders / files
      # for inp in self.gparams.input:
      #   if not os.path.exists(inp):
      #     error_msg = wx.MessageDialog(None,
      #                                  'INPUT PATH {} NOT FOUND!\n'
      #                                  'Would you like to look for it?'
      #                                  ''.format(inp),
      #                                  wx.YES_NO | wx.ICON_ERROR)
      #     if (error_msg.ShowModal() == wx.ID_YES):
      #       open_dlg = wx.DirDialog(self, "Choose the input folder:",
      #                               style=wx.DD_DEFAULT_STYLE)
      #       if open_dlg.ShowModal() == wx.ID_OK:
      #         path = os.path.abspath(open_dlg.GetPath())
      #         if path not in self.gparams.input:
      #           self.gparams.input.append(path)
      #       else:
      #         return
      #     else:
      #       return
      #
      # # Test / fix output folder
      # if not os.path.isdir(self.gparams.output):
      #   error_msg = wx.MessageDialog(None,
      #                                'OUTPUT FOLDER {} NOT FOUND!\n'
      #                                'Would you like to look for it?'
      #                                ''.format(self.gparams.output),
      #                                wx.YES_NO | wx.ICON_ERROR)
      #   if (error_msg.ShowModal() == wx.ID_YES):
      #     open_dlg = wx.DirDialog(self, "Choose the output folder:",
      #                             style=wx.DD_DEFAULT_STYLE)
      #     if open_dlg.ShowModal() == wx.ID_OK:
      #       self.gparams.output = os.path.abspath(open_dlg.GetPath())
      #       int_no = rec_init.int_base.split('/')[-1]
      #       cnv_no = rec_init.conv_base.split('/')[-1]
      #       rec_init.input_base = self.gparams.input
      #       rec_init.conv_base = os.path.join(self.gparams.output,
      #                                         'converted_pickles', cnv_no)
      #       rec_init.int_base = os.path.join(self.gparams.output,
      #                                        'integration', int_no)
      #       rec_init.fin_base = os.path.join(rec_init.int_base, 'final')
      #       rec_init.log_base = os.path.join(rec_init.int_base, 'logs')
      #       rec_init.obj_base = os.path.join(rec_init.int_base, 'image_objects')
      #       rec_init.viz_base = os.path.join(rec_init.int_base, 'visualization')
      #       rec_init.logfile = os.path.join(int_path, 'iot.log')
      #     else:
      #       return
      #   else:
      #     return

      # Re-populate input window with settings from read-in run (check that
      # nothing has been moved)
      rec_target_phil_file = os.path.join(rec_init.int_base, 'target.phil')
      with open(rec_target_phil_file, 'r') as pf:
        rec_target_phil = pf.read()
      self.target_phil = rec_target_phil


      self.update_input_window()

      # Re-open processing window with results of the run
      self.proc_window = frm.ProcWindow(self, -1, title='Image Processing',
                                        target_phil=rec_target_phil,
                                        phil=self.iota_phil)
      self.proc_window.recover(int_path=rec_init.int_base,
                               init=rec_init,
                               status=selected[0],
                               params=self.gparams)
      self.proc_window.Show(True)

  def onRun(self, e):
    # Run full processing

    if e.GetId() == self.tb_btn_test.GetId():  # Not testing right now
      self.init_settings()
      title = 'Test'
    else:
      self.init_settings()
      title = 'Image Processing'

    self.proc_window = frm.ProcWindow(self, -1, title=title,
                                      target_phil=self.target_phil,
                                      phil=self.iota_phil)
    init = InitAll(iver=misc.iota_version)
    self.proc_window.run(init)

    if self.proc_window.good_to_go:
      self.term_file = self.proc_window.tmp_abort_file
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
        if self.gparams.advanced.integrate_with == 'cctbx':
          phil_filepath = os.path.join(os.path.dirname(script_filepath),
                                       'cctbx.phil')
          self.gparams.cctbx.target = phil_filepath
        if self.gparams.advanced.integrate_with == 'dials':
          phil_filepath = os.path.join(os.path.dirname(script_filepath),
                                       'dials.phil')
          self.gparams.dials.target = phil_filepath

      # Generate text of params
      final_phil = self.iota_phil.format(python_object=self.gparams)

      test_params = final_phil.extract()

      with misc.Capturing() as txt_output:
        final_phil.show()
      txt_out = ''
      for one_output in txt_output:
        txt_out += one_output + '\n'

      # Save files
      with open(script_filepath, 'w') as param_file:
        param_file.write(txt_out)

      if self.target_phil is not None:
        with open(phil_filepath, 'w') as phil_file:
          phil_file.write(self.target_phil)

  def onLoadScript(self, e):
    '''
    Widget event either for Load Script menu item or toolbar button
    :param e: event object
    :return:
    '''
    load_dlg = wx.FileDialog(self,
                             message="Load script file",
                             defaultDir=os.curdir,
                             defaultFile="*.param",
                             wildcard="*.param",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if load_dlg.ShowModal() == wx.ID_OK:
      self.load_script(load_dlg.GetPaths()[0])

  def load_script(self, filepath):
    '''
    Clears settings and loads new settings from IOTA param file

    :param filepath: path to script file
    :return:
    '''

    self.reset_settings()
    # Extract params from file
    with open(filepath, 'r') as pf:
      phil_string = '\n'.join(pf.readlines())
    user_phil = ip.parse(phil_string)

    self.iota_phil = self.fix_old_phil(phil=user_phil)
    self.prefs_phil = self.iota_phil
    self.gparams = self.iota_phil.extract()

    # Pass on param PHIL to input window
    self.input_window.input_phil = self.iota_phil

    # Pass on target PHIL (if found) to input window
    if self.gparams.advanced.integrate_with == 'cctbx':
      target = self.gparams.cctbx.target
    elif self.gparams.advanced.integrate_with == 'dials':
      target = self.gparams.dials.target

    if target is not None:
      try:
        with open(target, 'r') as pf:
          self.target_phil = pf.read()
      except Exception:
        self.target_phil = None
    else:
      self.target_phil = None

    # Update input window with loaded parameters
    self.update_input_window()


  def reset_settings(self):
    ''' Clear all controls in input window '''

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
    self.gparams.n_processors = multiprocessing.cpu_count()
    self.update_input_window()

  def update_input_window(self):
    ''' Update input window with parameters in PHIL '''

    # Choice of backend
    idx = self.input_window.int_box.ctr.FindString(
      self.gparams.advanced.integrate_with)
    self.input_window.int_box.ctr.SetSelection(idx)

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

    # Number of processors
    self.input_window.opt_spc_nprocs.ctr.SetValue(self.gparams.n_processors)

    # PHILs
    self.input_window.input_phil = self.iota_phil
    self.input_window.target_phil = self.target_phil

  def fix_old_phil(self, phil):
    ''' Backwards compatibility: convert settings from old format to new '''

    temp_phil = inp.master_phil.fetch(source=phil)
    params = temp_phil.extract()

    # Renaming of imported images
    if not hasattr(params.image_conversion, 'rename_pickle'):
      prefix = params.image_conversion.rename_pickle_prefix
      if 'none' in str(prefix).lower():
        params.image_conversion.__inject__('rename_pickle',
                                           'keep_file_structure')
        params.image_conversion.rename_pickle_prefix = None
      elif 'auto' in str(prefix).lower():
        params.image_conversion.__inject__('rename_pickle', 'auto_filename')
        params.image_conversion.rename_pickle_prefix = None
      else:
        params.image_conversion.__inject__('rename_pickle', 'custom_filename')
        if hasattr(prefix, '__iter__'):
          if len(prefix) > 1:
            params.image_conversion.rename_pickle_prefix = \
              [i for i in prefix if "*" in i][0].replace('*', '')
          else:
            params.image_conversion.rename_pickle_prefix = prefix[0]

    if str(params.image_conversion.square_mode).lower() == 'none':
      params.image_conversion.square_mode = 'no_modification'

    if str(params.image_triage.type).lower() == 'none':
      params.image_triage.type = 'no_triage'

    if str(params.cctbx.grid_search.type).lower() == 'none':
      params.cctbx.grid_search.type = 'no_grid_search'

    if str(params.analysis.viz).lower() == 'none':
      params.analysis.viz = 'no_visualization'

    fixed_phil = temp_phil.format(python_object=params)

    return fixed_phil


  def onQuit(self, e):
    if self.term_file is not None:
      with open(self.term_file, 'w') as tf:
        tf.write('')
    self.Close()


# ------------------------------ Initialization  ----------------------------- #


class InitAll(object):
  """ Class to initialize current IOTA run in GUI

      iver = IOTA version (hard-coded)
      help_message = description (hard-coded)

  """

  def __init__(self, iver):
    from datetime import datetime
    self.iver = iver
    self.user_id = user
    self.now = "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())
    self.input_base = None
    self.conv_base = None
    self.obj_base = None
    self.int_base = None


  def make_input_list(self):
    """ Reads input directory or directory tree and makes lists of input images.
        Optional selection of a random subset
    """
    input_entries = [i for i in self.params.input if i != None]
    input_list = ginp.make_input_list(input_entries)

    # Pick a randomized subset of images
    if self.params.advanced.random_sample.flag_on and \
                    self.params.advanced.random_sample.number < len(input_list):
      inp_list = self.select_random_subset(input_list)
    else:
      inp_list = input_list

    return inp_list


  def select_random_subset(self, input_list):
    """ Selects random subset of input entries """
    import random

    random_inp_list = []
    if self.params.advanced.random_sample.number == 0:
      if len(input_list) <= 5:
        random_sample_number = len(input_list)
      elif len(input_list) <= 50:
        random_sample_number = 5
      else:
        random_sample_number = int(len(input_list) * 0.1)
    else:
      random_sample_number = self.params.advanced.random_sample.number

    for i in range(random_sample_number):
      random_number = random.randrange(0, len(input_list))
      if input_list[random_number] in random_inp_list:
        while input_list[random_number] in random_inp_list:
          random_number = random.randrange(0, len(input_list))
        random_inp_list.append(input_list[random_number])
      else:
        random_inp_list.append(input_list[random_number])

    return random_inp_list


  def make_int_object_list(self):
    """ Generates list of image objects from previous grid search """
    from libtbx import easy_pickle as ep

    if self.params.cctbx.selection.select_only.grid_search_path == None:
      int_dir = misc.set_base_dir('integration', True)
    else:
      int_dir = self.params.cctbx.selection.select_only.grid_search_path

    img_objects = []

    # Inspect integration folder for image objects
    for root, dirs, files in os.walk(int_dir):
      for filename in files:
        found_file = os.path.join(root, filename)
        if found_file.endswith(('int')):
          obj = ep.load(found_file)
          img_objects.append(obj)

    # Pick a randomized subset of images
    if self.params.advanced.random_sample.flag_on and \
                  self.params.advanced.random_sample.number < len(img_objects):
      gs_img_objects = self.select_random_subset(img_objects)
    else:
      gs_img_objects = img_objects

    return gs_img_objects


  def sanity_check(self):
    ''' Check for conditions necessary to starting the run
    @return: True if passed, False if failed
    '''

    # Check for existence of appropriate target files. If none are specified,
    # ask to generate defaults; if user says no, fail sanity check. If file is
    # specified but doesn't exist, show error message and fail sanity check
    if self.target_phil == None:
      if self.params.advanced.integrate_with == 'cctbx':
        write_def = wx.MessageDialog(None,
                                     'WARNING! No target file for CCTBX.XFEL. '
                                     'Generate defaults?','WARNING',
                                     wx.YES_NO | wx.NO_DEFAULT | wx.ICON_EXCLAMATION)
        if (write_def.ShowModal() == wx.ID_YES):
          self.target_phil, _ = inp.write_defaults(current_path=self.params.output,
                                                   write_param_file=False,
                                                   method='cctbx')
          return True
        else:
          return False

      elif self.params.advanced.integrate_with == 'dials':
          write_def = wx.MessageDialog(None,
                                       'WARNING! No target file for DIALS. '
                                       'Generate defaults?', 'WARNING',
                                       wx.YES_NO | wx.NO_DEFAULT | wx.ICON_EXCLAMATION)
          if (write_def.ShowModal() == wx.ID_YES):
            self.target_phil, _ = inp.write_defaults(current_path=self.params.output,
                                                     write_param_file=False,
                                                     method='dials')
            return True
          else:
            return False
    else:
      return True

  def run(self, gparams, target_phil=None, list_file=None):
    ''' Run initialization for IOTA GUI

        gparams = IOTA parameters from the GUI elements in PHIL format
        gtxt = text version of gparams
        list_file = if "Write Input List" button pressed, specifies name
                    of list file
    '''

    from iota.components.iota_init import parse_command_args
    self.args, self.phil_args = parse_command_args(self.iver, '').parse_known_args()
    self.params = gparams
    self.target_phil = target_phil

    # Call function to read input folder structure (or input file) and
    # generate list of image file paths
    if self.params.cctbx.selection.select_only.flag_on:
      self.gs_img_objects = self.make_int_object_list()
      self.input_list = [i.conv_img for i in self.gs_img_objects]
    else:
      self.input_list = self.make_input_list()

    # Check for data not found
    if len(self.input_list) == 0:
      wx.MessageBox('ERROR: Data Not Found!', 'ERROR', wx.OK | wx.ICON_ERROR)
      return False

    # If list-only option selected, output list only
    if list_file != None:
      with open(list_file, "w") as lf:
        for i, input_file in enumerate(self.input_list, 1):
          lf.write('{}\n'.format(input_file))
      return True

    # Run the sanity check procedure
    if not self.sanity_check():
      return False

    # If fewer images than requested processors are supplied, set the number of
    # processors to the number of images
    if self.params.n_processors > len(self.input_list):
      self.params.n_processors = len(self.input_list)

    # Generate base folder paths
    self.conv_base = misc.set_base_dir('converted_pickles',
                                       out_dir=self.params.output)
    self.int_base = misc.set_base_dir('integration', out_dir=self.params.output)
    self.obj_base = os.path.join(self.int_base, 'image_objects')
    self.fin_base = os.path.join(self.int_base, 'final')
    self.log_base = os.path.join(self.int_base, 'logs')
    self.viz_base = os.path.join(self.int_base, 'visualization')
    if str(self.params.advanced.temporary_output_folder).lower() in ('none',''):
      self.tmp_base = os.path.join(self.int_base, 'tmp')
    else:
      self.tmp_base = os.path.join(self.params.advanced.temporary_output_folder)

    # Generate base folders
    os.makedirs(self.int_base)
    os.makedirs(self.obj_base)
    os.makedirs(self.fin_base)
    os.makedirs(self.log_base)
    try:
      if not os.path.isdir(self.tmp_base):
        os.makedirs(self.tmp_base)
    except OSError:
      pass

    # Determine input base
    self.input_base = os.path.abspath(os.path.dirname(os.path.commonprefix(self.input_list)))

    # Initialize main log
    self.logfile = os.path.abspath(os.path.join(self.int_base, 'iota.log'))

    # Write target file and record its location in params
    local_target_file = os.path.join(self.int_base, 'target.phil')
    if type(self.target_phil) == list:
      self.target_phil = '\n'.join(self.target_phil)
    with open(local_target_file, 'w') as tf:
      tf.write(self.target_phil)

    if self.params.advanced.integrate_with == 'cctbx':
      self.params.cctbx.target = local_target_file
    elif self.params.advanced.integrate_with == 'dials':
      self.params.dials.target = local_target_file

    # Collect final params and convert to PHIL object
    final_phil = inp.master_phil.format(python_object=self.params)

    # Generate text of params
    with misc.Capturing() as txt_output:
      final_phil.show()
    self.txt_out = ''
    for one_output in txt_output:
      self.txt_out += one_output + '\n'

    # Log starting info
    misc.main_log(self.logfile, '{:=^80} \n'.format(' IOTA MAIN LOG '))
    misc.main_log(self.logfile, '{:-^80} \n'.format(' SETTINGS FOR THIS RUN '))
    misc.main_log(self.logfile, self.txt_out)

    # Log cctbx.xfel / DIALS settings
    misc.main_log(self.logfile, '{:-^80} \n\n'
                                ''.format(' TARGET FILE ({}) CONTENTS '
                                          ''.format(local_target_file)))
    misc.main_log(self.logfile, self.target_phil)

    return True
