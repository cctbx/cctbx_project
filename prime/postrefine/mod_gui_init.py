from __future__ import absolute_import, division, print_function
from six.moves import map

'''
Author      : Lyubimov, A.Y.
Created     : 05/01/2016
Last Changed: 10/21/2018
Description : PRIME GUI Initialization module
'''

import os
import wx

from libtbx import easy_run
import iotbx.phil as ip

from iota import prime_description, prime_license
from iota.components.iota_utils import Capturing
from iota.components.iota_ui_base import IOTABaseFrame
import iota.components.iota_ui_controls as ct
from prime.postrefine import mod_gui_frames as frm
from prime.postrefine import mod_gui_dialogs as dlg
from prime.postrefine.mod_input import master_phil

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

user = os.getlogin()
icons = os.path.join(os.path.dirname(os.path.abspath(ct.__file__)), 'icons/')

def str_split(string, delimiters=(' ', ','), maxsplit=0):
  import re
  rexp = '|'.join(map(re.escape, delimiters))
  return re.split(rexp, string, maxsplit)

# ------------------------------ PRIME Windows ------------------------------- #

class PRIMEWindow(IOTABaseFrame):

  def __init__(self, parent, id, title, prefix='prime'):
    IOTABaseFrame.__init__(self, parent, id, title, size=(800, 500))

    self.prime_filename = '{}.phil'.format(prefix)
    self.prime_phil = master_phil

    # Toolbar
    self.initialize_toolbar()
    self.tb_btn_quit = self.add_tool(label='Quit',
                                     bitmap=('actions', 'exit'),
                                     shortHelp='Exit PRIME')

    self.tb_btn_prefs = self.add_tool(label='Preferences',
                                      bitmap=('apps', 'advancedsettings'),
                                      shortHelp='PRIME GUI Settings')
    self.add_toolbar_separator()
    self.tb_btn_load = self.add_tool(label='Load Script',
                                     bitmap=('actions', 'open'),
                                     shortHelp='Load PRIME Script')
    self.tb_btn_save = self.add_tool(label='Save Script',
                                     bitmap=('actions', 'save'),
                                     shortHelp='Save PRIME Script')
    self.tb_btn_reset = self.add_tool(label='Reset',
                                                  bitmap=('actions', 'reload'),
                                                  shortHelp='Reset Settings')
    self.toolbar.AddSeparator()
    self.tb_btn_analysis = self.add_tool(label='Recover',
                                         bitmap=(
                                         'mimetypes', 'text-x-generic-2'),
                                         shortHelp='Recover PRIME run')
    self.tb_btn_run = self.add_tool(label='Run',
                                    bitmap=('actions', 'run'),
                                    shortHelp='Run PRIME')
    # These buttons will be disabled until input path is provided
    self.set_tool_state(self.tb_btn_run, enable=False)
    self.realize_toolbar()

    # Status bar
    self.sb = self.CreateStatusBar()
    self.sb.SetFieldsCount(3)
    self.sb.SetStatusWidths([320, 200, -2])

    # Menu bar
    menubar = wx.MenuBar()
    m_help = wx.Menu()
    m_file = wx.Menu()
    self.mb_load_script = m_file.Append(wx.ID_OPEN, '&Load Script...')
    self.mb_save_script = m_file.Append(wx.ID_SAVE, '&Save Script...')
    self.mb_about = m_help.Append(wx.ID_ANY, '&About')
    menubar.Append(m_file, '&File')
    menubar.Append(m_help, '&Help')
    self.SetMenuBar(menubar)

    # Place elements in main PRIME window
    main_box = wx.BoxSizer(wx.VERTICAL)

    # Instantiate windows
    self.input_window = frm.PRIMEInputWindow(self, phil=self.prime_phil)

    # Single input window
    main_box.Add(self.input_window, 1, flag=wx.ALL | wx.EXPAND, border=10)
    main_box.Add((-1, 20))

    # Menubar button bindings
    self.Bind(wx.EVT_MENU, self.OnAboutBox, self.mb_about)
    self.Bind(wx.EVT_MENU, self.onSaveScript, self.mb_save_script)
    self.Bind(wx.EVT_MENU, self.onLoadScript, self.mb_load_script)

    # Toolbar button bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onPreferences, self.tb_btn_prefs)
    self.Bind(wx.EVT_TOOL, self.onRecovery, self.tb_btn_analysis)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_run)
    self.Bind(wx.EVT_TOOL, self.onLoadScript, self.tb_btn_load)
    self.Bind(wx.EVT_TOOL, self.onSaveScript, self.tb_btn_save)
    self.Bind(wx.EVT_TOOL, self.onReset, self.tb_btn_reset)

    # Input window bindings
    self.Bind(wx.EVT_BUTTON, self.onAdvancedOptions, self.input_window.opt_btn)


    # Draw the main window sizer
    self.SetSizer(main_box)

  def OnAboutBox(self, e):
    ''' About dialog '''
    info = wx.AboutDialogInfo()
    info.SetName('PRIME')
    info.SetWebSite('http://cci.lbl.gov/xfel')
    info.SetLicense(prime_license)
    info.SetDescription(prime_description)
    info.AddDeveloper('Monarin Uervirojnangkoorn')
    info.AddDeveloper('Axel Brunger')
    wx.AboutBox(info)

  def onPreferences(self, e):
    prefs = dlg.PRIMEPreferences(self,
                                 phil=self.prime_phil,
                                 title='Advanced PRIME Options',
                                 style=wx.DEFAULT_DIALOG_STYLE |
                                       wx.STAY_ON_TOP |
                                       wx.RESIZE_BORDER)
    prefs.SetMinSize((600, -1))
    prefs.Fit()

    if prefs.ShowModal() == wx.ID_OK:
      self.prime_phil = self.prime_phil.fetch(prefs.pref_phil)
      self.input_window.regenerate_params(self.prime_phil)

  def onAdvancedOptions(self, e):
    self.input_window.update_settings()
    self.prime_phil = self.prime_phil.fetch(self.input_window.prime_phil)
    advanced = dlg.PRIMEAdvancedOptions(self,
                                        phil=self.prime_phil,
                                        title='Advanced PRIME Options',
                                        style=wx.DEFAULT_DIALOG_STYLE |
                                              wx.STAY_ON_TOP |
                                              wx.RESIZE_BORDER)
    advanced.SetMinSize((600, -1))
    advanced.Fit()

    if advanced.ShowModal() == wx.ID_OK:
      self.prime_phil = self.prime_phil.fetch(advanced.phil_script)
      self.prime_phil = self.prime_phil.fetch(advanced.new_prime_phil)
      self.input_window.regenerate_params(self.prime_phil)

    advanced.Destroy()


  def onInput(self, e):
    if self.input_window.inp_box.ctr.GetValue() != '':
      self.set_tool_state(self.tb_btn_run, enable=True)
    else:
      self.set_tool_state(self.tb_btn_run, enable=False)

  def init_settings(self):
    self.input_window.update_settings()
    self.pparams = self.input_window.pparams
    self.out_dir = self.input_window.out_dir

  def sanity_check(self):
    '''
    Goes through and checks that the key parameters are populated; pops
    up an error message if they are not
    :return: True if satisfied, False if not
    '''

    # Check to see that pixel size is specified (PRIME won't run without it)
    if self.pparams.pixel_size_mm is None:
      warning = wx.MessageDialog(None,
                                 caption='Warning!',
                                 message='Pixel size not specified!',
                                 style=wx.OK)
      warning.ShowModal()
      warning.Destroy()
      return False

    return True

  def onRecovery(self, e):
    # Find finished runs and display results
    p_folder = os.path.abspath('{}/prime'.format(os.curdir))

    if not os.path.isdir(p_folder):
      open_dlg = wx.DirDialog(self, "Choose the integration run:",
                              style=wx.DD_DEFAULT_STYLE)
      if open_dlg.ShowModal() == wx.ID_OK:
        p_folder = open_dlg.GetPath()
        open_dlg.Destroy()
      else:
        open_dlg.Destroy()
        return

    paths = [os.path.join(p_folder, p) for p in os.listdir(p_folder)]
    paths = [p for p in paths if (os.path.isdir(p) and
                                  os.path.basename(p).isdigit())]

    path_dlg = dlg.RecoveryDialog(self)
    path_dlg.insert_paths(paths)

    if path_dlg.ShowModal() == wx.ID_OK:
      selected = path_dlg.selected
      recovery = path_dlg.recovery_mode
      prime_path = selected[1]
      prime_status = selected[0]
      settings_file = os.path.join(prime_path, 'settings.phil')

      # If neither log-file nor stat file are found, terminate; otherwise
      # import settings
      if not os.path.isfile(settings_file):
        wx.MessageBox('Cannot Import This Run \n(No Files Found!)',
                      'Info', wx.OK | wx.ICON_ERROR)
        return
      else:
        self.reset_settings()
        with open(settings_file, 'r') as sf:
          phil_string = sf.read()
        read_phil = ip.parse(phil_string)
        self.prime_phil = master_phil.fetch(source=read_phil)
        self.pparams = self.prime_phil.extract()
        self.input_window.regenerate_params(self.prime_phil)
        self.update_input_window()

      # If any cycles (or full run) were completed, show results
      if prime_status == 'Unknown':
        return

      if recovery == 0:
        self.prime_run_window = frm.PRIMERunWindow(self, -1,
                                                   title='PRIME Output',
                                                   params=self.pparams,
                                                   prime_file=settings_file,
                                                   recover=True)
        self.prime_run_window.Show(True)
        self.prime_run_window.recover()

  def onRun(self, e):
    # Run full processing

    self.init_settings()
    if self.sanity_check():
      prime_phil = master_phil.format(python_object=self.pparams)

      with Capturing() as output:
        prime_phil.show()

      txt_out = ''
      for one_output in output:
        txt_out += one_output + '\n'

      source_dir = os.path.dirname(self.out_dir)
      prime_file = os.path.join(source_dir, self.prime_filename)
      out_file = os.path.join(self.out_dir, 'stdout.log')
      with open(prime_file, 'w') as pf:
        pf.write(txt_out)

      self.prime_run_window = frm.PRIMERunWindow(self, -1,
                                                 title='PRIME Output',
                                                 params=self.pparams,
                                                 prime_file=prime_file)
      self.prime_run_window.prev_pids = easy_run.fully_buffered('pgrep -u {} {}'
                                        ''.format(user, python)).stdout_lines
      self.prime_run_window.place_and_size(set_by='parent')
      self.prime_run_window.Show(True)

  def onSequence(self, e):
    pass

  def onSaveScript(self, e):
    self.init_settings()

    # Generate text of params
    final_phil = master_phil.format(python_object=self.pparams)
    with Capturing() as txt_output:
      final_phil.show()
    txt_out = ''
    for one_output in txt_output:
      txt_out += one_output + '\n'

    # Save param file
    save_dlg = wx.FileDialog(self,
                             message="Save PRIME Script",
                             defaultDir=os.curdir,
                             defaultFile="*.phil",
                             wildcard="*.phil",
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                             )
    if save_dlg.ShowModal() == wx.ID_OK:
      with open(save_dlg.GetPath(), 'w') as savefile:
        savefile.write(txt_out)

    save_dlg.Destroy()

  def onLoadScript(self, e):
    # Extract params from file
    load_dlg = wx.FileDialog(self,
                             message="Load script file",
                             defaultDir=os.curdir,
                             defaultFile="*.phil",
                             wildcard="*.phil",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if load_dlg.ShowModal() == wx.ID_OK:
      script = load_dlg.GetPaths()[0]
      out_dir = os.path.dirname(script)
      self.prime_filename=os.path.basename(script)
      self.load_script(out_dir=out_dir)
    load_dlg.Destroy()

  def onReset(self, e):
    self.reset_settings()

  def load_script(self, out_dir):
    ''' Loads PRIME script '''
    script = os.path.join(out_dir, self.prime_filename)
    with open(script, 'r') as sf:
      phil_string = sf.read()

    self.reset_settings()

    user_phil = ip.parse(phil_string)
    self.prime_phil = master_phil.fetch(source=user_phil)
    self.pparams = self.prime_phil.extract()
    self.input_window.pparams = self.pparams
    self.input_window.phil_string = phil_string
    self.update_input_window()


  def update_input_window(self):
    ''' Update input window with current (or default) params'''
    for input_item in self.pparams.data:
      if input_item is not None:
        self.input_window.inp_box.add_item(input_item)

    if self.pparams.run_no is not None:
      current_dir = os.path.dirname(self.pparams.run_no)
    else:
      current_dir = os.path.abspath(os.curdir)
    self.input_window.out_box.ctr.SetValue(str(current_dir))
    if str(self.input_window.out_box.ctr.GetValue).lower() == '':
      self.input_window.out_box.ctr.SetValue(self.out_dir)
    if str(self.pparams.title).lower() != 'none':
      self.input_window.project_title.ctr.SetValue(str(self.pparams.title))
    if str(self.pparams.hklisoin).lower() != 'none':
      self.input_window.inp_box.add_item(self.pparams.hklisoin)
    elif str(self.pparams.hklrefin).lower() != 'none':
      self.input_window.inp_box.add_item(self.pparams.hklrefin)
      self.input_window.opt_chk_useref.Enable()
      self.input_window.opt_chk_useref.SetValue(True)
    if str(self.pparams.n_residues).lower() == 'none':
      self.input_window.opt_spc_nres.ctr.SetValue(500)
    else:
      self.input_window.opt_spc_nres.ctr.SetValue(int(self.pparams.n_residues))
    self.input_window.opt_spc_nproc.ctr.SetValue(int(self.pparams.n_processors))

  def reset_settings(self):
    self.prime_phil = master_phil
    self.pparams = self.prime_phil.extract()
    self.input_window.inp_box.delete_all()
    self.input_window.out_box.reset_default()
    self.input_window.project_title.reset_default()
    self.input_window.opt_chk_useref.SetValue(False)
    self.input_window.opt_spc_nproc.reset_default()
    self.input_window.opt_spc_nres.reset_default()

    # Generate Python object and text of parameters
    with Capturing() as txt_output:
      master_phil.show()
    self.phil_string = ''
    for one_output in txt_output:
      self.phil_string += one_output + '\n'

  def onQuit(self, e):
    self.Destroy()
