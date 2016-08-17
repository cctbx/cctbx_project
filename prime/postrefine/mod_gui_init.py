from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 05/01/2016
Last Changed: 08/16/2016
Description : PRIME GUI Initialization module
'''

import os
import wx
from threading import Thread

from cctbx.uctbx import unit_cell
from libtbx import easy_run
from libtbx import easy_pickle as ep
import numpy as np

import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

import iota.components.iota_misc as misc
import iota.components.iota_controls as ct
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

class PRIMEWindow(wx.Frame):

  def __init__(self, parent, id, title, phil=None, prefix='prime'):
    wx.Frame.__init__(self, parent, id, title, size=(800, 500))

    self.prime_filename = '{}.phil'.format(prefix)

    # Toolbar
    self.toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.tb_btn_quit = self.toolbar.AddLabelTool(wx.ID_EXIT, label='Quit',
                                                 bitmap=wx.Bitmap('{}/32x32/exit.png'.format(icons)),
                                                 shortHelp='Quit',
                                                 longHelp='Exit PRIME')
    self.tb_btn_prefs = self.toolbar.AddLabelTool(wx.ID_ANY,
                                                  label='Preferences',
                                                  bitmap=wx.Bitmap('{}/32x32/config.png'.format(icons)),
                                                  shortHelp='Preferences',
                                                  longHelp='PRIME Preferences')
    self.toolbar.AddSeparator()
    self.tb_btn_run = self.toolbar.AddLabelTool(wx.ID_ANY, label='Run',
                                                bitmap=wx.Bitmap('{}/32x32/run.png'.format(icons)),
                                                shortHelp='Run',
                                                longHelp='Run scaling, merging and post-refinement with PRIME')
    # These buttons will be disabled until input path is provided
    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.Realize()

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
    self.input_window = PRIMEInputWindow(self, phil)
    #self.pparams = self.input_window.pparams

    # Single input window
    main_box.Add(self.input_window, flag=wx.ALL | wx.EXPAND, border=10)
    main_box.Add((-1, 20))

    # Menubar button bindings
    self.Bind(wx.EVT_MENU, self.OnAboutBox, self.mb_about)
    self.Bind(wx.EVT_MENU, self.onSaveScript, self.mb_save_script)
    self.Bind(wx.EVT_MENU, self.onLoadScript, self.mb_load_script)

    # Toolbar button bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onPreferences, self.tb_btn_prefs)
    self.Bind(wx.EVT_BUTTON, self.onInput,
              self.input_window.inp_box.btn_browse)
    self.Bind(wx.EVT_TEXT, self.onInput, self.input_window.inp_box.ctr)
    self.Bind(wx.EVT_BUTTON, self.onIsoRef,
              self.input_window.ref_box.btn_browse)
    self.Bind(wx.EVT_TEXT, self.onIsoRef, self.input_window.ref_box.ctr)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_run)
    # self.Bind(wx.EVT_BUTTON, self.onSequence,
    #           self.input_window.seq_box.btn_browse)
    # self.Bind(wx.EVT_TEXT, self.onSequence, self.input_window.seq_box.ctr)

    # Draw the main window sizer
    self.SetSizer(main_box)

  def OnAboutBox(self, e):
    ''' About dialog '''
    info = wx.AboutDialogInfo()
    info.SetName('PRIME')
    info.SetWebSite('http://cci.lbl.gov/xfel')
    info.SetLicense(misc.prime_license)
    info.SetDescription(misc.prime_description)
    info.AddDeveloper('Monarin Uervirojnangkoorn')
    info.AddDeveloper('Axel Brunger')
    wx.AboutBox(info)

  def onPreferences(self, e):
    self.pparams = self.input_window.pparams
    dlg = PRIMEPreferences(self)
    dlg.set_choices(method=self.pparams.queue.mode,
                          queue=self.pparams.queue.qname)

    if dlg.ShowModal() == wx.ID_OK:
      if dlg.method == 'multiprocessing':
        self.pparams.queue.mode = None
      else:
        self.pparams.queue.mode = dlg.method
      self.pparams.queue.qname = dlg.queue

  def onInput(self, e):
    if self.input_window.inp_box.ctr.GetValue() != '':
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), True)
    else:
      self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)

  def onIsoRef(self, e):
    if self.input_window.ref_box.ctr.GetValue() != '':
      self.input_window.opt_chk_useref.Enable()
    else:
      self.input_window.opt_chk_useref.Disable()

  def init_settings(self):
    self.pparams = self.input_window.pparams
    self.pparams.data = [self.input_window.inp_box.ctr.GetValue()]
    self.out_dir = self.input_window.out_box.ctr.GetValue()
    self.pparams.run_no = misc.set_base_dir(out_dir=self.out_dir)               # Need to change
    self.pparams.title = self.input_window.title_box.ctr.GetValue()
    if str(self.input_window.ref_box.ctr.GetValue()).lower() != '':
      self.pparams.hklisoin = self.input_window.ref_box.ctr.GetValue()
      if self.input_window.opt_chk_useref.GetValue():
        self.pparams.hklrefin = self.input_window.ref_box.ctr.GetValue()
    self.pparams.n_residues = self.input_window.opt_spc_nres.GetValue()
    self.pparams.n_processors = self.input_window.opt_spc_nproc.GetValue()

  def onRun(self, e):
    # Run full processing

    self.init_settings()
    prime_phil = master_phil.format(python_object=self.pparams)

    with misc.Capturing() as output:
      prime_phil.show()

    txt_out = ''
    for one_output in output:
      txt_out += one_output + '\n'

    prime_file = os.path.join(self.out_dir, self.prime_filename)
    out_file = os.path.join(self.out_dir, 'stdout.log')
    with open(prime_file, 'w') as pf:
      pf.write(txt_out)

    self.prime_run_window = PRIMERunWindow(self, -1,
                                           title='PRIME Output',
                                           params=self.pparams,
                                           prime_file=prime_file,
                                           out_file=out_file)
    self.prime_run_window.prev_pids = easy_run.fully_buffered('pgrep -u {} {}'
                                           ''.format(user, python)).stdout_lines
    self.prime_run_window.Show(True)

  def onSequence(self, e):
    pass

  def onSaveScript(self, e):
    self.init_settings()

    # Generate text of params
    final_phil = master_phil.format(python_object=self.pparams)
    with misc.Capturing() as txt_output:
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

  def load_script(self, out_dir):
    ''' Loads PRIME script '''
    import iotbx.phil as ip

    script = os.path.join(out_dir, self.prime_filename)
    user_phil = ip.parse(open(script).read())
    self.pparams = master_phil.fetch(sources=[user_phil]).extract()
    self.input_window.pparams = self.pparams

    self.input_window.inp_box.ctr.SetValue(str(self.pparams.data[0]))
    current_dir = os.path.dirname(self.pparams.run_no)
    self.input_window.out_box.ctr.SetValue(str(current_dir))
    if str(self.input_window.out_box.ctr.GetValue).lower() == '':
      self.input_window.out_box.ctr.SetValue(self.out_dir)
    if str(self.pparams.title).lower() != 'none':
      self.input_window.title_box.ctr.SetValue(str(self.pparams.title))
    if str(self.pparams.hklisoin).lower() != 'none':
      self.input_window.ref_box.ctr.SetValue(str(self.pparams.hklisoin))
    elif str(self.pparams.hklrefin).lower() != 'none':
      self.input_window.ref_box.ctr.SetValue(str(self.pparams.hklrefin))
      self.input_window.opt_chk_useref.SetValue(True)
    if str(self.pparams.n_residues).lower() == 'none':
      self.input_window.opt_spc_nres.SetValue(500)
    else:
      self.input_window.opt_spc_nres.SetValue(int(self.pparams.n_residues))
    self.input_window.opt_spc_nproc.SetValue(int(self.pparams.n_processors))

  def onQuit(self, e):
    self.Destroy()


class PRIMEInputWindow(wx.Panel):
  ''' Main PRIME Window panel '''

  def __init__(self, parent, phil=None):
    self.parent = parent
    super(PRIMEInputWindow, self).__init__(self.parent)

    #Generate default parameters
    if phil is not None:
      self.pparams = master_phil.format(python_object=phil).extract()
    else:
      self.pparams = master_phil.extract()

    main_box = wx.StaticBox(self, label='Main Settings')
    vbox = wx.StaticBoxSizer(main_box, wx.VERTICAL)

    # Input file box
    self.inp_box = ct.InputCtrl(self, label='Input: ',
                                label_size=(120, -1),
                                label_style='bold',
                                button=True)
    vbox.Add(self.inp_box,
             flag=wx.LEFT | wx.TOP | wx.RIGHT| wx.EXPAND,
             border=15)

    # Output file box
    self.out_box = ct.InputCtrl(self, label='Output: ',
                                label_size=(120, -1),
                                label_style='bold',
                                button=True)
    vbox.Add(self.out_box,
             flag=wx.LEFT | wx.TOP | wx.RIGHT| wx.EXPAND,
             border=15)

    # # Sequence file box
    # self.seq_box = ct.InputCtrl(self, label='Sequence: ',
    #                             label_size=(120, -1),
    #                             button=True)
    # vbox.Add(self.seq_box,
    #          flag=wx.LEFT | wx.TOP | wx.RIGHT| wx.EXPAND,
    #          border=15)

    # Isomorphous reference file box
    self.ref_box = ct.InputCtrl(self, label='Reference: ',
                                label_size=(120, -1),
                                button=True)
    vbox.Add(self.ref_box,
             flag=wx.LEFT | wx.TOP | wx.RIGHT| wx.EXPAND,
             border=15)

    # Title box
    self.title_box = ct.InputCtrl(self, label='Job Title: ',
                                  label_size=(120, -1))
    vbox.Add(self.title_box,
             flag=wx.LEFT | wx.TOP | wx.RIGHT| wx.EXPAND,
             border=15)

    # Options and number of processors
    opt_box = wx.FlexGridSizer(2, 4, 15, 10)
    self.opt_chk_useref = wx.CheckBox(self, label='Use reference in refinement')
    self.opt_chk_useref.SetValue(False)
    self.opt_chk_useref.Disable()
    self.opt_txt_nres = wx.StaticText(self, label='No. of Residues: ')
    self.opt_spc_nres = wx.SpinCtrl(self, value='500',
                                    max=(15000), size=(80,-1))
    self.opt_txt_nproc = wx.StaticText(self, label='No. of Processors: ')
    self.opt_spc_nproc = wx.SpinCtrl(self, value='8', size=(80, -1))

    self.opt_btn = wx.Button(self, label='Advanced Options...')

    opt_box.Add((120, -1))
    opt_box.Add(self.opt_chk_useref)
    opt_box.Add(self.opt_txt_nres, flag=wx.ALIGN_RIGHT)
    opt_box.Add(self.opt_spc_nres, flag=wx.ALIGN_RIGHT)
    opt_box.Add((120, -1))
    opt_box.Add(self.opt_btn)
    opt_box.Add(self.opt_txt_nproc, flag=wx.ALIGN_RIGHT)
    opt_box.Add(self.opt_spc_nproc, flag=wx.ALIGN_RIGHT)
    opt_box.AddGrowableCol(1, 1)
    vbox.Add(opt_box,
             flag=wx.LEFT | wx.TOP | wx.RIGHT| wx.EXPAND,
             border=15)

    vbox.Add((-1, 15))
    self.SetSizer(vbox)

    # Button bindings
    self.inp_box.btn_browse.Bind(wx.EVT_BUTTON, self.onInputBrowse)
    self.out_box.btn_browse.Bind(wx.EVT_BUTTON, self.onOutputBrowse)
    # self.seq_box.btn_browse.Bind(wx.EVT_BUTTON, self.onSequenceBrowse)
    self.ref_box.btn_browse.Bind(wx.EVT_BUTTON, self.onIsoRefBrowse)
    self.opt_btn.Bind(wx.EVT_BUTTON, self.onAdvancedOptions)

  def onInputBrowse(self, e):
    """ On clincking the Browse button: show the DirDialog and populate 'Input'
        box w/ selection """
    dlg = wx.DirDialog(self, "Choose the input directory:",
                       style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.inp_box.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()
    e.Skip()

  def onOutputBrowse(self, e):
    """ On clicking the Browse button: show the DirDialog and populate 'Output'
        box w/ selection """
    dlg = wx.DirDialog(self, "Choose the output directory:",
                       style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.out_box.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()
    e.Skip()

  # def onSequenceBrowse(self, e):
  #   dlg = wx.FileDialog(self,
  #                       message="Select Sequence file",
  #                       defaultDir=os.curdir,
  #                       defaultFile="*",
  #                       wildcard="*",
  #                       style=wx.OPEN | wx.FD_FILE_MUST_EXIST)
  #   if dlg.ShowModal() == wx.ID_OK:
  #     self.seq_box.ctr.SetValue(dlg.GetPaths()[0])
  #   dlg.Destroy()
  #   e.Skip()

  def onIsoRefBrowse(self, e):
    dlg = wx.FileDialog(self,
                        message="Select isomorphous reference file",
                        defaultDir=os.curdir,
                        defaultFile="*.mtz",
                        wildcard="*.mtz",
                        style=wx.OPEN | wx.FD_FILE_MUST_EXIST)
    if dlg.ShowModal() == wx.ID_OK:
      self.ref_box.ctr.SetValue(dlg.GetPaths()[0])
    dlg.Destroy()
    e.Skip()

  def onAdvancedOptions(self, e):
    advanced = PRIMEAdvancedOptions(self,
                               title='Advanced PRIME Options',
                               style=wx.DEFAULT_DIALOG_STYLE)
    advanced.Fit()

    # Set values to default parameters
    advanced.res.high.SetValue('{:4.2f}'.format(self.pparams.scale.d_max))
    advanced.res.low.SetValue('{:4.2f}'.format(self.pparams.scale.d_min))
    advanced.sg.spacegroup.SetValue(str(self.pparams.target_space_group))
    if str(self.pparams.target_unit_cell).lower() != 'none':
      uc = ' '.join(list(map(str, self.pparams.target_unit_cell.parameters())))
      advanced.uc.unit_cell.SetValue(uc)
    else:
      advanced.uc.unit_cell.SetValue(str(self.pparams.target_unit_cell))
    advanced.uc_override.SetValue(self.pparams.flag_override_unit_cell)
    advanced.anom.SetValue(self.pparams.target_anomalous_flag)
    advanced.cc.cc_cutoff.SetValue(str(self.pparams.frame_accept_min_cc))
    advanced.pix.pixel_size.SetValue(str(self.pparams.pixel_size_mm))
    advanced.cycles.ctr.SetValue(int(self.pparams.n_postref_cycle))

    if advanced.ShowModal() == wx.ID_OK:
      self.pparams.scale.d_max = float(advanced.res.high.GetValue())
      self.pparams.scale.d_min = float(advanced.res.low.GetValue())
      self.pparams.merge.d_max = float(advanced.res.high.GetValue())
      self.pparams.merge.d_min = float(advanced.res.low.GetValue())
      self.pparams.postref.scale.d_max = float(advanced.res.high.GetValue())
      self.pparams.postref.scale.d_min = float(advanced.res.low.GetValue())
      self.pparams.postref.crystal_orientation.d_max = \
        float(advanced.res.high.GetValue())
      self.pparams.postref.crystal_orientation.d_min = \
        float(advanced.res.low.GetValue())
      self.pparams.postref.reflecting_range.d_max = \
        float(advanced.res.high.GetValue())
      self.pparams.postref.reflecting_range.d_min = \
        float(advanced.res.low.GetValue())
      self.pparams.postref.unit_cell.d_max = float(advanced.res.high.GetValue())
      self.pparams.postref.unit_cell.d_min = float(advanced.res.low.GetValue())
      self.pparams.postref.allparams.d_max = float(advanced.res.high.GetValue())
      self.pparams.postref.allparams.d_min = float(advanced.res.low.GetValue())
      self.pparams.target_space_group = advanced.sg.spacegroup.GetValue()
      if advanced.uc.unit_cell.GetValue().lower() != 'none':
        uc = str_split(advanced.uc.unit_cell.GetValue())
        self.pparams.target_unit_cell = unit_cell(list(map(float, uc)))
      else:
        self.pparams.target_unit_cell = None
      self.pparams.flag_override_unit_cell = advanced.uc_override.GetValue()
      self.pparams.target_anomalous_flag = advanced.anom.GetValue()
      if advanced.cc.cc_cutoff.GetValue().lower() != 'none':
        self.pparams.frame_accept_min_cc = float(advanced.cc.cc_cutoff.GetValue())
      else:
        self.pparams.frame_accept_min_cc = None
      if advanced.pix.pixel_size.GetValue().lower() != 'none':
        self.pparams.pixel_size_mm = float(advanced.pix.pixel_size.GetValue())
      else:
        self.pparams.pixel_size_mm = None
      self.pparams.n_postref_cycle = int(advanced.cycles.ctr.GetValue())


    advanced.Destroy()
    e.Skip()


class PRIMEAdvancedOptions(wx.Dialog):
  ''' Advanced Options Dialog'''

  def __init__(self, *args, **kwargs):
    super(PRIMEAdvancedOptions, self).__init__(*args, **kwargs)

    main_sizer = wx.BoxSizer(wx.VERTICAL)

    main_box = wx.StaticBox(self, label='Advanced PRIME Options')
    vbox = wx.StaticBoxSizer(main_box, wx.VERTICAL)

    # Resolution
    self.res = ct.OptionCtrl(self,
                             label='Resolution: ',
                             label_size=(120, -1),
                             label_style='normal',
                             ctrl_size=wx.DefaultSize,
                             items=[('high', 50),
                                    ('low', 1.5)])
    vbox.Add(self.res, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    # Target space group
    self.sg = ct.OptionCtrl(self,
                            label='Space Group: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(100, -1),
                            items=[('spacegroup','P212121')])
    vbox.Add(self.sg, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    # Target unit cell
    self.uc = ct.OptionCtrl(self,
                            label='Unit Cell: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(300, -1),
                            items=[('unit_cell', '72 120 134 90 90 90')])
    self.uc_override = wx.CheckBox(self,
                                   label='Unit cell override')
    self.uc_override.SetValue(False)
    self.anom = wx.CheckBox(self, label='Anomalous')
    self.anom.SetValue(False)
    vbox.Add(self.uc, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)
    vbox.Add(self.uc_override, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)
    vbox.Add(self.anom, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    # CC cutoff
    self.cc = ct.OptionCtrl(self,
                            label='CC cutoff: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(100, -1),
                            items=[('cc_cutoff', 0.25)])
    vbox.Add(self.cc, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    # Pixel size
    self.pix = ct.OptionCtrl(self,
                             label='Pixel size: ',
                             label_size=(120, -1),
                             label_style='normal',
                             ctrl_size=(100, -1),
                             items=[('pixel_size', 0.172)])
    vbox.Add(self.pix, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    self.cycles = ct.SpinCtrl(self,
                              label='No. of Cycles:',
                              label_size=(120, -1),
                              label_style='normal',
                              ctrl_size=(60, -1))
    vbox.Add(self.cycles, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    # Dialog controls
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    main_sizer.Add(vbox, flag=wx.ALL, border=15)
    main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=15)

    self.SetSizer(main_sizer)

# -------------------------------- Threading --------------------------------- #

# Set up events for finishing one cycle and for finishing all cycles
tp_EVT_ALLDONE = wx.NewEventType()
EVT_ALLDONE = wx.PyEventBinder(tp_EVT_ALLDONE, 1)

class AllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class PRIMEThread(Thread):
  ''' Worker thread; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               prime_file,
               out_file,
               command=None):
    Thread.__init__(self)
    self.parent = parent
    self.prime_file = prime_file
    self.out_file = out_file
    self.command = command

  def run(self):
    if os.path.isfile(self.out_file):
      os.remove(self.out_file)
    if self.command is None:
      cmd = 'prime.run {}'.format(self.prime_file, self.out_file)
    else:
      cmd = self.command

    easy_run.fully_buffered(cmd, join_stdout_stderr=True)
    #evt = AllDone(tp_EVT_ALLDONE, -1)
    #wx.PostEvent(self.parent, evt)

# ----------------------------  Processing Window ---------------------------  #

class LogTab(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)

    self.log_sizer = wx.BoxSizer(wx.VERTICAL)
    self.log_window = wx.TextCtrl(self,
                                  style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_DONTWRAP)
    self.log_window.SetFont(wx.Font(9, wx.TELETYPE, wx.NORMAL, wx.NORMAL, False))
    self.log_sizer.Add(self.log_window, proportion=1, flag= wx.EXPAND | wx.ALL, border=10)
    self.SetSizer(self.log_sizer)

class RuntimeTab(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent)
    self.prime_sizer = wx.BoxSizer(wx.VERTICAL)
    self.prime_figure = Figure()
    self.prime_figure.patch.set_alpha(0)

    plt.rc('font', family='sans-serif', size=10)
    plt.rc('mathtext', default='regular')

    # Create nested GridSpec
    gsp = gridspec.GridSpec(2, 2, height_ratios=[2, 3])
    self.cc_axes = self.prime_figure.add_subplot(gsp[0])
    self.comp_axes = self.prime_figure.add_subplot(gsp[1])
    self.mult_axes = self.comp_axes.twinx()
    self.rej_table = self.prime_figure.add_subplot(gsp[2])

    gsub = gridspec.GridSpecFromSubplotSpec(3, 1,
                                            subplot_spec=gsp[3],
                                            hspace=0)
    self.bcc_axes = self.prime_figure.add_subplot(gsub[0])
    self.bcomp_axes = self.prime_figure.add_subplot(gsub[1],
                                                    sharex=self.bcc_axes)
    self.bmult_axes = self.prime_figure.add_subplot(gsub[2],
                                                    sharex=self.bcc_axes)
    self.draw_axes()

    # Set layout, create canvas, add to sizer
    self.prime_figure.set_tight_layout(True)
    self.canvas = FigureCanvas(self, -1, self.prime_figure)
    self.prime_sizer.Add(self.canvas, proportion=1, flag =wx.EXPAND)

    # Initialize main sizer
    self.SetSizer(self.prime_sizer)

  def draw_axes(self):
    self.rej_table.axis('off')
    self.cc_axes.set_title(r'$CC_{1/2}$', fontsize=12)
    self.cc_axes.set_xlabel('Cycle')
    self.cc_axes.set_ylabel(r'$CC_{1/2}$ (%)')
    self.comp_axes.set_title('Completeness / Multiplicity', fontsize=12)
    self.comp_axes.set_xlabel('Cycle')
    self.comp_axes.set_ylabel('Completeness (%)')
    self.mult_axes.set_ylabel('# of Observations')
    self.bcc_axes.yaxis.get_major_ticks()[0].label1.set_visible(False)
    self.bcc_axes.yaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.bcc_axes.set_ylabel(r'$CC_{1/2}$ (%)')
    plt.setp(self.bcc_axes.get_xticklabels(), visible=False)
    self.bcomp_axes.yaxis.get_major_ticks()[0].label1.set_visible(False)
    self.bcomp_axes.yaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.bcomp_axes.set_ylabel("Comp (%)")
    plt.setp(self.bcomp_axes.get_xticklabels(), visible=False)
    self.bmult_axes.yaxis.get_major_ticks()[0].label1.set_visible(False)
    self.bmult_axes.yaxis.get_major_ticks()[-1].label1.set_visible(False)
    self.bmult_axes.set_xlabel("Resolution ($\AA$)")
    self.bmult_axes.set_ylabel("# of Obs")

    self.prime_sizer.Layout()

  def draw_plots(self, info, total_cycles):

    # Plot mean CC1/2
    meanCC = info['total_cc12']
    cycles = range(len(meanCC))
    self.cc_axes.clear()
    self.cc_axes.plot(cycles, meanCC, 'o', c='#2b8cbe', ls='-', lw=3)
    self.cc_axes.set_xlim(0, total_cycles)

    # Plot mean completeness and multiplicity
    mean_comp = info['total_completeness']
    mean_mult = info['total_n_obs']
    cycles = range(len(mean_comp))
    self.comp_axes.clear()
    self.mult_axes.clear()
    self.comp_axes.set_xlim(0, total_cycles)
    self.mult_axes.set_xlim(0, total_cycles)
    self.comp_axes.plot(cycles, mean_comp, c='#f03b20', ls='-', lw=2)
    comp = self.comp_axes.scatter(cycles, mean_comp, marker='o', s=25,
                                  edgecolors='black', color='#f03b20')
    self.mult_axes.plot(cycles, mean_mult, c='#feb24c', ls='-', lw=2)
    mult = self.mult_axes.scatter(cycles, mean_mult, marker='o', s=25,
                                  edgecolors='black', color='#feb24c')
    labels = ['Completeness', 'Multiplicity']
    self.comp_axes.legend([comp, mult], labels, loc='upper right',
                          fontsize=9, fancybox=True)

    # Binned bar plots
    x = info['binned_resolution'][-1]
    bins = np.arange(len(x))
    xticks = bins[0::len(bins) // 6]
    xlabels = ["{:.2f}".format(i) for i in x]

    # plot binned stats
    self.bcc_axes.clear()
    self.bcc_axes.bar(bins, info['binned_cc12'][-1], color='#2b8cbe',
                      alpha=0.5, width=1, lw=0)
    self.bcc_axes.step(bins, info['binned_cc12'][-1], color='blue',
                       where='post')
    self.bcomp_axes.clear()
    self.bcomp_axes.bar(bins, info['binned_completeness'][-1],
                        alpha=0.5, color='#f03b20', width=1, lw=0)
    self.bcomp_axes.step(bins, info['binned_completeness'][-1], color='red',
                         where='post')
    self.bmult_axes.clear()
    self.bmult_axes.bar(bins, info['binned_n_obs'][-1],
                        alpha=0.5, color='#feb24c', width=1, lw=0)
    self.bmult_axes.step(bins, info['binned_n_obs'][-1], color='orange',
                         where='post')

    # Set x-axis tick labels
    self.bmult_axes.set_xticks(xticks)
    self.bmult_axes.set_xticklabels(xlabels)
    self.draw_axes()

    # Rejection table
    txt = 'No. good frames:           {}\n' \
          'No. bad CC frames:         {}\n' \
          'No. bad G frames:          {}\n' \
          'No. bad unit cell frames:  {}\n' \
          'No. bad gamma_e frames:    {}\n' \
          'No. bad SE frames:         {}\n' \
          'No. observations:          {}\n' \
          ''.format(info['n_frames_good'][-1],
                    info['n_frames_bad_cc'][-1],
                    info['n_frames_bad_G'][-1],
                    info['n_frames_bad_uc'][-1],
                    info['n_frames_bad_gamma_e'][-1],
                    info['n_frames_bad_SE'][-1],
                    info['n_observations'][-1])

    self.rej_table.clear()
    self.rej_table.axis('off')
    font = {'family': 'monospace',
            'color': 'darkblue',
            'weight': 'normal',
            'size': 13,
            'linespacing': 2.5
            }
    self.rej_table.text(0, 0.85, txt, fontdict=font,
                                  transform=self.rej_table.transAxes,
                                  va='top')

    # Redraw canvas
    self.canvas.draw()

class SummaryTab(wx.Panel):
  def __init__(self,
               parent,
               pparams,
               info):
    wx.Panel.__init__(self, parent)

    from prime.postrefine.mod_plotter import Plotter

    self.info = info
    self.pparams = pparams
    self.plot = Plotter(self.pparams, self.info)

    self.summary_sizer = wx.BoxSizer(wx.VERTICAL)

    sfont = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
    bfont = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.BOLD)
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
    self.summary_sizer.Add(run_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Table 1
    tb1_box = wx.StaticBox(self, label='Merging Statistics')
    tb1_box.SetFont(sfont)
    tb1_box_sizer = wx.StaticBoxSizer(tb1_box, wx.HORIZONTAL)

    self.tb1_labels, self.tb1_data = self.plot.table_one()
    self.tb1 = ct.TableCtrl(self,
                            rlabels=self.tb1_labels,
                            contents=self.tb1_data,
                            label_style='bold')

    # Buttons (placeholder for now)
    btn_box_sizer = wx.BoxSizer(wx.VERTICAL)
    self.btn_stats = ct.GradButton(self,
                                  bmp=wx.Bitmap('{}/24x24/line.png'.format(icons)),
                                  label=' Statistical charts', size=(250, -1))
    self.btn_table1 = ct.GradButton(self,
                                  bmp=wx.Bitmap('{}/24x24/txt.png'.format(icons)),
                                  label=' Output Table 1', size=(250, -1))
    btn_box_sizer.Add(self.btn_stats)
    self.Bind(wx.EVT_BUTTON, self.onPlotStats, self.btn_stats)
    btn_box_sizer.Add(self.btn_table1, flag=wx.TOP, border=5)
    self.Bind(wx.EVT_BUTTON, self.onWriteTableOne, self.btn_table1)

    tb1_box_sizer.Add(self.tb1, flag=wx.EXPAND | wx.ALL, border=10)
    tb1_box_sizer.AddStretchSpacer()
    tb1_box_sizer.Add(btn_box_sizer, flag=wx.ALIGN_RIGHT | wx.ALL, border=10)
    self.summary_sizer.Add(tb1_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    self.SetSizer(self.summary_sizer)

  def onPlotStats(self, e):
    self.plot.stat_charts()

  def onWriteTableOne(self, e):
    ''' Write Table 1 to a file '''

    save_dlg = wx.FileDialog(self,
                             message="Save Table 1",
                             defaultDir=os.curdir,
                             defaultFile="table_1.txt",
                             wildcard="*.txt",
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                             )
    if save_dlg.ShowModal() == wx.ID_OK:
      with open(save_dlg.GetPath(), 'a') as tb1_file:
       for i in range(len(self.tb1_data)):
          line = u'{:<25} {:<40}\n'.format(self.tb1_labels[i].decode('utf-8'),
                                           self.tb1_data[i][0]).encode('utf-8')
          tb1_file.write(line)

class PRIMERunWindow(wx.Frame):
  ''' New frame that will show processing info '''

  def __init__(self, parent, id, title, params,
               prime_file, out_file, mp_method='python', command=None):
    wx.Frame.__init__(self, parent, id, title, size=(800, 900),
                      style= wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER)

    self.logtext = ''
    self.pparams = params
    self.prime_file = prime_file
    self.out_file = os.path.join(self.pparams.run_no, 'log.txt')
    self.bookmark = 0
    self.prev_pids = []
    self.aborted = False
    self.command=command
    self.mp_method = mp_method
    self.current_cycle = -1

    self.main_panel = wx.Panel(self)
    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    # Toolbar
    self.prime_toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.tb_btn_abort = self.prime_toolbar.AddLabelTool(wx.ID_ANY, label='Abort',
                                                        bitmap=wx.Bitmap('{}/32x32/stop.png'.format(icons)),
                                                        shortHelp='Abort')
    self.prime_toolbar.Realize()

    # Status box
    self.status_panel = wx.Panel(self.main_panel)
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
    self.prime_panel = wx.Panel(self.main_panel)
    self.prime_nb = wx.Notebook(self.prime_panel, style=0)
    self.log_tab = LogTab(self.prime_nb)
    self.graph_tab = RuntimeTab(self.prime_nb)
    self.prime_nb.AddPage(self.log_tab, 'Log')
    self.prime_nb.AddPage(self.graph_tab, 'Charts')
    self.prime_nb.SetSelection(1)
    self.prime_sizer = wx.BoxSizer(wx.VERTICAL)
    self.prime_sizer.Add(self.prime_nb, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.prime_panel.SetSizer(self.prime_sizer)

    self.main_sizer.Add(self.status_panel, flag=wx.EXPAND | wx.ALL, border=3)
    self.main_sizer.Add(self.prime_panel, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.main_panel.SetSizer(self.main_sizer)

    #Processing status bar
    self.sb = self.CreateStatusBar()
    self.sb.SetFieldsCount(2)
    self.sb.SetStatusWidths([-1, -2])

    # Output gauge in status bar
    self.gauge_prime = wx.Gauge(self.sb, -1,
                                style=wx.GA_HORIZONTAL | wx.GA_SMOOTH)
    rect = self.sb.GetFieldRect(0)
    self.gauge_prime.SetPosition((rect.x + 2, rect.y + 2))
    self.gauge_prime.SetSize((rect.width - 4, rect.height - 4))
    self.gauge_prime.Hide()

    # Output polling timer
    self.timer = wx.Timer(self)

    # Event bindings
    self.Bind(EVT_ALLDONE, self.onFinishedProcess)
    self.sb.Bind(wx.EVT_SIZE, self.onStatusBarResize)
    self.Bind(wx.EVT_TIMER, self.onTimer, id=self.timer.GetId())

    # Button bindings
    self.Bind(wx.EVT_TOOL, self.onAbort, self.tb_btn_abort)

    self.run()

  def onStatusBarResize(self, e):
    rect = self.sb.GetFieldRect(0)
    self.gauge_prime.SetPosition((rect.x + 2, rect.y + 2))
    self.gauge_prime.SetSize((rect.width - 4, rect.height - 4))

  def onAbort(self, e):
    self.status_txt.SetForegroundColour('red')
    self.status_txt.SetLabel('Aborting...')
    self.prime_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)

    if self.mp_method == 'python':
      self.pids = easy_run.fully_buffered('pgrep -u {} {}'
                                          ''.format(user, python)).stdout_lines
      self.pids = [i for i in self.pids if i not in self.prev_pids]
      for i in self.pids:
        easy_run.fully_buffered('kill -9 {}'.format(i))
        print 'killing PID {}'.format(i)

    self.aborted = True

  def run(self):
    self.status_txt.SetForegroundColour('black')
    self.status_txt.SetLabel('Running...')
    self.gauge_prime.SetRange(self.pparams.n_postref_cycle)

    prime_process = PRIMEThread(self, self.prime_file, self.out_file,
                                command=self.command)
    prime_process.start()
    self.timer.Start(5000)


  def display_log(self):
    ''' Display PRIME stdout '''
    if os.path.isfile(self.out_file):
      with open(self.out_file, 'r') as out:
        out.seek(self.bookmark)
        output = out.readlines()
        self.bookmark = out.tell()

      ins_pt = self.log_tab.log_window.GetInsertionPoint()
      for i in output:
        self.log_tab.log_window.AppendText(i)
        self.log_tab.log_window.SetInsertionPoint(ins_pt)


  def plot_runtime_results(self):
    ''' Plot results for each cycle upon cycle completion '''

    stat_file = os.path.join(self.pparams.run_no, 'pickle.stat')
    if os.path.isfile(stat_file):
      info = ep.load(stat_file)
    else:
      info = {}

    if 'binned_resolution' in info:
      self.graph_tab.draw_plots(info, self.pparams.n_postref_cycle)
      self.current_cycle = len(info['total_cc12']) - 1

  def plot_final_results(self):
    ''' Plot final results '''

    stat_file = os.path.join(self.pparams.run_no, 'pickle.stat')
    if os.path.isfile(stat_file):
      info = ep.load(stat_file)
    else:
      info = {}

    self.summary_tab = SummaryTab(self.prime_nb,
                                  self.pparams,
                                  info)

    self.summary_tab.title_txt.SetLabel(self.pparams.title)
    self.summary_tab.folder_txt.SetLabel(self.pparams.run_no)


    # Display summary
    self.prime_nb.AddPage(self.summary_tab, 'Analysis')
    self.prime_nb.SetSelection(2)


  def onTimer(self, e):

    # Inspect output and update gauge
    self.gauge_prime.Show()
    if self.current_cycle == -1:
      self.sb.SetStatusText(('Merging...'), 1)
      self.gauge_prime.SetValue(0)
    else:
      self.sb.SetStatusText('Macrocycle {} of {} completed...' \
                            ''.format(self.current_cycle,
                                      self.pparams.n_postref_cycle), 1)
      self.gauge_prime.SetValue(self.current_cycle)

    # Plot runtime results
    self.plot_runtime_results()

    # Update log
    self.display_log()

    # Sense end of cycle
    if self.current_cycle >= self.pparams.n_postref_cycle:
      self.final_step()


  def onFinishedProcess(self, e):
    self.final_step()


  def final_step(self):
    font = self.status_txt.GetFont()
    font.SetWeight(wx.BOLD)

    # Check for aborted run
    if self.aborted:
      self.sb.SetStatusText(('Aborted by user'), 1)
      self.status_txt.SetForegroundColour('orange')
      self.status_txt.SetLabel('ABORTED BY USER')
      self.status_txt.SetFont(font)
    else:
      self.sb.SetStatusText(('Run finished'), 1)
      self.status_txt.SetForegroundColour('blue')
      self.status_txt.SetLabel('DONE')
      self.status_txt.SetFont(font)

      # Show final results
      self.plot_runtime_results()
      self.plot_final_results()

    # Finish up
    self.display_log()
    self.gauge_prime.Hide()
    self.prime_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)
    self.timer.Stop()

# ---------------------------------  Dialogs --------------------------------  #

class PRIMEPreferences(wx.Dialog):
  def __init__(self, *args, **kwargs):
    super(PRIMEPreferences, self).__init__(*args, **kwargs)

    self.method = None
    self.queue = None
    main_sizer = wx.BoxSizer(wx.VERTICAL)

    main_box = wx.StaticBox(self, label='IOTA Preferences')
    vbox = wx.StaticBoxSizer(main_box, wx.VERTICAL)

    self.SetSizer(main_sizer)

    q_choices = ['psanaq', 'psnehq', 'psfehq'] + ['custom']
    self.queues = ct.ChoiceCtrl(self,
                                label='Queue:',
                                label_size=(120, -1),
                                label_style='bold',
                                ctrl_size=wx.DefaultSize,
                                choices=q_choices)
    vbox.Add(self.queues, flag=wx.ALL, border=10)

    self.custom_queue = ct.OptionCtrl(self,
                                      items=[('cqueue', '')],
                                      label='Custom Queue:',
                                      label_size=(120, -1),
                                      label_style='normal',
                                      ctrl_size=(150, -1))
    self.custom_queue.Disable()
    vbox.Add(self.custom_queue, flag=wx.ALL, border=10)

    mp_choices = ['multiprocessing', 'bsub']
    self.mp_methods = ct.ChoiceCtrl(self,
                                    label='Method:',
                                    label_size=(120, -1),
                                    label_style='bold',
                                    ctrl_size=wx.DefaultSize,
                                    choices=mp_choices)
    vbox.Add(self.mp_methods, flag=wx.ALL, border=10)

    main_sizer.Add(vbox, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    self.Bind(wx.EVT_CHOICE, self.onQueue, self.queues.ctr)
    self.Bind(wx.EVT_CHOICE, self.onMethod, self.mp_methods.ctr)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def set_choices(self, method, queue):
    # Set queue to default value
    if queue is None:
      queue = 'None'
    inp_queue = self.queues.ctr.FindString(queue)
    if inp_queue != wx.NOT_FOUND:
      self.queues.ctr.SetSelection(inp_queue)
    else:
      self.custom_queue.Enable()
      self.custom_queue.cqueue.SetValue(queue)

    # Set method to default value
    print "method: ", method
    inp_method = self.mp_methods.ctr.FindString(str(method))
    if inp_method != wx.NOT_FOUND:
      self.mp_methods.ctr.SetSelection(inp_method)

    self.check_method()

  def onMethod(self, e):
    self.check_method()

  def check_method(self):
    choice = self.mp_methods.ctr.GetString(self.mp_methods.ctr.GetSelection())
    if choice == 'multiprocessing':
      self.queues.Disable()
      self.custom_queue.Disable()
    else:
      self.queues.Enable()
      queue = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
      if queue == 'custom':
        self.custom_queue.Enable()

  def onQueue(self, e):
    choice = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
    if choice == 'custom':
      self.custom_queue.Enable()
    else:
      self.custom_queue.Disable()

  def onOK(self, e):
    self.method = self.mp_methods.ctr.GetString(self.mp_methods.ctr.GetSelection())
    queue_selection = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
    if queue_selection == 'custom':
      if self.custom_queue.cqueue.GetValue() == '':
        wx.MessageBox('Please choose or enter a queue', wx.OK)
      else:
        self.queue = self.custom_queue.cqueue.GetValue()
        e.Skip()
    else:
      self.queue = queue_selection
      e.Skip()
